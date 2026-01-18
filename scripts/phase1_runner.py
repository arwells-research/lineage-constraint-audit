#!/usr/bin/env python3
"""
phase1_runner.py — Phase-1/2/3 GO/NO-GO test + Phase-4 localization report
for Σ₂-style "state sufficiency failure" on GEO GSE126954.

Phase-4 additions:
  --report_by {none,timebin,fate}
    * timebin: per embryo.time.bin subgroup stats + permutation p-values
    * fate: per cell.type subgroup stats + permutation p-values (top N by count)

  Group stats are computed using the SAME null permutations already run for global p-values.
  No new modeling knobs are introduced.

Key outputs:
  - v_obs: fraction of centers whose neighborhood is validity-positive (>=2 founders)
  - A_obs: mean G-test statistic across validity-positive centers only
  - Null permutations: shuffle founder within strata (timebin or timebin+batch)
  - Null reference: 95th percentile of the stratified permutation null
  - p-values: P(null >= observed) reported for diagnostics only

Controls:
  - control_mode can shuffle fate within strata (negative control)

GO/NO-GO:
  Stage-1A: Viability gates
    - Requires at least one validity-positive neighborhood
    - Requires >=100 validity-positive centers
    - Requires successful null permutations

  Stage-1B: Detectability GO if:
    (v_obs > v_null_95) OR (A_obs > A_null_95)

  Stage-1C: Strength grading is descriptive only and always reported.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors


def founder_group(lineage: object) -> Optional[str]:
    if not isinstance(lineage, str):
        return None
    s = lineage.strip()
    if not s or s.upper() == "NA":
        return None
    if s.startswith("AB"):
        return "AB"
    if s.startswith("MS"):
        return "MS"
    if s.startswith("C"):
        return "C"
    if s.startswith("D"):
        return "D"
    if s.startswith("E"):
        return "E"
    return None


def compute_g_stat(y: np.ndarray, g: np.ndarray) -> float:
    """G-test (likelihood ratio) for independence between categorical y and g."""
    uy, yinv = np.unique(y, return_inverse=True)
    ug, ginv = np.unique(g, return_inverse=True)
    ny = uy.size
    ng = ug.size
    obs = np.zeros((ny, ng), dtype=np.int32)
    np.add.at(obs, (yinv, ginv), 1)

    n = int(obs.sum())
    if n <= 0:
        return 0.0

    row = obs.sum(axis=1).astype(np.float64)
    col = obs.sum(axis=0).astype(np.float64)
    exp = (row[:, None] * col[None, :]) / float(n)

    mask = obs > 0
    obs_f = obs[mask].astype(np.float64)
    exp_f = np.maximum(exp[mask], 1e-300)
    return 2.0 * float(np.sum(obs_f * np.log(obs_f / exp_f)))


def l2_normalize_rows(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    norms = np.sqrt(np.sum(x * x, axis=1, keepdims=True))
    norms = np.maximum(norms, eps)
    return x / norms


@dataclass
class PhaseCohort:
    df: pd.DataFrame
    col_to_new: Dict[int, int]
    fate: np.ndarray
    founder: np.ndarray
    timebin: np.ndarray
    batch: np.ndarray
    cell_ids: np.ndarray
    fate_names: np.ndarray
    timebin_names: np.ndarray
    batch_names: np.ndarray


def load_and_filter_annotations(path_gz: str) -> PhaseCohort:
    df = pd.read_csv(path_gz, compression="gzip")
    if len(df) <= 0:
        raise RuntimeError("Annotation file appears empty.")

    # 1-based MatrixMarket column index equals original annotation row order
    df = df.copy()
    df["col_index_1based"] = np.arange(1, len(df) + 1, dtype=np.int32)

    qc_col = "passed_initial_QC_or_later_whitelisted"
    if qc_col in df.columns:
        df = df[df[qc_col] == True].copy()

    required = ("cell.type", "lineage", "embryo.time.bin", "cell")
    for col in required:
        if col not in df.columns:
            raise RuntimeError(f"Missing required column: {col}")

    # batch is optional for Phase-3/4, but we still try to use it if present
    has_batch = "batch" in df.columns
    if not has_batch:
        df["batch"] = "NO_BATCH_FIELD"

    df = df[df["cell.type"].notna()].copy()
    df = df[df["lineage"].notna()].copy()

    df["lineage"] = df["lineage"].astype(str)
    df = df[~df["lineage"].str.contains("/", regex=False)].copy()

    df["founder"] = df["lineage"].map(founder_group)
    df = df[df["founder"].notna()].copy()

    fate_codes, fate_uniques = pd.factorize(df["cell.type"], sort=True)
    founder_codes, founder_uniques = pd.factorize(df["founder"], sort=True)
    time_codes, time_uniques = pd.factorize(df["embryo.time.bin"], sort=True)
    batch_codes, batch_uniques = pd.factorize(df["batch"].astype(str), sort=True)

    df["fate_code"] = fate_codes.astype(np.int32)
    df["founder_code"] = founder_codes.astype(np.int32)
    df["time_code"] = time_codes.astype(np.int32)
    df["batch_code"] = batch_codes.astype(np.int32)

    col_indices = df["col_index_1based"].to_numpy(dtype=np.int32)
    col_to_new = {int(c): i for i, c in enumerate(col_indices)}

    cohort = PhaseCohort(
        df=df.reset_index(drop=True),
        col_to_new=col_to_new,
        fate=df["fate_code"].to_numpy(dtype=np.int32),
        founder=df["founder_code"].to_numpy(dtype=np.int32),
        timebin=df["time_code"].to_numpy(dtype=np.int32),
        batch=df["batch_code"].to_numpy(dtype=np.int32),
        cell_ids=df["cell"].to_numpy(dtype=str),
        fate_names=np.asarray(fate_uniques, dtype=str),
        timebin_names=np.asarray(time_uniques, dtype=str),
        batch_names=np.asarray(batch_uniques, dtype=str),
    )

    print("=== Phase cohort ===")
    print(f"Annotation rows (after filters): {len(cohort.df)}")
    print("Founder counts:")
    print(cohort.df["founder"].value_counts().to_string())
    print("Top fate counts:")
    print(cohort.df["cell.type"].value_counts().head(15).to_string())
    print("Time-bin counts:")
    print(cohort.df["embryo.time.bin"].value_counts().sort_index().to_string())
    print("Batch counts:")
    print(cohort.df["batch"].value_counts().head(20).to_string())
    print()
    return cohort


def stream_read_mtx_subset_cells(
    mtx_gz: str,
    n_genes: int,
    n_cells_total: int,
    col_to_new: Dict[int, int],
    max_nnz: Optional[int] = None,
) -> sp.csr_matrix:
    n_filtered = len(col_to_new)
    rows: List[int] = []
    cols: List[int] = []
    data: List[int] = []

    kept = 0
    with gzip.open(mtx_gz, "rt") as f:
        header1 = f.readline()
        if not header1.startswith("%%MatrixMarket"):
            raise RuntimeError("Unexpected MatrixMarket header.")

        dim_line = f.readline()
        while dim_line.startswith("%"):
            dim_line = f.readline()
        parts = dim_line.strip().split()
        if len(parts) != 3:
            raise RuntimeError("Unexpected MatrixMarket dimension line.")
        nrows_file, ncols_file, _nnz_file = map(int, parts)
        if nrows_file != n_genes or ncols_file != n_cells_total:
            raise RuntimeError(
                f"Matrix dims mismatch: file has {nrows_file}x{ncols_file}, expected {n_genes}x{n_cells_total}"
            )

        for line in f:
            if not line or line.startswith("%"):
                continue
            r_s, c_s, v_s = line.strip().split()
            c = int(c_s)
            new_i = col_to_new.get(c)
            if new_i is None:
                continue
            r = int(r_s) - 1
            v = int(v_s)
            rows.append(new_i)
            cols.append(r)
            data.append(v)
            kept += 1
            if max_nnz is not None and kept >= max_nnz:
                break

    if kept == 0:
        raise RuntimeError("No nonzero entries retained for filtered cohort; mapping may be wrong.")

    X = sp.coo_matrix(
        (np.asarray(data, dtype=np.float32),
         (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32))),
        shape=(n_filtered, n_genes),
    ).tocsr()

    print("=== Expression matrix (subset) ===")
    print(f"Filtered cells:   {n_filtered}")
    print(f"Genes:            {n_genes}")
    print(f"Kept nnz:         {X.nnz}")
    print(f"Sparsity:         {1.0 - (X.nnz / float(n_filtered * n_genes)):.6f}")
    print()
    return X


def normalize_log1p_cpm(X: sp.csr_matrix, target_sum: float = 1e4) -> sp.csr_matrix:
    X = X.tocsr(copy=True)
    lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    lib = np.maximum(lib, 1.0)
    scale = (target_sum / lib).astype(np.float32)
    X = X.multiply(scale[:, None])
    X.data = np.log1p(X.data)
    return X


def build_knn(Z: np.ndarray, k: int) -> np.ndarray:
    Z = l2_normalize_rows(Z)
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine")
    nn.fit(Z)
    return nn.kneighbors(return_distance=False).astype(np.int32)


def choose_centers(n: int, max_centers: int, seed: int) -> np.ndarray:
    if max_centers <= 0 or max_centers >= n:
        return np.arange(n, dtype=np.int32)
    rng = np.random.default_rng(seed)
    return rng.choice(n, size=max_centers, replace=False).astype(np.int32)


def _shuffle_within_groups(labels: np.ndarray, group_keys: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = labels.copy()
    for g in np.unique(group_keys):
        idx = np.where(group_keys == g)[0]
        if idx.size <= 1:
            continue
        out[idx] = rng.permutation(out[idx])
    return out


def make_group_keys(timebin: np.ndarray, batch: np.ndarray, mode: str) -> np.ndarray:
    if mode == "timebin":
        return timebin
    if mode == "timebin_batch":
        # Stable, future-proof factorization via MultiIndex
        mi = pd.MultiIndex.from_arrays([timebin, batch], names=["time", "batch"])
        keys = mi.factorize(sort=True)[0].astype(np.int32)
        return keys
    raise ValueError(f"Unknown null_mode: {mode}")


def neighborhood_validity_and_T(
    neighbor_idx: np.ndarray,
    centers: np.ndarray,
    fate: np.ndarray,
    founder: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    m = centers.size
    valid = np.zeros(m, dtype=bool)
    T = np.zeros(m, dtype=np.float64)

    for j, i in enumerate(centers):
        nb = neighbor_idx[i]
        y = fate[nb]
        g = founder[nb]
        if np.unique(y).size < 2:
            continue
        if np.unique(g).size < 2:
            continue
        valid[j] = True
        T[j] = compute_g_stat(y, g)

    return valid, T


def summarize_arr(name: str, x: np.ndarray) -> None:
    print(
        f"{name}: min={float(np.min(x)):.6f}  median={float(np.median(x)):.6f}  "
        f"p90={float(np.quantile(x, 0.90)):.6f}  "
        f"p95={float(np.quantile(x, 0.95)):.6f}  "
        f"max={float(np.max(x)):.6f}"
    )


def _p_ge(null: np.ndarray, obs: float) -> float:
    # +1 correction
    B = int(null.size)
    return float((1 + np.sum(null >= obs)) / (1 + B))


def phase4_report(
    *,
    report_by: str,
    report_min_valid: int,
    report_top_fates: int,
    cohort: PhaseCohort,
    centers: np.ndarray,
    neighbor_idx: np.ndarray,
    fate_obs: np.ndarray,
    founder_obs: np.ndarray,
    timebin_obs: np.ndarray,
    batch_obs: np.ndarray,
    null_mode: str,
    perms_founder: List[np.ndarray],   # list of shuffled founder arrays, length B
) -> None:
    if report_by == "none":
        return

    # Map each center to its group label (in the observed arrays)
    if report_by == "timebin":
        grp_codes = timebin_obs[centers]
        grp_names = cohort.timebin_names
        title = "=== Phase-4 report: by embryo.time.bin (centers) ==="
        # keep all bins
        keep_codes = np.unique(grp_codes)
    elif report_by == "fate":
        grp_codes = fate_obs[centers]
        grp_names = cohort.fate_names
        title = "=== Phase-4 report: by cell.type (centers; top-N) ==="
        # only top-N fates by count among centers
        vc = pd.Series(grp_codes).value_counts()
        keep_codes = vc.index.to_numpy(dtype=np.int32)[: int(report_top_fates)]
    else:
        raise ValueError(f"Unknown report_by: {report_by}")

    print(title)

    rows = []
    B = len(perms_founder)

    # Precompute observed validity/T for all centers once, then slice per group.
    valid_all, T_all = neighborhood_validity_and_T(neighbor_idx, centers, fate_obs, founder_obs)

    for c in keep_codes:
        mask = (grp_codes == c)
        if not np.any(mask):
            continue
        centers_g = centers[mask]
        valid_g = valid_all[mask]
        T_g = T_all[mask]

        n_centers = int(mask.sum())
        n_valid = int(valid_g.sum())
        v_obs_g = float(n_valid / n_centers) if n_centers > 0 else 0.0
        if n_valid == 0:
            # if no valid centers, report but no p-values
            rows.append(
                dict(
                    group=str(grp_names[int(c)]) if int(c) < len(grp_names) else str(c),
                    n_centers=n_centers,
                    n_valid=n_valid,
                    v_obs=v_obs_g,
                    A_obs=float("nan"),
                    p_v=float("nan"),
                    p_A=float("nan"),
                )
            )
            continue

        A_obs_g = float(np.mean(T_g[valid_g]))

        # Null arrays for this group, computed from the same permutations
        v_null_g = np.zeros(B, dtype=np.float64)
        A_null_g = np.zeros(B, dtype=np.float64)

        for b, f_sh in enumerate(perms_founder):
            valid_b, T_b = neighborhood_validity_and_T(neighbor_idx, centers_g, fate_obs, f_sh)
            v_null_g[b] = float(np.mean(valid_b))
            A_null_g[b] = float(np.mean(T_b[valid_b])) if np.any(valid_b) else 0.0

        p_v_g = _p_ge(v_null_g, v_obs_g)
        p_A_g = _p_ge(A_null_g, A_obs_g)

        # apply min-valid filter *after* computing, so you still see counts
        rows.append(
            dict(
                group=str(grp_names[int(c)]) if int(c) < len(grp_names) else str(c),
                n_centers=n_centers,
                n_valid=n_valid,
                v_obs=v_obs_g,
                A_obs=A_obs_g,
                p_v=p_v_g,
                p_A=p_A_g,
            )
        )

    out = pd.DataFrame(rows)

    # filter tiny-valid groups for readability, but keep them if user wants
    if report_min_valid > 0:
        out = out[(out["n_valid"].fillna(0).astype(int) >= int(report_min_valid)) | (out["n_valid"].isna())]

    # sort by smallest p_A then largest A_obs
    if "p_A" in out.columns:
        out = out.sort_values(by=["p_A", "A_obs"], ascending=[True, False], na_position="last")

    # Print in a stable plain-text table
    if len(out) == 0:
        print("(no groups met reporting criteria)")
        print()
        return

    with pd.option_context("display.max_rows", 200, "display.max_colwidth", 120):
        print(out.to_string(index=False, justify="left", formatters={
            "v_obs": lambda x: f"{x:.6f}" if pd.notna(x) else "nan",
            "A_obs": lambda x: f"{x:.6f}" if pd.notna(x) else "nan",
            "p_v":  lambda x: f"{x:.6g}" if pd.notna(x) else "nan",
            "p_A":  lambda x: f"{x:.6g}" if pd.notna(x) else "nan",
        }))
    print()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot", default="GSE126954_cell_annotation.csv.gz")
    ap.add_argument("--mtx", default="GSE126954_gene_by_cell_count_matrix.txt.gz")
    ap.add_argument("--n_genes", type=int, default=20222)
    ap.add_argument("--n_cells_total", type=int, default=89701)
    ap.add_argument("--svd_components", type=int, default=50)
    ap.add_argument("--knn_k", type=int, default=20)
    ap.add_argument("--max_centers", type=int, default=5000, help="0 means all")
    ap.add_argument("--global_perms", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max_nnz", type=int, default=0, help="DEBUG cap nnz read (0 = no cap)")

    ap.add_argument(
        "--null_mode",
        choices=["timebin", "timebin_batch"],
        default="timebin",
        help="Stratification for founder shuffle in the null.",
    )
    ap.add_argument(
        "--control_mode",
        choices=["none", "fate_timebin", "fate_timebin_batch"],
        default="none",
        help="Negative control: shuffle fate within strata (founder fixed).",
    )

    # Phase-4 reporting
    ap.add_argument(
        "--report_by",
        choices=["none", "timebin", "fate"],
        default="none",
        help="Phase-4 localization report by group.",
    )
    ap.add_argument("--report_min_valid", type=int, default=30, help="Min valid centers per group to show.")
    ap.add_argument("--report_top_fates", type=int, default=20, help="Top-N fates (by center count) for fate report.")

    args = ap.parse_args()

    if not os.path.exists(args.annot):
        print(f"ERROR: annotation not found: {args.annot}", file=sys.stderr)
        return 2
    if not os.path.exists(args.mtx):
        print(f"ERROR: matrix not found: {args.mtx}", file=sys.stderr)
        return 2

    cohort = load_and_filter_annotations(args.annot)

    max_nnz = None if args.max_nnz <= 0 else int(args.max_nnz)
    X_raw = stream_read_mtx_subset_cells(
        mtx_gz=args.mtx,
        n_genes=int(args.n_genes),
        n_cells_total=int(args.n_cells_total),
        col_to_new=cohort.col_to_new,
        max_nnz=max_nnz,
    )
    X = normalize_log1p_cpm(X_raw, target_sum=1e4)

    print("=== Embedding ===")
    print(f"TruncatedSVD components: {args.svd_components}")
    svd = TruncatedSVD(n_components=int(args.svd_components), random_state=int(args.seed))
    Z = svd.fit_transform(X)
    explained = float(np.sum(svd.explained_variance_ratio_))
    print(f"Explained variance ratio (sum): {explained:.4f}")
    print()

    print("=== Neighborhoods ===")
    print(f"k (excluding self): {args.knn_k}")
    neighbor_idx = build_knn(Z, k=int(args.knn_k))
    n = neighbor_idx.shape[0]
    centers = choose_centers(n, max_centers=int(args.max_centers), seed=int(args.seed))
    print(f"Filtered cells: {n}")
    print(f"Centers tested: {centers.size}")
    print()

    # Determine which fate/founder arrays are used in "observed"
    if args.control_mode == "none":
        fate_obs = cohort.fate
        founder_obs = cohort.founder
        control_desc = "none"
    else:
        if args.control_mode == "fate_timebin":
            gkeys = make_group_keys(cohort.timebin, cohort.batch, mode="timebin")
            fate_obs = _shuffle_within_groups(cohort.fate, gkeys, seed=int(args.seed) + 777)
            founder_obs = cohort.founder
            control_desc = "shuffle fate within embryo.time.bin (founder fixed)"
        elif args.control_mode == "fate_timebin_batch":
            gkeys = make_group_keys(cohort.timebin, cohort.batch, mode="timebin_batch")
            fate_obs = _shuffle_within_groups(cohort.fate, gkeys, seed=int(args.seed) + 777)
            founder_obs = cohort.founder
            control_desc = "shuffle fate within (embryo.time.bin, batch) (founder fixed)"
        else:
            raise RuntimeError("Unreachable control_mode")

    time_obs = cohort.timebin
    batch_obs = cohort.batch

    print("=== Config ===")
    print(f"null_mode:    {args.null_mode}")
    print(f"control_mode: {args.control_mode} ({control_desc})")
    print()

    # Observed
    valid_obs, T_obs = neighborhood_validity_and_T(neighbor_idx, centers, fate_obs, founder_obs)
    n_valid_centers = int(valid_obs.sum())
    n_centers = int(valid_obs.size)    
    v_obs = float(np.mean(valid_obs))
    if np.any(valid_obs):
        A_obs = float(np.mean(T_obs[valid_obs]))
        T_valid_obs = T_obs[valid_obs]
    else:
        A_obs = float("nan")
        T_valid_obs = np.array([], dtype=np.float64)

    print("=== Observed ===")
    print(f"Validity rate v_obs: {v_obs:.6f}  (valid centers: {n_valid_centers}/{n_centers})")
    if np.any(valid_obs):
        print(f"Conditional dependence A_obs = mean(T_i | valid): {A_obs:.6f}")
        summarize_arr("T_i (valid only)", T_valid_obs)
    else:
        print("No valid neighborhoods under the validity rule; NO-GO by design.")
        print()
        print("=== GO / NO-GO ===")
        print("Decision: NO-GO (viability V1 failed: no valid neighborhoods)")
        print()
        return 0
    print()

    # Stage-1A viability gate V2 (normative; see docs/audit_logic.md)
    if n_valid_centers < 100:
        print("=== GO / NO-GO ===")
        print("Decision: NO-GO (viability V2 failed: n_valid_centers < 100)")
        print(f"n_valid_centers={n_valid_centers} (required >= 100)")
        print()
        return 0

    # Null: shuffle founder within selected strata
    print("=== Global null (shuffle founder within strata) ===")
    B = int(args.global_perms)
    v_null = np.zeros(B, dtype=np.float64)
    A_null = np.zeros(B, dtype=np.float64)

    gkeys_null = make_group_keys(cohort.timebin, cohort.batch, mode=args.null_mode)

    # Store the shuffled founders for Phase-4 reuse (cheap: int32 vectors)
    perms_founder: List[np.ndarray] = []

    for b in range(B):
        f_sh = _shuffle_within_groups(cohort.founder, gkeys_null, seed=int(args.seed) + 1000 + b)
        perms_founder.append(f_sh)

        valid_b, T_b = neighborhood_validity_and_T(neighbor_idx, centers, fate_obs, f_sh)
        v_null[b] = float(np.mean(valid_b))
        A_null[b] = float(np.mean(T_b[valid_b])) if np.any(valid_b) else 0.0

        if (b + 1) % max(1, B // 10) == 0:
            print(f"  perm {b+1:>4}/{B}  v_null={v_null[b]:.6f}  A_null={A_null[b]:.6f}")

    p_v = _p_ge(v_null, v_obs)
    p_A = _p_ge(A_null, A_obs)

    print()
    print("=== Global p-values ===")
    print(f"p_v (v_null >= v_obs): {p_v:.6g}")
    print(f"p_A (A_null >= A_obs): {p_A:.6g}")
    summarize_arr("v_null", v_null)
    summarize_arr("A_null", A_null)
    print()

    # GO / NO-GO
    alpha = 0.01
    go = (p_A <= alpha) or (p_v <= alpha)
    print("=== GO / NO-GO ===")
    print(f"Rule: GO if (p_A <= {alpha}) OR (p_v <= {alpha})")
    print(f"Decision: {'GO' if go else 'NO-GO'}")
    print()

    # Phase-4 report (optional)
    phase4_report(
        report_by=str(args.report_by),
        report_min_valid=int(args.report_min_valid),
        report_top_fates=int(args.report_top_fates),
        cohort=cohort,
        centers=centers,
        neighbor_idx=neighbor_idx,
        fate_obs=fate_obs,
        founder_obs=founder_obs,
        timebin_obs=time_obs,
        batch_obs=batch_obs,
        null_mode=str(args.null_mode),
        perms_founder=perms_founder,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
