#!/usr/bin/env python3
"""
phaseD_founder_dispersion_gse126954.py — Phase-D: Founder dispersion / transition admissibility (GSE126954)

Phase-D question (worm / GSE126954):
  Do different founders occupy distinct admissible developmental "transition regions",
  beyond contribution imbalance?

Operationalization (worm):
  - Unit: validity-positive neighborhood centers (same validity rule as Phase-1 runner)
      validity := neighborhood has >=2 founders AND >=2 fates
  - Development axis: embryo.time.bin
  - Dispersion metric: Shannon entropy (bits) of embryo.time.bin within neighborhood
  - Condition by founder: founder label of the center cell
  - Overall test statistic:
      For each founder f with n_units >= min_units_per_founder:
         m_f = median(time_entropy_bits | center-founder=f)
      S_obs = Var_f(m_f)

Null model:
  - Shuffle founder labels within strata (timebin or timebin+batch), preserving geometry and time/batch structure
  - Recompute valid centers and time entropy per valid neighborhood
  - Recompute S_null

Viability gates (Phase-D):
  - D1: n_valid_centers >= min_units_total
  - D2: at least 2 founders have n_units >= min_units_per_founder
  - D3: perms completed; enough finite null draws

Inputs:
  --annot: data/raw/GSE126954_cell_annotation.csv.gz
  --mtx:   data/raw/GSE126954_gene_by_cell_count_matrix.txt.gz

Outputs (out_dir):
  phaseD_summary.txt
  unit_dispersion_table.tsv
  dispersion_by_founder.tsv
  null_dist_overall.tsv
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


def l2_normalize_rows(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    norms = np.sqrt(np.sum(x * x, axis=1, keepdims=True))
    norms = np.maximum(norms, eps)
    return x / norms


def _entropy_bits_from_counts(counts: np.ndarray) -> float:
    total = float(np.sum(counts))
    if total <= 0.0:
        return 0.0
    p = counts.astype(np.float64) / total
    p = p[p > 0.0]
    if p.size == 0:
        return 0.0
    return float(-np.sum(p * np.log2(p)))


def _p_ge(null: np.ndarray, obs: float) -> float:
    B = int(null.size)
    return float((1 + np.sum(null >= obs)) / (1 + B))


@dataclass
class Cohort:
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
    founder_names: np.ndarray


def load_and_filter_annotations(path_gz: str) -> Cohort:
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

    # batch optional
    if "batch" not in df.columns:
        df["batch"] = "NO_BATCH_FIELD"

    df = df[df["cell.type"].notna()].copy()
    df = df[df["lineage"].notna()].copy()

    df["lineage"] = df["lineage"].astype(str)
    # exclude ambiguous lineage strings
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

    cohort = Cohort(
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
        founder_names=np.asarray(founder_uniques, dtype=str),
    )

    print("=== Cohort (filters applied) ===")
    print(f"rows={len(cohort.df)}")
    print("Founder counts:")
    print(cohort.df["founder"].value_counts().to_string())
    print("Time-bin counts:")
    print(cohort.df["embryo.time.bin"].value_counts().sort_index().to_string())
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
    print(f"filtered_cells={n_filtered}")
    print(f"genes={n_genes}")
    print(f"nnz={X.nnz}")
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


def make_group_keys(timebin: np.ndarray, batch: np.ndarray, mode: str) -> np.ndarray:
    if mode == "timebin":
        return timebin
    if mode == "timebin_batch":
        mi = pd.MultiIndex.from_arrays([timebin, batch], names=["time", "batch"])
        keys = mi.factorize(sort=True)[0].astype(np.int32)
        return keys
    raise ValueError(f"Unknown null_mode: {mode}")


def shuffle_within_groups(labels: np.ndarray, group_keys: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = labels.copy()
    for g in np.unique(group_keys):
        idx = np.where(group_keys == g)[0]
        if idx.size <= 1:
            continue
        out[idx] = rng.permutation(out[idx])
    return out


def validity_mask_for_centers(
    neighbor_idx: np.ndarray,
    centers: np.ndarray,
    fate: np.ndarray,
    founder: np.ndarray,
) -> np.ndarray:
    m = centers.size
    valid = np.zeros(m, dtype=bool)
    for j, i in enumerate(centers):
        nb = neighbor_idx[i]
        y = fate[nb]
        g = founder[nb]
        if np.unique(y).size < 2:
            continue
        if np.unique(g).size < 2:
            continue
        valid[j] = True
    return valid


def time_entropy_bits_for_centers(
    neighbor_idx: np.ndarray,
    centers: np.ndarray,
    timebin: np.ndarray,
) -> np.ndarray:
    m = centers.size
    H = np.zeros(m, dtype=np.float64)
    for j, i in enumerate(centers):
        nb = neighbor_idx[i]
        t = timebin[nb]
        vc = np.bincount(t.astype(np.int32))
        H[j] = _entropy_bits_from_counts(vc)
    return H


def compute_founder_medians_and_S(
    founder_center: np.ndarray,
    time_entropy_bits: np.ndarray,
    valid: np.ndarray,
    founder_names: np.ndarray,
    min_units_per_founder: int,
) -> Tuple[pd.DataFrame, float, int]:
    rows = []
    eligible_meds = []

    # Only valid units participate
    fc = founder_center[valid]
    H = time_entropy_bits[valid]

    for code in np.unique(fc):
        code = int(code)
        mask = (fc == code)
        n = int(np.sum(mask))
        med = float(np.median(H[mask])) if n > 0 else float("nan")
        mean = float(np.mean(H[mask])) if n > 0 else float("nan")
        name = str(founder_names[code]) if code < len(founder_names) else str(code)
        rows.append(dict(founder=name, founder_code=code, n_units=n, median_time_entropy_bits=med, mean_time_entropy_bits=mean))
        if n >= int(min_units_per_founder):
            eligible_meds.append(med)

    summ = pd.DataFrame(rows).sort_values(by=["n_units", "founder"], ascending=[False, True])
    n_eligible = int(len(eligible_meds))
    if n_eligible >= 2:
        S = float(np.var(np.asarray(eligible_meds, dtype=np.float64), ddof=0))
    else:
        S = float("nan")
    return summ, S, n_eligible


def write_summary(path: str, lines: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for ln in lines:
            f.write(ln.rstrip() + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot", required=True)
    ap.add_argument("--mtx", required=True)
    ap.add_argument("--out_dir", required=True)

    ap.add_argument("--n_genes", type=int, default=20222)
    ap.add_argument("--n_cells_total", type=int, default=89701)
    ap.add_argument("--svd_components", type=int, default=50)
    ap.add_argument("--knn_k", type=int, default=20)
    ap.add_argument("--max_centers", type=int, default=5000)
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)

    ap.add_argument("--null_mode", choices=["timebin", "timebin_batch"], default="timebin_batch")

    # Phase-D viability gates
    ap.add_argument("--min_units_total", type=int, default=100)
    ap.add_argument("--min_units_per_founder", type=int, default=30)

    args = ap.parse_args()

    if not os.path.exists(args.annot):
        print(f"ERROR: annot not found: {args.annot}", file=sys.stderr)
        return 2
    if not os.path.exists(args.mtx):
        print(f"ERROR: mtx not found: {args.mtx}", file=sys.stderr)
        return 2

    os.makedirs(args.out_dir, exist_ok=True)

    cohort = load_and_filter_annotations(args.annot)

    X_raw = stream_read_mtx_subset_cells(
        mtx_gz=args.mtx,
        n_genes=int(args.n_genes),
        n_cells_total=int(args.n_cells_total),
        col_to_new=cohort.col_to_new,
        max_nnz=None,
    )
    X = normalize_log1p_cpm(X_raw, target_sum=1e4)

    print("=== Embedding ===")
    svd = TruncatedSVD(n_components=int(args.svd_components), random_state=int(args.seed))
    Z = svd.fit_transform(X)
    print(f"explained_var_sum={float(np.sum(svd.explained_variance_ratio_)):.4f}")
    print()

    print("=== Neighborhoods ===")
    neighbor_idx = build_knn(Z, k=int(args.knn_k))
    n = int(neighbor_idx.shape[0])
    centers = choose_centers(n, max_centers=int(args.max_centers), seed=int(args.seed))
    print(f"filtered_cells={n}")
    print(f"centers_tested={centers.size}")
    print()

    # Observed founder labels (center), plus validity, plus time entropy
    founder_center_obs = cohort.founder[centers]
    valid_obs = validity_mask_for_centers(neighbor_idx, centers, cohort.fate, cohort.founder)
    H_obs_all = time_entropy_bits_for_centers(neighbor_idx, centers, cohort.timebin)

    n_valid = int(np.sum(valid_obs))
    v_obs = float(n_valid / float(valid_obs.size)) if valid_obs.size > 0 else 0.0

    # Viability gates
    viability_ok = True
    reasons: List[str] = []

    if n_valid < int(args.min_units_total):
        viability_ok = False
        reasons.append(f"D1 failed: n_valid_centers < {int(args.min_units_total)} (got {n_valid})")

    summ_obs, S_obs, n_founders_eligible = compute_founder_medians_and_S(
        founder_center=founder_center_obs,
        time_entropy_bits=H_obs_all,
        valid=valid_obs,
        founder_names=cohort.founder_names,
        min_units_per_founder=int(args.min_units_per_founder),
    )

    if n_founders_eligible < 2:
        viability_ok = False
        reasons.append(
            f"D2 failed: need >=2 founders with n_units >= {int(args.min_units_per_founder)} (got {n_founders_eligible})"
        )

    B = int(args.perms)
    if B < 1:
        viability_ok = False
        reasons.append("D3 failed: perms < 1")

    # Null
    S_null = np.full(B, np.nan, dtype=np.float64)
    gkeys = make_group_keys(cohort.timebin, cohort.batch, mode=str(args.null_mode))

    if viability_ok:
        for b in range(B):
            f_sh = shuffle_within_groups(cohort.founder, gkeys, seed=int(args.seed) + 1000 + b)
            founder_center_b = f_sh[centers]
            valid_b = validity_mask_for_centers(neighbor_idx, centers, cohort.fate, f_sh)
            H_b_all = H_obs_all  # time entropy depends only on timebin, which is unchanged by founder shuffle

            _, S_b, n_eligible_b = compute_founder_medians_and_S(
                founder_center=founder_center_b,
                time_entropy_bits=H_b_all,
                valid=valid_b,
                founder_names=cohort.founder_names,
                min_units_per_founder=int(args.min_units_per_founder),
            )
            if np.isfinite(S_b) and (n_eligible_b >= 2):
                S_null[b] = S_b

        S_clean = S_null[np.isfinite(S_null)]
        if S_clean.size < max(10, int(0.8 * B)):
            viability_ok = False
            reasons.append(f"D3 failed: too many degenerate null perms (finite {S_clean.size}/{B})")
            S_med = float("nan")
            S_95 = float("nan")
            p_S = float("nan")
        else:
            S_med = float(np.median(S_clean))
            S_95 = float(np.quantile(S_clean, 0.95))
            p_S = _p_ge(S_clean, float(S_obs)) if np.isfinite(S_obs) else float("nan")
    else:
        S_med = float("nan")
        S_95 = float("nan")
        p_S = float("nan")

    # Decision
    if not viability_ok:
        decision = "NO-GO (viability)"
    else:
        detect_go = bool(np.isfinite(S_obs) and np.isfinite(S_95) and (float(S_obs) > float(S_95)))
        decision = "GO" if detect_go else "NO-GO (detectability)"

    uplift = (float(S_obs) / float(S_med)) if (np.isfinite(S_obs) and np.isfinite(S_med) and S_med > 0) else float("nan")

    # Write outputs
    unit_tbl = pd.DataFrame(
        dict(
            center_idx=np.arange(centers.size, dtype=np.int32),
            center_cell_id=cohort.cell_ids[centers],
            center_founder_code=founder_center_obs.astype(np.int32),
            center_founder=np.asarray([cohort.founder_names[int(c)] for c in founder_center_obs], dtype=str),
            valid=valid_obs.astype(bool),
            time_entropy_bits=H_obs_all.astype(np.float64),
        )
    )
    unit_path = os.path.join(args.out_dir, "unit_dispersion_table.tsv")
    unit_tbl.to_csv(unit_path, sep="\t", index=False)

    founder_path = os.path.join(args.out_dir, "dispersion_by_founder.tsv")
    summ_obs.to_csv(founder_path, sep="\t", index=False)

    null_path = os.path.join(args.out_dir, "null_dist_overall.tsv")
    pd.DataFrame({"S_null": S_null}).to_csv(null_path, sep="\t", index=False)

    summary_lines: List[str] = []
    summary_lines.append("=== Phase-D (GSE126954) — Founder dispersion / transition admissibility ===")
    summary_lines.append(f"annot={args.annot}")
    summary_lines.append(f"mtx={args.mtx}")
    summary_lines.append(f"out_dir={args.out_dir}")
    summary_lines.append("")
    summary_lines.append("--- Config ---")
    summary_lines.append(f"svd_components={int(args.svd_components)}")
    summary_lines.append(f"knn_k={int(args.knn_k)}")
    summary_lines.append(f"max_centers={int(args.max_centers)}")
    summary_lines.append(f"null_mode={str(args.null_mode)}")
    summary_lines.append(f"perms={int(args.perms)}")
    summary_lines.append("")
    summary_lines.append("--- Validity (same as Phase-1) ---")
    summary_lines.append("valid := neighborhood has >=2 founders AND >=2 fates")
    summary_lines.append(f"v_obs={v_obs:.6f}  (n_valid={n_valid}/{int(centers.size)})")
    summary_lines.append("")
    summary_lines.append("--- Statistic ---")
    summary_lines.append("dispersion metric = entropy_bits(embryo.time.bin within neighborhood)")
    summary_lines.append("S = Var_f(median_time_entropy_bits_by_founder) over eligible founders")
    summary_lines.append(f"min_units_per_founder={int(args.min_units_per_founder)}")
    summary_lines.append(f"S_obs={S_obs if np.isfinite(S_obs) else 'nan'}")
    summary_lines.append("")
    summary_lines.append("--- Null ---")
    summary_lines.append("null = shuffle founder within strata (timebin or timebin+batch); recompute valid centers and S")
    summary_lines.append(f"S_null_median={S_med if np.isfinite(S_med) else 'nan'}")
    summary_lines.append(f"S_null_95={S_95 if np.isfinite(S_95) else 'nan'}")
    summary_lines.append(f"uplift_S=S_obs/median(S_null)={uplift if np.isfinite(uplift) else 'nan'}")
    summary_lines.append(f"p_S(P(null>=obs))={p_S if np.isfinite(p_S) else 'nan'}")
    summary_lines.append("")
    summary_lines.append("--- Viability ---")
    summary_lines.append(f"min_units_total={int(args.min_units_total)}")
    summary_lines.append(f"viability_ok={str(viability_ok).lower()}")
    if reasons:
        summary_lines.append("viability_reasons:")
        for r in reasons:
            summary_lines.append(f"  - {r}")
    summary_lines.append("")
    summary_lines.append("--- Decision ---")
    summary_lines.append(f"Decision={decision}")
    summary_lines.append("")
    summary_lines.append("Wrote:")
    summary_lines.append(f"  - {os.path.basename(unit_path)}")
    summary_lines.append(f"  - {os.path.basename(founder_path)}")
    summary_lines.append(f"  - {os.path.basename(null_path)}")

    summary_path = os.path.join(args.out_dir, "phaseD_summary.txt")
    write_summary(summary_path, summary_lines)

    print("=== Phase-D summary ===")
    print(f"v_obs={v_obs:.6f} (n_valid={n_valid}/{int(centers.size)})")
    print(f"viability_ok={str(viability_ok).lower()}")
    if reasons:
        for r in reasons:
            print(f"  {r}")
    print(f"S_obs={S_obs if np.isfinite(S_obs) else 'nan'}")
    print(f"S_null_median={S_med if np.isfinite(S_med) else 'nan'}")
    print(f"S_null_95={S_95 if np.isfinite(S_95) else 'nan'}")
    print(f"p_S={p_S if np.isfinite(p_S) else 'nan'}")
    print(f"Decision={decision}")
    print(f"Wrote: {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
