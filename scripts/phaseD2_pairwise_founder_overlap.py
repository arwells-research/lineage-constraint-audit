#!/usr/bin/env python3
"""
phaseD2_pairwise_founder_overlap.py — Phase-D2: Pairwise founder “shared admissibility space”
for GSE126954-style (worm) neighborhood audits.

Purpose
-------
Test whether certain founder pairs co-occupy local (valid) expression neighborhoods
more than expected under stratified founder-shuffle nulls.

This is a Phase-D extension that stays within audit discipline:
- same embedding + kNN geometry
- same validity rule as Phase-1 runner: neighborhood has >=2 founders AND >=2 fates
- no learned distances, no mechanistic claims
- null = shuffle founder labels within strata (timebin or timebin+batch)

Definitions
-----------
Neighborhood: the returned neighbor index list includes self because we build kNN with n_neighbors=k+1.
Validity (same as Phase-1 runner): neighborhood has >=2 founders AND >=2 fates.

Pairwise overlap weights (both are reported; the statistical test uses PRODUCT by default):
- Binary co-presence:
    w_bin(a,b) = 1 if count_a>0 and count_b>0 else 0
- Product weight:
    w_prod(a,b) = (count_a * count_b) / (K^2) where K = neighborhood size (k+1)

Aggregate (conditional on validity):
    O_ab = mean_j w(a,b) over valid centers j

Null
----
For each permutation:
- shuffle founder labels within strata (timebin or timebin+batch)
- recompute validity (because validity depends on founder mixing)
- recompute O_ab

Outputs
-------
out_dir/
  phaseD2_summary.txt
  pairwise_overlap_obs.tsv
  pairwise_overlap_null_summary.tsv
  pairwise_overlap_pvals.tsv
  null_dist_global.tsv
  (optional) tree_consistency.tsv   # only if --tree_dist_tsv is provided

Notes
-----
- This script is worm/GSE126954-oriented (uses embryo.time.bin + optional batch).
- If you want to apply to LARRY, you would need a different strata key and founder labels.
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


# ----------------------------
# Worm founder extraction
# ----------------------------

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


def _p_ge(null: np.ndarray, obs: float) -> float:
    B = int(null.size)
    return float((1 + np.sum(null >= obs)) / (1 + B))


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

    if "batch" not in df.columns:
        df["batch"] = "NO_BATCH_FIELD"

    df = df[df["cell.type"].notna()].copy()
    df = df[df["lineage"].notna()].copy()

    df["lineage"] = df["lineage"].astype(str)
    # Exclude ambiguous lineage strings
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
        (
            np.asarray(data, dtype=np.float32),
            (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32)),
        ),
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


# ----------------------------
# Pairwise overlap core
# ----------------------------

def _pairs(F: int) -> List[Tuple[int, int]]:
    out = []
    for a in range(F):
        for b in range(a + 1, F):
            out.append((a, b))
    return out


def _validity_for_centers(
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


def _pairwise_overlap_for_centers(
    neighbor_idx: np.ndarray,
    centers: np.ndarray,
    founder: np.ndarray,
    valid: np.ndarray,
    F: int,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Returns:
      O_bin: mean binary co-presence per pair (conditional on validity)
      O_prod: mean product-weight per pair (conditional on validity)
      n_valid: number of valid centers used
    """
    pairs = _pairs(F)
    P = len(pairs)
    sum_bin = np.zeros(P, dtype=np.float64)
    sum_prod = np.zeros(P, dtype=np.float64)

    # Neighborhood size includes self (k+1)
    K = int(neighbor_idx.shape[1])
    denom = float(K * K)

    n_valid = 0
    for j, i in enumerate(centers):
        if not bool(valid[j]):
            continue
        nb = neighbor_idx[i]
        g = founder[nb].astype(np.int32)
        counts = np.bincount(g, minlength=F)
        present = counts > 0
        for p_idx, (a, b) in enumerate(pairs):
            if present[a] and present[b]:
                sum_bin[p_idx] += 1.0
                sum_prod[p_idx] += (float(counts[a]) * float(counts[b])) / denom
        n_valid += 1

    if n_valid <= 0:
        return np.full(P, 0.0), np.full(P, 0.0), 0

    O_bin = sum_bin / float(n_valid)
    O_prod = sum_prod / float(n_valid)
    return O_bin, O_prod, n_valid


def _format_pair_name(a: int, b: int, founder_names: np.ndarray) -> str:
    an = str(founder_names[a]) if a < len(founder_names) else str(a)
    bn = str(founder_names[b]) if b < len(founder_names) else str(b)
    return f"{an}–{bn}"


def _write_text(path: str, lines: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for ln in lines:
            f.write(ln.rstrip() + "\n")


def _try_tree_consistency(
    out_dir: str,
    founder_names: np.ndarray,
    pair_tbl: pd.DataFrame,
    tree_dist_tsv: str,
) -> Optional[str]:
    """
    Optional: if user provides a founder-pair distance TSV, compute correlation vs overlap uplift.
    Expected input columns:
      founder_a, founder_b, dist
    Names must match cohort founder labels (AB, MS, C, D, E).

    Writes:
      tree_consistency.tsv
    """
    if not tree_dist_tsv:
        return None
    if not os.path.exists(tree_dist_tsv):
        raise RuntimeError(f"--tree_dist_tsv not found: {tree_dist_tsv}")

    dist_df = pd.read_csv(tree_dist_tsv, sep="\t")
    for col in ("founder_a", "founder_b", "dist"):
        if col not in dist_df.columns:
            raise RuntimeError(f"--tree_dist_tsv missing column: {col}")

    # normalize pair key
    def key(a: str, b: str) -> str:
        aa = str(a)
        bb = str(b)
        return "–".join(sorted([aa, bb]))

    dist_df = dist_df.copy()
    dist_df["pair"] = [key(a, b) for a, b in zip(dist_df["founder_a"], dist_df["founder_b"])]
    dist_map = dict(zip(dist_df["pair"], dist_df["dist"]))

    # pair_tbl has pair_name like "AB–MS"
    x = []
    y = []
    rows = []
    for _, r in pair_tbl.iterrows():
        pair_name = str(r["pair"])
        # dist keys are sorted, our pair names are already in that form (a–b)
        d = dist_map.get(pair_name)
        if d is None:
            continue
        uplift = r.get("uplift_prod", np.nan)
        if not np.isfinite(uplift) or uplift <= 0:
            continue
        log_uplift = float(np.log(uplift))
        x.append(float(d))
        y.append(log_uplift)
        rows.append({"pair": pair_name, "dist": float(d), "log_uplift_prod": log_uplift})

    out_path = os.path.join(out_dir, "tree_consistency.tsv")
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)

    return out_path


# ----------------------------
# Main
# ----------------------------

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

    # Viability (keep consistent with your earlier Stage-1A style floors)
    ap.add_argument("--min_units_total", type=int, default=100)

    # Optional tree distance table (if you want the AB–MS vs AB–E style check)
    ap.add_argument("--tree_dist_tsv", default="", help="Optional TSV with founder_a founder_b dist columns.")

    args = ap.parse_args()

    if not os.path.exists(args.annot):
        print(f"ERROR: annot not found: {args.annot}", file=sys.stderr)
        return 2
    if not os.path.exists(args.mtx):
        print(f"ERROR: mtx not found: {args.mtx}", file=sys.stderr)
        return 2

    os.makedirs(args.out_dir, exist_ok=True)

    cohort = load_and_filter_annotations(args.annot)
    F = int(len(cohort.founder_names))
    if F < 2:
        print("ERROR: fewer than 2 founders after filtering.", file=sys.stderr)
        return 2

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
    print(f"neighborhood_size_including_self={int(neighbor_idx.shape[1])}")
    print()

    # Observed validity and overlaps
    valid_obs = _validity_for_centers(neighbor_idx, centers, cohort.fate, cohort.founder)
    n_valid_obs = int(np.sum(valid_obs))
    v_obs = float(n_valid_obs / float(valid_obs.size)) if valid_obs.size > 0 else 0.0

    print("=== Observed validity ===")
    print(f"v_obs={v_obs:.6f}  (n_valid={n_valid_obs}/{int(centers.size)})")
    print()

    viability_ok = True
    reasons: List[str] = []
    if n_valid_obs < int(args.min_units_total):
        viability_ok = False
        reasons.append(f"D2-V1 failed: n_valid_centers < {int(args.min_units_total)} (got {n_valid_obs})")
    if int(args.perms) < 1:
        viability_ok = False
        reasons.append("D2-V2 failed: perms < 1")

    pairs = _pairs(F)
    pair_names = [_format_pair_name(a, b, cohort.founder_names) for a, b in pairs]
    P = len(pairs)

    O_bin_obs, O_prod_obs, n_valid_used_obs = _pairwise_overlap_for_centers(
        neighbor_idx=neighbor_idx,
        centers=centers,
        founder=cohort.founder,
        valid=valid_obs,
        F=F,
    )
    assert n_valid_used_obs == n_valid_obs

    # Null
    B = int(args.perms)
    O_bin_null = np.full((B, P), np.nan, dtype=np.float64)
    O_prod_null = np.full((B, P), np.nan, dtype=np.float64)

    gkeys = make_group_keys(cohort.timebin, cohort.batch, mode=str(args.null_mode))

    global_score_obs = float(np.mean(np.log((O_prod_obs + 1e-15) / (np.median(O_prod_obs) + 1e-15)))) if np.any(O_prod_obs > 0) else 0.0
    global_score_null = np.full(B, np.nan, dtype=np.float64)

    if viability_ok:
        print("=== Null (shuffle founders within strata) ===")
        for b in range(B):
            f_sh = shuffle_within_groups(cohort.founder, gkeys, seed=int(args.seed) + 1000 + b)
            valid_b = _validity_for_centers(neighbor_idx, centers, cohort.fate, f_sh)
            O_bin_b, O_prod_b, n_valid_b = _pairwise_overlap_for_centers(neighbor_idx, centers, f_sh, valid_b, F)

            O_bin_null[b, :] = O_bin_b
            O_prod_null[b, :] = O_prod_b

            # A simple global scalar for sanity + optional summary:
            # mean log-odds-like uplift of product overlaps vs their own median (per run).
            med = float(np.median(O_prod_b)) if np.isfinite(np.median(O_prod_b)) else 0.0
            if med > 0:
                global_score_null[b] = float(np.mean(np.log((O_prod_b + 1e-15) / (med + 1e-15))))
            else:
                global_score_null[b] = 0.0

            if (b + 1) % max(1, B // 10) == 0:
                print(f"  perm {b+1:>4}/{B}  n_valid={int(np.sum(valid_b))}")

        print()
    else:
        print("=== Null skipped (viability NO-GO) ===")
        print()

    # Summaries per pair
    obs_tbl = pd.DataFrame(
        {
            "pair": pair_names,
            "O_bin_obs": O_bin_obs,
            "O_prod_obs": O_prod_obs,
        }
    )

    null_summ_rows = []
    pval_rows = []
    for p_idx, name in enumerate(pair_names):
        if not viability_ok:
            null_summ_rows.append(
                dict(pair=name, O_bin_null_median=np.nan, O_bin_null_95=np.nan, O_prod_null_median=np.nan, O_prod_null_95=np.nan)
            )
            pval_rows.append(
                dict(pair=name, p_bin=np.nan, p_prod=np.nan, uplift_bin=np.nan, uplift_prod=np.nan)
            )
            continue

        bin_clean = O_bin_null[:, p_idx]
        prod_clean = O_prod_null[:, p_idx]

        bin_med = float(np.median(bin_clean))
        bin_95 = float(np.quantile(bin_clean, 0.95))
        prod_med = float(np.median(prod_clean))
        prod_95 = float(np.quantile(prod_clean, 0.95))

        p_bin = _p_ge(bin_clean, float(O_bin_obs[p_idx]))
        p_prod = _p_ge(prod_clean, float(O_prod_obs[p_idx]))

        uplift_bin = (float(O_bin_obs[p_idx]) / bin_med) if bin_med > 0 else np.nan
        uplift_prod = (float(O_prod_obs[p_idx]) / prod_med) if prod_med > 0 else np.nan

        null_summ_rows.append(
            dict(
                pair=name,
                O_bin_null_median=bin_med,
                O_bin_null_95=bin_95,
                O_prod_null_median=prod_med,
                O_prod_null_95=prod_95,
            )
        )
        pval_rows.append(
            dict(pair=name, p_bin=p_bin, p_prod=p_prod, uplift_bin=uplift_bin, uplift_prod=uplift_prod)
        )

    null_summ_tbl = pd.DataFrame(null_summ_rows)
    pvals_tbl = pd.DataFrame(pval_rows)

    # Decide (detectability) using PRODUCT overlaps:
    # "Any-pair exceeds null95" is too permissive; instead use a conservative global criterion:
    # - Require at least one pair with p_prod <= 0.01 AND uplift_prod >= 1.2
    # This matches the audit style: binary detectability + minimal effect size.
    decision = "NO-GO (viability)"
    if viability_ok:
        hit = pvals_tbl[
            (pvals_tbl["p_prod"].astype(float) <= 0.01) &
            (pvals_tbl["uplift_prod"].astype(float) >= 1.2)
        ]
        decision = "GO" if len(hit) > 0 else "NO-GO (detectability)"

    # Global sanity stats
    if viability_ok:
        gclean = global_score_null[np.isfinite(global_score_null)]
        p_global = _p_ge(gclean, global_score_obs) if gclean.size > 0 else np.nan
        g_med = float(np.median(gclean)) if gclean.size > 0 else np.nan
        g_95 = float(np.quantile(gclean, 0.95)) if gclean.size > 0 else np.nan
    else:
        p_global = np.nan
        g_med = np.nan
        g_95 = np.nan

    # Write outputs
    obs_path = os.path.join(args.out_dir, "pairwise_overlap_obs.tsv")
    null_path = os.path.join(args.out_dir, "pairwise_overlap_null_summary.tsv")
    pval_path = os.path.join(args.out_dir, "pairwise_overlap_pvals.tsv")
    global_path = os.path.join(args.out_dir, "null_dist_global.tsv")

    obs_tbl.to_csv(obs_path, sep="\t", index=False)
    null_summ_tbl.to_csv(null_path, sep="\t", index=False)
    pvals_tbl.to_csv(pval_path, sep="\t", index=False)
    pd.DataFrame({"global_score_null": global_score_null}).to_csv(global_path, sep="\t", index=False)

    # Optional tree consistency
    tree_out = None
    if args.tree_dist_tsv:
        try:
            merged = obs_tbl.merge(null_summ_tbl, on="pair", how="left").merge(pvals_tbl, on="pair", how="left")
            tree_out = _try_tree_consistency(args.out_dir, cohort.founder_names, merged, args.tree_dist_tsv)
        except Exception as e:
            # Keep audit strict: do not fail the whole run for optional analysis
            print(f"WARNING: tree consistency skipped due to error: {e}", file=sys.stderr)
            tree_out = None

    # Summary file
    lines: List[str] = []
    lines.append("=== Phase-D2 — Pairwise founder overlap (shared admissibility space) ===")
    lines.append(f"annot={args.annot}")
    lines.append(f"mtx={args.mtx}")
    lines.append(f"out_dir={args.out_dir}")
    lines.append("")
    lines.append("--- Config ---")
    lines.append(f"svd_components={int(args.svd_components)}")
    lines.append(f"knn_k={int(args.knn_k)}")
    lines.append(f"max_centers={int(args.max_centers)}")
    lines.append(f"null_mode={str(args.null_mode)}")
    lines.append(f"perms={int(args.perms)}")
    lines.append("")
    lines.append("--- Validity (same as Phase-1) ---")
    lines.append("valid := neighborhood has >=2 founders AND >=2 fates")
    lines.append(f"v_obs={v_obs:.6f}  (n_valid={n_valid_obs}/{int(centers.size)})")
    lines.append("")
    lines.append("--- Decision rule (Phase-D2 detectability) ---")
    lines.append("Test channel: PRODUCT overlap (count_a * count_b / K^2)")
    lines.append("GO if exists at least one pair with (p_prod <= 0.01) AND (uplift_prod >= 1.2).")
    lines.append("")
    lines.append("--- Global sanity (not a gate) ---")
    lines.append(f"global_score_obs={global_score_obs}")
    lines.append(f"global_score_null_median={g_med}")
    lines.append(f"global_score_null_95={g_95}")
    lines.append(f"p_global={p_global}")
    lines.append("")
    lines.append("--- Viability ---")
    lines.append(f"min_units_total={int(args.min_units_total)}")
    lines.append(f"viability_ok={str(viability_ok).lower()}")
    if reasons:
        lines.append("viability_reasons:")
        for r in reasons:
            lines.append(f"  - {r}")
    lines.append("")
    lines.append(f"Decision={decision}")
    lines.append("")
    lines.append("Wrote:")
    lines.append(f"  - {os.path.basename(obs_path)}")
    lines.append(f"  - {os.path.basename(null_path)}")
    lines.append(f"  - {os.path.basename(pval_path)}")
    lines.append(f"  - {os.path.basename(global_path)}")
    if tree_out:
        lines.append(f"  - {os.path.basename(tree_out)}")

    summ_path = os.path.join(args.out_dir, "phaseD2_summary.txt")
    _write_text(summ_path, lines)

    print("=== Phase-D2 summary ===")
    print(f"v_obs={v_obs:.6f} (n_valid={n_valid_obs}/{int(centers.size)})")
    print(f"viability_ok={str(viability_ok).lower()}")
    if reasons:
        for r in reasons:
            print(f"  {r}")
    print(f"Decision={decision}")
    print(f"Wrote: {summ_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
    