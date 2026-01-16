#!/usr/bin/env python3
"""
Phase-5A: Founder contribution asymmetry decomposition.

This is NOT a new hypothesis test. It decomposes the already-established
Phase-4 Σ₂ signal inside the frozen locus:
  - centers: cell.type == Hypodermis
  - time bins: 210-270,270-330,330-390  (default)
  - neighborhoods from the same embedding + kNN as Phase-1/4

For each valid center i:
  - build contingency O[y, g] where y=cell.type (fate), g=founder
  - compute G-test statistic T_i
  - compute per-founder contribution T_i(g) so that sum_g T_i(g) = T_i

Outputs a table of:
  founder, contrib_sum, contrib_share, n_valid_centers_with_founder_present, mean_neighbor_frac
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
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


def stream_read_mtx_subset_cells(
    mtx_gz: str,
    n_genes: int,
    n_cells_total: int,
    col_to_new: Dict[int, int],
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

    if kept == 0:
        raise RuntimeError("No nonzero entries retained for filtered cohort; mapping may be wrong.")

    X = sp.coo_matrix(
        (np.asarray(data, dtype=np.float32),
         (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32))),
        shape=(n_filtered, n_genes),
    ).tocsr()
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


def load_and_filter_annotations(path_gz: str) -> Tuple[pd.DataFrame, Dict[int, int], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(path_gz, compression="gzip")
    if len(df) <= 0:
        raise RuntimeError("Annotation file appears empty.")

    df = df.copy()
    df["col_index_1based"] = np.arange(1, len(df) + 1, dtype=np.int32)

    qc_col = "passed_initial_QC_or_later_whitelisted"
    if qc_col in df.columns:
        df = df[df[qc_col] == True].copy()

    df = df[df["cell.type"].notna()].copy()
    df = df[df["lineage"].notna()].copy()
    df["lineage"] = df["lineage"].astype(str)
    df = df[~df["lineage"].str.contains("/", regex=False)].copy()

    df["founder"] = df["lineage"].map(founder_group)
    df = df[df["founder"].notna()].copy()

    if "batch" not in df.columns:
        df["batch"] = "NO_BATCH_FIELD"

    fate_codes, _fate_uniques = pd.factorize(df["cell.type"], sort=True)
    founder_codes, _founder_uniques = pd.factorize(df["founder"], sort=True)
    time_codes, _time_uniques = pd.factorize(df["embryo.time.bin"], sort=True)
    batch_codes, _batch_uniques = pd.factorize(df["batch"].astype(str), sort=True)

    df["fate_code"] = fate_codes.astype(np.int32)
    df["founder_code"] = founder_codes.astype(np.int32)
    df["time_code"] = time_codes.astype(np.int32)
    df["batch_code"] = batch_codes.astype(np.int32)

    col_indices = df["col_index_1based"].to_numpy(dtype=np.int32)
    col_to_new = {int(c): i for i, c in enumerate(col_indices)}

    df = df.reset_index(drop=True)
    fate = df["fate_code"].to_numpy(dtype=np.int32)
    founder = df["founder_code"].to_numpy(dtype=np.int32)
    timebin = df["time_code"].to_numpy(dtype=np.int32)
    batch = df["batch_code"].to_numpy(dtype=np.int32)
    return df, col_to_new, fate, founder, timebin, batch


def g_test_and_founder_contrib(y: np.ndarray, g: np.ndarray) -> Tuple[float, Dict[int, float]]:
    """
    Returns:
      T: float, G-test statistic
      contrib: dict founder_code -> contribution sum so that sum(contrib)=T
    """
    uy, yinv = np.unique(y, return_inverse=True)
    ug, ginv = np.unique(g, return_inverse=True)
    ny = uy.size
    ng = ug.size

    obs = np.zeros((ny, ng), dtype=np.int32)
    np.add.at(obs, (yinv, ginv), 1)

    n = int(obs.sum())
    if n <= 0:
        return 0.0, {}

    row = obs.sum(axis=1).astype(np.float64)
    col = obs.sum(axis=0).astype(np.float64)
    exp = (row[:, None] * col[None, :]) / float(n)

    contrib: Dict[int, float] = {}
    T = 0.0
    for j_g in range(ng):
        og = obs[:, j_g]
        eg = exp[:, j_g]
        mask = og > 0
        if not np.any(mask):
            contrib[int(ug[j_g])] = 0.0
            continue
        og_f = og[mask].astype(np.float64)
        eg_f = np.maximum(eg[mask], 1e-300)
        c = 2.0 * float(np.sum(og_f * np.log(og_f / eg_f)))
        contrib[int(ug[j_g])] = c
        T += c

    return float(T), contrib


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot", default="GSE126954_cell_annotation.csv.gz")
    ap.add_argument("--mtx", default="GSE126954_gene_by_cell_count_matrix.txt.gz")
    ap.add_argument("--n_genes", type=int, default=20222)
    ap.add_argument("--n_cells_total", type=int, default=89701)

    ap.add_argument("--svd_components", type=int, default=50)
    ap.add_argument("--knn_k", type=int, default=20)
    ap.add_argument("--max_centers", type=int, default=5000, help="0 means all")
    ap.add_argument("--seed", type=int, default=0)

    ap.add_argument("--center_fate", default="Hypodermis")
    ap.add_argument("--timebins", default="210-270,270-330,330-390")
    ap.add_argument("--min_valid_centers", type=int, default=30, help="Just a sanity print gate; no filtering of computation.")

    args = ap.parse_args()

    if not os.path.exists(args.annot):
        print(f"ERROR: annotation not found: {args.annot}", file=sys.stderr)
        return 2
    if not os.path.exists(args.mtx):
        print(f"ERROR: matrix not found: {args.mtx}", file=sys.stderr)
        return 2

    df, col_to_new, fate, founder, timebin, batch = load_and_filter_annotations(args.annot)

    # Names for readability
    founder_names = pd.Index(df["founder"].astype(str)).unique().tolist()
    # But factorization order might differ; reconstruct mapping from codes -> labels:
    code_to_founder = df.drop_duplicates("founder_code")[["founder_code", "founder"]].set_index("founder_code")["founder"].to_dict()
    code_to_fate = df.drop_duplicates("fate_code")[["fate_code", "cell.type"]].set_index("fate_code")["cell.type"].to_dict()
    code_to_time = df.drop_duplicates("time_code")[["time_code", "embryo.time.bin"]].set_index("time_code")["embryo.time.bin"].to_dict()

    X_raw = stream_read_mtx_subset_cells(
        mtx_gz=args.mtx,
        n_genes=int(args.n_genes),
        n_cells_total=int(args.n_cells_total),
        col_to_new=col_to_new,
    )
    X = normalize_log1p_cpm(X_raw, target_sum=1e4)

    svd = TruncatedSVD(n_components=int(args.svd_components), random_state=int(args.seed))
    Z = svd.fit_transform(X)

    neighbor_idx = build_knn(Z, k=int(args.knn_k))
    n = neighbor_idx.shape[0]
    centers = choose_centers(n, max_centers=int(args.max_centers), seed=int(args.seed))

    # Freeze Phase-4 locus for centers
    timebins_keep = [t.strip() for t in str(args.timebins).split(",") if t.strip()]
    time_keep_codes = set(k for k, v in code_to_time.items() if v in timebins_keep)

    center_fate = str(args.center_fate)
    center_fate_codes = set(k for k, v in code_to_fate.items() if v == center_fate)
    if not center_fate_codes:
        print(f"ERROR: center_fate={center_fate} not found in annotations (post-filter).", file=sys.stderr)
        return 2

    # Restrict centers by fate and timebin
    centers_mask = np.zeros(centers.size, dtype=bool)
    for j, i in enumerate(centers):
        if (int(fate[i]) in center_fate_codes) and (int(timebin[i]) in time_keep_codes):
            centers_mask[j] = True
    centers_frozen = centers[centers_mask]

    print("=== Phase-5A frozen cohort (centers) ===")
    print(f"center_fate: {center_fate}")
    print(f"timebins: {', '.join(timebins_keep)}")
    print(f"centers sampled: {centers.size}")
    print(f"centers in frozen cohort: {centers_frozen.size}")
    print()

    if centers_frozen.size == 0:
        print("NO-GO: no centers in frozen cohort under current sampling.")
        return 0

    # Compute per-founder contribution aggregates over VALID centers only
    contrib_sum: Dict[int, float] = {}
    present_count: Dict[int, int] = {}
    neighbor_frac_sum: Dict[int, float] = {}

    n_valid = 0
    T_total = 0.0

    for i in centers_frozen:
        nb = neighbor_idx[i]
        y = fate[nb]
        g = founder[nb]

        # validity rule (same as Phase-1/4)
        if np.unique(y).size < 2:
            continue
        if np.unique(g).size < 2:
            continue

        n_valid += 1
        T_i, c_i = g_test_and_founder_contrib(y, g)
        T_total += T_i

        # neighbor founder fractions (descriptive)
        ug, counts = np.unique(g, return_counts=True)
        frac = counts.astype(np.float64) / float(counts.sum())

        for code, cval in c_i.items():
            contrib_sum[code] = contrib_sum.get(code, 0.0) + float(cval)
            present_count[code] = present_count.get(code, 0) + (1 if code in ug else 0)

        for code_u, frac_u in zip(ug.tolist(), frac.tolist()):
            neighbor_frac_sum[int(code_u)] = neighbor_frac_sum.get(int(code_u), 0.0) + float(frac_u)

    print("=== Phase-5A validity ===")
    print(f"valid centers in frozen cohort: {n_valid}/{centers_frozen.size}")
    if n_valid < int(args.min_valid_centers):
        print(f"WARNING: n_valid < {int(args.min_valid_centers)}; founder-share estimates may be noisy.")
    print()

    if n_valid == 0 or T_total <= 0.0:
        print("NO-GO: no valid centers (or zero total statistic) in frozen cohort.")
        return 0

    # Build output table
    all_codes = sorted(contrib_sum.keys())
    rows = []
    for code in all_codes:
        name = str(code_to_founder.get(code, f"founder_code_{code}"))
        csum = float(contrib_sum.get(code, 0.0))
        share = csum / float(T_total) if T_total > 0 else 0.0
        n_pres = int(present_count.get(code, 0))
        mean_frac = float(neighbor_frac_sum.get(code, 0.0)) / float(n_valid) if n_valid > 0 else 0.0
        rows.append(
            dict(
                founder=name,
                founder_code=int(code),
                contrib_sum=csum,
                contrib_share=share,
                n_valid_centers_with_founder_present=n_pres,
                mean_neighbor_frac=mean_frac,
            )
        )

    out = pd.DataFrame(rows).sort_values(by=["contrib_share", "contrib_sum"], ascending=[False, False])

    print("=== Phase-5A founder contribution decomposition ===")
    print(f"T_total (sum over valid centers): {T_total:.6f}")
    with pd.option_context("display.max_rows", 50, "display.max_colwidth", 120):
        print(out.to_string(index=False, justify="left", formatters={
            "contrib_sum": lambda x: f"{x:.6f}",
            "contrib_share": lambda x: f"{x:.6f}",
            "mean_neighbor_frac": lambda x: f"{x:.6f}",
        }))
    print()

    # Simple keep/discard gate suggestion (you decide whether to apply)
    top1 = float(out.iloc[0]["contrib_share"]) if len(out) >= 1 else 0.0
    top2 = float(out.iloc[1]["contrib_share"]) if len(out) >= 2 else 0.0
    print("=== Phase-5A keep/discard gate (suggested) ===")
    print("KEEP if (top1 >= 0.70) OR (top1 + top2 >= 0.85)")
    print(f"top1={top1:.6f}  top1+top2={(top1+top2):.6f}")
    print(f"Decision: {'KEEP' if (top1 >= 0.70 or (top1 + top2) >= 0.85) else 'DISCARD'}")
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())