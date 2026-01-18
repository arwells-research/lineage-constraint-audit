#!/usr/bin/env python3
"""
LARRY Phase-5A: Founder contribution asymmetry decomposition (clone-internal).

This is a decomposition step, not a new hypothesis test.

Given a Phase-1 joined table with per-cell fields:
  - Cell barcode (or --obs_col)
  - clone_id + clone_valid
  - Starting population (founder)
  - Cell type annotation (fate)

For each valid clone:
  - build contingency O[y, g] where y=fate, g=founder
  - compute G-test statistic T_clone
  - compute per-founder contribution T_clone(g) so that sum_g T_clone(g) = T_clone

Aggregate across clones:
  - founder, contrib_sum, contrib_share, n_valid_clones_with_founder_present, mean_within_clone_founder_frac

Validity rule per clone (conservative, to match “asymmetry” intent):
  - at least 2 unique fates AND at least 2 unique founders inside the clone

Outputs:
  - out_dir/founder_contrib_table.tsv
  - out_dir/clone_gtest_table.tsv
  - out_dir/summary.txt
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def _parse_bool_series(x: pd.Series) -> pd.Series:
    """
    Accepts True/False, 1/0, "true"/"false", "yes"/"no".
    Unknown -> False.
    """
    if x.dtype == bool:
        return x.fillna(False)
    s = x.astype(str).str.strip().str.lower()
    return s.isin(["1", "true", "t", "yes", "y"])


def g_test_and_founder_contrib(y: np.ndarray, g: np.ndarray) -> Tuple[float, Dict[str, float]]:
    """
    y: fate labels
    g: founder labels

    Returns:
      T: float, G-test statistic
      contrib: dict founder_label -> contribution sum so that sum(contrib)=T
    """
    # Factorize to compact integer codes, but keep label mapping
    y_codes, y_uniques = pd.factorize(pd.Series(y), sort=True)
    g_codes, g_uniques = pd.factorize(pd.Series(g), sort=True)

    ny = int(len(y_uniques))
    ng = int(len(g_uniques))

    obs = np.zeros((ny, ng), dtype=np.int64)
    np.add.at(obs, (y_codes, g_codes), 1)

    n = int(obs.sum())
    if n <= 0:
        return 0.0, {}

    row = obs.sum(axis=1).astype(np.float64)
    col = obs.sum(axis=0).astype(np.float64)
    exp = (row[:, None] * col[None, :]) / float(n)

    contrib: Dict[str, float] = {}
    T = 0.0

    for j in range(ng):
        og = obs[:, j]
        eg = exp[:, j]
        mask = og > 0
        if not np.any(mask):
            contrib[str(g_uniques[j])] = 0.0
            continue
        og_f = og[mask].astype(np.float64)

        # Guard tiny expected counts
        eg_f = np.maximum(eg[mask], 1e-300)

        # 2 * sum( O * ln(O/E) )
        c = 2.0 * float(np.sum(og_f * np.log(og_f / eg_f)))
        contrib[str(g_uniques[j])] = c
        T += c

    return float(T), contrib


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_tsv", required=True, help="Phase-1 joined TSV (cells x fields).")
    ap.add_argument("--out_dir", required=True, help="Output directory.")

    ap.add_argument("--obs_col", default="Cell barcode", help="Observation id column.")
    ap.add_argument("--clone_col", default="clone_id", help="Clone id column.")
    ap.add_argument("--clone_valid_col", default="clone_valid", help="Clone-valid flag column.")
    ap.add_argument("--fate_col", default="Cell type annotation", help="Fate column.")
    ap.add_argument("--founder_col", default="Starting population", help="Founder column.")

    ap.add_argument("--min_cells_per_clone", type=int, default=2, help="Skip clones smaller than this.")
    ap.add_argument("--min_unique_fates", type=int, default=2, help="Per-clone validity: unique fates >= this.")
    ap.add_argument("--min_unique_founders", type=int, default=2, help="Per-clone validity: unique founders >= this.")

    args = ap.parse_args()

    if not os.path.exists(args.input_tsv):
        print(f"ERROR: input_tsv not found: {args.input_tsv}", file=sys.stderr)
        return 2

    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str)
    required = [args.obs_col, args.clone_col, args.clone_valid_col, args.fate_col, args.founder_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: missing required columns {missing}. Columns={df.columns.tolist()[:50]}", file=sys.stderr)
        return 2

    # Normalize fields
    df = df.copy()
    df[args.obs_col] = df[args.obs_col].astype(str).str.strip()
    df[args.clone_col] = df[args.clone_col].astype(str).str.strip()
    df[args.fate_col] = df[args.fate_col].astype(str).str.strip()
    df[args.founder_col] = df[args.founder_col].astype(str).str.strip()

    clone_valid = _parse_bool_series(df[args.clone_valid_col])
    df = df[clone_valid].copy()

    # Exclude NA/blank clone ids defensively
    bad_clone = df[args.clone_col].isna() | (df[args.clone_col].str.upper() == "NA") | (df[args.clone_col] == "")
    df = df[~bad_clone].copy()

    # Exclude missing fate/founder (strict)
    df = df[(df[args.fate_col] != "") & (df[args.fate_col].str.upper() != "NA")].copy()
    df = df[(df[args.founder_col] != "") & (df[args.founder_col].str.upper() != "NA")].copy()

    if df.empty:
        print("NO-GO: after filtering, no valid rows remain.", file=sys.stderr)
        return 2

    grouped = df.groupby(args.clone_col, sort=False)

    contrib_sum: Dict[str, float] = {}
    present_count: Dict[str, int] = {}
    frac_sum: Dict[str, float] = {}
    T_total = 0.0
    n_valid_clones = 0

    clone_rows: List[Dict[str, object]] = []

    for clone_id, sub in grouped:
        n_cells = int(len(sub))
        if n_cells < int(args.min_cells_per_clone):
            continue

        y = sub[args.fate_col].to_numpy(dtype=str)
        g = sub[args.founder_col].to_numpy(dtype=str)

        n_fates = int(pd.unique(y).size)
        n_founders = int(pd.unique(g).size)

        valid = (n_fates >= int(args.min_unique_fates)) and (n_founders >= int(args.min_unique_founders))

        T_i = 0.0
        if valid:
            T_i, contrib_i = g_test_and_founder_contrib(y, g)
            if T_i > 0.0:
                n_valid_clones += 1
                T_total += T_i

                founders_u, counts_u = np.unique(g, return_counts=True)
                frac_u = counts_u.astype(np.float64) / float(counts_u.sum())

                for f in founders_u.tolist():
                    present_count[f] = present_count.get(f, 0) + 1
                for f, fr in zip(founders_u.tolist(), frac_u.tolist()):
                    frac_sum[f] = frac_sum.get(f, 0.0) + float(fr)

                for f, cval in contrib_i.items():
                    contrib_sum[f] = contrib_sum.get(f, 0.0) + float(cval)

        major_fate = ""
        major_fate_count = 0
        vc = sub[args.fate_col].value_counts()
        if len(vc) > 0:
            major_fate = str(vc.index[0])
            major_fate_count = int(vc.iloc[0])

        major_founder = ""
        major_founder_count = 0
        vcg = sub[args.founder_col].value_counts()
        if len(vcg) > 0:
            major_founder = str(vcg.index[0])
            major_founder_count = int(vcg.iloc[0])

        clone_rows.append(
            dict(
                clone_id=str(clone_id),
                clone_size=n_cells,
                n_fates=n_fates,
                n_founders=n_founders,
                valid=bool(valid),
                T_gtest=float(T_i),
                major_fate=major_fate,
                major_fate_count=major_fate_count,
                major_founder=major_founder,
                major_founder_count=major_founder_count,
            )
        )

    out_clone = pd.DataFrame(clone_rows).sort_values(
        by=["valid", "T_gtest", "clone_size"],
        ascending=[False, False, False],
    )

    out_clone_path = os.path.join(args.out_dir, "clone_gtest_table.tsv")
    out_clone.to_csv(out_clone_path, sep="\t", index=False)

    founder_rows: List[Dict[str, object]] = []
    all_founders = sorted(set(list(contrib_sum.keys()) + list(present_count.keys()) + list(frac_sum.keys())))

    for fndr in all_founders:
        csum = float(contrib_sum.get(fndr, 0.0))
        share = (csum / float(T_total)) if T_total > 0.0 else 0.0
        n_pres = int(present_count.get(fndr, 0))
        mean_frac = (float(frac_sum.get(fndr, 0.0)) / float(n_pres)) if n_pres > 0 else 0.0
        founder_rows.append(
            dict(
                founder=str(fndr),
                contrib_sum=csum,
                contrib_share=share,
                n_valid_clones_with_founder_present=n_pres,
                mean_within_clone_founder_frac=mean_frac,
            )
        )

    out_founder = pd.DataFrame(founder_rows).sort_values(
        by=["contrib_share", "contrib_sum"],
        ascending=[False, False],
    )

    out_founder_path = os.path.join(args.out_dir, "founder_contrib_table.tsv")
    out_founder.to_csv(out_founder_path, sep="\t", index=False)

    summary_path = os.path.join(args.out_dir, "summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("# LARRY Phase-5A founder contribution summary\n")
        f.write(f"input_tsv={args.input_tsv}\n")
        f.write(f"rows_after_valid_filter={len(df)}\n")
        f.write(f"clones_total={df[args.clone_col].nunique()}\n")
        f.write(f"valid_clones={n_valid_clones}\n")
        f.write(f"T_total={T_total:.6f}\n")
        f.write("\n")

        top1 = float(out_founder.iloc[0]["contrib_share"]) if len(out_founder) >= 1 else 0.0
        top2 = float(out_founder.iloc[1]["contrib_share"]) if len(out_founder) >= 2 else 0.0
        f.write("## Suggested keep/discard gate (optional)\n")
        f.write("KEEP if (top1 >= 0.70) OR (top1 + top2 >= 0.85)\n")
        f.write(f"top1={top1:.6f}\n")
        f.write(f"top1_plus_top2={(top1 + top2):.6f}\n")
        f.write(f"Decision={'KEEP' if (top1 >= 0.70 or (top1 + top2) >= 0.85) else 'DISCARD'}\n")
        f.write("\n")

    print("=== LARRY Phase-5A founder contribution decomposition ===")
    print(f"input_tsv: {args.input_tsv}")
    print(f"rows_after_valid_filter: {len(df)}")
    print(f"clones_total: {df[args.clone_col].nunique()}")
    print(f"valid_clones: {n_valid_clones}")
    print(f"T_total: {T_total:.6f}")
    print(f"Wrote: {out_founder_path}")
    print(f"Wrote: {out_clone_path}")
    print(f"Wrote: {summary_path}")
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
