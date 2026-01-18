#!/usr/bin/env python3
"""
phase1_tsv_audit.py

Generic Phase-1 audit for a joined TSV containing per-observation rows with
explicit clone IDs.

Required columns:
  - obs column (default: "obs_id", override with --obs_col)
  - clone column (default: "clone_id", override with --clone_col)
Optional:
  - clone_valid column (default: "clone_valid", override with --clone_valid_col)
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Tuple

import pandas as pd


def _coerce_bool_series(s: pd.Series) -> pd.Series:
    m = s.astype(str).str.strip().str.lower()
    true_set = {"1", "true", "t", "yes", "y"}
    false_set = {"0", "false", "f", "no", "n"}
    out = pd.Series([True] * len(m), index=m.index, dtype="boolean")
    out[m.isin(true_set)] = True
    out[m.isin(false_set)] = False
    out[out.isna()] = True
    return out.astype(bool)


def _ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def _basic_stats(
    df: pd.DataFrame,
    obs_col: str,
    clone_col: str,
    clone_valid_col: str,
) -> Dict[str, object]:
    obs = df[obs_col].astype(str).str.strip()
    clone = df[clone_col].astype(str).str.strip()

    if clone_valid_col in df.columns:
        clone_valid = _coerce_bool_series(df[clone_valid_col])
    else:
        clone_valid = pd.Series([True] * len(df), index=df.index)

    clone_missing = clone.eq("") | clone.str.lower().isin({"na", "nan", "none", "null"})
    clone_present = (~clone_missing) & clone_valid

    n_rows = int(len(df))
    n_obs_unique = int(obs.nunique(dropna=False))

    n_clone_present = int(clone_present.sum())
    n_clone_missing = int((~clone_present).sum())

    coverage_rate = n_clone_present / max(1, n_rows)

    df_ok = pd.DataFrame({obs_col: obs, clone_col: clone})[clone_present].copy()
    clone_sizes = df_ok.groupby(clone_col, sort=False)[obs_col].nunique()

    n_clones = int(clone_sizes.shape[0])
    max_clone = int(clone_sizes.max()) if n_clones > 0 else 0
    med_clone = float(clone_sizes.median()) if n_clones > 0 else 0.0

    return dict(
        rows=n_rows,
        obs_unique=n_obs_unique,
        clone_present=n_clone_present,
        clone_missing_or_invalid=n_clone_missing,
        coverage_rate=coverage_rate,
        clones=n_clones,
        clone_size_max=max_clone,
        clone_size_median=med_clone,
    )


def _clone_tables(
    df: pd.DataFrame,
    obs_col: str,
    clone_col: str,
    clone_valid_col: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    obs = df[obs_col].astype(str).str.strip()
    clone = df[clone_col].astype(str).str.strip()

    if clone_valid_col in df.columns:
        clone_valid = _coerce_bool_series(df[clone_valid_col])
    else:
        clone_valid = pd.Series([True] * len(df), index=df.index)

    clone_missing = clone.eq("") | clone.str.lower().isin({"na", "nan", "none", "null"})
    clone_present = (~clone_missing) & clone_valid

    df_ok = pd.DataFrame({obs_col: obs, clone_col: clone})[clone_present].copy()

    clone_sizes = (
        df_ok.groupby(clone_col, sort=False)[obs_col]
        .nunique()
        .reset_index()
        .rename(columns={obs_col: "n_obs"})
        .sort_values(["n_obs", clone_col], ascending=[False, True])
        .reset_index(drop=True)
    )

    hist = (
        clone_sizes.groupby("n_obs", sort=True)[clone_col]
        .count()
        .reset_index()
        .rename(columns={clone_col: "n_clones"})
        .sort_values("n_obs", ascending=True)
        .reset_index(drop=True)
    )

    return clone_sizes, hist


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_tsv", required=True, help="Joined Phase-1 TSV.")
    ap.add_argument("--out_dir", required=True, help="Output directory for audit artifacts.")
    ap.add_argument("--sep", default="\t", help="Delimiter for input TSV (default: tab).")

    ap.add_argument(
        "--obs_col",
        default="obs_id",
        help='Observation ID column (default: "obs_id"). For LARRY use: "Cell barcode".',
    )
    ap.add_argument(
        "--clone_col",
        default="clone_id",
        help='Clone ID column (default: "clone_id").',
    )
    ap.add_argument(
        "--clone_valid_col",
        default="clone_valid",
        help='Optional clone validity flag column (default: "clone_valid").',
    )

    ap.add_argument(
        "--min_coverage",
        type=float,
        default=0.20,
        help="GO threshold on clone assignment coverage rate (default: 0.20).",
    )
    args = ap.parse_args()

    _ensure_dir(args.out_dir)

    df = pd.read_csv(args.input_tsv, sep=args.sep, dtype=str)

    missing = [c for c in (args.obs_col, args.clone_col) if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns {missing}. Columns={list(df.columns)}")

    stats = _basic_stats(df, args.obs_col, args.clone_col, args.clone_valid_col)
    go = float(stats["coverage_rate"]) >= float(args.min_coverage)

    summary_path = os.path.join(args.out_dir, "audit_summary.txt")
    lines = []
    lines.append("# Phase-1 TSV audit summary")
    lines.append(f"input_tsv={args.input_tsv}")
    lines.append(f"obs_col={args.obs_col}")
    lines.append(f"clone_col={args.clone_col}")
    lines.append(f"clone_valid_col={args.clone_valid_col if args.clone_valid_col in df.columns else '(absent)'}")
    lines.append(f"rows={stats['rows']}")
    lines.append(f"obs_unique={stats['obs_unique']}")
    lines.append(f"clone_present={stats['clone_present']}")
    lines.append(f"clone_missing_or_invalid={stats['clone_missing_or_invalid']}")
    lines.append(f"coverage_rate={stats['coverage_rate']:.6f}")
    lines.append(f"clones={stats['clones']}")
    lines.append(f"clone_size_max={stats['clone_size_max']}")
    lines.append(f"clone_size_median={stats['clone_size_median']}")
    lines.append(f"min_coverage={args.min_coverage:.6f}")
    lines.append(f"GO={'true' if go else 'false'}")
    lines.append("")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    clone_sizes, hist = _clone_tables(df, args.obs_col, args.clone_col, args.clone_valid_col)

    clone_sizes_path = os.path.join(args.out_dir, "clone_size_table.tsv")
    clone_sizes.to_csv(clone_sizes_path, sep="\t", index=False)

    hist_path = os.path.join(args.out_dir, "clone_size_hist.tsv")
    hist.to_csv(hist_path, sep="\t", index=False)

    print(open(summary_path, "r", encoding="utf-8").read().rstrip())
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {clone_sizes_path}")
    print(f"Wrote: {hist_path}")

    return 0 if go else 2


if __name__ == "__main__":
    raise SystemExit(main())
