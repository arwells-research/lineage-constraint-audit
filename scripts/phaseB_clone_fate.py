#!/usr/bin/env python3
import argparse
import os
import math
from typing import Dict, List, Tuple

import pandas as pd


def shannon_entropy_bits(counts: List[int]) -> float:
    """
    Shannon entropy in bits for a discrete distribution given integer counts.
    """
    total = sum(counts)
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts:
        if c <= 0:
            continue
        p = c / total
        h -= p * math.log(p, 2)
    return float(h)


def write_text(path: str, s: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(s)


def make_hist(values: pd.Series, bin_edges: List[float]) -> pd.DataFrame:
    """
    Return a simple histogram table with fixed bin edges.
    Output columns: bin_left, bin_right, count
    """
    # pandas cut: right-open by default if right=False
    cats = pd.cut(values, bins=bin_edges, right=False, include_lowest=True)
    vc = cats.value_counts().sort_index()
    rows = []
    for interval, cnt in vc.items():
        rows.append(
            {
                "bin_left": float(interval.left),
                "bin_right": float(interval.right),
                "count": int(cnt),
            }
        )
    return pd.DataFrame(rows)


def run_phase_b(input_tsv: str, out_dir: str) -> int:
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(input_tsv, sep="\t", dtype=str, low_memory=False)

    required = ["clone_id", "clone_valid", "Cell type annotation"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns {missing}. Columns={list(df.columns)[:30]}")

    # Normalize booleans / NA
    df["clone_id"] = df["clone_id"].astype(str).str.strip()
    df["clone_valid"] = df["clone_valid"].astype(str).str.strip().str.lower()
    df["Cell type annotation"] = df["Cell type annotation"].astype(str).str.strip()

    # Keep only valid clone assignments
    valid = df[(df["clone_valid"].isin(["true", "1", "yes"])) & (df["clone_id"] != "") & (df["clone_id"] != "NA")]
    if len(valid) == 0:
        write_text(os.path.join(out_dir, "phaseB_summary.txt"), "ERROR: no valid clone rows found.\n")
        return 2

    # Group clone -> cell-type counts
    # table: index=(clone_id, cell_type) -> count
    g = (
        valid.groupby(["clone_id", "Cell type annotation"], dropna=False)
        .size()
        .reset_index(name="n")
    )

    # Per-clone totals
    totals = g.groupby("clone_id")["n"].sum().rename("clone_size").reset_index()

    # Per-clone type stats
    ntypes = g.groupby("clone_id")["Cell type annotation"].nunique().rename("n_cell_types").reset_index()

    # Majority type per clone
    # sort within clone by count desc, then by type name for stability
    g_sorted = g.sort_values(["clone_id", "n", "Cell type annotation"], ascending=[True, False, True])
    major = g_sorted.groupby("clone_id").head(1).rename(columns={"Cell type annotation": "major_type", "n": "major_count"})
    major = major[["clone_id", "major_type", "major_count"]].reset_index(drop=True)

    # Entropy per clone
    # collect counts list per clone
    ent_rows: List[Tuple[str, float]] = []
    for cid, sub in g.groupby("clone_id"):
        counts = sub["n"].astype(int).tolist()
        ent_rows.append((cid, shannon_entropy_bits(counts)))
    ent = pd.DataFrame(ent_rows, columns=["clone_id", "entropy_bits"])

    # Merge
    out = totals.merge(ntypes, on="clone_id", how="left").merge(major, on="clone_id", how="left").merge(ent, on="clone_id", how="left")
    out["clone_size"] = out["clone_size"].astype(int)
    out["major_count"] = out["major_count"].astype(int)
    out["n_cell_types"] = out["n_cell_types"].astype(int)
    out["purity"] = out["major_count"] / out["clone_size"]

    # Fate class bins
    def fate_class(n_types: int) -> str:
        if n_types <= 1:
            return "single-type"
        if n_types == 2:
            return "mixed-2"
        return "mixed-3+"

    out["fate_class"] = out["n_cell_types"].apply(fate_class)

    # Write main table
    path_summary = os.path.join(out_dir, "clone_fate_summary.tsv")
    out.sort_values(["clone_size", "purity"], ascending=[False, True]).to_csv(path_summary, sep="\t", index=False)

    # Histograms
    purity_bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0000001]
    ent_bins = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 8.0]

    purity_hist = make_hist(out["purity"].astype(float), purity_bins)
    purity_hist.to_csv(os.path.join(out_dir, "clone_purity_hist.tsv"), sep="\t", index=False)

    ent_hist = make_hist(out["entropy_bits"].astype(float), ent_bins)
    ent_hist.to_csv(os.path.join(out_dir, "clone_entropy_hist.tsv"), sep="\t", index=False)

    # Phase-B summary text
    n_clones = int(out.shape[0])
    n_cells = int(valid.shape[0])
    single = int((out["n_cell_types"] == 1).sum())
    mixed2 = int((out["n_cell_types"] == 2).sum())
    mixed3p = int((out["n_cell_types"] >= 3).sum())
    purity_med = float(out["purity"].median())
    purity_p10 = float(out["purity"].quantile(0.10))
    purity_p90 = float(out["purity"].quantile(0.90))
    ent_med = float(out["entropy_bits"].median())
    ent_p90 = float(out["entropy_bits"].quantile(0.90))

    summary_lines = [
        "# Phase-B fate summary",
        f"input_tsv={input_tsv}",
        f"rows_valid_clone={n_cells}",
        f"clones={n_clones}",
        f"single_type_clones={single}",
        f"mixed_2_clones={mixed2}",
        f"mixed_3p_clones={mixed3p}",
        f"purity_median={purity_med:.6f}",
        f"purity_p10={purity_p10:.6f}",
        f"purity_p90={purity_p90:.6f}",
        f"entropy_bits_median={ent_med:.6f}",
        f"entropy_bits_p90={ent_p90:.6f}",
        "",
        f"Wrote: {path_summary}",
        f"Wrote: {os.path.join(out_dir, 'clone_purity_hist.tsv')}",
        f"Wrote: {os.path.join(out_dir, 'clone_entropy_hist.tsv')}",
        "",
    ]
    write_text(os.path.join(out_dir, "phaseB_summary.txt"), "\n".join(summary_lines))

    print("\n".join(summary_lines))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_tsv", required=True, help="Phase-1 joined TSV (metadata + clone_id + clone_valid).")
    ap.add_argument("--out_dir", required=True, help="Output directory for Phase-B fate artifacts.")
    args = ap.parse_args()
    return run_phase_b(args.input_tsv, args.out_dir)


if __name__ == "__main__":
    raise SystemExit(main())
