#!/usr/bin/env python3
"""
viability_audit.py

Phase-1 viability audit for GEO GSE126954 (Packer et al.).

Purpose:
  - Quantify usable sample sizes under the *locked* Phase-1 exclusion rules
  - Provide fate (cell.type) and lineage founder-group distributions
  - Optionally sanity-check overlap between annotation cell IDs and the count-matrix columns

Phase-1 locked exclusion rules (conservative):
  1) Keep QC-whitelisted rows (if column exists)
  2) Drop rows with NA cell.type
  3) Drop rows with NA lineage
  4) Drop rows where lineage contains "/" (ambiguous lineage assignment)
  5) Parse founder group from lineage as one of {AB, MS, E, C, D}; drop if not parseable

Notes:
  - This script does NOT perform any modeling, clustering, PCA, or neighborhood tests.
  - It is strictly a viability/accounting audit to protect the GO/NO-GO gate.
"""

import gzip
import re
import sys
from typing import Optional

import pandas as pd


ANNOT_GZ_DEFAULT = "GSE126954_cell_annotation.csv.gz"
COUNT_MATRIX_GZ_DEFAULT = "GSE126954_gene_by_cell_count_matrix.txt.gz"


FOUNDERS = ("AB", "MS", "E", "C", "D")


def founder_group(lineage: object) -> Optional[str]:
    """
    Map a lineage string to a founder group among {AB, MS, E, C, D}.

    Rules:
      - If lineage starts with "AB" -> "AB"
      - Else if starts with "MS" -> "MS"
      - Else if starts with "E"  -> "E"
      - Else if starts with "C"  -> "C"
      - Else if starts with "D"  -> "D"
      - Else -> None

    We ignore any deeper prefix letters beyond the founder.
    """
    if not isinstance(lineage, str):
        return None
    s = lineage.strip()
    if not s or s.upper() == "NA":
        return None

    # Fast-path known founders
    if s.startswith("AB"):
        return "AB"
    if s.startswith("MS"):
        return "MS"
    if s.startswith("E"):
        return "E"
    if s.startswith("C"):
        return "C"
    if s.startswith("D"):
        return "D"

    # As a fallback, try first 2 letters or first letter
    # (kept conservative: only return if it matches known founders)
    m2 = re.match(r"^([A-Za-z]{2})", s)
    if m2 and m2.group(1) in FOUNDERS:
        return m2.group(1)

    m1 = re.match(r"^([A-Za-z])", s)
    if m1 and m1.group(1) in FOUNDERS:
        return m1.group(1)

    return None


def read_count_matrix_header_cell_ids(path_gz: str) -> Optional[set]:
    """
    Lightweight overlap check: read only the first line of the gene-by-cell matrix
    to extract the cell IDs (column headers). Returns a set of cell IDs, or None if
    the file is missing or unreadable.

    This function assumes a tab-delimited matrix with first column as gene ID.
    """
    try:
        with gzip.open(path_gz, "rt") as f:
            header_line = f.readline()
        if not header_line:
            return None
        header = header_line.rstrip("\n").split("\t")
        if len(header) < 2:
            return None
        # Common format: first column header is empty or a gene-id label; cell IDs follow.
        cell_ids = header[1:]
        return set(cell_ids)
    except FileNotFoundError:
        return None
    except Exception:
        return None


def pct(n: int, denom: int) -> float:
    return (100.0 * n / denom) if denom else 0.0


def main() -> int:
    annot_gz = ANNOT_GZ_DEFAULT
    count_gz = COUNT_MATRIX_GZ_DEFAULT

    # Allow optional argv overrides:
    #   python3 viability_audit.py [cell_annotation.csv.gz] [gene_by_cell_count_matrix.txt.gz]
    if len(sys.argv) >= 2:
        annot_gz = sys.argv[1]
    if len(sys.argv) >= 3:
        count_gz = sys.argv[2]

    df = pd.read_csv(annot_gz, compression="gzip")
    n0 = len(df)

    print("=== Viability audit: GSE126954_cell_annotation ===")
    print(f"Annotation file: {annot_gz}")
    print(f"Total rows:                          {n0:>10} ({pct(n0, n0):6.2f}%)")

    # QC whitelist
    qc_col = "passed_initial_QC_or_later_whitelisted"
    if qc_col in df.columns:
        df_qc = df[df[qc_col] == True].copy()
        n_qc = len(df_qc)
        print(f"QC-whitelisted rows:                 {n_qc:>10} ({pct(n_qc, n0):6.2f}%)")
    else:
        df_qc = df.copy()
        print("QC column not found; skipping QC filter")

    # Drop NA fate label
    fate_col = "cell.type"
    if fate_col not in df_qc.columns:
        print(f"\nERROR: required column missing: {fate_col}")
        return 2
    df_fate = df_qc[df_qc[fate_col].notna()].copy()
    n_fate = len(df_fate)
    print(f"After drop NA cell.type:             {n_fate:>10} ({pct(n_fate, n0):6.2f}%)")

    # Drop NA lineage
    lin_col = "lineage"
    if lin_col not in df_fate.columns:
        print(f"\nERROR: required column missing: {lin_col}")
        return 2
    df_lin = df_fate[df_fate[lin_col].notna()].copy()
    n_lin = len(df_lin)
    print(f"After drop NA lineage:               {n_lin:>10} ({pct(n_lin, n0):6.2f}%)")

    # Drop ambiguous lineage with "/"
    df_lin[lin_col] = df_lin[lin_col].astype(str)
    df_unamb = df_lin[~df_lin[lin_col].str.contains("/", regex=False)].copy()
    n_unamb = len(df_unamb)
    print(f"After drop ambiguous lineage ('/'):  {n_unamb:>10} ({pct(n_unamb, n0):6.2f}%)")

    # Parse founder group (AB/MS/E/C/D)
    df_unamb["founder"] = df_unamb[lin_col].map(founder_group)
    df_good = df_unamb[df_unamb["founder"].notna()].copy()
    n_good = len(df_good)
    print(f"After parse founder group:           {n_good:>10} ({pct(n_good, n0):6.2f}%)")

    # Sanity: how many were dropped due to founder parsing?
    dropped_founder = n_unamb - n_good
    print(f"Drop count (unparseable founder):    {dropped_founder:>10} ({pct(dropped_founder, n0):6.2f}%)")

    # Fate distribution
    print("\nTop cell.type counts (post-filters):")
    print(df_good[fate_col].value_counts().head(25).to_string())

    # Founder distribution
    print("\nFounder counts (post-filters):")
    print(df_good["founder"].value_counts().to_string())

    # Optional: time-bin distribution (useful for later stratification, still not modeling)
    timebin_col = "embryo.time.bin"
    if timebin_col in df_good.columns:
        print("\nembryo.time.bin counts (post-filters):")
        print(df_good[timebin_col].value_counts().sort_index().to_string())

    # Optional: overlap with expression matrix columns
    cell_col = "cell"
    if cell_col not in df_good.columns:
        print(f"\nWARNING: cell ID column missing: {cell_col} (overlap check skipped)")
        return 0

    cell_ids = read_count_matrix_header_cell_ids(count_gz)
    if cell_ids is None:
        print("\nExpression matrix overlap check: skipped (count matrix header not available).")
        print(f"Expected at: {count_gz}")
        return 0

    overlap = int(df_good[cell_col].isin(cell_ids).sum())
    print("\nExpression matrix overlap check:")
    print(f"Count-matrix file:                   {count_gz}")
    print(f"Cells in annotation (post-filters):  {n_good}")
    print(f"Cells also present in matrix cols:   {overlap}")
    print(f"Overlap fraction:                    {overlap / n_good if n_good else 0.0:.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())