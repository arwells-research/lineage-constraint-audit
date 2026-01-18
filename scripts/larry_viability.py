import argparse
import pandas as pd

def load_larry_clone_map(path: str) -> pd.DataFrame:
    # tab-delimited, no header: cell_barcode \t clone_id
    df = pd.read_csv(path, sep="\t", header=None, names=["cell_barcode", "clone_id"])
    df["clone_valid"] = df["clone_id"].astype("Int64").fillna(0) > 0
    return df

def load_larry_metadata(path: str) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t")
    # normalize expected join key
    if "Cell barcode" in meta.columns:
        meta = meta.rename(columns={"Cell barcode": "cell_barcode"})
    return meta

def run_viability(metadata_path: str, clone_map_path: str) -> int:
    import pandas as pd

    # Load metadata (tab-delimited)
    md = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if "Cell barcode" not in md.columns:
        raise SystemExit(
            f"ERROR: metadata missing required column 'Cell barcode'. Columns={list(md.columns)[:20]}"
        )

    md_bar = md["Cell barcode"].astype(str).str.strip()
    md_set = set(md_bar.unique())
    metadata_unique = len(md_set)

    # Load clone map (tsv: barcode \t clone_id), allow blank/missing clone_id rows
    cm = pd.read_csv(
        clone_map_path,
        sep="\t",
        header=None,
        names=["barcode", "clone_id"],
        dtype=str,
        keep_default_na=False,  # don't auto-convert "NA" to NaN; we will treat "NA" explicitly
    )
    cm["barcode"] = cm["barcode"].astype(str).str.strip()
    cm["clone_id"] = cm["clone_id"].astype(str).str.strip()

    cm_set = set(cm["barcode"].unique())
    clone_map_unique = len(cm_set)

    # Barcode coverage: does clone_map contain every metadata barcode?
    inter = md_set.intersection(cm_set)
    barcode_overlap_barcodes = len(inter)
    barcode_coverage_rate = barcode_overlap_barcodes / max(1, metadata_unique)

    # Clone assignment: among metadata barcodes, how many have a NON-missing clone_id?
    cm2 = cm.copy()
    cm2["clone_id_norm"] = cm2["clone_id"].astype(str).str.strip()
    cm2.loc[
        (cm2["clone_id_norm"] == "") | (cm2["clone_id_norm"].str.upper() == "NA"),
        "clone_id_norm",
    ] = ""

    # Choose the first non-empty clone_id per barcode (if multiple rows exist)
    cm2 = cm2.sort_values(["barcode", "clone_id_norm"])
    cm2 = cm2.groupby("barcode", as_index=False).agg({"clone_id_norm": "max"})

    # Left-join metadata barcodes to clone map
    md_keys = pd.DataFrame({"barcode": sorted(md_set)})
    j = md_keys.merge(cm2, on="barcode", how="left")

    clone_assigned_barcodes = int((j["clone_id_norm"].fillna("").astype(str).str.strip() != "").sum())
    clone_assigned_rate = clone_assigned_barcodes / max(1, metadata_unique)

    # --- Report (two-phase gate) ---
    # Phase-A1: joinability (deterministic barcode overlap)
    # Phase-A2: clone coverage (non-missing clone_id rate) — allowed to be < 1 if missingness is representable

    print("# LARRY viability report")

    print(f"metadata_rows={len(md)}")
    print(f"metadata_unique_barcodes={metadata_unique}")

    print(f"clone_map_rows={len(cm)}")
    print(f"clone_map_unique_barcodes={clone_map_unique}")

    print(f"barcode_overlap_barcodes={barcode_overlap_barcodes}")
    print(f"barcode_coverage_rate={barcode_coverage_rate:.6f}")

    print(f"clone_assigned_barcodes={clone_assigned_barcodes}")
    print(f"clone_assigned_rate={clone_assigned_rate:.6f}")

    # Explicit decision logic (make the intended GO gate unambiguous)
    # GO is primarily about *joinability* (A1). A2 is a quality/coverage signal.
    GO_A1 = barcode_coverage_rate >= 0.95
    GO_A2 = True  # missing clone_id is allowed if representable; keep as informational

    print("## Gate interpretation")
    print(f"GO_A1_joinable={'true' if GO_A1 else 'false'}  (requires deterministic join key coverage)")
    print(f"GO_A2_clone_coverage={'true' if GO_A2 else 'false'}  (missing clone_id allowed; represent via clone_valid=false)")
    print(f"GO={'true' if GO_A1 else 'false'}")

    return 0 if GO_A1 else 2

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clone_map", required=True, help="larry_barcode_to_clone.tsv")
    ap.add_argument("--metadata", required=True, help="GSM4185642_*_metadata.txt.gz")
    ap.add_argument("--out_join", default="results/larry/join_check.tsv")
    ap.add_argument("--out_report", default="results/larry/viability_report.txt")
    args = ap.parse_args()

    return run_viability(args.metadata, args.clone_map)

if __name__ == "__main__":
    raise SystemExit(main())
