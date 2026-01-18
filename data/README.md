# data/README.md

This repository intentionally does **not** commit raw GEO matrices, images,
or extracted working data.

Reproducibility is ensured via **scripted download**, **local extraction**,
and (where applicable) **checksum verification**.

Normative audit logic and decision gates are specified in `docs/audit_logic.md`.

---

## Completed audit dataset

### GSE126954 — *C. elegans* embryogenesis scRNA-seq

This dataset supports the **completed lineage detectability audit**.

Phases applied:
  - **Phase-D:** founder-conditioned developmental dispersion (evaluated; NO-GO)
  - **Phase-D2:** pairwise founder neighborhood overlap (evaluated; NO-GO)

Phase-E (founder-conditioned admissibility) is implemented and validated
only for the **LARRY dataset** and is not applied to GSE126954.

Only the supplementary processed matrices required by the audit pipeline
are used; no FASTQs or raw sequencing data are required.

### Required files

Place the following GEO supplementary files in:

- `data/raw/`

1. `GSE126954_cell_annotation.csv.gz`
2. `GSE126954_gene_by_cell_count_matrix.txt.gz`

Optional (not required for this audit pipeline as written):

- `GSE126954_gene_annotation.csv.gz`
- `GSE126954_RAW.tar`

---

### One-command download

From repo root:

    bash data/download_gse126954.sh

This will:

- create `data/raw/`
- download the two required files
- (optionally) verify checksums if `data/GSE126954.sha256` is populated

---

### Checksum verification (optional)

This repo does not require checksum tools. If you want verification, install
`sha256sum` (often via `coreutils`) and add your own checksum step locally.

---

### Quick smoke-check (optional)

Confirm annotation header:

    zcat data/raw/GSE126954_cell_annotation.csv.gz | head -n 2

Confirm MatrixMarket header:

    zcat data/raw/GSE126954_gene_by_cell_count_matrix.txt.gz | head -n 2

Expected first two lines:

- `%%MatrixMarket matrix coordinate integer general`
- `<n_genes> <n_cells> <nnz>`

---

### Next step: run Phase-1 (completed audit)

Canonical command (from repo root):

    python3 scripts/phase1_runner.py \
      --annot data/raw/GSE126954_cell_annotation.csv.gz \
      --mtx   data/raw/GSE126954_gene_by_cell_count_matrix.txt.gz \
      --knn_k 20 --svd_components 50 --max_centers 5000 \
      --global_perms 1000 --seed 0



