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

---

## Active extension dataset (no results asserted)

### GSE153424 — Mouse brain scRNA-seq + Visium spatial transcriptomics

This dataset is included for an **active extension** of the audit framework
testing whether lineage / clonal detectability collapses under **explicit
spatial control**.

At present, GSE153424 is used **only for feasibility and viability checks**.
No biological or statistical claims from this dataset are asserted.

---

### Data policy (GSE153424)

The following directories are **local-only** and must remain gitignored:

- `data/gse153424/raw/`   — downloaded GEO / SRA archives
- `data/gse153424/work/`  — extracted matrices, spatial files, derived tables
- `data/gse153424/results/` — intermediate or exploratory outputs

Only **small, human-authored files** (scripts, READMEs, checksums) may be
committed under `data/gse153424/`.

---

### Expected local layout (GSE153424)

After local download and extraction, a typical working tree will look like:

    data/gse153424/
      raw/
        GSE153424_family.soft.gz
        GSE153424_family.xml.tgz
        scRNA/*.tar.gz
        spatial/*.tar.gz
      work/
        extracted/
          scRNA/<GSM...>/
          spatial/<GSM...>/
        derived/
        logs/

Exact contents vary by GSM and modality.

---

### Current status

- Dataset organization scripts exist (see repo root).
- Clone / lineage label availability is **not yet confirmed**.
- Proceeding beyond feasibility requires explicit, mappable clone labels
  aligned to scRNA cells and/or Visium spots.

See `SESSION_HANDOFF.md` for the authoritative phase plan and GO / NO-GO gates
for this dataset.

---

## Important reminders

- **Do not commit raw or extracted data.**
- **Do not infer results from GSE153424 until feasibility is explicitly GO.**
- **Do not modify audit logic to accommodate a dataset.**

Negative outcomes (NO-GO) are valid and expected under this framework.
