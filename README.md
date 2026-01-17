# Lineage Constraint Audit — Local Detectability Under Stratified Nulls

> **Repository status:**  
> This repository contains (i) a completed lineage-detectability audit in *C. elegans* embryogenesis (GSE126954), and  
> (ii) an active extension testing spatial vs lineage contributions using mouse scRNA-seq and Visium spatial data (GSE153424).  
>  
> Only results explicitly labeled as complete should be interpreted as validated findings.  
> The normative definitions, null construction, and decision criteria are specified in `docs/audit_logic.md`.

---

## Completed audit summary (C. elegans)

**Result (GSE126954):** After controlling for developmental time and batch, lineage information remains locally detectable in expression space, with constraint strength **~30–60% above stratified null expectations** in **validity-positive neighborhoods (~3–11% of local contexts)**.

Validity-positive neighborhoods are those in which multiple founder lineages are locally represented; constraint strength is evaluated **conditionally** on this validity criterion (see `docs/audit_logic.md` for formal definitions of validity rate `v`, neighborhood statistic `Tᵢ`, and aggregate strength `A`).

This repository implements a rigorous statistical audit asking a narrowly defined question:

> **Does lineage information remain locally detectable in single-cell expression space after controlling for developmental time and batch effects?**

Rather than proposing a predictive model or biological mechanism, this project evaluates **detectability and constraint structure** under explicit null hypotheses, with predefined GO / NO-GO decision gates and negative controls.

---

## What this project is

- A **statistical audit**, not a predictive or mechanistic model
- A **constraint-first analysis**, asking what structural regularities *must* exist if lineage labels remain informative
- A **falsification-oriented pipeline**, designed to accept negative outcomes as informative
- A **local neighborhood analysis**, testing detectability where biological constraints are most plausible

The core output is not a classifier, trajectory, or gene list, but a set of **validated statistical statements** about when and where lineage-associated structure persists in expression space.

---

## What this project is not

- ❌ Not a claim that lineage *causes* observed expression structure
- ❌ Not a developmental model or fate-prediction framework
- ❌ Not a gene regulatory network inference
- ❌ Not a global correlation or embedding-wide association test

No mechanistic interpretation is asserted beyond detectability under controlled nulls.

---

## Methodological overview

The analysis proceeds through a sequence of **audit stages**, each with explicit decision gates (specified normatively in `docs/audit_logic.md`):

- **Audit Stage 1: Global validity & strength**
  Does lineage information remain locally detectable at all?

- **Audit Stage 2: Robustness & sensitivity checks**
  Does the signal persist across neighborhood size, sampling strategy, and random seeds?

- **Audit Stage 3: Negative controls**
  Does the signal vanish when lineage or fate labels are deliberately destroyed?

- **Audit Stage 4: Contextual stratification**
  Is detectability structured by developmental time or cell type?

- **Audit Stage 5A: Contribution decomposition**
  Are detected constraints dominated by a single founder lineage or distributed across founders?

All null models are **stratified by developmental time and batch**, and all GO / NO-GO criteria are defined *a priori*.

---

## Interpretation boundaries

This project establishes **statistical detectability**, not biological causation.

Specifically:

- Lineage-associated structure may reflect ancestry, spatial proximity, or their entanglement in *C. elegans* embryogenesis.
- No claim is made that lineage directly regulates transcriptional state.
- No mechanistic or regulatory interpretation is inferred from the observed constraints.
- Absence of signal under negative controls is treated as a successful validation, not a failure.

Results should be interpreted as statements about **residual lineage-associated information that persists in local expression neighborhoods after accounting for developmental time and batch**.

---

## Known limitations (completed audit)

- Lineage and spatial organization are entangled in *C. elegans*; this analysis does not disentangle ancestry-driven from spatially mediated structure.
- Findings are based on a single large-scale embryonic scRNA-seq dataset.
- The audit does not identify specific genes or regulatory mechanisms driving constraints.
- The local-neighborhood framework may miss global coordination effects operating at larger spatial or temporal scales.

These limitations motivate extensions rather than undermining the completed audit.

---

## Active extension: spatial decoupling audit (mouse brain)

An active extension of this audit framework is underway using mouse scRNA-seq and Visium spatial transcriptomics data (GSE153424).

The goal of this extension is **not** to assert new biological conclusions, but to test whether lineage-associated detectability observed in *C. elegans* collapses under null models that explicitly control for spatial proximity.

At present:

- Dataset organization and feasibility checks are in progress
- No claims or results from the mouse dataset are asserted
- Proceeding past feasibility requires explicit, mappable clone/lineage labels

See `SESSION_HANDOFF.md` for the authoritative phase plan, decision gates, and implementation constraints.

---

## Repository structure

    .
    ├── scripts/        # Core audit pipeline (Stages 1–5A; dataset adapters included)
    ├── data/           # Download scripts + checksums (raw files exist locally but are not committed)
    │   ├── download_gse126954.sh
    │   ├── GSE126954.sha256
    │   └── raw/        # Local-only (gitignored): GEO supplementary files
    ├── docs/           # Normative definitions, null models, and decision criteria
    ├── results/        # Computed statistics, tables, and figures (tracked; large blobs ignored)
    └── README.md

---

## Visual summary

Placeholder for key audit visualization.

Example figure:
- Distribution of constraint strength in validity-positive neighborhoods
- Comparison against stratified null distributions

See `results/figures/` for completed audit visualizations.

---

## Datasets

### Completed audit
- *C. elegans* embryogenesis scRNA-seq
- GEO accession: **GSE126954**
- Includes lineage annotations, developmental timing, batch metadata, and UMI count matrices

### Active extension
- Mouse brain scRNA-seq and Visium spatial transcriptomics
- GEO accession: **GSE153424**
- Included for spatial vs lineage decoupling audit; results pending feasibility

No additional preprocessing beyond documented filters is performed.

---

## Who should use this repository?

This audit framework is designed for researchers who want to:

- Test whether lineage or clonal information remains detectable in single-cell datasets
- Implement rigorous stratified null models for local neighborhood statistics
- Distinguish **detectability claims** from **mechanistic claims**
- Build constraint-first analyses with explicit decision gates and negative controls

This repository is **not** intended for trajectory inference, fate prediction, or regulatory network reconstruction.

---

## Philosophy

This project is built on the principle that **negative results are informative**, provided the null hypotheses are explicit and the decision criteria are enforced.

The goal is not to demonstrate that lineage “matters,” but to test whether it remains **detectable at all** once major confounders are removed.

If the answer had been “no,” that outcome would have been reported with equal weight.

---

## Citation

If you use this audit framework or adapt the methodology, please cite this repository and the original dataset for the completed audit:

Packer JS, et al. (2019).
*A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution.*
Science, 365(6459).

A `CITATION.cff` file is provided for convenience.