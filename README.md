# Lineage Constraint Audit in Single-Cell Systems

**A constraint-first statistical audit of lineage detectability in single-cell expression space.**

---

## Summary

**Result:** After controlling for developmental time and batch, lineage information remains locally detectable in expression space, with constraint strength **30–60% above stratified null expectations** in **validity-positive neighborhoods (~3–11% of local contexts)**.

Validity-positive neighborhoods are those in which multiple founder lineages are locally represented; constraint strength is evaluated **conditionally** on this validity criterion.

This repository implements a rigorous statistical audit asking a narrowly defined question:

> **Does lineage information remain locally detectable in single-cell expression space after controlling for developmental time and batch effects?**

Rather than proposing a predictive model or biological mechanism, this project evaluates **detectability and constraint structure** under explicit null hypotheses, with predefined GO / NO-GO decision gates and negative controls.

---

## What this project is

- A **statistical audit**, not a predictive or mechanistic model
- A **constraint-first analysis**, asking what structural regularities *must* exist if lineage matters
- A **falsification-oriented pipeline**, designed to accept negative outcomes as informative
- A **local neighborhood analysis**, testing detectability where biological constraints are most plausible

The core output is not a classifier, trajectory, or gene list, but a set of **validated statistical statements** about when and where lineage-associated structure persists in expression space.

---

## What this project is not

- Not a claim that lineage causes observed expression structure
- Not a developmental model or fate-prediction framework
- Not a gene regulatory network inference
- Not a global correlation or embedding-wide association test

No mechanistic interpretation is asserted beyond detectability under controlled nulls.

---

## Methodological overview

The analysis proceeds through a sequence of **audit stages**, each with explicit decision gates:

- **Audit Stage 1: Global validity and strength**  
  Does lineage information remain locally detectable at all?

- **Audit Stage 2: Robustness and sensitivity checks**  
  Does the signal persist across neighborhood size, sampling strategy, and random seeds?

- **Audit Stage 3: Negative controls**  
  Does the signal vanish when lineage or fate labels are deliberately destroyed?

- **Audit Stage 4: Contextual stratification**  
  Is detectability structured by developmental time or cell type?

- **Audit Stage 5A: Contribution decomposition**  
  Are detected constraints dominated by a single founder lineage or distributed across founders?

All null models are stratified by developmental time and batch, and all GO / NO-GO criteria are defined a priori.

## Phase-D / Phase-D2 — Founder admissibility diagnostics

Beyond detectability, the audit evaluates whether founders occupy distinct
admissible developmental regions or preferentially share local neighborhoods.

These diagnostics are descriptive only and do not alter Stage-1/2 decisions.

Dataset | Founders | Phase-D (dispersion) | Phase-D2 (pairwise overlap)
------- | -------- | -------------------- | ---------------------------
GSE126954 (worm) | 5 | **NO-GO** — no founder-conditioned time dispersion beyond stratified null | **NO-GO** — founder pair overlap matches null expectations
LARRY (mouse) | 2 | **GO** — founders occupy distinct time-dispersion structure (p ≈ 0.001) | **GO** — founders partition local admissibility space

**Interpretation:**  
Founder-conditioned admissibility is supported in the LARRY system, but not
resolved in the worm system under the present constraints. Negative outcomes are
retained intentionally as part of the audit’s falsification-first philosophy.

---

### Phase-E (Founder-conditioned admissibility) — LARRY

On LARRY (mouse hematopoiesis), Phase-E evaluates whether founder identity constrains
clone-internal developmental admissibility beyond time-stratified null expectations.

Score:
S = mean_c ( w_c · purity_c )

where:
  w_c = 1 − H(Time|clone=c) / log₂(n_times)
  purity_c = max_f P(founder=f | clone=c)

Null shuffles founder labels within time strata.

**Canonical Phase-E diagnostic panel (LARRY):**

| Mode | Description | S_obs | S_null_median | p | Decision |
|------|-------------|-------|---------------|----|----------|
| baseline | Standard Phase-E | **0.3479** | 0.3152 | 0.001 | **GO** |
| sanity_time_within_clone | Preserve time histograms | 0.3479 | 0.3046 | 0.001 | GO |
| sanity_founder_within_clone | Preserve founder histograms | 0.3479 | 0.3046 | 0.001 | GO |
| destruct_time_across | Break clone↔time coupling | 0.1447 | 0.1275 | 0.001 | GO (magnitude reduced) |
| destruct_time_resample | Resample time per clone | 0.1443 | 0.1273 | 0.001 | GO (magnitude reduced) |
| destruct_founder_across | Break clone↔founder coupling | **0.2988** | 0.3002 | 0.916 | **NO-GO** |
| diag_perm_w | Break w↔purity alignment | 0.3309 | 0.3029 | 0.001 | GO |

**Interpretation:**

- Founder-conditioned admissibility is **robust in LARRY.**
- Destroying founder identity across clones eliminates the effect (**NO-GO**).
- Destroying time structure changes the score magnitude but does not eliminate detectability.
- These are structural diagnostics only and do not imply mechanism.

---

## Viability boundaries

This audit distinguishes between:

- whether lineage detectability exists, and
- whether the data supports estimating it.

For datasets with sparse founder mixing (e.g., early embryogenesis), a minimum number
of validity-positive neighborhoods is required before constraint strength can be
meaningfully evaluated.

**Empirical rule:**

At least **100 validity-positive neighborhood centers** are required to pass the
Stage-1 viability gate.  

Given typical validity rates (~3–4%), this implies:

```
max_centers ≥ 5000
```

Runs with fewer sampled centers are marked **Viability NO-GO by design**, rather than
being treated as evidence against detectability.

---

## Verified audit pillars

This framework has been validated across two independent lineage-tracing systems:

Dataset | Organism | Mapping Quality | Status
------- | -------- | --------------- | ------
GSE126954 | *C. elegans* embryogenesis | Lineage provided | GO
LARRY (Weinreb et al., 2020) | Mouse hematopoiesis | Deterministic barcode-to-clone table | GO

---

## Interpretation boundaries

This project establishes statistical detectability, not biological causation.

Specifically:

- Lineage-associated structure may reflect ancestry, spatial proximity, or their entanglement.
- No claim is made that lineage directly regulates transcriptional state.
- Absence of signal under negative controls is treated as successful validation, not failure.

Results should be interpreted as statements about residual lineage-associated information that persists in local expression neighborhoods after accounting for developmental time and batch.

---

## Known limitations

- Lineage and spatial organization may be entangled.
- Findings depend on dataset quality and joinability.
- The audit does not infer gene regulatory mechanisms.
- The local-neighborhood framework may miss global coordination effects.

These limitations motivate future extensions rather than undermine the present results.

---

## Repository structure

    .
    ├── scripts/        # Core audit pipeline (Audit Stages 1–5A)
    ├── data/           # Download + checksums (raw files exist locally but are not committed)
    ├── docs/           # Detailed methodology, null models, and decision criteria
    ├── results/        # Computed statistics, tables, and figures
    └── README.md

---

## Dataset sources

This audit has been run on:

- GSE126954 — *C. elegans* embryogenesis at single-cell resolution (Packer et al., 2019)
- GSM4185642 — LARRY hematopoietic lineage tracing dataset (Weinreb et al., 2020)

No additional preprocessing beyond documented filters is performed.

---

## Who should use this repository?

This framework is designed for researchers who want to:

- Test whether lineage information remains detectable in their own single-cell datasets
- Implement rigorous stratified null models for local neighborhood statistics
- Distinguish detectability claims from mechanistic claims
- Build constraint-first analyses with explicit decision gates and negative controls

This repository is not intended for trajectory inference, fate prediction, or regulatory network reconstruction.

---

## Philosophy

This project is built on the principle that negative results are informative, provided the null hypotheses are explicit and the decision criteria are enforced.

The goal is not to demonstrate that lineage matters, but to test whether it remains detectable at all once major confounders are removed.

If the answer had been “no,” that outcome would have been reported with equal weight.

---

## Citation

If you use this audit framework or adapt the methodology, please cite:

Packer JS et al. (2019).  
*A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution.* Science.

Weinreb C et al. (2020).  
*Lineage tracing on transcriptional landscapes links state to fate.* Science.


A `CITATION.cff` file is provided for convenience.
