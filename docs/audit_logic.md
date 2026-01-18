# Audit Logic and Decision Framework

This document specifies the statistical logic, null hypotheses, and decision criteria
underlying the **Lineage Constraint Audit**.

It is a **normative specification**: changes to logic, thresholds, or criteria require
explicit justification.

The goal is to test **detectability** of lineage-associated structure under controlled
nulls—not to infer mechanism or causation.

---

## Core question

Does lineage information remain **locally detectable** in single-cell expression space
after controlling for:

- developmental time
- batch effects

Detectability is evaluated without assuming lineage causality or regulatory mechanism.

---

## Definitions

### Local neighborhood

A local neighborhood is defined as the k nearest neighbors of a cell in a reduced
expression embedding.

All statistics are computed at the neighborhood level.

---

### Validity

A neighborhood is **validity-positive** if multiple founder lineages are locally
represented.

Validity is a *structural prerequisite* for evaluating lineage constraint.
Neighborhoods failing validity are excluded by construction.

Validity rate (v) is defined as:

- v = fraction of tested neighborhood centers that are validity-positive

---

### Constraint strength

Constraint strength (Tᵢ) quantifies the deviation of observed lineage composition
from expectation under a stratified null model.

Strength is evaluated **conditionally on validity**.

Aggregate constraint strength (A) is defined as:

- A = mean(Tᵢ | validity-positive neighborhoods)

---

## Statistical significance testing

All statements that an observed statistic “exceeds null expectation” mean:

- Observed statistic > 95th percentile of the stratified permutation null distribution
- Equivalent to a one-sided permutation test with α = 0.05
- Null distributions generated from N ≥ 1000 permutations unless otherwise stated

Both validity rate (v) and constraint strength (A) are evaluated independently
against their respective null distributions.

---

## Permutation null construction

For each null iteration:

1. Group cells by (developmental time bin × batch) strata
2. Randomly shuffle lineage labels within each stratum
3. Preserve expression-space neighborhoods and embeddings
4. Recompute:
   - validity rate
   - constraint strength statistics

This preserves:

- neighborhood geometry in expression space
- time–batch structure
- cell-type distributions insofar as they correlate with time/batch

This destroys:

- lineage–expression coupling (the target signal)

---

## Audit stages

## Audit Stage 1 — Detectability and Measurement Viability

### Purpose

Establish whether lineage information is locally detectable in expression space,
given the defined validity criterion and stratified null.

This stage separates:

1) whether the statistic is meaningful  
2) whether lineage detectability is present  
3) how strong the signal is  

---

### Stage-1A: Viability Gates (binary)

These ensure the measurement can be interpreted.

A run is viability-GO if all are satisfied:

Gate | Definition | Purpose
---- | ---------- | -------
V1 | At least one neighborhood meets the validity rule | Ensures founder mixing exists at all
V2 | n_valid_centers ≥ 100 | Ensures statistical estimation is meaningful
   | (empirical boundary: if v_obs ≈ 0.03–0.04, require max_centers ≥ 5000) |
V3 | Stratified null permutations completed successfully | Ensures comparison against null is valid

If any viability gate fails → **NO-GO by design.**

No detectability or strength statements are made.

---

### Stage-1B: Detectability Test (binary)

Performed only if Stage-1A is satisfied.

Lineage is considered detectability-GO if either holds:

1) Validity detectability  
 v_obs > v_null_95  

2) Strength detectability  
 A_obs > A_null_95  

If neither holds → **detectability NO-GO.**

Subsequent audit stages are not performed.

---

### Stage-1C: Strength Grading (reported, not gated)

These quantify the strength of the detected constraint.

They are descriptive only and never determine GO/NO-GO.

Metric | Meaning
------ | -------
Uplift ratio = A_obs / median(A_null) | Effect size relative to null
Validity rate v_obs | Fraction of neighborhoods where founder mixing exists
Share concentration (if computed) | Whether a single founder dominates constraint strength

These values must always be reported, regardless of GO/NO-GO outcome.

---

### Interpretation

Outcome | Meaning
------- | -------
Viability NO-GO | Cannot evaluate detectability; data does not support the audit stage
Detectability NO-GO | No lineage-associated structure detectable beyond stratified null
Detectability GO | Lineage-associated structure is statistically detectable
Strength grading | Quantifies effect magnitude, not its validity

---

### Design Rationale

Detectability is a binary scientific question.  
Viability is a data-quality prerequisite.  
Strength is a continuous measure and should never stifle valid science.

This ensures that:

- borderline effects are preserved  
- only non-interpretable analyses are rejected  
- subsequent robustness and control stages are justified by detectability  

---

### Audit Stage 2: Robustness & sensitivity

Evaluates whether detectability is robust to reasonable parameter variation.

**Parameter variations tested:**

- Neighborhood size: k ∈ predefined range
- Random seeds: multiple independent draws
- Center sampling strategies:
  - uniform random
  - stratified by time
  - density-weighted (if applicable)

**Robustness criterion (GO):**

- Stage 1 GO criteria satisfied in ≥80% of parameter combinations
- Direction of effect is consistent across all tests:
  - v_obs > v_null
  - A_obs > A_null

**NO-GO outcome:**
If signal is parameter-dependent or inconsistent in direction, conclude that
detectability is not robust.

---

### Audit Stage 3: Negative controls

Deliberately destroys lineage-associated structure while preserving time/batch strata.

Examples include:

- shuffling lineage labels
- shuffling fate labels
- combined shuffles within strata

**Expected outcome:**

- loss of detectability
- failure of Stage 1 GO criteria

Successful failure is treated as validation of the audit framework.

---

### Audit Stage 4: Contextual stratification

Evaluates detectability and strength within subsets defined by:

- developmental time bins
- cell-type categories

This stage is descriptive and does not alter global GO / NO-GO decisions.

---

### Audit Stage 5A: Signal decomposition by founder lineage

Decomposes aggregate constraint strength by founder lineage contribution.

Used to assess whether detected signal is:

- dominated by a single lineage
- or distributed across multiple founders

---

## Audit Stage D — Founder-conditioned admissibility dispersion
This stage is **descriptive only** and does not alter Stage-1 or Stage-2 GO / NO-GO decisions.

**Purpose:**  
Test whether different founder lineages occupy distinct admissible developmental
transition regions, rather than merely contributing unequally.

This extends Stage-1 detectability by examining founder-specific admissibility structure.

**Method:**  
Within validity-positive neighborhoods, compute developmental time dispersion
(entropy or variance) conditioned on founder identity.  
Compare against the same stratified founder-shuffle null used in Stage-1.

**Interpretation:**

Dataset | Founders | Outcome
------- | -------- | -------
**LARRY (mouse)** | 2 | **GO** — founders show distinct dispersion structure (p ≈ 0.001)
**GSE126954 (worm)** | 5 | **NO-GO (detectability)** — no dispersion signal beyond stratified null

A NO-GO at this stage indicates that founder-specific dispersion cannot be
distinguished from developmental time or batch effects alone.

---

## Audit Stage D2 — Pairwise founder neighborhood overlap
This stage is **descriptive only** and does not alter Stage-1 or Stage-2 GO / NO-GO decisions.

**Purpose:**  
Test whether certain founder pairs preferentially share admissibility space in
local expression neighborhoods.

**Method:**  
Measure the frequency with which founder pairs co-occur within validity-positive
neighborhoods.  
Evaluate against the same stratified founder-shuffle null as Stage-1.

**Interpretation:**

Dataset | Outcome
------- | -------
**LARRY (mouse)** | **GO** — founder pairs partition local admissibility space
**GSE126954 (worm)** | **NO-GO (detectability)** — founder pair overlap matches stratified null

A NO-GO implies that founder adjacency on the lineage tree does not yield
measurable overlap beyond time and batch controls.

---

## Example decision pathways

**Scenario A: Strong signal**

- Stage 1: GO (validity and/or strength exceed null)
- Stage 2: GO (robust across parameters)
- Stage 3: GO (negative controls fail as expected)

Conclusion:
Lineage-associated structure is robustly detectable.

---

**Scenario B: Null result**

- Stage 1: NO-GO

Conclusion:
No evidence for local lineage detectability under these controls.
Subsequent stages are not performed.

---

**Scenario C: Fragile signal**

- Stage 1: GO
- Stage 2: NO-GO (signal lost under parameter variation)

Conclusion:
Evidence insufficient; signal may reflect methodological sensitivity.

---

**Scenario D: Failed negative control**

- Stage 1–2: GO
- Stage 3: NO-GO (negative controls retain signal)

Conclusion:
Signal likely spurious; audit framework or assumptions require revision.

---

## Interpretation boundaries

This audit establishes:

- statistical detectability
- residual information content
- constraint structure under null control

It does not establish:

- causation
- regulatory mechanism
- developmental necessity

Lineage-associated structure may reflect ancestry, spatial proximity, or their
entanglement.

---

## Reporting requirements

All audit reports must include:

1. Stage 1 results:
   - v_obs and v_null distribution
   - A_obs and A_null distribution
2. Explicit GO / NO-GO decision with cited criteria
3. Stage 2 robustness summary across parameter sweeps
4. Stage 3 negative control outcomes (expected vs. observed)
5. Quantitative effect sizes (uplift over null, not only p-values)
6. Failure modes:
   - which criterion failed
   - why it failed

Selective reporting of only positive outcomes violates audit principles.

---

## Frequently misinterpreted

**Q: Does “GO” mean lineage causes expression differences?**  
A: No. GO means lineage labels carry detectable information in local neighborhoods.
This may reflect ancestry, spatial organization, or other lineage-correlated factors.

**Q: Does “NO-GO” mean lineage is biologically irrelevant?**  
A: No. NO-GO means lineage is not detectable in local expression neighborhoods under
these specific controls.

**Q: Why not use a classifier to test detectability?**  
A: Classifiers optimize global accuracy and may exploit spurious correlations.
This audit tests local constraint structure under explicit nulls.

**Q: Can this be used to identify fate-determining genes?**  
A: No. This audit does not identify genes, mechanisms, or causal relationships.

---

## Decision philosophy

Negative results are informative.

All decision criteria are defined prior to evaluation.
Post-hoc rationalization is explicitly disallowed.

Phase-D, Phase-D2, and Phase-E provide additional structural diagnostics
and do not constitute independent evidence for or against lineage detectability.