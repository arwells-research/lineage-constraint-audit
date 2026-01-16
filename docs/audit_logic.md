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

### Audit Stage 1: Global validity & strength

Tests whether lineage information is detectable at all.

**GO criteria** (at least one must be satisfied):

1. Validity rate:
   - v_obs > v_null_95
2. Constraint strength:
   - A_obs > A_null_95

**Quantitative thresholds:**

- Minimum effect size:
  - A_obs / median(A_null) ≥ 1.2 (≥20% uplift)
- Minimum validity rate:
  - v_obs ≥ 0.03 (≥3% of neighborhoods)

**NO-GO outcome:**
If neither criterion is satisfied, conclude that lineage information is not locally
detectable under these controls. Subsequent stages are not performed.

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