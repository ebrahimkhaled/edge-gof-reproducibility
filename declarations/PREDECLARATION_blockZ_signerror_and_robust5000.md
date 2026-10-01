# Pre-declaration: block Z, the two closing experiments the robustness referee asked for (written before either is run)

Date: 2026-09-20, after a specialist robust-statistics review of the manuscript. Two gaps in the
evidence were identified, both cheap to close, and both are closed here under one document frozen before
either experiment runs.

- **Part A, the sign error.** The paper's headline corruption, C1, multiplies a covariate by four. That
  exaggerates a prediction that was already directionally right: because the outcome is drawn at the
  uncorrupted covariate, it contradicts the exaggerated prediction only about three times in ten. C1 is
  therefore a mild member of its own family. The classical worst case is a **sign error**, and it is not
  run anywhere in the study.
- **Part B, the robust fit at the paper's own boundary.** Block 9d compares maximum likelihood with a
  robust fit at `n = 1000` and at most twenty-five corrupted records. Section 6.5's boundary, where the
  directed test fails, is at `n = 5000` with twenty-five and fifty. The paper points forward from the
  boundary to a study that never reaches it.

## Part A: C1b, a sign error in one covariate

### A1. Design

`c1b_gen()` reproduces `b9_gen()` of block 9 exactly --- the same draw order, the same truths, the same
seeds --- with one change: the corrupted records have their covariate replaced by `-4x` rather than
`4x`. The outcome is still drawn from the truth at the **original** covariate, so a corrupted record
carries a high-risk patient's outcome at a low-risk recorded position, or the reverse: a bad leverage
point in the sense of Rousseeuw and Leroy.

The code asserts on every run that `c1b_gen()` with the corruption set to `C1` reproduces `b9_gen()`'s
output bit for bit; a failure stops the block. This is the identity gate for Part A.

Scenarios: the logistic truth, `n = 1000` and `n = 5000`, `k/n` in `0.001, 0.002, 0.005, 0.01`, plus the
matched clean scenarios --- ten scenarios, 1000 replicates each, on fresh seed bases
`700000000 + cell_id * 10000` that no block of this study has used. Every test of the battery is
computed by the battery's own code.

### A2. What would count against the paper

The paper claims that pooling protects a grouped test against corrupted covariates up to about ten
corrupted records, and that the protection is a property of the partition rather than of the corruption.
Part A goes against that claim if **the directed test at the rule G exceeds 0.10 at ten corrupted
records in 1000 or at ten in 5000** under the sign error --- the two scenarios where the paper's own C1
readings are 0.097 and 0.060. If it does, the scope sentence of Section 6.5 must be restated for the
mild member only, in the abstract as well as in Section 6, and the sign error reported as the case that
breaks it.

### A3. Written prediction

The mechanism is indifferent to the sign of the error: both extreme groups are affected, the corrupted
records still land in the extreme groups, and the group denominator is unchanged. What does change is
the probability that the outcome contradicts the recorded prediction, which rises from about 0.3 to
close to 0.8, so the per-record shift should be roughly three times larger. I therefore expect the
directed test to move **more** under C1b than under C1 --- perhaps 0.15 to 0.30 at ten corrupted records
in 1000, against 0.097 --- and the record-level tests to move further still. If that happens the
protection claim survives only in its partition form and must be re-scoped by corruption severity, which
is the honest outcome and is why the experiment is worth running.

## Part B: 9dX, the robust fit at n = 5000

### B1. Design

Block 9d's machinery unchanged --- the same two estimators on the same data sets, the same replacement
of `bt_fit`, the same battery --- with the configuration grid extended to `n = 5000` under the logistic
truth at `k = 25, 50, 100`, and the matched clean configuration at `k = 0`. Four configurations, two
estimators, eight cells, 1000 replicates each, seed bases `650000000 + config_id * 10000`.

The estimator is `robustbase::glmrob(method = "Mqle")` with `weights.on.x` at its package default,
`"none"`. That is the **Huber-type** quasi-likelihood estimator of Cantoni and Ronchetti at tuning
constant `c = 1.345`, not the Mallows-type one; the manuscript's present description of it as
Mallows-type is an error found by the same review and corrected in the same revision.

### B2. What would count against the paper

Section 6.6 concludes that a robust fit does not reduce any test's false alarms. Part B goes against
that conclusion if **the robust arm reduces the directed test's false-alarm rate at `n = 5000` by more
than 0.05 at any of `k = 25, 50, 100`**, or if it brings any record-level test back within
`0.05 +/- 3` nominal standard errors at `k = 25`. Either would mean the conclusion was an artefact of
the smaller sample, and Section 6.6 would have to say so.

Level is read first: a cell whose `k = 0` rejection rate under either estimator lies outside
`0.029` to `0.071` has its contaminated readings reported but not used for the conclusion, which is
block 9d's own gate D9d.1 applied unchanged.

### B3. Written prediction

The robust fit restores the slope, and the corrupted record's linear predictor is formed at its
corrupted covariate whatever the slope is, so restoring the slope makes that record **more** extreme,
not less. At `n = 1000` and `k = 25` this already shows: the Hosmer--Lemeshow statistic rises from
0.114 to 0.171 under the robust fit. I therefore expect the robust arm at `n = 5000` to be no better
than maximum likelihood at every `k`, and probably worse, with the gap widening in `k`.

## Rules for both parts

Rejection at `p <= 0.05`. A replicate with no p-value counts as no rejection and is reported. The
Monte Carlo standard error of a 5% rate at `B = 1000` is 0.007, and readings are compared against three
of those. Nothing here changes the default test or any existing declaration; both parts are additions
that can only add rows to the record.

## Outputs

`simulations/battery/C1b/` and `simulations/battery/9dX/`, one per-replicate file per cell, with an
analysis read against this document.
