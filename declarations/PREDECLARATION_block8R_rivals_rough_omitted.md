# Pre-declaration: block 8R, the three outside rivals on omitted terms and rough misfit (written before any 8R scenario is run)

Date: 2026-09-19. Author's decision of the same day (option 1 of the three offered).

## 0. Why

Blocks 8 and 8L put le Cessie's smoothed test, the Liu projection test and BAGofT on twenty
scenarios from four departure families (symmetric tails, asymmetric links, the Stukel family, one
off-index scenario). None came from the two families where the directed test is weakest -- omitted
terms and rough misfit, where its declared hypothesis H3 went against it -- and where these rivals
may be strongest: the Liu projection test was designed on quadratic and interaction departures. This
block fills that gap, so every rival is judged on the same families as the directed test, the
unfavourable ones included.

## 1. Scenarios, fixed here by rule and not by any test's power

All scenarios are the battery's own (`battery_cells()`), families 4 (omitted terms) and 5 (rough
misfit) of `rule_A_membership.csv`.

- **le Cessie:** every alternative of the two families with n <= 1000 -- 48 omitted-term and 8 rough
  scenarios -- and their 11 matched nulls; replicates 1-500.
- **Liu projection test:** every one of those alternatives at n = 500 -- the interaction ladder
  `binint` (5 effect sizes), `contint` (4), the quadratic ladder `quad` (7) and the rough shapes
  `rough_osc2`, `rough_osc4`, `rough_sawtooth` -- 19 alternatives, and their 4 matched nulls;
  replicates 1-500.
- **BAGofT:** at n = 500, the three rough shapes and, from each omitted-term ladder, its median rung
  (the upper of the two middle rungs when the count is even): `binint_0.3`, `contint_0.5`,
  `quad_0.05`; 6 alternatives, and the same 4 nulls; replicates 1-200.

n = 500 for the two resampling tests is set by cost measured in block 8: at n = 1000 they take a
median of 291 and 546 seconds a data set, at n = 500 about 20 and 377.

## 2. Implementations

- **le Cessie:** exactly block 8L's, `run.all.gof(fit, tests = "le-Cessie")` of ebrahim.gof >= 2.6.0,
  checked against its definition (the corrected (I - H)' R (I - H)) before the run.
- **Liu projection test:** exactly block 8's, `proj_pvalue(y, X, B = 250)` of `_proj_test.R`, X the
  model matrix of the working fit.
- **BAGofT:** BAGofT 1.0.0 at its package defaults, with one line changed. The package's partition
  function `parRF` takes the covariates with `Train.data[, -which(names(Train.data) == Rsp)]`; with a
  single covariate R drops that to a vector and the call stops with an error, so the package cannot
  test any model with one covariate -- the quadratic and rough scenarios here, and their nulls. The
  line is given `drop = FALSE` (the fix of the DeepGOF study's `parRF_dropfix.R`) and nothing else
  changes. With two or more covariates the change is a no-op; the self-test checks that bit for bit on
  two-covariate data before the run, and checks that the unpatched function does fail on a
  one-covariate model. The patched function is used in every BAGofT call of this block.
- **Data:** nothing is simulated afresh. Every data set is regenerated from the seed its block stored,
  with block 8's identity gate (EDGE-poly3 unit at G = 10 equal to the stored value to 1e-8; seed, n and
  events equal). Each rival continues the replicate's random-number stream, as in block 8.
- **Order:** le Cessie, then the projection test, then BAGofT; 20 workers; the block starts only after
  the DeepGOF study's confirmatory run has finished and no R process has run for ten minutes.

## 3. Rules

As blocks 8 and 8L. The two resampling tests have p-values on a grid and reject at p < alpha; le
Cessie and the battery's tests at p <= alpha. No p-value counts as no rejection, and its rate is
reported. Size is checked against 0.05 +/- 3 nominal standard errors: 0.021-0.079 at B = 500, and
0.004-0.096 at B = 200. Power is size-adjusted with the critical value read from the matched null of
the same test on its own replicates, EDGE likewise on the same replicates. A rival leads or ties with
the margin used throughout the paper (EDGE size-adjusted power > rival - 0.01). Every rival is paired
with EDGE-poly3 (unit, rule G) on the same replicates by an exact McNemar test at 0.05, Holm-corrected
across all pairs of block 8R.

## 4. Claims

- **R1 (level).** Each rival's rejection rate in each of its null scenarios lies within the band of
  section 3. *Against:* any null scenario outside it; the power in that null's alternatives is then
  reported with the failure beside it.
- **R2 (reported, no direction predicted).** For each rival and each of the two families, the number of
  scenarios in which the directed test leads or ties, the paired comparisons, and the mean
  size-adjusted power of both. No direction is predicted: the directed test is expected to be weak
  here, but whether a smoother with a wide default bandwidth, or a cumulative projection statistic,
  can follow a high-frequency oscillation is not known in advance.
- **R3 (reported).** Seconds per data set of each rival; the number of data sets on which any rival
  returns no p-value.

## 5. Outputs

`simulations/battery/8R/`: one per-replicate file per scenario and test, `_summary.csv`,
`_paired.csv`, the identity-gate record, and an analysis restricted to block 8R.
