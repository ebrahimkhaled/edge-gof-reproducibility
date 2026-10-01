# Pre-declaration: block EXT, the directed test against the usual external-validation tests

Date: 2026-09-30. Written and hashed before any replicate of this block is computed.

## Why
In external validation a published model with frozen coefficients is checked on new patients. The usual
tests there are Cox's recalibration test (intercept and slope), Spiegelhalter's z and the GiViTI belt in
external mode. EDGE paper 3 claims that in this setting its bound on a corrupted record is exact
(Proposition 1). This block measures (1) level, (2) power against five kinds of miscalibration, and (3)
false alarms under corrupted covariates, for the directed test in external mode and those tests.

## Design (data generator identical to paper_deepgof/theory/external/run_external.R)
Five N(0, 1) covariates; published model eta0 = -0.5 + 0.5 x1 + 0.8 x2 + 0.6 x3 + 0.4 x4 + 0.3 x5, p = expit(eta0).
Truth eta0 + delta:

| id | family | delta | n |
|---|---|---|---|
| 1-3 | null | 0 | 250, 500, 1000 |
| 4-5 | large | C in {.2, .4} | 500 |
| 6-7 | slope | (C - 1) eta0, C in {.8, .6} | 500 |
| 8-9 | ushape | C (x1^2 - 1), C in {.4, .8} | 500 |
| 10-11 | thresh | C max(x2 - 1, 0), C in {2, 4} | 500 |
| 12-13 | inter | C x1 x3, C in {.6, 1} | 500 |

Seeds for cells 1-13: set.seed(7000000 + id * 10000 + rep), exactly as run_external.R, so replicate r is
the same data set in both studies.

Corruption under a correct frozen model (cells 14-21), n = 1000: the outcome is drawn from p at the true
covariates; k random records have x2 recorded as 4 x2 (exaggeration, ids 14-17) or -4 x2 (sign error,
ids 18-21), k in {1, 2, 5, 10}, and their recorded prediction is expit(eta0) at the recorded x2. Seeds
set.seed(7600000 + id * 10000 + rep).

1000 replicates in every cell.

## Tests (all on the recorded predictions p and outcomes y; nothing refitted)
- EDGE-ext default: edge_external(y, p, G = "auto") (cubic basis plus constant, G = max(10, ceil(n/25))).
- EDGE-ext G10: edge_external(y, p, G = 10).
- Cox recalibration LR test, chi-square(2); Spiegelhalter's z, two-sided; GiViTI belt, devel = "external";
  Hosmer-Lemeshow external form, chi-square(10); Stukel's two terms on the offset eta0, LR chi-square(2).
  Code identical to run_external.R.
A replicate with no p-value counts as a non-rejection.

## Analysis, fixed now
- Level: rejection at 5% in cells 1-3; a test holds if within 0.05 +/- 3 SE (0.029 to 0.071).
- Power: size-adjusted at the 5% quantile of each test's p-values in cell 2 (null, n = 500), per cell.
- Protection: false-alarm rate in cells 14-21. **Reading rule:** the directed test (ten groups) is said to
  hold where its rate is at most 0.10; the same bar is applied to every test and reported for all.
- Expectation written in advance (not a pass/fail): Cox, Spiegelhalter and GiViTI lead on calibration-in-
  the-large and slope; the directed test is competitive on the u-shape and threshold; nothing along the
  predicted risk detects the interaction well.

## Outputs
simulations/battery/EXT/<cell>_pvalues.csv.gz, _summary.csv, _progress.log; runner simulations/run_M_blockEXT.R.
