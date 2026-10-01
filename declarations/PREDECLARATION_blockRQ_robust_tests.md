# Pre-declaration: block RQ, robust versions of the record-level test under corrupted covariates

Date: 2026-09-30. Written and hashed before any replicate of this block is computed.

## Why

Both referee rounds on EDGE paper 3 asked for the fair robust comparison. The paper so far fitted the model
robustly and put that fit inside non-robust tests (block 9d); a referee's point is that the competitor for
"a directed calibration test under a few corrupted records" is a robust *test* of the same alternative. The
robust test of Stukel's alternative is the quasi-deviance test of Cantoni and Ronchetti (2001) for the two
constructed covariates, computed with `robustbase`.

## Tests (every test on the same data sets)

1. **RQD-Huber**: `glmrob(method = "Mqle", weights.on.x = "none", tcc = 1.345)` null fit `y ~ x + d`; Stukel's
   constructed variables z1 = 0.5 eta^2 1{eta >= 0}, z2 = -0.5 eta^2 1{eta < 0} from the null fit's linear
   predictor; alternative fit `y ~ x + d + z1 + z2`, same settings; `anova(null, alt, test = "QD")`, p-value
   from its column `Pr(>chisq)`. Bounds the Pearson residual only.
2. **RQD-Mallows**: the same with `weights.on.x = "hat"`, which also down-weights high-leverage covariates.
   This is the variant a robustness referee would call the fair competitor for corrupted covariates.
3. **EDGE-default**: `ebrahim.gof::edge.gof(fit, G = "auto", basis = "poly3")` on the maximum-likelihood fit.
4. **EDGE-G10**: the same with `G = 10`.

## Data

The generator of block C1b (`_blockC1b_contam.R::c1b_gen`), which reproduces block 9's exactly for C1.
x ~ U(-3, 3), d ~ Bernoulli(0.5), eta = 0.6 x + 0.5 d.

| cell | truth | n | corruption | k |
|---|---|---|---|---|
| RQ_logit_clean_n1000 | logistic | 1000 | none | 0 |
| RQ_logit_C1_k{1,2,5,10}_n1000 | logistic | 1000 | x -> 4x | 1, 2, 5, 10 |
| RQ_logit_C1b_k{1,2,5,10}_n1000 | logistic | 1000 | x -> -4x | 1, 2, 5, 10 |
| RQ_cloglog_clean_n1000 | complementary log-log | 1000 | none | 0 |

1000 replicates a cell; seed `970000000 + 100000 * cell_index + rep`, disjoint from every earlier block.
A replicate on which a test returns no p-value counts as a non-rejection.

## Analysis, fixed now

- False-alarm rate at the nominal 5% (p <= 0.05) for every test in every logistic cell.
- Size-adjusted power in the cloglog cell, each test at the 5% quantile of its own p-values in the clean
  logistic cell.
- **Reading rule.** RQD-Mallows is said to *match the directed test's protection* if its false-alarm rate
  is at or below 0.10 at k = 10 under both corruptions. Whatever the outcome, the result replaces the
  paper's statement that robust tests were not studied, in Section 6.5, and is reported as measured.
- If either robust test fails to hold its level in the clean logistic cell (outside 0.05 +/- 3 standard
  errors, i.e. above 0.071 or below 0.029), that is reported, and its contaminated rates are still printed
  but read as mixing the level failure with the contamination.

## Outputs

`simulations/battery/RQ/<cell>_pvalues.csv.gz` (per replicate), `battery/RQ/_summary.csv`,
`battery/RQ/_progress.log`; runner `simulations/run_M_blockRQ.R`.
