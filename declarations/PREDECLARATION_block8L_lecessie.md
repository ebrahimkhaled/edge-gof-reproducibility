# Pre-declaration: block 8L, the le Cessie-van Houwelingen test (written before any 8L cell is run)

Date: 2026-09-19. Author's decision of the same day.

## 0. Why it is added, and when

The le Cessie-van Houwelingen test smooths the individual residuals in covariate space
(le Cessie and van Houwelingen, 1991, 1995). It is the archetype of the ungrouped route that paper 2
contrasts with grouping -- its Introduction names it -- and on the Liu et al. (2024) grid it was among
the strongest tests in this family's own benchmark work. It was not in the battery. The reason given
for leaving it out of the earlier EDGE manuscript (that such tests are resampling-calibrated) does not
hold for it: the implementation used here refers its statistic to a closed-form moment reference.

It is therefore added now, after blocks 0-9d have reported, in the same way block 8 added the
projection test and BAGofT: this document fixes the design, the rules and the claims before any 8L
cell is computed, and 8L sees exactly the data sets the other tests saw.

## 1. The implementation

`gof_lecessie()` of **ebrahim.gof >= 2.6.0**, the version with the corrected moment reference
`M = (I - H)' R (I - H)` (the transpose fix of 2026-07-29, measured then as the difference between a
null size of 0.0705 and 0.0555 at n = 200). Uniform kernel `R = max(1 - d / h, 0)` on the scaled
covariate distances, bandwidth `h` the mean pairwise distance (the package default), and the
package's own reference distribution for the p-value. It uses no random numbers.

**Not used:** the vendor copy `_thesis_paper/logistic-gof-benchmark/R/vendor/lecessie1995.r`, which
builds `(I - H) R (I - H)` and is anti-conservative. It is deposited unmodified for the benchmark's
reproducibility and must not be substituted here.

The test builds an n x n kernel matrix, so it is O(n^2) in memory. Every cell below has n <= 1000.

## 2. Design

**Part A -- the block 8 cells.** The 30 cells of block 8 (n <= 1000; ten matched nulls and twenty
alternatives), replicates 1-500 of each, the same replicates the projection test used, so the two
and the battery's tests can be paired inside a cell.

**Part B -- corrupted covariates.** The block 9 cells at n = 1000 under the logistic truth: the clean
cell and corruption C1 (x multiplied by 4, outcome drawn at the original x) at k = 1, 2, 5 and 10,
replicates 1-1000 of each. The fitted model is correct for every uncorrupted record, so every
rejection is a false alarm.

**Data.** Nothing is simulated afresh. Every data set is regenerated from the seed the original block
stored for it, and block 8's identity gate applies: a data set is used only if the directed test's
p-value recomputed on it matches the stored one to 1e-8. A failure stops the cell and is reported.

## 3. Rules

1. Rejection at p <= 0.05. The reference is closed form, so there is no Monte Carlo grid and no atom.
2. A data set on which the test gives no p-value counts as no rejection, and the rate is reported.
3. Size is checked against a band of three nominal standard errors around 0.05, the standard error
   being sqrt(0.05 x 0.95 / B): 0.021 to 0.079 at B = 500 (Part A), and 0.029 to 0.071 at B = 1000
   (Part B).
4. Power is size-adjusted with the critical value read from the matched null, as in block 8.
5. Inside a cell, le Cessie is paired with EDGE-poly3 (unit form, rule G) on the same replicates by an
   exact McNemar test at 0.05, Holm-corrected within block 8L.

## 4. Claims, and what would count against them

- **L1 (level).** In each of the ten Part A null cells, le Cessie's rejection rate lies within three
  nominal standard errors of 0.05. *Against:* any null cell outside the band. If L1 fails in a cell,
  that cell's power comparisons are read with that failure stated beside them.

- **L2 (power, reported).** Size-adjusted power against the twenty Part A alternatives, and the paired
  comparison with EDGE-poly3. No direction is predicted: the test is added because it is strong.

- **L3 (corrupted covariates).** In Part B, le Cessie's false-alarm rate exceeds 0.10 at k = 5 and
  at k = 10. *Against:* a rate at or below 0.10 at either. The prediction comes from the mechanism of
  block 9: a covariate multiplied by four is isolated in covariate space, so the kernel finds few
  neighbours to smooth it with, and its standardised residual is large because its fitted risk is
  extreme. If L3 fails, le Cessie is reported as robust to this corruption alongside the grouped
  tests, in the same place and type.

- **L4 (reported).** The clean Part B cell's rate, the rate at k = 1 and 2, and the median seconds per
  data set.

## 5. Outputs

`simulations/battery/8L/`, one per-replicate file per cell, the identity-gate record, and an analysis
restricted to block 8L.
