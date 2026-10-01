# Pre-declaration: block 9R, four remedies for the extreme groups (written before any 9R scenario is run)

Date: 2026-09-20. Author's request of the same day: a variant in which the analyst states how much
corruption they are willing to tolerate at each extreme, to be used when the data are suspected of
carrying more than about ten corrupted records, or one per cent, since below that the default test is
already at its level.

## 0. The mechanism this attacks

With `m_c` corrupted records among `m` in a group, the grouped residual is about
`m_c / sqrt((m - m_c) * v_bar)`: the corrupted records add up in the numerator while the denominator
grows only with the square root of the clean variance they are pooled against. At the rule G every
group holds about twenty-five records, so a corrupted record that lands in the top group is diluted by
twenty-four clean ones at most. Every remedy below either enlarges that denominator or removes the
corrupted records from the statistic; none can do both, and each pays in power exactly where the
extremes carry signal.

## 1. The variants (all EDGE-poly3, unit form)

- **V0, the default.** Equal-frequency grouping at the rule `G = max(10, ceil(n/25))`. The published test.
- **V1, EDGE-FR(alpha).** The analyst's tolerance: the `floor(alpha*n)` records with the lowest fitted
  risk and the `floor(alpha*n)` with the highest are dropped, and the remainder is grouped at the rule
  G computed from the retained count. Run at `alpha = 0.01` (the default the author proposes) and at
  `alpha = 0.02`.
- **V2, outer groups dropped(alpha).** All records are grouped at the rule G; the `ceil(alpha*G)`
  lowest and highest groups are dropped from the residual vector and the basis. `alpha = 0.01`.
- **V3, tail-pooled partition(alpha).** The lowest `floor(alpha*n)` records form one group and the
  highest `floor(alpha*n)` another; the middle is split equal-frequency into groups of about
  twenty-five. `alpha = 0.01`.
- **V4, coarse partition.** Equal-frequency grouping at `G = 10`, the lever the paper already measures.

In every variant the estimation adjustment keeps the full-sample information `X'WX` of the fit that
produced the predictions, and the statistic is referred to the same closed-form weighted chi-squared
null; only the partition, or the set of records entering it, changes.

## 2. Scenarios (all regenerated from the seeds their blocks stored)

- **Level:** block 9's clean logistic scenarios at n = 1000 and 5000, 1000 replicates.
- **False alarms:** block 9's C1 scenarios under the logistic truth -- n = 1000 at k = 1, 2, 5, 10 and
  n = 5000 at k = 5, 10, 25, 50 -- 1000 replicates, the rows of the paper's contamination table.
- **Power, no corruption:** block 9's probit and complementary log-log scenarios at k = 0, n = 1000 and
  5000, 1000 replicates.
- **Power on the paper's own ground:** the battery scenarios `cauchit_n1000`, `t4_n1000`,
  `loglog_n1000`, `stk_short_n1000`, `stk_long_n1000`, `stk_asym_n1000`, replicates 1-1000.
- **Power under corruption:** block 9's probit and complementary log-log C1 scenarios at k = 5,
  n = 1000, 1000 replicates.

## 3. Rules

Rejection at p <= 0.05. A replicate with no p-value counts as no rejection and is reported. Size is
read against 0.05 +/- 3 nominal standard errors (0.029 to 0.071 at B = 1000). Power is reported raw,
not size-adjusted, and beside it the size of the same variant in the matched clean scenario, so that a
variant that buys power by running hot cannot be read as more powerful.

## 4. What would make a remedy worth keeping (fixed here)

A variant is worth reporting as an option when all three hold:

1. **It holds its level:** the clean scenarios at n = 1000 and 5000 lie inside the band.
2. **It protects where the default fails:** false alarms at or below 0.10 at n = 5000 with k = 25 and
   with k = 50, where the default reads 0.151 and 0.740.
3. **It keeps the power the paper sells:** mean power over the six battery scenarios of section 2 at
   most 0.05 below the default's.

A variant that fails (3) is reported with the loss stated, not dropped. Nothing here changes the
default test: EDGE-FR is an option for data already suspected of corruption, and the paper's scope
claim for the default is unchanged.

## 5. Outputs

`simulations/battery/9R/`: one per-replicate file per scenario with a column per variant, `_summary.csv`,
and an analysis restricted to block 9R.
