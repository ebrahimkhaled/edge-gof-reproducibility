# Pre-declaration: block 9R2, the hybrid partition -- coarse tails, fine middle (written before any 9R2 scenario is run)

Date: 2026-09-20. Author's proposal of the same day, after block 9R: keep the two extreme groups as large as
they would be under a ten-group partition, and group the middle of the risk scale at the rule of about
twenty-five records a group. It is the one combination block 9R did not run: `POOL1` pooled the tails but
made them one per cent of the sample, which is far too small to dilute anything, and `G10` enlarged every
group, including the middle ones that carry most of the shape.

## 0. The mechanism, and what the hybrid is supposed to do

With `m_c` corrupted records among `m` in a group the grouped residual is about
`m_c / sqrt((m - m_c) * v_bar)`. Only the size of the group a corrupted record lands in decides whether it
is diluted, so protection is bought at the extremes alone; resolution, in contrast, is what the poly3 basis
uses to describe a bend, and a bend can sit anywhere on the risk scale. The hybrid buys the first without
paying the second everywhere -- if the departures the paper tests bend in the middle. If they bend in the
tails, the hybrid will cost as much power as the coarse partition and there is nothing to gain.

## 1. The variants (all EDGE-poly3, unit form)

- **V0, the default.** Equal frequency at the rule `G = max(10, ceil(n/25))`. The published test, and the
  identity gate: its p-value must reproduce the stored `EDGE.poly3.u.Grule` to 1e-8 on every replicate.
- **HYB05, HYB10, HYB20.** The `floor(f*n)` records of lowest fitted risk form one group and the
  `floor(f*n)` of highest fitted risk form another; the remaining `(1-2f)n` records are split equal
  frequency into groups of about twenty-five. Run at `f = 0.05, 0.10, 0.20`. At `f = 0.10` each tail group
  holds exactly what the top group of a ten-group partition holds, which is the comparison the author asked
  for.
- **G10, the coarse partition.** Equal frequency at `G = 10`, re-computed here as an anchor; it must
  reproduce block 9R's `G10` column exactly, since both are computed on the same regenerated data.

In every variant the estimation adjustment keeps the full-sample information `X'WX` of the fit that produced
the predictions, and the statistic is referred to the same closed-form weighted chi-squared null; only the
partition changes.

## 2. Scenarios

The twenty-two scenarios of block 9R, unchanged and regenerated from the same stored seeds: the clean
logistic scenarios at n = 1000 and 5000; the C1 corrupted-covariate scenarios under the logistic truth at
n = 1000 with k = 1, 2, 5, 10 and at n = 5000 with k = 5, 10, 25, 50; the clean probit and complementary
log-log scenarios at both sample sizes; the six battery scenarios `cauchit_n1000`, `t4_n1000`,
`loglog_n1000`, `stk_short_n1000`, `stk_long_n1000`, `stk_asym_n1000`; and the probit and complementary
log-log C1 scenarios at k = 5, n = 1000. One thousand replicates each.

## 3. Rules

Rejection at p <= 0.05. A replicate with no p-value counts as no rejection and is reported. Size is read
against 0.05 +/- 3 nominal standard errors (0.029 to 0.071 at B = 1000). Power is reported raw, beside the
size of the same variant in the matched clean scenario.

## 4. What would make the hybrid worth keeping (fixed here)

The same three criteria as block 9R, so that the answer is comparable:

1. **It holds its level:** the clean logistic scenarios at n = 1000 and 5000 lie inside the band.
2. **It protects where the default fails:** false alarms at or below 0.10 at n = 5000 with k = 25 and with
   k = 50, where the default reads 0.151 and 0.740 and the coarse partition reads 0.065 and 0.107.
3. **It keeps more power than the coarse partition:** mean power over the six battery scenarios strictly
   above the coarse partition's 0.490, and within 0.05 of the default's 0.508.

A variant that fails any of these is reported with the loss stated, not dropped. Nothing here changes the
default test.

## 5. Written prediction (recorded before the run)

Protection should track the tail-group size and nothing else, so HYB10 should read about what the coarse
partition reads at n = 5000 -- near 0.065 at k = 25 and near 0.107 at k = 50 -- HYB20 a little better and
HYB05 clearly worse, since fifty corrupted records in a 250-record group are a fifth of it. Power is the
open question. The departures tested here are link and tail misspecifications, which bend the calibration
curve most where the fitted risk is extreme, so the two coarse groups sit exactly where the signal is
strongest and we expect the hybrid to recover only part of the coarse partition's loss: somewhere between a
third and two thirds of the 0.018 mean gap, not all of it. If the hybrid recovers essentially all of it, the
signal is in the middle of the risk scale and the coarse tails cost nothing; if it recovers none, the
partition cannot be refined in the middle alone and the question is closed.

## 6. Outputs

`simulations/battery/9R2/`: one per-replicate file per scenario with a column per variant, and an analysis
restricted to block 9R2, read beside block 9R's table.
