# Pre-declaration: block 8S, the three rivals at n = 200 on the link and tail scenarios (written before any 8S scenario is run)

Date: 2026-09-20, at the author's request of the same day. Block 8 put the Liu projection test and
BAGofT on the link and tail scenarios at n = 500 and n = 1000, and block 8R put all three rivals on the
omitted-term and rough-misfit scenarios at n = 200, 500 and 1000. The one gap is here: **no rival has
been run on the link and tail family at n = 200**, which is the family the directed test is built for and
the sample size at which a clinical model is most often fitted. It is also the cheapest arm left, because
both resampling tests scale steeply with n.

## 1. Scenarios

The six link and tail alternatives at n = 200 and their matched null, all existing battery cells,
regenerated from the seeds the battery stored:

| cell | block | role |
|---|---|---|
| `cauchit_n200` | 3 | alternative |
| `t4_n200` | 3 | alternative |
| `loglog_n200` | 3 | alternative |
| `stk_short_n200` | 3 | alternative |
| `stk_long_n200` | 3 | alternative |
| `stk_asym_n200` | 3 | alternative |
| `null_link_n200` | 1a | null, matched to all six |

## 2. Tests, replicates and rules

The three rivals of blocks 8 and 8R, at their published defaults, with the replicate counts those blocks
used: **le Cessie's smoothed test** (500 replicates), the **Liu projection test** with B = 250 model-based
bootstrap draws (500), and **BAGofT** (200), the last with the same one-line repair declared in block 8R
(the covariate subset given `drop = FALSE`, nothing else changed). The directed test is recomputed on
every data set and must reproduce the battery's stored p-value to 1e-8; a replicate that fails that gate
is dropped and counted, as in blocks 8L, 8R and 9R.

The two resampling tests return Monte Carlo p-values on a grid and therefore reject at `p < 0.05`; le
Cessie's test and the directed test reject at `p <= 0.05`. Power is size-adjusted using each test's own
critical value read from `null_link_n200`, the matched null, exactly as the main study does. Comparisons
are paired inside the scenario by exact McNemar on the discordant pairs, Holm-corrected within block 8S.

At n = 200 the rule `G = max(10, ceil(n/25))` is ten groups, so the directed test has **no partition
refinement available here**: this arm tests it without the lever the paper sells.

## 3. How these scenarios enter the paper (fixed here, before the numbers exist)

The census counts of Section 5 -- the scenarios in which the directed test leads or ties each rival --
are **recomputed over the enlarged set**, and the previous counts are stated beside the new ones so that
the effect of adding this arm is visible rather than absorbed. Six comparisons are added against each
rival. No scenario already run is removed, re-weighted or re-read.

## 4. What would count against the paper (fixed here)

The paper claims the directed test is the most powerful of the pooled tests. This arm goes against that
claim if, **against any single rival, the directed test loses in two or more of the six scenarios with a
Holm-corrected paired p-value below 0.05**. Should that happen, the claim is qualified by sample size in
the abstract, in Section 5 and in the discussion -- not softened elsewhere and left standing there.

A rival that fails its size check at `null_link_n200` has its power reported size-adjusted, as everywhere
else, and the failure is reported.

## 5. Written prediction (recorded before the run)

At n = 500 the directed test led or tied the projection test in five of these six scenarios and BAGofT in
all six, the two clear wins being `stk_short` (0.46 against 0.25, and 0.42 against 0.10) and `stk_asym`
against BAGofT (0.92 against 0.54); `t4` was the one cell where a rival was nominally ahead, by 0.010 and
0.015, neither significant. At n = 200 every test will be near its level in the weak scenarios, so most
of the six comparisons should be ties decided by a handful of discordant pairs. I expect the directed
test to lead or tie **four to six of six against each rival**, with `t4` and `stk_long` the likeliest
losses, and I do not expect two Holm-significant losses against any one rival. The honest risk in this
arm is the opposite of a win: with ten groups and 200 records, several cells may be too flat to separate
any pair of tests, in which case the arm adds scenarios without adding information, and that is what the
report will say.

## 6. Outputs

`simulations/battery/8S/`: one per-replicate file per scenario and test, `_summary.csv` and
`_paired.csv` written by the block's own summary, and an analysis read against this document.
