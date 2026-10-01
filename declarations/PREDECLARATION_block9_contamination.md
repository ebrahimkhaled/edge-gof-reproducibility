# Pre-declaration: block 9, a few corrupted records (written before any block 9 cell is run)

Date: 2026-09-16. Author's idea of the same day: Stukel's test uses every record at full weight, EDGE averages inside
risk groups, so a handful of corrupted records should disturb the ungrouped tests much more. Two exploratory runs that
afternoon (`scratchpad/contam_check.R`, `contam_check2.R`, n = 1000, 300 replicates) supported it and are named here as
**exploratory only**: their numbers motivate the design and are never quoted as evidence. Block 9 repeats the question
with fresh seeds, more replicates, a wider grid and the rules below, all fixed before it runs.

Block 9 belongs to Paper A (individual patient data). It does not touch blocks 0-8, whose pre-declaration (A-D,
E0-E13) is unchanged.

---

## 1. Design

Base design of the battery: x ~ U(-3, 3), d ~ Bernoulli(0.5), eta = 0.6x + 0.5d, fitted as y ~ x + d.

**Truths.** (a) logistic, so the model is correct for every uncorrupted record; (b) probit at n = 5000 and (c) cloglog
at n = 1000, both real misfits with moderate clean power, to see whether contamination costs a test its power.

**Corruptions.** k records out of n, with k/n in {0.001, 0.002, 0.005, 0.01}:
- **C1 corrupted covariate**: k records chosen at random have x multiplied by 4; their outcome is still drawn from the
  truth at the original x. The record's prediction is wrong; its outcome is not.
- **C2 missed events**: the k records with the highest true risk have y set to 0. The prediction is right; the recorded
  outcome is wrong.
- **C3 unit error**: k records chosen at random have x divided by 10 (a units mix-up), outcome from the original x.
- **C4 duplicated high-risk records** (at k/n = 0.005 only): the k highest-risk records are copied over k randomly
  chosen rows.
- **clean**: k = 0, one cell per truth and n.

**Sample sizes.** n = 1000 and 5000. **Replicates.** B = 1000 per cell. Cells: 3 truths x 3 corruptions x 4 rates x 2 n,
plus C4 (3 x 2) and clean (3 x 2) = 84 cells, about 13-23 CPU-hours.

**Tests.** Stukel joint score, Stukel one-parameter score, Stukel likelihood-ratio refit; GiViTI (default threshold);
the cubic calibration LR; EDGE-poly3 and EDGE-sym, each unit and score weighted, at G = 10 and at the rule G; HL and
HL_F at both G. Each replicate also stores the fitted coefficients and the number of corrupted records that land in
the highest-risk group.

**Seeds.** Per replicate, `set.seed(seed_base + rep)`, `seed_base = 300000000 + cell_id * 10000`; disjoint from every
earlier run.

## 2. Rules

1. Under truth (a) the model is correct for the uncorrupted records, so a rejection is a **false alarm**, not a size in
   the strict sense. It is reported as such, with the nominal standard error sqrt(0.05 x 0.95 / B) as in E13.1.
2. Under truths (b) and (c), power is compared with the clean cell of the same truth and n (both raw and the change).
3. A test that gives no p-value counts as no rejection, and its rate is reported.
4. Comparisons between EDGE (unit) and each ungrouped test are paired inside the cell, with an exact McNemar test at
   0.05, Holm-corrected within block 9.
5. Every claim below is judged on the C1 and C3 cells (corrupted predictions); C2 and C4 are reported without a claim.

## 3. Claims, and what would count against them

- **B9.1 (false alarms).** With corrupted covariates at 0.1% to 0.5% and a correct bulk model,
  EDGE-poly3 (unit) and EDGE-sym (unit) stay within 3 nominal standard errors of 0.05 in every such cell, while at
  least two of Stukel's joint score, Stukel's one-parameter score, GiViTI and the cubic LR exceed 0.10 in at least half
  of them. *Against:* either half fails.
- **B9.2 (missed events).** Reported without a claim: the expectation is that every directed test rises, grouped or
  not, because a wrong outcome at the top of the risk scale is exactly the pattern these tests look for.
- **B9.3 (power retention), revised on 2026-09-16 before any block 9 cell ran.** At 0.1% corrupted covariates, EDGE's
  unit forms lose at most 0.05 of their clean power against the probit and cloglog truths. *Against:* a larger loss.
  At 0.5% and 1% the loss is measured and reported without a claim: the exploratory follow-up of the same day
  suggests it is substantial, because the corrupted predictions and the genuine tail misfit live in the same extreme
  groups, so the same averaging that ignores the first also dulls the second. The paper reports this trade rather
  than hiding it.
- **B9.4 (the corruption is not influential for the fit).** In the C1 and C3 cells the mean fitted coefficients differ
  from the clean cells by less than 10%. Reported; if it fails, the cell is described as influential and B9.1 is judged
  without it.
- **B9.5 (which weighting).** Reported: the same quantities for the score form, which the exploratory run suggests
  behaves differently from the unit form.

- **B9.6 (the mechanism is leverage, not any covariate error), added 2026-09-16 before the run.** The unit-error
  corruption (C3), which makes a record's prediction milder rather than extreme, moves no test at any rate. *Against:*
  any test whose rejection rate changes by more than 3 nominal standard errors in the C3 cells.

## 4. Reading of the result, fixed in advance

If B9.1 and B9.3 hold, the paper states: among directed calibration tests, only the grouped unit-weighted ones ignore a
handful of corrupted predictions, and they keep their power against real misfit. The reason is that the influence of
one record on an ungrouped score grows with the square of its linear predictor, while inside a group it is bounded by
that group's standard deviation. The paper will also say, in the same place, that not reacting to corrupted records is
the right behaviour for deciding whether to update a model, and that checking data quality is a separate task that
needs its own tools.

If B9.1 fails, the robustness claim is dropped and the exploratory finding is reported as not confirmed.

**Revision note (2026-09-16, before any block 9 cell was run).** After the design above was written and hashed, the
second exploratory run (n = 1000 and 5000, rates 0.1% to 1%, 300 replicates) reported. It showed three things that the
claims now reflect: a single corrupted record in 1000 already lifts the ungrouped tests (Stukel's joint score from
0.070 to 0.163) while EDGE's unit forms do not move; at 1% every test reacts, EDGE included; and under a real probit
misfit with 0.5% corrupted records EDGE's power falls well below its clean value. The rate grid gained 0.1%, B9.1 was
narrowed to 0.1-0.5%, B9.3 became a claim at 0.1% and a measurement above it, and B9.6 was added. The exploratory
numbers remain exploratory and are not quoted as evidence.

## 5. Outputs

`simulations/battery/9/` with one per-replicate file per cell and a summary in the battery's format; the analysis is an
extension of `analyse_M_battery.R` restricted to block 9.
