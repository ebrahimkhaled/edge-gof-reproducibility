# Pre-declaration: the EDGE restructure battery (written BEFORE any of its cells is run)

Date: 2026-09-13 (night). Author decisions recorded in `ACADEMIC_TRACKER.md` (same date). This file fixes, before
the full re-run, (A) how the weighting of each EDGE basis will be chosen, (B) the tests and their exact forms,
(C) the analysis rules, and (D) the hypotheses with what would count against them. Section E (the cell list) is
completed from `MAP_battery.md` before launch; nothing in A–D may change after the first result is read. Any cell
or analysis added after results are seen is labelled *post hoc* wherever it is reported.

---

## A. Choosing unit or sqrt(V_g) weighting (the author asked for an empirical choice)

**Candidates, per basis.** For each of EDGE-poly3 (default), EDGE-sym (new, column eta|eta| at the group-mean
logit) and EDGE-stk:
- *unit form*: the current construction, S = r'P_Z r, weighted chi-square null (Satterthwaite/Imhof);
- *score form*: columns multiplied by sqrt(V_g), statistic u'I^-1 u with I the post-fit information of the
  weighted columns, plain chi-square on the column rank.

**Eligibility (validity gate).** A variant is eligible only if, in every null cell of the battery, its realised size
at alpha = 0.05 lies within three Monte Carlo standard errors of 0.05 or is conservative (below), and at
alpha = 0.01 it does not exceed 0.01 by more than three standard errors. A variant that fails the gate in any cell is
reported and not chosen.

**Primary criterion.** Macro-averaged size-adjusted power at alpha = 0.05: average over the alternative cells within
each departure family first, then over families with equal weight. Families: (1) symmetric-tail links (probit,
cauchit, t4); (2) asymmetric links (cloglog, loglog, Stukel asymmetric); (3) Stukel symmetric-family links (heavy,
light); (4) omitted terms (quadratic, binary interaction, continuous interaction); (5) rough misfit (oscillation,
sawtooth); (6) off-index (crossover). Each family is averaged over its designs (base, high AUC, low event rate,
where run) and sample sizes.

**Decision.** Chosen separately for each basis:
1. If one eligible variant has a macro-average higher by at least 0.014 (the paper's tie margin), it is chosen, unless
   it loses more than 0.10 against the other variant in any single family average, in which case the other variant is
   chosen.
2. Otherwise (difference below 0.014) the unit form is kept, for continuity with the published definition.

Both variants of every basis are reported in Supporting Information whatever is chosen.

## B. Tests and their exact forms

- **EDGE**: poly3 (headline default), sym (new named basis, headline on symmetric tails), stk and poly2 (sensitivity);
  each in the unit and the score form; partition G = max(10, round(n/25)) (the paper's rule) as the headline and
  G = 10 in Supporting Information.
- **GiViTI**: `givitiR::givitiCalibrationTest(..., devel = "internal")`, default threshold 0.95 (headline rival), and
  threshold 0.50 (sensitivity row).
- **Stukel**:
  - *joint score* u'I^-1 u on 2 df with the one-df fallback when a half-column is identically zero (the corrected
    package statistic; headline "Stukel's score test"); the fallback frequency is reported;
  - *likelihood-ratio refit* (augmented glm; non-convergence or an error is recorded as "no p-value");
  - *one-parameter symmetric score* (eta|eta|, ungrouped; named as the ungrouped counterpart of EDGE-sym);
  - *marginal sum* (the LogisticDx / ebrahim.gof 2.6.0 statistic): size section and one power column only.
- **Partition and covariate-space family**: HL (G = 10 and the rule), HL_F, Pigeon–Heyse, Tsiatis, Xie,
  Pulkstenis–Robinson, as in the current manuscript.
- **Projection test**: boundary (off-index and interaction) cells only, as now.

## C. Analysis rules

1. Size first: every test's realised size at 0.01, 0.05, 0.10 in every null cell, with Monte Carlo standard errors.
2. Power is size-adjusted (5% quantile of the matched null at the same design, n and G); raw power is also reported.
3. Tie margin 0.014; differences are paired (same replicates) and tested with McNemar's exact test at the nominal
   level for tests that hold size, Holm-corrected within each claim family below.
4. "Cannot return a p-value" is counted for EVERY test, including samples with fewer events than groups for the
   grouped tests, and reported both ways (no sample is silently excluded).
5. Per-replicate p-values are stored for every cell.
6. Seeds: one stream per cell, offset by design, n and G; no replicate is shared between cells.

## D. Hypotheses and what would count against them

- **H1 (GiViTI).** On symmetric-tail links (family 1, all designs run), EDGE-poly3 exceeds default GiViTI by at least
  0.10 in family-averaged size-adjusted power. *Against:* average gain below 0.10, or any design where GiViTI leads by
  more than 0.014 (to be reported as a limit, not hidden).
- **H2 (Stukel).** On family 1, EDGE-sym is no worse than Stukel's joint score by more than 0.014 on average.
  *Against:* average deficit beyond 0.014.
- **H3 (bows).** On family 2, EDGE-poly3 is within 0.05 of default GiViTI on average. *Against:* a deficit above 0.05.
- **H4 (partition family).** On families 2–4 at the base design, EDGE-poly3 leads or ties the four omnibus partition
  tests in at least as many detectable cells as in the current manuscript (19 of 22). *Against:* fewer.
- **H5 (grouping cost).** At the base design EDGE-sym is within 0.02 of Stukel's one-parameter ungrouped score test;
  under a skewed covariate the cost is larger. *Against:* a cost above 0.02 at the base design.
- **H6 (computability).** Reported both ways, no hypothesis: the share of samples in which each test returns no
  p-value, including events < groups.

## E. Cell list

Written 2026-09-13 (late night), after the author's decisions on `MAP_battery.md` and BEFORE any battery cell was run.
Sections A-D are unchanged: their sha256 (`aad03b19...`) still matches, and a byte copy is kept as
`PREDECLARATION_restructure_battery_FROZEN_AD.md`. E only fixes what A-D left open; where it reads A-D, it says so.

### E0. Conventions that apply to every block

1. **"Stukel" links in A and D mean Stukel's exact h-family** (Stukel 1988, Sec. 2): logit(mu) = h_alpha(eta), with
   h = (exp(a1 eta) - 1)/a1 for eta >= 0 and a1 > 0, -log(1 - a1 eta)/a1 for a1 < 0, and the mirror rule with a2 for
   eta < 0. Shape values are those of Hosmer, Hosmer, le Cessie and Lemeshow (1997), Table V: **long (heavy) tails
   (-1, -1); short (light) tails (1, 1); asymmetric long-short (-1, 1)**. The linear predictor is the design's eta.
2. **The capped generators of the earlier versions are kept** (author decision), under names that describe what they
   generate: `plateau_upper` = old `stukel_heavy` inv_stukel(eta, -1, -1); `plateau_lower` = old `stukel_light`
   (1, 1); `plateau_both` = old `stukel_asym` (-1, 1). They keep their old seeds and cell ids. **Placement is fixed
   now, whatever the results:** all three rows appear in the rebuilt Table 4 directly under the three exact Stukel
   rows, as a group "Risk plateau (capped generators used in earlier versions)", and the full n ladder of both groups
   is in Supporting Information. They are not part of the A macro-average and not part of H1-H6. No row of either
   group may be moved, dropped or renamed after results are read.
3. **Paired arms share replicates.** Within a cell, every test, both G arms (G = 10 and G = max(10, round(n/25))),
   both weighting forms and all four Stukel forms are computed on the same data set. C6's "no replicate is shared
   between cells" applies across cells only.
4. **Seeds.** Every replicate uses `set.seed(seed_base + rep)` (per-replicate seeding, as `ek_run_cell()`), so data do
   not depend on the worker count. Cells of the July grids (blocks 1a, 3, 4) reuse their old seed bases and cell ids
   (`20260705`, `20260708`, `20260710`, `20260711`, `20260720`, each `+ cell_id*1e5`), so their data are byte-identical
   to the stored runs. New cells use `seed_base = 100000000 + block*1e7 + cell_id*1e4` (no overlap for B <= 10,000,
   within integer range, disjoint from the old bases). Runs A-L used worker-dependent cluster streams and cannot be
   regenerated; block 5 re-runs them with new per-replicate seeds, and their tables will move within Monte Carlo error.
5. **"No p-value".** A test that errors, declines (for example a non-converged refit) or returns NaN counts as
   **no rejection** in every power and size figure, and its rate is reported in a `declined` column (C4). GiViTI
   warnings are suppressed, never turned into NA; only an error is "no p-value".
6. **Covariate-space and partition rivals** (HL_w, Pigeon-Heyse, Tsiatis, Xie, Pulkstenis-Robinson) run at G = 10 only,
   and are **not run at n >= 10,000** (reported as "not run"). HL and HL_F run at G = 10 and at the rule G. The
   projection test is not re-run: in the omit_2cov and boundary cells its stored per-replicate p-values are paired
   with the new tests through the reused seeds.
7. **Sample sizes of new cells** come from a test-agnostic rule, fixed now: from a population draw (N = 400,000) of the
   design, compute D = E[(p_true - pi*)^2 / (pi*(1 - pi*))], where pi* is the limit of the fitted working logistic model
   (fit with the true probabilities as weights), and n_env = 7.85/D (the 80% power size of an oracle one-degree-of-
   freedom test at alpha = 0.05). Then **n_lo = 0.6 n_env and n_hi = n_env, rounded to two significant figures and
   kept within [300, 20,000]**. A cell whose n_env exceeds 40,000 is run once at n = 20,000 as an undetectable control.
   `battery_nplan.R` computes these numbers; they are appended below as E7 before any battery cell runs. No test's
   power enters this rule.
8. **Workers.** Block 0 times 200 replicates of one n = 5000 cell with 12, 16 and 20 workers and uses the fastest.

### E1. Designs

x ~ U(-3, 3), d ~ Bernoulli(0.5), eta = c0 + s(0.6x + 0.5d), always fitted as y ~ x + d, unless a block says otherwise.
- **base**: s = 1, c0 = 0 (the paper's base design).
- **auc**: s = 2, c0 = 0 (AUC 0.86-0.96 depending on the link; `MAP_battery.md` Sec. 2).
- **e12**: s = 1, c0 giving exactly 12% events for that link: logit -2.64, probit -2.01, t4 -2.18, cauchit -3.15,
  loglog -1.60, cloglog -2.76. The matched null is the logistic model at c0 = -2.64.
- **skew**: as base, but x = 1.5 * standardised chi-square(4) draws (as `run_L3_closing.R`).
Links: probit pnorm; cauchit pcauchy; t4 pt(., 4); loglog exp(-exp(-eta)); cloglog 1 - exp(-exp(eta)).

### E2. Tests computed in every replicate

- EDGE-poly3, EDGE-poly2, EDGE-stk, EDGE-sym: unit form and score form, each at G = 10 and at the rule G (one shared
  grouping per G). EDGE-ao (unit) in block 5 run K2 cells only.
- Stukel: joint score (`Stk.joint`, with the one-df fallback flag), LR refit (`Stk.LR`), one-parameter symmetric score
  (`Stk.sym1`), marginal sum (`Stk.marg`).
- GiViTI internal at thres 0.95 (`GiViTI`) and 0.50 (`GiViTI.t50`); degree-3 calibration LR (`Cubic.LR`, sensitivity).
- HL (G = 10 and rule), HL_F (`ef.gof`, G = 10 and rule), and the E0.6 rivals.
- Flags per replicate: n, events, min(events, non-events) < G for each G, a Stukel half-column identically zero, the
  `stuk_diag` refit-failure flag (Figure 8 comes from the same replicates).
- Implementation: EDGE is computed in the harness and must agree with `ebrahim.gof` 2.8.0 `def.gof(..., weights =)`
  to 1e-8 in block 0; the joint score must agree with `gof_stukel(form = "joint")` to 1e-8.

### E3. Cells by block (B = replicates per cell)

| Block | Cells | B |
|---|---|---|
| 0 validation | (i) old-seed identity: `grid_null` quad n = 1000; `grid_power_broad` link cloglog n = 1000; crossover n = 1000 (EDGE poly2/poly3/stk unit at G = 10, HL, HL_F must equal the stored p-values to 1e-6, or the cause must be found and written down before launch); (ii) harness against package (E2); (iii) workers (E0.8) | 200 |
| 1a size grid | `dgp_null` x {quad, binint, contint, link, rough} x n {200, 500, 1000, 2000, 5000}, old seeds and cell ids; link and quad at n = 1000 also G {6, 8, 12, 14, 20} inside the same replicates | 10,000 |
| 1b new nulls | logistic null for every (design, n) used in block 2 (base, auc, e12); skew null at every n of the block 4 skew cells; `gen_sparse_link(-4.9)` at n {100, 150, 200, 300, 500, 1000, 2000, 5000} and `gen_sparse_link(-1)` at the same n (shared with Figure 8); run G sparse null at n {500, 2000, 10,000}; run A null at p* = 20, n = 1000; quad null at n = 100 with G {5, 10} | 2000; 5000 for the `gen_sparse_link` cells; 10,000 at n = 100 |
| 2 symmetric tails | {probit, cauchit, t4, loglog} x {base, auc, e12} x {n_lo, n_hi} from E7 (controls at 20,000 where E0.7 says so) | 2000 |
| 3 Table 4 battery | all 117 `grid_power_broad.R` cells with their old seeds and cell ids (this includes probit, cloglog and the three plateau generators); NEW: exact Stukel long, short and asymmetric x n {200, 500, 1000, 2000, 5000}; cauchit, t4 and loglog x the same n; the 8 June census scenarios of `_master_add.R` (logx, int_binbin, skew, corr, omit_x2x3, omit_2int, joint, omit_2cov; generators as `MAP_battery.md` Sec. 2.2) at n = 1000 | 5000 (10,000 at n = 200) |
| 4 boundary and grouping | `grid_edge_loses.R` 6 cells with old seeds; omit_2cov n {500, 1000} with the `grid_proj_power.R` seeds (paired with the stored projection p-values); **skew** x {probit, cauchit} at the same two n as their block 2 base cells (H5; the covariate changes, n does not); sparse cloglog `1 - exp(-exp(-4.9 + x))` on `gen_sparse_link` covariates at n {200, 300, 500} (events < groups: size and power reported unconditionally and separately for samples with events < G and >= G) | 5000 (grid cells); 500 (omit_2cov); 2000 (others) |
| 5 rest of the paper | run A (32 cells, `gen_pstar`); run F (75 cells, n = 20,000 de-duplicated to 9); run G (6 scenarios x {500, 2000, 10,000}, both G arms in one replicate); run K2 (14 cells, plus EDGE-ao); run H (36 cells, external battery plus an EDGE-sym external column, per-replicate p-values stored this time); Figure 8 from the block 1b `gen_sparse_link` cells | as in the original runs |
| 6 real data | `_realdata_power.R` GLOW 4 cells at n = 500 with the new tests (new per-replicate seeds); one-off fits of beetle, LBW, vaso, UIS (model 19) and Diabetes-130 with the new tests | 2000 / real data |
| 7 long tail | run D probit ladder to n = 50,000 (cells already in blocks 2-4 dropped); run J (EDGE only); `null_calibration_checks.R` unchanged; projection not re-run | as in the original runs |

### E4. Run order, outputs and stopping rules

- Order: 0; 2 with its 1b nulls; 1a and the rest of 1b; 3 (the n = 1000 slice first); 4; 5; 6; 7.
- Outputs: one `battery/<block>/<cell>_pvalues.csv.gz` per cell (one row per replicate: seed, n, events, flags, every
  p-value), a per-block summary, and `battery/<block>/_progress.log`. A cell whose output exists is skipped, so the
  driver can resume after an interruption.
- `battery/RULE_weighting.txt` holds Section A verbatim; the driver prints its sha256 and refuses to start block 2 if
  the file is missing.
- Stop before launch if block 0 (i) or (ii) fails. A test that fails size is reported and not used for comparative
  claims (C1); that is not a reason to stop.

### E7. Sample sizes from `battery_nplan.R`

Appended 2026-09-14, before any battery cell was run at its real B. Source: `simulations/battery/nplan.csv`
(population N = 400,000, seed 20260913). Only smoke tests had run (block 0 at B = 20, two cells per block at B = 10);
their outputs check that the code runs and matches the stored runs and the package, and no power or size figure from
them was read.

| Design | Link | Event rate | AUC | D | n_env | n used in block 2 |
|---|---|---|---|---|---|---|
| base | probit | .564 | .867 | 5.03e-4 | 15,594 | 9,400; 16,000 |
| base | cauchit | .547 | .767 | 1.97e-3 | 3,985 | 2,400; 4,000 |
| base | t4 | .559 | .837 | 1.18e-4 | 66,267 | 20,000 (control) |
| base | loglog | .444 | .846 | 7.77e-3 | 1,010 | 610; 1,000 |
| auc | probit | .570 | .961 | 1.33e-3 | 5,888 | 3,500; 5,900 |
| auc | cauchit | .558 | .857 | 1.03e-2 | 762 | 460; 760 |
| auc | t4 | .568 | .937 | 1.00e-3 | 7,837 | 4,700; 7,800 |
| auc | loglog | .493 | .943 | 1.24e-2 | 636 | 380; 640 |
| e12 | probit | .120 | .872 | 1.08e-3 | 7,245 | 4,300; 7,200 |
| e12 | cauchit | .120 | .622 | 7.31e-4 | 10,740 | 6,400; 11,000 |
| e12 | t4 | .121 | .811 | 5.40e-4 | 14,531 | 8,700; 15,000 |
| e12 | loglog | .120 | .892 | 4.60e-3 | 1,708 | 1,000; 1,700 |

The E1 intercepts give event rates of 0.1200-0.1208 for every link, within the 12% +/- 0.5% check (logit null
0.1203; e12 cloglog 0.1203, n_env 71,501, not a block 2 cell). The block 4 skew cells use the base n above (probit
9,400 and 16,000; cauchit 2,400 and 4,000).

**Stated exception to C6, inherited from the July grids.** Reusing the July seed bases (E0.4) keeps their data
byte-identical, but those grids already shared seed ranges across families: `grid_null`, `grid_power_broad` and
`grid_edge_loses` cells with the same id use bases 3, 5 or 2 apart (`battery/seed_overlaps.csv`, 41 pairs). Such pairs
use common random numbers shifted by a few replicates. This does not bias any size or power estimate; it only
correlates Monte Carlo error between those cells, as in the published runs. All new cells are disjoint from each
other and from the old ranges (checked by `battery_selftest.R`).

**Size of the run.** Block 0 plus 436 cells, 1,707,808 replicates; projected 42.4 CPU-hours, about 5.7 wall-hours at
the parallel efficiency observed on this machine (`battery/_estimate_blocks.csv`). The worker count is chosen by the
real block 0 (E0.8).

### E8. Pre-launch changes after the independent review

Written 2026-09-14, after two independent reviews of the battery code (`BATTERY_review_conformance.md`,
`BATTERY_review_numerics.md`) and of the package (`PKG280_verify_*.md`), and before any cell was run at its real B.
No power or size figure had been read. Sections A-D are unchanged; these items only close gaps the reviews found.

1. **Samples with no event or no non-event.** The working model does not exist, so every test gets "no p-value"
   (`flag.degenerate`), counted as no rejection and in `declined` (E0.5). Before this rule, four statistics returned
   p = 0 on such samples and would have counted as rejections (size inflated by up to 0.15 at n = 100).
2. **Information guard.** Stukel's joint and one-parameter score statistics and both EDGE score statistics drop a
   column whose post-fit information is below 1e-10 of its unadjusted information, and give no p-value when no column
   is left (`flag.info_guard`, rate reported). `ebrahim.gof` 2.8.0 applies the same guard, so harness and package agree.
3. **EDGE basis columns (harness = package).** A column with sum |z| <= 1e-8 is dropped, as in 2.7.0; the kept columns
   are scaled to unit length before any solve. The statistics do not depend on that scale, so every result that ran
   before is unchanged, and the singular stop on a tiny Stukel half-column is gone. A scale-free drop rule (relative
   length below 1e-6) was tried on 2026-09-14 and withdrawn, because it changed valid results.
4. **Cells added** (E0.4 new seeds unless stated):
   - crossover logistic nulls at n = 1000 and 2000 (B 5000; truth = the working model's population limit, eta = 0.62x),
     so Section A family 6 has a matched null;
   - one logistic null per census design at n = 1000 (8) and omit_2cov at n = 500 (B 5000; truth = the working model's
     population limit on the same covariates);
   - all 24 `grid_proj_power.R` cells in block 4 with their old seeds (B 500), paired with the stored projection
     p-values (B and E0.6 promised the interaction cells; only omit_2cov had been included);
   - `null_rough_n10000` (B 1000); the osc4 alternatives of runs G and K2 are now matched to the rough nulls of the same
     design (y ~ x) instead of the base-design null;
   - the Diabetes-130 G sweep (10, 50, 200, 1000, 2000) on both halves, and a beetle cloglog refit reporting the tests
     that do not assume a logit.
   Totals: block 0 plus 471 cells, 1,774,809 replicates, about 46 CPU-hours projected.
5. **Real data and tied fitted risks.** Equal-frequency grouping splits tied fitted risks by row position (in the
   package as in the harness), and the beetle rows list events first within each dose, so grouped p-values depend on the
   row order. Every real one-off fit is therefore evaluated on random row orders (seed 20260914; 1000 orders, 100 for
   each Diabetes-130 half) and each test is reported as the median p-value with its 5%-95% range, beside the value in
   the stored row order.
6. **Block 0 (ii)** compares harness and installed package on every replicate, including six constructed samples with no
   event or no non-event, where both must give no p-value.
7. **E7 corrections.** With the 24 projection cells, the inherited July seed overlaps are 91 pairs with offsets 2, 3, 5,
   6, 12 and 15 (`grid_proj_power.R` included); none involves a new cell. E7's parallel efficiency (7.5x at 20 workers)
   was observed in runs A, F and H, not yet by the battery code.
8. **Reporting clarifications.** C4 "both ways" means the unconditional rejection rate with the declined rate, and the
   rejection rate among samples with a p-value. Pulkstenis-Robinson is left out of averages over cells without a
   categorical covariate, where it is undefined. The rule G is round(n/25) (the manuscript's floor will be corrected).
   Run A computes G = 10 and 20 in one replicate (16 cells). `null_calibration_checks.R` is not re-run. The G
   sensitivity is run for size only. The Stukel LR refit declines only on an error or non-convergence (the old harness
   also declined on a separation warning), so runs A, F, G and K2 can move for that reason as well as from new seeds.
   Each alternative uses the null named in `battery/cells.csv`.

### E9. Block 8: two recent slow rivals on the same data (author decision, declared before launch)

Written 2026-09-14, before any battery cell was run at its real B. Block 8 runs after block 7 and does not change
blocks 0-7.

- **Tests.** (i) The projection test of Liu, Li, Chen, Haerdle and Liang (2024, Statistics and Computing 34:175),
  `proj_pvalue(y, X, B = 250)` from `simulations/_proj_test.R`, as in `grid_proj_power.R`. (ii) BAGofT (Zhang, Ding and
  Yang 2023, JASA), `BAGofT(testModel = testGlmBi(formula = <the cell's fitted formula>, link = "logit"), data = ...)`
  at the package defaults (BAGofT 1.0.0: nsplits = 100, ne = floor(5 sqrt(n)), nsim = 100, parFun = parRF()). The
  earlier `nsplits = 1` benchmark is handicapped and is not used.
- **Cells (30), all with n <= 1000.**
  - block 3: `stk_long`, `stk_short`, `stk_asym`, `cauchit`, `t4`, `loglog` at n = 500 and 1000, and their nulls
    `1a null_link_n500` and `null_link_n1000`;
  - block 2: every cell with n <= 1000 (`loglog_base_n610`, `loglog_base_n1000`, `cauchit_auc_n460`,
    `cauchit_auc_n760`, `loglog_auc_n380`, `loglog_auc_n640`, `loglog_e12_n1000`) and their seven nulls;
  - block 4: `crossover_n1000` and `null_crossover_n1000` (off-index, where the projection test is expected to win).
- **Replicates.** The projection test on replicates 1-500 of each cell, BAGofT on replicates 1-200. Each data set is
  regenerated from its stored seed (`seed_base + rep`), and the rival's own random numbers continue that replicate's
  stream. A replicate is kept only if the regenerated data reproduce the stored `EDGE.poly3.u.G10` p-value to 1e-8;
  failures are counted and reported. An error or no p-value counts as no rejection (E0.5).
- **Analysis.** Raw and size-adjusted power (against the same replicates of the matched null), paired with EDGE-poly3,
  EDGE-sym, Stukel's joint score test and GiViTI on the shared replicates (McNemar at the nominal level for tests that
  hold size, Holm-corrected within this block). Median computing time per data set is reported for each rival and n.
  Results go to Supporting Information with one paragraph in the main text. Block 8 is not part of the Section A
  weighting decision or of H1-H6.
- **Scope limit, stated in advance.** Neither rival is run above n = 1000: the projection test builds an n x n kernel
  and refits the model in every bootstrap replicate (cost grows like n^3), and BAGofT takes about 7.5-12.5 minutes per
  data set at n = 1000. The symmetric-tail cells of block 2 at n = 1,700-20,000 therefore have no slow-rival column; the
  paper says so and reports the measured times instead.
- **Cost.** About 30 CPU-equivalent wall-hours for the projection test and 35-40 for BAGofT, run after the main battery.

### E10. Final pre-launch record (corrections to E8 after the last verification)

Written 2026-09-14, after the package correction and its independent verification (`PKG280_verify_correction.md`,
`BATTERY_verify_prelaunch.md`) and before any cell was run at its real B. It corrects E8's wording and records three
small definitions; nothing in A-D changes.

1. **Fitted risks in the Stukel score statistics.** Stukel's joint score and the one-parameter symmetric score use the
   raw fitted risks of the working model, as `gof_stukel` does (the EDGE statistics keep the package's clamp to
   [1e-6, 1 - 1e-6]). This moves p-values only on samples with a fitted risk outside that range (sparse cells); on the
   verification draws no decision at 0.01, 0.05 or 0.10 changed.
2. **Guard flag.** `flag.info_guard` is 1 whenever the guard of E8.2 drops a column, including a column whose computed
   post-fit information is zero or negative.
3. **E8.4 totals.** Block 0 plus 471 cells; 1,782,009 rows including 7,200 row-order rows of block 6; about 46 CPU-hours.
   Block 0's timing measured a 12.3x speed-up at 20 workers on a 3.6-second test, so the run should take 4-6.5 hours.
4. **E8.5.** Row order k uses `set.seed(20260914 + k)`; order 0 is the stored order. Medians and 5%-95% ranges are over
   the random orders only.
5. **E8.6.** On samples with no event or no non-event, the harness gives no p-value for every test. The package does the
   same for `def.gof` and the Stukel forms; its HL, HL_F, Tsiatis and Xie rows still return (non-rejecting) p-values.
   Block 0 (ii) therefore compares harness and package on the EDGE and Stukel statistics there.
6. **E8.7.** The inherited July seed overlaps are 91 pairs with offsets 2, 3, 5, 6, 9, 10, 12 and 15; none involves a
   new cell.
7. **Block 8 (E9)** is not part of `run_M_all.R`; it gets its own script and step after block 7.

### E11. Block 8 analysis details (written after its independent review, before any block 8 result exists)

Written 2026-09-14 while blocks 1b-7 were still running. No battery summary and no block 8 output had been read.
Sources: `BLOCK8_build_report.md`, `BLOCK8_review.md`.

1. **Decision rule for the two Monte Carlo rivals.** The projection test's p-value lies on a 1/250 grid and BAGofT's
   on a 1/100 grid (`mean(meanPv > simulated)` over nsim = 100), so "p <= alpha" would give BAGofT a level of 0.059
   at 0.05. Both rivals therefore reject at **p < alpha** (as the earlier benchmarks did); EDGE, Stukel and GiViTI keep
   p <= alpha. Size-adjusted power uses the matched null's critical value as in the battery, with the null's own
   rejection rate at that value reported.
2. **Comparison family for the Holm correction.** Per cell and rival: EDGE-poly3 and EDGE-sym at the rule G in the
   weighting form that Section A selects, Stukel's joint score test, and GiViTI (4 comparisons). The G = 10 and
   other-form rows are reported beside them, outside the family. A comparison enters the family only if both tests
   hold size in that cell's null (size at 0.05 at most 0.05 + 3 MCSE; the rival's size from its block 8 null
   replicates, the comparator's from the full battery null).
3. **Identity gate** also requires the seed, n and number of events to match the stored replicate. Each rival's
   random numbers start right after that replicate's data, so the projection p-value does not depend on whether
   BAGofT ran.
4. **Order and start.** Block 8 starts only after `battery/launch.log` records that blocks 0-7 finished; the projection
   test runs first, then BAGofT.
5. **Cost, corrected from E9.** Measured timings give about 8-11 wall-hours for the projection test and about 40-55 for
   BAGofT at 20 workers. The paper reports the measured median seconds per data set, noting that they were measured
   with 20 workers busy.

### E12. How Section A and H1-H6 are computed (written while blocks 3-7 run, before any summary value is read)

Written 2026-09-14. Only the column names of `_summary.csv` had been looked at. Sections A-D are unchanged; this fixes the
readings they leave open. The analysis script implements exactly this.

1. **Numbers used.** Rows with `subset = all` at alpha = 0.05 (0.01 as well for the size gate). Size = `rejection` in a
   null cell (declined = no rejection); power = `size_adj_power` of the alternative against its matched null.
2. **Partition arm.** Section A and H1-H3, H5 use the rule-G arm (`.Grule`) of EDGE, HL and HL_F, as B makes it the
   headline; the same computations at G = 10 go to Supporting Information. HL_w, Pigeon-Heyse, Tsiatis, Xie and PR exist
   at G = 10 only. H4 uses G = 10 throughout, because it repeats the current manuscript's census, and reports the rule-G
   version beside it.
3. **Size gate of A.** "Every null cell of the battery" = every null cell, in any block, whose summary contains the
   variant at the rule-G arm. With some 120 null cells and two levels, a correctly sized test fails at least one 3-MCSE
   check by chance with probability of roughly 25-30%; the gate is applied as declared, and every failing cell is
   reported with its z-value. If exactly one variant is eligible it is chosen; if neither is, the unit form is kept.
4. **The six families of A** (alternative cells, blocks 2-4):
   1. symmetric tails: block 2 probit, cauchit and t4 cells (all designs and n, the t4 control included) and block 3
      `link_probit_n*`, `cauchit_n*`, `t4_n*`;
   2. asymmetric links: block 2 loglog cells and block 3 `link_cloglog_n*`, `loglog_n*`, `stk_asym_n*`;
   3. Stukel symmetric family: block 3 `stk_long_n*`, `stk_short_n*`;
   4. omitted terms: block 3 `quad_*`, `binint_*`, `contint_*` (every severity and n);
   5. rough misfit: block 3 `rough_osc2_n*`, `rough_osc4_n*`, `rough_sawtooth_n*` and block 4 `osc4_n*`, `sawtooth_n*`
      (the bump is not named in A and is left out);
   6. off-index: block 4 `crossover_n1000`, `crossover_n2000`.
   Not in A: the plateau cells (E0.2), the 8 census cells, the 24 projection-grid cells of block 4 (B 500, kept for the
   pairing), the skew and sparse cells of block 4, and blocks 5-8. A family average is the plain mean over its cells.
5. **Hypotheses**, with EDGE-poly3 and EDGE-sym in the weighting form that A selects (both forms also reported):
   - H1: over family 1 as defined above, mean(EDGE-poly3 minus GiViTI) >= 0.10. Against: a smaller mean, or a design
     (block 2 base, auc, e12; block 3 cells count as base) whose mean difference is below -0.014.
   - H2: over family 1, mean(EDGE-sym minus Stukel joint) >= -0.014.
   - H3: over family 2, mean(EDGE-poly3 minus GiViTI) >= -0.05.
   - H4: block 3 cells at n = 1000 of families 2-4 plus the 8 census cells, the census rule of
     `make_headline_recount.R`: detectable = some compared test (EDGE-poly3, HL_F, HL, HL_w, PH, Tsiatis, Xie) reaches
     size-adjusted power 0.15; saturated = every one of EDGE-poly3, HL_F, HL, HL_w, PH reaches 0.97; EDGE-poly3 leads or
     ties when it exceeds the best of HL_F, HL, HL_w, PH minus 0.014. Holds if at least 19 detectable, non-saturated
     cells are led or tied (the share is reported too).
   - H5: block 2 base-design probit and cauchit cells: mean(Stukel one-parameter score minus EDGE-sym) <= 0.02; the same
     mean over the block 4 skew cells at the same n is reported as the skewed-covariate cost.
   - H6: declined rate of every test in every cell, overall and within events < G and events >= G.
6. **Paired tests (C3).** In each cell of a hypothesis, an exact McNemar test at 0.05 on the two tests' raw decisions over
   the shared replicates, when both tests hold size in that cell's null (size at 0.05 at most 0.05 + 3 MCSE, C1);
   Holm-corrected within each hypothesis. Cells where a test fails size are listed and not used as support.

### E13. Readings left open by E12, fixed after the analysis review and before any result is read

Written 2026-09-16. Blocks 0-7 finished on 2026-09-14 at 13:13 (every step exit 0, no replicate error); block 8 started
on 2026-09-16. No summary value or per-replicate value of the battery has been read. Source: the independent review of
`analyse_M_battery.R` (`battery/_review/analysis/review/`).

1. **Monte Carlo standard error for size.** The size gate of A (E12.3), "holds size" in C1, E11.2 and E12.6 all use the
   nominal standard error sqrt(alpha (1 - alpha) / B), with B the replicates of the null cell: the standard error of the
   realised size when the true size equals alpha. The standard error of the realised size is reported beside it but
   decides nothing.
2. **Missing values.** A family mean, a macro-average or a hypothesis mean is computed only when every cell in it has a
   value for every test involved. Otherwise the analysis gives no decision for that item and lists the cells; it never
   averages over fewer cells or families. A size-adjusted power whose matched null gave no p-value in at least
   1 - alpha of its replicates counts as no value.
3. **Size failures inside a hypothesis.** The verdict follows E12.5 over all cells. Each hypothesis row also lists the
   cells where either test fails size in its matched null and gives the mean without those cells. If the two verdicts
   differ, the hypothesis is reported as unresolved because of size (C1).
4. **H4 paired test.** EDGE-poly3 against the best of HL_F, HL, HL_w and PH in the same cell, over the detectable and
   non-saturated cells only; Holm over those cells.
5. **Count.** The size gate covers 135 null cells (131 without the four run H external nulls, where the unit and score
   forms are not computed), not "some 120"; the chance that a correctly sized form fails at least one check by chance
   is therefore about 30%.

### E14. A test that cannot run in a cell is "not applicable", not "no value"

Written 2026-09-16, after the verified analysis refused the real root and named the cause. Disclosure first: the
provisional analysis I ran on the same day (my own implementation, reported to the author) already used the convention
below, because it is the convention of `make_headline_recount.R`, and it produced the H4 count 17 of 19 at G = 10. This
addendum is therefore a disclosed clarification, not a blind decision.

**The case.** In the census cell `census_int_binbin_n1000` the model has two binary covariates and nothing continuous.
Pigeon-Heyse, Tsiatis and Xie cannot be computed there at all: they need distinct fitted risks or a covariate space to
cluster, so they give no p-value in the alternative and in its matched null. Under E13.2 that makes H4 "no value",
which would discard a declared hypothesis because two covariate-space rivals do not apply to a categorical design.

**The rule.** A test is **not applicable** in a cell when it gives no p-value in at least 1 - alpha of the replicates of
both the alternative and its matched null. A not-applicable test is left out of that cell's comparisons (the census
maximum, the best partition test, the detectability and saturation checks), and the cell records which tests were left
out. This is what `make_headline_recount.R` does with `na.rm = TRUE`, so the rebuilt census is comparable with the
published one. E13.2 keeps its force for every test that is defined in the cell but lacks a value: that still gives no
value and no decision.

**Limits.** Not applicable is decided per cell from the replicate rates, never by hand. If EDGE-poly3 itself were not
applicable in a cell, the cell has no value and H4 has no value. The tests left out are named in `h4_census.csv` and in
the paper's table, so a reader can see that the census of that cell rests on fewer rivals.
