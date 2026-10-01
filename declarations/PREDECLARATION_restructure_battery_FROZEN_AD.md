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

To be completed from `MAP_battery.md` before launch (scenario x design x n x G x B), reproducing every current table
and figure and adding the symmetric-tail cells (probit, cauchit, t4, loglog at the base design, AUC about 0.95 and a
12% event rate) and the skewed-covariate cells used for H5.
