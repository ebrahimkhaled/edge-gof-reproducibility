# Pre-declaration: block 9b, are the slow rivals robust to a corrupted record? (written before any block 9b cell is run)

Date: 2026-09-16. Author decision of the same day. Motivation, named here as **exploratory only**:
`simulations/robust_influence_rivals.R` moved one record's covariate along a grid on a single data set (n = 500). BAGofT
never rejected, the projection test of Liu et al. (2024) held until the record's linear predictor passed about 8 and then
rejected, and EDGE-poly3's unit form held across the grid. That was one data set with one call per point, and BAGofT's
p-value is itself random, so it cannot support a claim. Block 9b replicates the question with independent data sets.

It belongs to Paper A's robustness section, beside block 9, and changes nothing in blocks 0-9.

---

## 1. Design

**Data.** Base design: x ~ U(-3, 3), d ~ Bernoulli(0.5), eta = 0.6x + 0.5d, logistic truth, fitted as y ~ x + d, n = 500
(the largest n at which BAGofT is affordable).

**Cells (3).** In each replicate, k = 3 records chosen at random have their covariate corrupted while their outcome is
still drawn from the truth at the original x:
- `clean`: no corruption;
- `x4`: the three covariates multiplied by 4 (linear predictor of a corrupted record up to about 7.7, the region where
  the exploratory curve showed Stukel and GiViTI already rejecting and the projection test still holding);
- `x8`: the three covariates multiplied by 8 (up to about 15, the region where the projection test broke).

**Replicates.** B = 100 per cell. At this B the nominal standard error of a 5% rate is 0.022, so the 3-standard-error
bound is 0.115; a rate of about 0.2 or more is clearly separated from it. This resolution is chosen because BAGofT takes
15-19 minutes per data set here.

**Tests, all on the same data set.** The projection test (`proj_pvalue(y, X, B = 250)`, X the model matrix), BAGofT at the
package defaults (BAGofT 1.0.0: nsplits = 100, nsim = 100, ne = floor(5 sqrt(n)), parFun = parRF()), EDGE-poly3 and
EDGE-sym unit forms at the rule G, Stukel's joint score, one-parameter score and likelihood-ratio refit, GiViTI (default
threshold), the cubic calibration LR, and HL at G = 10. Seconds per data set are stored for the two slow rivals.

**Seeds.** Per replicate, `set.seed(seed_base + rep)`, `seed_base = 400000000 + cell_id * 10000`; the rivals' own random
numbers continue that stream (as E9).

**When.** After block 8 finishes and after block 9, on the battery's 20 workers: about 5 hours, almost all BAGofT.

## 2. Rules

1. The bulk model is correct in every cell, so a rejection is a **false alarm**. The projection test and BAGofT reject at
   p < 0.05 (their p-values lie on grids, E11.1); every other test at p <= 0.05.
2. A test that gives no p-value counts as no rejection, and its rate is reported.
3. Within each cell, each rival is paired with EDGE-poly3's unit form on the same replicates (exact McNemar at 0.05,
   Holm-corrected within block 9b).
4. Every per-replicate p-value and time is stored.

## 3. Claims, and what would count against them

- **B9b.1 (BAGofT is robust).** BAGofT's false-alarm rate is at most 0.115 in all three cells. *Against:* above 0.115 in
  any cell.
- **B9b.2 (the projection test is robust, then breaks).** Its rate is at most 0.115 in `clean` and `x4`, and above 0.115
  in `x8`. *Against:* either half fails.
- **B9b.3 (EDGE-poly3 unit is robust).** Its rate is at most 0.115 in all three cells. *Against:* above 0.115 in any cell.
- **B9b.4 (reported, no claim).** The rates of Stukel's three tests, GiViTI, the cubic LR, EDGE-sym's unit form and HL;
  the fitted coefficients; the seconds per data set of the two slow rivals.

## 4. Reading fixed in advance

Whatever the outcome, the paper reports all three cells for every test. If B9b.1 and B9b.3 hold, robustness is shared by
EDGE-poly3, HL and BAGofT, and the paper's claim is the combination: among tests with competitive power on
calibration-shape misfit, EDGE-poly3's unit form is the robust one, while the equally robust rivals lack power (HL as G
grows) or speed and power (BAGofT). If B9b.3 fails, the robustness claim for EDGE is reduced to what block 9 confirms.
If B9b.2 fails in the other direction (the projection test robust in `x8` too), the projection test joins the robust
group and the paper says so.

## 5. Outputs

`simulations/battery/9b/` with one per-replicate file per cell, a summary, and the paired tests.
