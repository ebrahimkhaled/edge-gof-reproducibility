# Pre-declaration: the Paper B pilot (written before any pilot cell is run)

Date: 2026-09-16. Author decision of the same day: Paper A (the pre-registered comparison on individual patient data,
battery blocks 0-8) and Paper B (calibration testing when the individual data or the fitted model are not available)
are two papers with separate questions. Run H (external mode) leaves Paper A and becomes Paper B material.

This pilot decides whether Paper B has a finding worth a paper. It is small on purpose: about one day of computing
after block 8 finishes. Nothing in it may change after the first pilot result is read; anything added later is
labelled *post hoc*.

---

## 0. Why a pilot is needed (the state of the evidence)

The battery and the earlier runs already show three things that limit what Paper B can claim.

1. **A grouped directed test cannot out-power its ungrouped twin on individual data.** The score form of a grouped
   basis is the Rao score test for adding that basis as a step covariate (proved, and checked to 1e-9 in the package
   tests). Grouping discards within-group information, so the smooth ungrouped test is at least as good. The battery
   agrees: EDGE-sym and Stukel's one-parameter score test differ by 0.005 on average over 32 symmetric-tail cells.
2. **Two ideas that sounded like niches are not.** On the local "bump" alternative, EDGE-sym scores 0.790 and Stukel's
   one-parameter test 0.793 at n = 500 (0.977 against 0.978 at n = 1000). With a skewed (chi-square) covariate every
   score-type test keeps its size (0.045-0.055); only the refit-based tests drift at small n (Stukel refit 0.081,
   cubic LR 0.075 at n = 100).
3. **The one power lead we have seen** is a maximum over several directed probes with an exact multivariate-normal
   null: cloglog n = 1000 max3 0.864 against Stukel's joint 0.803; the best implementable test in 19 of 36 population
   grid cells (max-probe 19, EDGE 8, GiViTI 8, Stukel 1). It works grouped (score weighting) and ungrouped.

So the pilot asks where a grouped test is the only one that can run, and whether the max-over-directions form is a
real gain.

## 1. Niches, tests and rivals

**N1. Calibration table only.** The analyst has, per risk group: the number of patients, the observed events and the
expected events (a published decile table, an audit dashboard, a federated summary). No individual data.
- EDGE external mode on the table (unit and score form), with the group variance either exact or approximated by
  n_g pbar_g (1 - pbar_g) when it is not reported.
- Rivals that can also run on a table: Hosmer-Lemeshow on the table; a grouped recalibration model (binomial GLM of
  the observed counts on the group logit, testing intercept and slope); the same model with Stukel's symmetric column
  added at the group logit; a grouped Spiegelhalter-type z.
- Stukel's own test cannot run here: it needs individual predictions.

**N2. Predictions only (external validation).** Individual predictions from a model fitted elsewhere, no design matrix.
- EDGE external mode (bases: slope2, poly4, stk4, sym3).
- Rivals: calibration intercept and slope likelihood-ratio test; Spiegelhalter z; GiViTI external; Hosmer-Lemeshow;
  and a Stukel-type recalibration test (logit(y) on logit(p) with Stukel's two columns added), which is the fair
  ungrouped competitor.

**N3. Maximum over directions.** On individual data, the grouped score-form maximum over the symmetric probe and the
two Aranda-Ordaz half-probes, with the exact multivariate-normal null, against Stukel's joint score test, EDGE-sym,
GiViTI and the cubic LR.

**N4. Clinically relevant miscalibration at large n (feasibility only).** At n >= 20,000 every test rejects. A
relevance test inverts the weighted chi-square null with a non-centrality set by a margin on the probability scale
(for example an average absolute calibration error of 1 percentage point). The pilot checks that the inversion is
calibrated in one null cell and one cell with a miscalibration below the margin; no comparison is claimed.

## 2. Data-generating models

Development and validation populations, x ~ U(-3,3), d ~ Bernoulli(0.5), eta = 0.6x + 0.5d.
- **True model**: the development model is the logistic model with those coefficients (so predictions carry no
  estimation error; estimation error is added in one sensitivity cell where the development sample has n = 2,000).
- **Miscalibration shapes in the validation population** (as run H, plus one new local shape):
  - `temp`: logit p_true = logit(p_model) / s, s = 0.85 and 0.7 (predictions too extreme);
  - `asym`: the same but on the upper half only;
  - `probit_comp`: a probit truth matched to the model;
  - `window`: p_true = p_model + delta * w(p_model), with w a smooth bump on [0.05, 0.20] and delta chosen so that the
    average absolute calibration error over that window is 3 percentage points (a decision-threshold departure).
- **Null**: p_true = p_model.
- n = 1,000 and 5,000; B = 2,000 replicates; G = 10 for the table niche (and G = max(10, round(n/25)) reported).

## 3. Rules

1. Size is checked first, with the nominal standard error sqrt(alpha (1 - alpha) / B) and a 3 standard error bound, as
   in the main pre-declaration (E13.1). A test that fails size is reported and not used for a claim.
2. Power is size-adjusted against the matched null of the same cell and n.
3. A test that cannot run in a niche is recorded as "not applicable", never as zero power.
4. Every replicate's p-values are stored; seeds are per replicate, `set.seed(seed_base + rep)`, with fresh bases.
5. Comparisons are paired (same replicates), with exact McNemar tests and Holm correction within each niche.

## 4. Go / no-go, decided before the pilot runs

Paper B goes ahead if **at least one** of these holds.

- **N1**: EDGE on the table beats the best table rival by at least 0.05 size-adjusted power in at least half of the
  alternative cells, while holding size; or it holds size where a table rival does not.
- **N3**: the maximum over directions beats Stukel's joint score test by at least 0.05 in at least half of the link
  cells, while holding size.
- **N2**: EDGE external beats the best prediction-only rival by at least 0.05 in at least half of the cells.

If none holds, Paper B stops. The table-only use then becomes one short section of Paper A ("testing a published
calibration table"), and the maximum over directions is offered to the ensemble paper, where combining bases is the
subject.

If the go comes from N3 alone, the finding is a better *directed* test rather than a grouped one, and it is written
that way: the gain comes from the maximum over directions, not from grouping.

## 5. Outputs

`simulations/pilotB/` with one per-replicate file per cell, a summary in the same shape as the battery's, a size
table, the niche comparisons with their McNemar tests, and a one-page verdict against Section 4.
