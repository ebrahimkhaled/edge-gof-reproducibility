# Pre-declaration: block 9d, does robust estimation solve the problem instead? (written before any 9d cell is run)

Date: 2026-09-19. Written after blocks 9 and 9c reported and before any block 9d data exists.

## 0. The question, and why it has to be answered with data

Block 9 showed that a handful of corrupted covariates raises the false-alarm rate of every ungrouped
goodness-of-fit test, while a grouped directed test is far less disturbed. The obvious objection is
that this is the wrong remedy: if the data contain corrupted records, fit the model with a robust
estimator and the problem goes away, whatever test is used afterwards.

The objection deserves an answer rather than an argument, because the answer is not obvious in
either direction. A robust fit down-weights the corrupted records, which should stop them reaching
the test; but every test in this study refers its statistic to a null distribution derived for
**maximum likelihood** fitting, and in particular the estimation adjustment
`Omega = I - U(X'WX)^{-1}U'` is the ML adjustment. Applying it after a robust fit is not justified,
and might break the tests' level on its own.

Block 9d therefore separates the two effects with a 2 x 2: contaminated or not, ML or robust.

## 1. Design

**Data.** The battery base design, as in block 9: x ~ U(-3, 3), d ~ Bernoulli(0.5),
eta = 0.6x + 0.5d, fitted as y ~ x + d, n = 1000.

**Truths.** (a) logistic, so the fitted model is correct for every uncorrupted record and a rejection
is a false alarm; (b) probit, a real misfit with moderate clean power, so the cost of robust fitting
in power can be measured.

**Corruption.** C1 of block 9, unchanged: k records chosen at random have x multiplied by 4, their
outcome still drawn from the truth at the original x. k in {0, 5, 10, 25}. k = 0 is the control.

**Estimators.** Each replicate is fitted twice on the same data:
- **ML**: `stats::glm`, as everywhere else in this study;
- **robust**: `robustbase::glmrob` with `method = "Mqle"`, the Mallows-type quasi-likelihood
  estimator of Cantoni and Ronchetti (2001), at the package default tuning constant `tcc = 1.345`.
  `method = "BY"` (Bianco-Yohai) is computed where it converges and reported without a claim.

**Cells.** 2 truths x 4 values of k x 2 estimators = **16 cells**, B = 1000 replicates each.

**Tests.** As in block 9 and on both fits: EDGE-poly3 and EDGE-sym, unit and score, at G = 10 and at
the rule G; Stukel's joint score, one-parameter score and likelihood-ratio refit; GiViTI; the cubic
calibration LR; HL and HL_F at both G. Each replicate stores the fitted coefficients of both fits,
the number of corrupted records in the highest-risk group, and whether the robust fit converged.

**Seeds.** `set.seed(seed_base + rep)`, `seed_base = 600000000 + cell_id * 10000`. Disjoint from
blocks 0-8, from block 9 (3e8), 9b (4e8) and 9c (5e8). The two fits of a replicate see the same data.

## 2. Rules

1. Every test rejects at p <= 0.05. Nominal standard error sqrt(0.05 x 0.95 / 1000) = 0.0069; the
   3-SE band is 0.029 to 0.071.
2. A test that gives no p-value counts as no rejection, and its rate is reported. A replicate whose
   robust fit does not converge is counted and excluded from that cell's robust column only; the
   count is reported with the cell.
3. Under truth (a) a rejection is a false alarm; under truth (b) the rejection rate is power and is
   compared with the same truth's k = 0 cell under the same estimator.
4. The ML columns of the k = 0 and k in {5, 10, 25} cells must reproduce block 9's corresponding
   cells to within 3 nominal SE. They are drawn on different seeds, so this is a consistency check,
   not an identity check; a larger discrepancy stops the block and is investigated.

## 3. Claims, and what would count against them

- **D9d.1 (the tests survive robust fitting at all).** Under the logistic truth with k = 0, every
  test's rejection rate under the robust fit is within 3 nominal SE of 0.05. *Against:* any test
  outside the band. **This claim gates the others**: if the reference distributions do not hold
  after a robust fit, the remaining comparisons measure the mismatch rather than the contamination,
  and the block reports that instead.

- **D9d.2 (robust fitting repairs the ungrouped tests).** Under the logistic truth at k = 10, the
  robust fit brings Stukel's joint score, Stukel's one-parameter score, GiViTI and the cubic LR each
  within 3 nominal SE of 0.05, from the 0.40 to 0.63 they reach under ML. *Against:* any of the four
  outside the band.

- **D9d.3 (what robust fitting costs, reported).** Against the probit truth with k = 0, the change in
  each test's power from the ML fit to the robust fit. No threshold.

- **D9d.4 (EDGE under a robust fit, reported).** The same quantities for the grouped directed test.
  No threshold; the grouping and the robust estimator address the same problem by different means,
  and whether they compose or interfere is a measurement, not a prediction.

- **D9d.5 (reported).** The convergence rate of `glmrob`, the mean fitted coefficients of both fits,
  and the Bianco-Yohai estimator's results where it converges.

## 4. Reading fixed in advance

- If **D9d.1 fails**, the block reports that the tests' references are not valid after robust fitting
  and no comparison of contamination remedies is drawn from it. This outcome would itself be worth
  reporting: it would mean the obvious remedy cannot simply be bolted onto the existing tests.
- If **D9d.1 and D9d.2 hold**, robust fitting is a genuine alternative remedy. The paper says so,
  states the cost measured in D9d.3, and the grouping argument is then one of two available routes
  rather than the only one. Two differences remain and are stated: robust fitting requires refitting,
  which is unavailable when a model is validated on frozen predictions, the setting of the paper's
  clinical section; and it changes the model being tested, so the question answered is no longer
  "is this model calibrated here".
- If **D9d.1 holds and D9d.2 fails**, robust fitting does not repair the ungrouped tests, and this is
  the strongest available answer to the objection.

Whichever occurs, block 9d is reported in full in the paper's contamination section, in the same
type as the results that favour the test.

## 5. Outputs

`simulations/battery/9d/`, one per-replicate file per cell, and an analysis restricted to block 9d.
