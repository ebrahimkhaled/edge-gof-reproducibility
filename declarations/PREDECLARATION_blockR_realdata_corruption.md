# Pre-declaration: block R, corrupted predictors on the real cohort (written before any block R scenario is run)

Date: 2026-09-20. The paper's headline property --- that a grouped directed test is unmoved by a handful of
corrupted predictors while record-level tests are not --- is measured only on simulated designs. Section 7
analyses a real cohort but never exhibits that property there. This block closes the gap on the cohort the
paper already uses, and it is declared before it is run because its outcome will be quoted in Section 7.

## 1. Why the demonstration is semi-synthetic, and what that costs

On the real outcome the fitted model is genuinely miscalibrated: every test in Section 7 rejects it at
$p<10^{-11}$, so adding corrupted records to it could not show a false alarm, only a shift in a statistic
that is already far past its threshold. A false alarm needs a model that is right. We therefore keep the
real cohort's covariates and the real frozen predictions, and draw the outcome from those predictions, so
that the model is correct by construction on real covariate structure, real risk distribution and real
sample size. Only the outcome is synthetic, and the demonstration is labelled as such wherever it is quoted.

## 2. Design

- **Data.** The validation half of the Diabetes 130-US Hospitals cohort as Section 7 builds it: the same
  exclusion, the same fixed split seed, the same main-effects model fitted on the development half, and the
  same frozen predictions $\hat\pi$ on the validation half.
- **Outcome.** In each replicate, $y_i \sim \text{Bernoulli}(\hat\pi_i)$ independently. The frozen model is
  then the true model for that replicate, and every rejection is a false alarm.
- **Corruption.** $k$ records are drawn at random and their `num_lab_procedures` is multiplied by four; the
  prediction of those records is recomputed from the frozen model at the corrupted covariate, while their
  outcome stays as drawn from the original prediction. This is mechanism C1 of
  Section~4.6, applied to a real covariate. Run at $k = 0, 10, 50, 100$.
- **Replicates.** $500$ for every $k$.
- **Tests.** All in external mode, since the predictions are frozen and nothing is estimated on this half:
  the directed test at the rule $G$ and at $G=10$; the Hosmer--Lemeshow sum at $G=10$; Stukel's joint score
  on the frozen linear predictor; and the cubic calibration likelihood-ratio test on the frozen linear
  predictor. Rejection at $p \le 0.05$.

## 3. What would count against the paper (fixed here)

The paper claims the directed test is protected up to about ten corrupted records and that record-level
tests are not. On this cohort that claim requires all three:

1. **Level.** At $k=0$ every test lies inside $0.05 \pm 3$ nominal standard errors ($0.021$ to $0.079$ at
   $B=500$).
2. **Protection.** At $k=10$ the directed test at the rule $G$ reads at or below $0.10$.
3. **Separation.** At $k=10$ at least one record-level test reads at least $0.20$, that is, at least four
   times the nominal level.

If (1) fails the block is reported as uninformative and nothing is quoted from it. If (2) fails, Section 7
reports that the protection measured on the simulated designs does not reproduce on this cohort, and the
scope claim of Section 6.5 is qualified accordingly. If (3) fails, the separation is reported as absent on
real covariates, which would weaken the paper's motivation and must be said in the discussion.

## 4. Written prediction (recorded before the run)

`num_lab_procedures` runs to about $130$ with a fitted coefficient near zero, so multiplying it by four
moves a record's linear predictor by much less than the $\pm7$ of the simulated design: the corruption is
milder, and the separation should therefore be smaller than Table~\ref{tab:contamination}'s. I expect the
directed test to hold its level at $k=10$ and $k=50$ and the record-level tests to rise, but I would not be
surprised if $k=10$ is too small to move anything at $n\approx50{,}000$ and the separation appears only at
$k=50$ or $k=100$. If nothing separates at any $k$, that is the honest result and Section 7 says the
mechanism does not bite on this covariate at this scale.

## 5. Outputs

`simulations/battery/R/`: one per-replicate file per $k$ with a column per test, and a summary read against
this document.
