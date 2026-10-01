# Pre-declaration: block R2, the same demonstration on a covariate the model actually uses (written before it is run)

Date: 2026-09-21, immediately after block R returned a null result. This document states what block R
found, why it found nothing, and what block R2 changes. It is written before block R2 is computed.

## 1. What block R found

Block R (declaration `e80aaef5`) corrupted `num_lab_procedures` in the validation half of the Diabetes
130-US Hospitals cohort, multiplying it by four in $k$ records, with the outcome drawn from the frozen
model. At $k = 0, 10, 50, 100$ every test stayed at its level:

| $k$ | EDGE (rule $G$) | EDGE ($G{=}10$) | HL ($G{=}10$) | Stukel joint | cubic LR |
|---|---|---|---|---|---|
| 0 | 0.044 | 0.064 | 0.064 | 0.048 | 0.072 |
| 10 | 0.054 | 0.046 | 0.030 | 0.060 | 0.054 |
| 50 | 0.042 | 0.040 | 0.032 | 0.052 | 0.042 |
| 100 | 0.036 | 0.034 | 0.056 | 0.058 | 0.042 |

Declared criterion 1 (level) holds. Criterion 2 (the directed test at or below 0.10 at $k=10$) holds.
**Criterion 3 fails: no record-level test reaches 0.20 at $k=10$, or anywhere.** Under the reading
fixed in advance, the separation is reported as absent on this covariate and the reason must be given.

## 2. Why it found nothing, and why that is worth reporting

The reason is measurable and it follows from the paper's own mechanism. Leverage is distance in the
\emph{fitted linear predictor}, not in the covariate. The fitted coefficient on `num_lab_procedures`
is $-0.0004$, so multiplying that covariate by four moves a record's linear predictor by $0.08$ at the
ninetieth percentile of its distribution, against a spread of $0.42$ across the cohort. The corrupted
records are not leverage points at all; they are ordinary records with a wrong number in a column the
model barely reads.

That is a finding rather than a failed experiment, and Section 7 will report it as one: **a gross
error in a covariate the model gives no weight is harmless to every test, grouped or not.** It is also
the practical rule a clinical analyst wants, because it says where to look first when data quality is
in doubt --- at the covariates with the largest coefficients, not the ones with the widest range.

## 3. What block R2 changes

One thing: the covariate. `number_inpatient` carries the model's largest coefficient, $0.271$, and a
fourfold error in it moves the linear predictor by $1.63$ at the ninetieth percentile --- close to four
standard deviations of the clean spread, which is a leverage point in the sense the paper studies.
Everything else is block R's design unchanged: the same cohort, the same split seed, the same frozen
model, the outcome drawn from the frozen predictions so that the model is correct by construction, the
same $k = 0, 10, 50, 100$, the same $500$ replicates, the same five tests, rejection at $p \le 0.05$.

This is a covariate chosen after seeing block R's result. It is chosen by a stated rule that does not
involve any test's behaviour --- the largest fitted coefficient --- and both blocks are reported, the
null one first.

## 4. What would count against the paper (fixed here)

The paper claims a grouped test keeps its level under gross covariate errors while record-level tests
do not. On this cohort that claim requires all three:

1. **Level.** At $k=0$ every test lies inside $0.05 \pm 3$ nominal standard errors, $0.021$ to $0.079$.
2. **Protection.** At $k=10$ the directed test at the rule $G$ reads at or below $0.10$.
3. **Separation.** At $k=10$ at least one record-level test reads at least $0.20$.

If (2) fails, the protection does not reproduce on real covariates and Section 7 says so, in the
abstract as well. If (3) fails again, then on this cohort no covariate produces the separation at ten
corrupted records, and the paper reports that the simulated designs are more adversarial than this
cohort allows --- which would qualify the motivation and belongs in the discussion.

## 5. Written prediction

A shift of $1.63$ in the linear predictor is about a fifth of the shift the simulated design produces
($\pm 7$), so I expect the separation to appear but to be smaller and to need more records than the
simulations do: little at $k=10$, something visible at $k=50$, and a clear gap at $k=100$. I expect
the directed test to hold its level throughout, because $n$ is fifty thousand and its groups hold
twenty-five records whatever happens at the extremes, and I expect the ranking to be the paper's:
record-level tests first and furthest.

## 6. Outputs

`simulations/battery/R2/`, one per-replicate file per $k$, read against this document beside block R's.
