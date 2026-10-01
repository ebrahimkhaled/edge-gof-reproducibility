# Pre-declaration: block TW, which property protects a test from corrupted records?

Date: 2026-10-01. Written and hashed before any replicate was computed. Prompted by the SiM-style Associate Editor
review of 2026-10-01 (M1-M4, M5.3, M6.6), which argued that the protection may come from bounded weights rather
than from grouping, since Spiegelhalter's z (record-level, bounded weights) stayed protected in block EXT.

## Arms (all on the same data sets; the model y ~ covariates is refitted by maximum likelihood)
Grouped:
1. EDGE.default  edge.gof(fit, G = "auto", basis = "poly3"), Satterthwaite
2. EDGE.imhof    the same with method = "imhof" (reviewer M5.3)
3. EDGE.G10      edge.gof(fit, G = 10, basis = "poly3")
4. HL.G10        Hosmer-Lemeshow, ten groups
Ungrouped, probability-scale (bounded) directions:
5. TWIN.score    Rao score test for adding the cubic orthogonal polynomial in the fitted risk, poly(pi_hat, 3), to the
                 model: u = Z'(y - pi_hat), I = Z'WZ - Z'WX (X'WX)^-1 X'WZ, chi-square(3). The ungrouped twin of EDGE.
6. EDGE.Gn       edge.gof(fit, G = n): one record a group, each residual standardised by its own variance (n = 1000 only).
7. SPZ           Spiegelhalter's z on the fitted risks, two-sided.
Ungrouped, logit-scale (unbounded) directions:
8. Stk.joint     Stukel's joint score test.  9. Cubic.LR  cubic calibration LR.  10. GiViTI  belt, internal, thres 0.95.
Diagnostics-then-test baseline (reviewer M3):
11-14. PRG.Stk, PRG.Cubic, PRG.EDGE.default, PRG.EDGE.G10: fit, drop every record with Pregibon's
       DeltaX2_i = rP_i^2 / (1 - h_i) > 4 or DeltaBeta_i = rP_i^2 h_i / (1 - h_i)^2 > 1, refit, then test.
       The number dropped is recorded.

## Part A: the corruption design of blocks 9 and C1b (c1b_gen: x ~ U(-3,3), d ~ Bern(0.5), eta = 0.6x + 0.5d)
n = 1000: clean logit; C1 (x4) k = 1, 2, 5, 10; C1b (x-4) k = 1, 2, 5, 10; cloglog truth clean; C2 (the ten
highest-risk outcomes set to 0, b9_gen) k = 10; localised tail misfit: true risk multiplied by 0.6 for the records
above the 98th percentile of true risk, model correct elsewhere.
n = 5000: clean logit; C1 k = 25, 50; C1b k = 10, 25; probit truth clean; localised tail misfit at 0.8.
Part B (reviewer M6.6): n = 2000, twelve covariates (eight N(0,1), four Bernoulli(0.3)), eta = -1 + sum b_j x_j with
b = (0.5, -0.4, 0.3, 0.3, -0.2, 0.2, 0.1, -0.1, 0.4, -0.3, 0.2, 0.2): clean logit; cloglog truth; omitted
interaction 0.5 x1 x2; C1 on x1 k = 10; C1b on x1 k = 10; C1 on x1 k = 25.
1000 replicates a cell; seed 980000000 + 10000 * cell + rep.

## Reading, fixed now
- False alarms: rejection at 5% in the corrupted cells; a test "holds" at <= 0.10, as in the paper.
- Size: clean cells, band 0.05 +/- 3 SE.
- Power: raw and size-adjusted at the clean cell of the same n and design.
- Mechanism verdict: if TWIN.score and SPZ hold where Stukel, cubic and GiViTI fail, bounded probability-scale
  directions protect without grouping, and the paper's mechanism must be restated. If TWIN.score fails where EDGE
  holds, grouping is (part of) the mechanism. EDGE.Gn separates standardisation by each record's own variance.
- Either outcome is reported.

## Outputs
simulations/battery/TW/<cell>_pvalues.csv.gz, _progress.log; runner simulations/run_M_blockTW.R.
