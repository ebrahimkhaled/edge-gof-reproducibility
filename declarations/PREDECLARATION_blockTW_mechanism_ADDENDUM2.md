# Addendum 2 to PREDECLARATION_blockTW_mechanism.md: the weighted grouped tests of Hosmer and Hjort (2002)

Date: 2026-10-01, before any Hosmer-Hjort replicate was computed. The referee asked for Hosmer and Hjort's weighted
grouped tests (Statistics in Medicine 21:2723-2738), the closest earlier directed grouped tests, to be compared
rather than argued away.

Implementation: simulations/_hosmer_hjort.R, hh_test(fit, g = 10), written from the paper and validated before this
addendum (level at n = 1000 and 5000 within 0.05 +/- 0.02 for the four general-purpose statistics).

Arms added, on the SAME data sets as block TW (same generator, cells and seeds; the runner checks that it regenerates
TW's data by reproducing TW's EDGE.G10 p-values exactly in every replicate):
- HH.HLnp: the weighted Hosmer-Lemeshow statistic with the paper's general-purpose weight z = pi_hat log pi_hat,
  chi-square(g - 2);
- HH.X2np: the same weight with the full covariance, chi-square(rank);
- HH.HL1 and HH.X2_1: the unweighted versions (w = 1), reported for completeness.
The arms that need a guessed omitted covariate (HLop, X2op) are not run: no single covariate is right across cells.

Reading: as declared for block TW (false alarms at 5%, holds at <= 0.10; size and raw and size-adjusted power).
Output: simulations/battery/TW_HH/<cell>_pvalues.csv.gz.
