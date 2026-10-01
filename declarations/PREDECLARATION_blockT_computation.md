# Design note: block T, what each test costs (written before any timing was taken)

Date: 2026-09-19. Author's request of the same day: computation time is one of the paper's strengths and
must be measured, not asserted. Until now the paper's "milliseconds for the directed test" was never
timed, and the resampling tests' seconds were recorded while 20 data sets ran in parallel.

## Design

- **Machine state:** one core (BLAS pinned to one thread), nothing else running: the block runs after the
  DeepGOF study's confirmatory run has finished and no R process has run for ten minutes, and before
  block 8R. AMD Ryzen 9 3900X, 32 GB, R 4.4.
- **Data:** the study's base design (x ~ U(-3, 3), d ~ Bernoulli(0.5), logit truth 0.6x + 0.5d), fitted
  as y ~ x + d; one data set per n, seed 20260920 + n. Every test receives the same fitted model; the
  time of the model fit is not counted, the time of any refit inside a test is.
- **Tests and sample sizes** (the grid fixed here, before any timing):
  - EDGE (`edge.gof`, cubic basis, unit form, rule G), Hosmer-Lemeshow (G = 10), Stukel's joint score,
    the cubic calibration LR (`cubic.calib.gof`) and the GiViTI calibration test
    (`givitiR::givitiCalibrationTest`, called directly): n = 500, 1000, 2000, 5000, 20000, 100000;
    11 repetitions, median reported.
  - le Cessie (`run.all.gof(fit, tests = "le-Cessie")`): n = 500, 1000, 2000 (3 repetitions) and 5000,
    10000 (1); peak memory recorded. n = 20000 needs an estimated 16 GB or more and is not run.
  - The Liu projection test (`proj_pvalue`, B = 250): n = 500 (3 repetitions) and 1000 (1).
  - BAGofT 1.0.0 at its defaults: n = 500, 1000 and 2000, 1 repetition each.
- **Reported, no claim:** median seconds per data set for each test and n, le Cessie's peak memory, and
  which test and n were not run and why.

Outputs: `simulations/battery/T/timing.csv`, `timing_log.txt`.
