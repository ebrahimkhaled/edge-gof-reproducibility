# Pre-declaration: block EXT_SP, the cumulative calibration test of Sadatsafavi and Petkau on block EXT

Date: 2026-11-19. Written and hashed before any replicate of this block is computed (the only computation done
before it: installing the R package cumulcalib 0.2.0 from CRAN and running its help-page example once).

## Why
The novelty search for the Computational Statistics version (PRECEDENTS.md) found that Sadatsafavi and Petkau
(2024, Statistics in Medicine 43:3524-3538, doi 10.1002/sim.10138) is the closest recent precedent for EDGE's
external mode: a calibration test for frozen predicted risks with no grouping, no tuning constant and no resampling,
built on standardized partial sums of prediction errors and their Brownian-motion and Brownian-bridge limits.
Statistics in Medicine desk-rejected the EDGE paper for novelty, and the paper currently says that this test was
not studied under corrupted records. This block puts it into the paper's external-validation study (block EXT,
Online Resource 1, Section S7) on the same data sets, so every comparison with EDGE is paired.

## The arm
SP-BB: cumulcalib::cumulcalib(y, p, method = "BB")$pval, cumulcalib 0.2.0. This is the package's default main
method (method = c("BB", "BM") makes "BB" the main method): the two-part Brownian-bridge test, which combines a
mean component (S_n, two-sided normal p-value) and a bridge component (sup |S(t) - t S_n|, Kolmogorov p-value) by
Fisher's method on chi-square(4). Defaults otherwise: ordered = FALSE, ties = "group", n_sim = 0 (analytic
p-values). The predictions are continuous, so no ties are expected; any tie message is logged.
Recorded alongside, descriptive only (no reading depends on them): the two BB component p-values and the one-part
Brownian-motion test p-value (method "BM").
A replicate with no p-value counts as a non-rejection, as in block EXT.

## The cells
All 21 cells of block EXT, unchanged: cells 1-3 correct model (n = 250, 500, 1000); cells 4-13 the ten
alternatives at n = 500 (shift 0.2, 0.4; slope 0.8, 0.6; U-shape 0.4, 0.8; threshold 2, 4; interaction 0.6, 1);
cells 14-21 corrupted records under a correct frozen model at n = 1000 (x2 recorded as 4 x2 or -4 x2 in
k = 1, 2, 5, 10 records). Same generator as run_M_blockEXT.R, same seeds (set.seed(7000000 + 10000 id + rep) for
ids 1-13, set.seed(7600000 + 10000 id + rep) for ids 14-21, Mersenne-Twister / Inversion / Rejection set
explicitly at every replicate), same 1000 replicates per cell.

Check that the data sets are block EXT's: for every replicate the regenerated data must reproduce block EXT's
stored number of events exactly and its stored Spiegelhalter p-value (all.equal at R's default tolerance, as
Addendum 3 to block TW did; the stored values passed through a csv file). If any replicate fails, that cell stops
and nothing from it is used.

EDGE's p-values for the comparison are block EXT's stored EDGE_G10 (ten groups) and EDGE_auto values, so the
pairing is exact; Cox's stored p-values are used the same way.

## Readings, fixed now
(a) Level. Rejection at 5% in cells 1-3. SP-BB holds if within 0.05 +/- 3 SE (0.029 to 0.071), the band of block EXT.

(b) Power. Size-adjusted at the 5% quantile (type 1) of SP-BB's p-values in cell 2 (null, n = 500), per
alternative cell 4-13, exactly as analyse_blockEXT.R does for the other tests. Paired differences per cell:
EDGE at ten groups minus SP-BB, and Cox minus SP-BB, each the mean over replicates of the difference of the two
size-adjusted rejection indicators, with Monte Carlo SE sd / sqrt(1000) as in analyse_blockEXT_paired.R and a 95%
interval of +/- 1.96 SE. Also reported by family (the mean of its two cells) and overall (the mean of the ten
cells), with SE sqrt(sum of the cell SE^2) / (number of cells), since the cells use independent replicates.
Mean size-adjusted power over the ten alternatives is reported for SP-BB next to block EXT's EDGE (ten groups and
default), Cox and GiViTI values.

(c) Protection. Raw false-alarm rate at 5% in cells 14-21. The bar is the paper's: a test holds in a cell where
its rate is at most 0.10. The number of the eight cells SP-BB holds is reported next to EDGE's (8 of 8 at ten groups).

## What counts against EDGE, fixed now
- Against EDGE (the paper must say SP-BB is at least as good on this study): SP-BB holds in all eight corrupted
  cells AND its mean size-adjusted power over the ten alternatives is at least EDGE-at-ten-groups' (overall paired
  difference EDGE minus SP-BB at or below zero).
- Strongly against EDGE (SP-BB dominates): as above, with the 95% interval of the overall paired difference
  wholly below zero.
- Neutral (EDGE's protection is shared, its power edge not shown): SP-BB holds in all eight corrupted cells and
  the overall paired difference is positive but its 95% interval contains zero.
- For EDGE: SP-BB fails the 0.10 bar in at least one corrupted cell, or the overall paired difference is positive
  with its 95% interval wholly above zero.
- Per family, a 95% interval wholly on one side of zero is reported as a lead for that test on that family; it does
  not change the overall reading.
- If SP-BB fails the level band in a null cell, that is reported, and its size-adjusted power is still the reading.

Expectation written in advance (not a pass/fail): SP-BB has a mean component close to Spiegelhalter-type and
Cox-type checks of calibration in the large, so it should be strong on the shift and slope cells and weak on the
interaction; its behaviour under the sign-error corruption at k = 10 is the open question.

## Timing
cumulcalib(y, p, method = "BB") on one core at n = 1000, 5000, 20000, 100000, data from the correct-model generator
of cell 1 (seed 7000000 + 10000 + 1), median and minimum of 5 runs (after one warm-up run), with edge_external at
ten groups timed the same way for comparison. Output: simulations/battery/T/timing_sp.csv.

## Outputs
Runner simulations/run_M_blockEXT_SP.R (every chunk of replicates written to disk as it finishes; restart resumes
from disk); results simulations/battery/EXT_SP/<cell>_pvalues.csv.gz, _progress.log; analysis
simulations/analyse_blockEXT_SP.R writing EXT_SP/_summary.csv and EXT_SP/_paired.csv.
