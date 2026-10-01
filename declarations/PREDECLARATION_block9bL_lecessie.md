# Pre-declaration: block 9bL, le Cessie on the data sets of block 9b (written before any 9bL data set is tested)

Date: 2026-09-19. Author's request of the same day: le Cessie is missing from the table of block 9b.

## 0. What is already known, and therefore why no claim is made

Block 8L (sha256 0c1a1731) found le Cessie's test NOT moved by the x4 corruption at n = 1000 (its claim L3 failed:
0.055, 0.052, 0.052 and 0.077 at k = 1, 2, 5 and 10). Block 9b (sha256 72d703a9) put three corrupted records in
n = 500 at x4 and at x8, and le Cessie was not among its tests. This addendum adds it on the same 300 data sets.
Because the x4 behaviour is already known, a prediction written now would not be a test of anything, so this
block makes **no claim**: every rate is reported, in the same place and type as the tests of block 9b.

## 1. Implementation

Exactly block 8L's: `run.all.gof(fit, tests = "le-Cessie")` of ebrahim.gof >= 2.6.0 (corrected moment reference
(I - H)' R (I - H), checked against its definition before the run; the vendor form must not match), package default
bandwidth, closed-form reference, no random numbers.

## 2. Data

Nothing is simulated afresh. The 300 data sets of block 9b -- cells clean, x4 and x8, replicates 1-100 each --
are regenerated from the seeds block 9b stored (`set.seed(400000000 + cell_id * 10000 + rep)`, `b9b_gen`), and a
data set is used only if EDGE-poly3 unit at G = 10 recomputed on it equals the stored value to 1e-8 and seed, n and
events agree (block 8's identity gate). The fitted model is `y ~ x + d`, logistic, as in block 9b.

## 3. Rules and outputs

1. Rejection at p <= 0.05; a data set with no p-value counts as no rejection, and that rate is reported.
2. Reported: the rejection rate in each of the three cells, beside EDGE-poly3 (unit, rule G) on the same replicates,
   with an exact McNemar test of the pair in each cell (Holm across the three); the median seconds per data set.
3. The run is serial (one core): the machine's other cores are in use, and 300 data sets need about a minute.

Outputs: `simulations/battery/9bL/<cell>_lecessie_pvalues.csv.gz`, `_summary.csv`, the identity-gate record.
