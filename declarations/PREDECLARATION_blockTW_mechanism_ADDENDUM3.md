# Addendum 3 to the block-TW declaration: the twin on frozen predictions

Date: 2026-10-01, before any replicate of this addendum ran. Prompted by the final pre-submission review: Spiegelhalter's
z, an ungrouped test with bounded weights, held in the external-validation study (block EXT), while the ungrouped twin
failed after a refit (block TW). Which ingredient protects, then: pooling, bounded weights, or the absence of a refit?

Arms, on block EXT's corrupted data sets (cells 14-21: five N(0,1) covariates, frozen model, C1 x4 and C1b x(-4) on
x2, k = 1, 2, 5, 10, n = 1000; seeds 7600000 + 10000 id + rep, regenerated and checked against EXT's Spiegelhalter
p-values, which must agree exactly) and on its correct-model cell 3 (n = 1000, seed 7000000 + 10000 * 3 + rep):
- TWIN.ext: the score statistic for the cubic orthogonal polynomial in the frozen prediction, u = Z'(y - p) with
  variance Z'VZ (nothing estimated), chi-square(3);
- SPZ.ext (recomputed, the check), and EDGE in external mode at ten groups and at the default partition (recomputed).
1000 replicates. Reading as in block EXT (holds at <= 0.10).
Expected in advance, stated so it can fail: the twin on frozen predictions holds, because no refit removes
information from its direction; if so, the refit, not the grouping alone, is what makes the bounded twin fragile.
