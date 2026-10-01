# Pre-declaration: block 9R3, does the tolerance rule hold away from one per cent? (written before any 9R3 scenario is run)

Date: 2026-09-20. Block 9R2 measured what each setting of EDGE-FR costs in power, but every corrupted
scenario in the study stops at one per cent of the sample, so the protection of the wider settings has never
been seen. The author's question is exactly that: what happens when the tolerance `c` is not one per cent.
This block answers it on fresh data, at corruption rates the study has not used.

## 0. The rule under test

EDGE-FR(c) groups the outer `5c` of the risk scale at each end into one group and leaves the middle at the
rule of about twenty-five records a group. The proposal is that the tail group should hold about five times
the corruption the analyst tolerates, so that `c` alone sets the partition:

| tolerance `c` | tail fraction `f = 5c` | the variant of block 9R2 |
|---|---|---|
| 1% | 5%  | HYB05 |
| 2% | 10% | HYB10 |
| 4% | 20% | HYB20 |

## 1. Scenarios (fresh seeds, disjoint from every other block)

Seed bases `390000000 + j * 10000`, `j = 1, ..., 6`, which no block of this study has used. Logistic truth,
corruption C1 (a corrupted covariate, the mechanism the default test fails on), 1000 replicates each:

- `fr_clean_n1000`, `fr_clean_n5000` -- no corruption, the level check on fresh seeds.
- `fr_C1_k020_n1000` (2%), `fr_C1_k040_n1000` (4%).
- `fr_C1_k100_n5000` (2%), `fr_C1_k200_n5000` (4%).

Variants: `V0` (the default), `HYB05`, `HYB10`, `HYB20`, `G10`, exactly as block 9R2 defines them.

## 2. Rules

Rejection at p <= 0.05. A replicate with no p-value counts as no rejection and is reported. The level band is
0.05 +/- 3 nominal standard errors (0.029 to 0.071 at B = 1000). The Monte Carlo standard error of a rate
near 0.10 is 0.009, and a reading is called above or below the bar only when it clears the bar by more than
two of those.

## 3. What would confirm the rule, and what would refute it (fixed here)

The rule is confirmed when both hold:

1. **It protects at its own tolerance.** HYB10 reads at or below 0.10 at the 2% scenarios, and HYB20 reads
   at or below 0.10 at the 4% scenarios, at both sample sizes.
2. **The knob is real, not slack.** A setting used beyond its tolerance fails: HYB05 reads above 0.10 at the
   2% scenarios, and HYB10 reads above 0.10 at the 4% scenarios. If every setting protects everywhere, the
   tolerance is not doing the work the rule claims for it and the parameter should be dropped in favour of
   one fixed partition.

Level must hold for every variant on the two clean scenarios; a variant that misses the band there is
reported as miscalibrated and its protection is not quoted.

If (1) holds and (2) fails, the honest reading is that a single wide partition is enough and EDGE-FR needs no
tolerance argument. If (1) fails, the rule is wrong and only the settings that were measured -- one per cent
-- may be recommended.

## 4. Written prediction (recorded before the run)

`5c` was chosen so that the corrupted records are about a fifth of their group, which is what HYB05 achieved
at one per cent (0.111 at k = 50 of 5000). The residual of a group scales as `m_c / sqrt((m - m_c) v_bar)`
and holds that ratio fixed as `c` grows, so protection should hold: HYB10 near 0.11 at 2% and HYB20 near
0.11 at 4%, both at the bar rather than comfortably under it. The second half is the one I am less sure of.
At 4% of n = 5000, two hundred corrupted records carry enough weight to bend the fitted model itself, not
only the group they sit in, and the whole partition family may fail together; if that happens, every variant
including HYB20 reads high and the honest conclusion is that the tolerance cannot be pushed that far.

## 5. Outputs

`simulations/battery/9R3/`: one per-replicate file per scenario with a column per variant, read beside the
table of block 9R2.
