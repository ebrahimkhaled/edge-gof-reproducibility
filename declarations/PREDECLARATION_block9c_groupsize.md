# Pre-declaration: block 9c, what governs EDGE's tolerance of corrupted predictions (written before any 9c cell is run)

Date: 2026-09-19. Written after block 9 finished and B9.1 **failed**, and before any block 9c data exists.

## 0. Why this block exists, stated plainly

Block 9's claim B9.1 required EDGE's unit forms to stay within 3 nominal standard errors of 0.05 in every C1
cell at 0.1-0.5% corrupted covariates. It failed in one cell of six: n = 5000 at 0.5% (k = 25), where
EDGE-poly3 reached 0.151 and EDGE-sym 0.099. That verdict stands and is reported as a failure; nothing in this
document revises it.

What block 9 also showed, in data already stored, is that the same k at the same n behaves completely
differently under a different partition:

| n | k | G = 10 (group size n/10) | G = rule (group size 25) |
|---|---|---|---|
| 5000 | 25 | 0.065 | 0.151 |
| 5000 | 50 | 0.107 | 0.740 |
| 1000 | 10 | 0.085 | 0.097 |

Same data, same corrupted records, same k: only the number of groups differs. So the driver is not the
contamination **rate** k/n. Block 9c asks what it is instead. The answer is not assumed here -- two
incompatible answers are written down below, and the design separates them.

**Discipline note.** The threshold "k at most 10" was proposed by the author after seeing block 9. Adopting it
as a claim would be choosing a boundary from the data, which this project has refused before (the Stukel
generator placement, 2026-09-16). Block 9c instead derives a mechanism, states what each possible outcome
would mean, and tests it on fresh seeds disjoint from every earlier run.

## 1. The mechanism, and the two hypotheses it leaves open

A C1 corruption multiplies x by 4, so the record's linear predictor reaches about +-7 and its predicted risk
goes to 0 or 1. It therefore sorts into an extreme risk group. Its outcome was drawn at the ORIGINAL x, so
relative to that extreme prediction it is close to a coin flip. In the group residual

    r_g = (O_g - E_g) / sqrt(V_g),   V_g = sum_i phat_i (1 - phat_i)

such a record contributes about +-1 to the numerator and almost nothing to the denominator, because
phat(1-phat) -> 0 as phat -> 1. With m_c corrupted records among m, and vbar the typical clean variance,

    r_g  ~  m_c / sqrt((m - m_c) vbar)  =  f sqrt(m) / sqrt((1 - f) vbar),    f = m_c / m.

The corrupted records are in the numerator; the **clean groupmates are the whole denominator**. This yields two
readings, and they disagree:

- **H-f (fraction only).** What matters is f, the fraction of a group that is corrupted. Cells with equal f
  behave alike whatever n, G and k are.
- **H-fm (fraction and group size).** The displacement grows like f sqrt(m) / sqrt(1-f), so at equal f a
  larger group is **worse**, not neutral.

These are separated below. A third possibility is recorded in advance: because EDGE's statistic is a
projection of the whole residual vector and its null reference depends on G, neither may describe the
rejection rate well, in which case the reported quantity is the measured surface and no functional form is
claimed.

## 2. Design

**Data.** The battery base design, logistic truth: x ~ U(-3, 3), d ~ Bernoulli(0.5), eta = 0.6x + 0.5d, fitted
as y ~ x + d. The bulk model is therefore CORRECT, and every rejection is a false alarm.

**Corruption.** C1 of block 9, unchanged: k records chosen at random have x multiplied by 4, their outcome
still drawn from the truth at the original x.

**Cells (10).** n in {1000, 5000} x k in {0, 5, 10, 25, 50}. k = 0 is the clean control.

**Partitions.** Each replicate is tested at every G below in the same fit, so G is manipulated with the data
held fixed -- the comparison across G is paired to the replicate:
- n = 1000: G = 10, 20, 25, 40 (rule). Group sizes 100, 50, 40, 25.
- n = 5000: G = 10, 25, 50, 100, 200 (rule). Group sizes 500, 200, 100, 50, 25.

This gives f = k / (n/G) from 0.01 to 2.0, with the SAME f reached by different (n, G, k) triples, which is
what separates H-f from H-fm.

**Replicates.** B = 1000 per cell. **Tests.** EDGE-poly3 and EDGE-sym, unit and score, at every G above; and
Stukel's joint score, which uses no partition, as the fixed point of comparison.

**Seeds.** set.seed(seed_base + rep), seed_base = 500000000 + cell_id * 10000. Disjoint from blocks 0-8, from
block 9 (300000000) and from block 9b (400000000).

**Stored per replicate.** Every p-value above, the fitted coefficients, and for each G the number of corrupted
records that land in the highest-risk group.

## 3. Rules

1. Every test rejects at p <= 0.05. No Monte Carlo rival takes part, so there is no grid rule here.
2. A test that gives no p-value counts as no rejection, and its rate is reported.
3. The nominal standard error is sqrt(0.05 x 0.95 / 1000) = 0.0069; "agree" below means within 3 of these,
   i.e. 0.021, and comparisons across G inside one cell are paired on the replicate.
4. Claims are judged on EDGE-poly3's unit form. The other three forms are reported without a claim.
5. The clean cells (k = 0) are a control: if EDGE does not hold its size there at some G, that G is excluded
   from the claims below and the exclusion is reported.

## 4. Claims, and what would count against them

- **C9c.1 (the partition governs, at fixed data).** At n = 5000 with k = 25 held fixed, EDGE-poly3's unit
  false-alarm rate at G = 200 exceeds its rate at G = 10 by at least 0.05. *Against:* a gap below 0.05.
  (Block 9's stored values were 0.151 and 0.065, a gap of 0.086, on different seeds.)

- **C9c.2 (the rate k/n does not govern).** Among the cells sharing k/n = 0.005 -- (n = 1000, k = 5) and
  (n = 5000, k = 25) -- the rule-G rates differ by more than 3 nominal SE. *Against:* they agree.

- **C9c.3 (monotone in the partition).** At each n and each k > 0, EDGE-poly3's unit rate is non-decreasing in
  G, allowing a fall of at most 3 nominal SE between neighbouring G. *Against:* any larger fall.

- **C9c.4 (which of H-f and H-fm).** Decided, not assumed, on the pairs of cells that share f but differ in
  group size m:
  - f = 0.25: (n = 1000, k = 25, G = 10, m = 100) against (n = 5000, k = 25, G = 50, m = 100) -- same f, same m,
    different n: a check that n itself does nothing.
  - f = 0.5: (n = 1000, k = 25, G = 20, m = 50) against (n = 5000, k = 50, G = 50, m = 100) -- same f, m
    doubled.
  - f = 1.0: (n = 1000, k = 25, G = 40, m = 25) against (n = 5000, k = 50, G = 100, m = 50) -- same f, m
    doubled.
  **Reading fixed in advance.** If the members of each pair agree within 3 nominal SE, H-f is supported and the
  paper says tolerance is governed by the corrupted fraction of a group. If the larger-m member is higher by
  more than 3 nominal SE in both of the last two pairs, H-fm is supported and the paper says tolerance falls as
  f sqrt(m). If the pattern is mixed, neither is claimed and the measured surface is reported as it stands.

- **C9c.5 (the control).** In the k = 0 cells EDGE-poly3's unit rate is within 3 nominal SE of 0.05 at every G
  used. *Against:* any G where it is not -- that G is then dropped from C9c.1-C9c.4, and this is reported.

- **C9c.6 (reported, no claim).** The same quantities for EDGE-sym, for both score forms, and for Stukel's
  joint score; the mean number of corrupted records in the top group at each G; and the fitted coefficients.

## 5. What this can and cannot rescue

Block 9's B9.1 failed and stays failed. If C9c.1 and C9c.3 hold, the paper does not re-assert the old claim;
it states a different and more useful one: **EDGE's tolerance of corrupted predictions is set by how many of
them share a group, which the analyst controls through G, and the same lever that buys power (section 6.7)
costs contamination tolerance.** That trade-off is reported whichever way C9c.4 falls.

If C9c.1 fails, the group-size explanation of block 9's failure is wrong, and the paper reports B9.1's failure
with no mechanism attached.

## 6. Outputs

`simulations/battery/9c/`, one per-replicate file per cell, and an analysis restricted to block 9c.
