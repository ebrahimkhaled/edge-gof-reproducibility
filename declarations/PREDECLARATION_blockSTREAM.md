# Pre-declaration: block STREAM, EDGE updated as patients arrive

Date: 2026-10-01. Written and hashed before any replicate was computed.

## Question
`edge.stream()` (ebrahim.gof 2.9.0) keeps four sums per risk group and updates them record by record, with
cut points fixed in advance from reference predictions. Does the streamed test keep its level when the
risk distribution of later patients drifts (so groups are no longer of equal size), what do repeated looks
cost, and how long does an update take?

## Design (frozen model, external mode, cubic basis plus constant, ten groups)
- Reference predictions: 10,000 draws of logit p ~ N(-1.5, 1); cut points = their deciles.
- Stream: batches of 100 patients, outcomes y ~ Bernoulli(p_true).
- Cells, 2000 replicates each, seed 9100000 + 10000 * cell + rep:
  1. null, no drift: logit p ~ N(-1.5, 1), p_true = p; read at n = 10,000 and n = 100,000.
  2. null, drift: logit p ~ N(-0.5, 1.3) (groups become unequal), p_true = p; read at n = 10,000 and 100,000.
  3. repeated looks: cell-1 data; test after each of ten blocks of 1000 patients (n = 1000, ..., 10,000);
     report P(at least one rejection at 5%) and the same with each look at 0.5% (Bonferroni over ten looks).
  4. alternative: overfitted model, p_true = plogis(0.8 * logit p), logit p ~ N(-1.5, 1); power at n = 2000 and
     n = 10,000.
- Timing: median time of `update()` for a batch of 100 and of `summary()`, with 1,000,000 records already
  streamed; 200 repetitions.

## Reading rule
Level holds if within 0.05 +/- 3 SE (0.035 to 0.065 at 2000 replicates). Reported as measured, whatever it is.

## Outputs
simulations/battery/STREAM/_summary.csv; runner simulations/run_M_blockSTREAM.R.
