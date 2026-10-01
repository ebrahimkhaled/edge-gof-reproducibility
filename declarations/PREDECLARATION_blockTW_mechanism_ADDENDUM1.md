# Addendum 1 to PREDECLARATION_blockTW_mechanism.md (sha256 80b5a83f...)

Date: 2026-10-01, after the smoke run (3 replicates a cell, a code check) and before any replicate of the full run.

The smoke run showed that the declared diagnostics rule (drop a record when DeltaX2 > 4 or DeltaBeta > 1) removes
about 3.5% of the records of a clean data set (35 of 1000) and that every test then rejected a correct model in all
three smoke replicates. DeltaX2 > 4 is a threshold for inspecting records, not for deleting them, so the declared rule
is kept and reported as declared, and one variant is added:

- PRG2.*: drop only records with Pregibon's DeltaBeta_i = rP_i^2 h_i / (1 - h_i)^2 > 1 (the influence criterion),
  refit, then Stukel's joint score, the cubic LR, EDGE at the default partition and at ten groups. The numbers dropped
  are recorded as PRG2.dropped and PRG2.dropped_corrupt.

Nothing else changes: cells, seeds, the other arms and the reading rules are as declared. Both rules are reported.
