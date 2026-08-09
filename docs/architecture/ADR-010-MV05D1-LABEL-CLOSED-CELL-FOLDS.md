# ADR-010: MV5-D1 label-closed cell-fold coordinates

## Status

Accepted 2026-08-07. MV5-D1 is complete. Production persistence homology and
all biological outcomes remain unauthorized.

## Context

ADR-009 accepted 450 runtime-complete SCT caches and required a separate review
before fold execution. The existing panel helper ranked variances from training
samples but formed its feature universe from all matrices, allowing held-out
feature availability to affect eligibility. Strictly training-derived panels
also exposed small cross-study schema gaps in held-out matrices.

## Decision

1. Derive the 500-gene feature universe, variance ranks, centers, scales, and
   30-PC model strictly from training studies in each LOSO fold.
2. Project both training and held-out cells through the one frozen model.
3. Map an absent held-out selected feature to the training mean, or zero after
   training standardization. Do not alter the panel using held-out availability.
4. Restrict this missing-feature rule to within-sample cell topology. It does
   not authorize gene topology, integration, or cross-sample coordinate
   baselines.
5. Accept 75 private fold-seed caches containing 6,750 typed 384-cell × 30-PC
   views as the complete MV5-D1 bundle.
6. Preserve separate H0/H1, all-active-level, exact/error-controlled landscapes.
7. Stop before PH and treat production cell-PH cost as unmeasured. The 8.5097-h
   known-components value is a lower bound, not a complete feasibility result.

## Evidence

- 529 source tests after the missing-feature rule and projection correction;
- four corrected real-data pilot folds, including both earlier incompatibilities;
- one independent repeat with zero maximum coordinate difference and a
  byte-identical serialized cache;
- 75/75 production entries and 6,750/6,750 independently validated views;
- 2.3756 worker-hours, 1.95 GiB maximum RSS, and 0.895 GB storage; and
- zero PH, landscape, distance, clustering, integration, gene-view, or outcome
  jobs.

## Consequences

The cells-as-observations coordinate layer is now complete and leakage-audited.
The fold caches can be the sole input to a later cell-primary PH stage. G-MV5
remains open, and no biological claim can be evaluated until predictions and
matrices are frozen under the later endpoint firewall.
