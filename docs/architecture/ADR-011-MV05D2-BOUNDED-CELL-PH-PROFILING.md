# ADR-011: MV5-D2 bounded cell-PH profiling

## Status

Accepted 2026-08-07. The bounded PH profile is complete. Full 6,750-view PH
production and every downstream landscape or biological stage remain
unauthorized pending a separate sprint decision.

## Context

MV5-D1 created 6,750 leakage-safe cell-coordinate views but left production PH
cost unmeasured. Prior PH timings came from different preprocessing paths and
could not establish correctness or resource use for the new fold-fitted
coordinates. The next decision also had to preserve the dissertation-aligned
all-level landscape definition without constructing landscapes prematurely.

## Decision

1. Use a deterministic 30-job pilot containing held-out and training views from
   all 15 folds, six jobs per seed, all nine mapped folds, and six unmapped
   held-out controls.
2. Compute complete Euclidean Vietoris-Rips H0/H1 through the typed corrected
   entry point; retain one explicit essential H0 interval.
3. Require the 383 finite H0 deaths to equal the full-view Euclidean MST edges.
4. Require five independent reruns to be byte-identical and five reduced
   Ripserr/GUDHI comparisons to agree in both H0 and H1 after explicit
   essential-class normalization.
5. Keep private diagrams immutable, identity-bound, resumable, and ignored by
   Git; publish only public-safe manifests, hashes, counts, and resources.
6. Accept the observed median, P90, and maximum projections as bounded planning
   scenarios, while treating them as projections rather than completed work.
7. Stop before full PH, landscapes, distances, clustering, integration, gene
   views, labels, or outcomes.

## Evidence

- 30/30 pilot jobs completed with zero failures;
- 11,520 H0 intervals and 10,322 H1 intervals retained;
- 30/30 full-view H0 MST checks passed with zero recorded maximum error;
- 5/5 repeated records and serialized files were byte-identical;
- 10/10 reduced Ripserr/GUDHI H0/H1 comparisons passed with zero recorded
  maximum error after removing GUDHI's capped essential H0 record;
- 110.554 worker-seconds, 5.433 seconds maximum per job, 239,341,568 bytes
  maximum process-tree RSS, and 907,330 bytes of pilot cache storage; and
- a projected 6.752/7.135/10.187 worker-hours and
  205/231/254 MB for median/P90/maximum assumptions.

Combining the measured normalization and coordinate stages with the previous
landscape projection gives 15.262, 15.645, and 18.697 total worker-hours. All
three remain below the 21.6-hour planning cap, with the conservative observed-
maximum scenario retaining 2.903 hours of margin.

## Consequences

The production cell-PH term is no longer unmeasured, and the SCT cell-primary
pipeline has a complete bounded feasibility projection. The next eligible
sprint is a separately authorized immutable full cell-PH cache stage. It must
still stop before landscapes. Gene topology, integration, clustering, outcome
labels, and G-MV5 remain closed.

