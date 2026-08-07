# ADR-003: MV-04 topological-distance production path

## Status

Accepted on 2026-08-05 for immutable corrected-pilot distance generation. This
does not activate a package-wide scientific default or authorize clustering,
fusion, or biological interpretation.

## Context

MV-03 produced 56 unique, scientifically eligible first-seed diagrams across
eight cohort-by-representation-by-view strata. The approved landscape contract
is dissertation-aligned: use every active consecutive landscape level, remove
essential intervals before construction, integrate the full squared L2
difference exactly or with declared error control, retain H0 and H1 separately,
and treat their unscaled Euclidean combination as secondary.

The R breakpoint implementation is a correctness oracle, not a production
engine. Persim 0.3.8 provides useful exact critical-pair construction, but its
built-in norm is invalid for tested differences that cross zero. ADR-002
therefore allowed only project-controlled exact integration over those critical
pairs, subject to eligible-diagram batching and validation.

## Decision

1. Use `persim_0.3.8_corrected_critical_pairs_batch_v1` with
   `full_l2_exact_critical_pairs_v1` for the immutable MV-04 pilot bundle.
2. Build each diagram-and-dimension landscape once per batch and integrate each
   signed piecewise-linear difference with the analytical segment identity.
3. Retain separate H0 and H1 matrices. Preserve the raw Euclidean combination
   only as a secondary dissertation-aligned output; do not call it balanced.
4. Provide fit-scope-bound median off-diagonal normalization separately for
   later leakage-safe comparisons. Fitted scales must carry their sample IDs,
   scope, input hash, and cache key, and degenerate components must fail.
5. Retain bottleneck only for matrix groups completed under the declared cap
   and technically exclude unfinished large gene-H1 groups. Retain
   Wasserstein-p=1 only as a completed representative feasibility panel and
   technically exclude full pilot matrices based on measured cost.
6. Keep Rust rejected. Eligible-diagram profiling identifies exact critical-pair
   construction and object growth as the dominant, geometry-dependent hotspot,
   but one production repetition does not establish uncertainty or a material
   end-to-end benefit, and existing-library/algorithmic options and cross-platform
   maintenance remain unresolved.

## Activation boundary

The decision authorizes reuse of the immutable MV-04 technical distance bundle
as input to the separately gated MV-05 benchmark. It does not authorize:

- a public default change;
- H0/H1 fusion weights or cell/gene multiview fusion;
- clustering selection or label-informed transformations;
- biological claims; or
- a Rust implementation.

## Consequences

The exact primary calculation is feasible at the frozen pilot scale but is not
yet an efficient large-study implementation. Gene-view H1 construction is the
clear optimization target. Future work should first reduce construction/object
overhead or evaluate a mathematically equivalent established implementation,
then repeat end-to-end profiles before reopening the Rust gate.

Detailed measurements and gate evidence are in
`docs/audits/MV04_TOPOLOGICAL_DISTANCE_VALIDATION_2026-08-05.md`.
