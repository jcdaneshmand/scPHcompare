# ADR-024: use training-only PAM with immutable held-out assignment

- **Date:** 2026-08-10
- **Status:** accepted; MV5-N real admission passed
- **Decision scope:** corrected cell-topology sample clustering

## Context

The original project clustered complete sample-distance matrices and sometimes
used the known number of classes. The corrected production artifacts instead
implement leave-one-study-out transforms and currently contain query-to-training
distances only. A scientifically valid clustering extension must preserve this
inductive boundary and must not use tissue, study, or approach outcomes to pick
`k` or a method.

## Decision

Fit PAM independently on each fold/seed training matrix. Select one `k` per
fold, representation, and distance method from `2:min(10,n-1)` using five-seed
partition ARI stability and the frozen one-SE rule. Canonicalize cluster labels
by sorted training-member signatures. Assign held-out samples to the nearest
frozen training medoid, resolving exact ties by canonical medoid sample ID.

Retain average linkage as a sensitivity at the PAM-selected `k`. Assign a
held-out sample to the cluster with minimum mean distance to frozen training
members, with canonical member-signature tie-breaking.

Use H0 and H1 separately under the exact all-consecutive-active-level landscape
definition. Retain the raw H0/H1 Euclidean composite only as descriptive and
apply the same clustering machinery to matched energy and pseudobulk baselines.

Instantiate and validate all pair identities, publish deterministic
global/group/chunk counts and identity-set digests, and discard non-admitted
wide rows until production is authorized. Profile only an outcome-closed,
size-selected bounded admission before authorizing full production.

## Rejected alternatives

- Complete 90-sample transductive clustering: leaks held-out samples into model
  construction and conflicts with the accepted LOSO benchmark.
- Oracle class-count `k`: outcome-informed and historical sensitivity only.
- Direct k-means on distances: treats dissimilarities as feature coordinates.
- Ward on arbitrary distances: requires an unproven Euclidean geometry.
- Current spectral implementation: lacks a frozen affinity and out-of-sample
  extension, so remains ineligible pending a separate gate.
- Fixed landscape level cap or uniform grid: conflicts with the corrected
  dissertation-aligned all-active-level exact definition.

## Consequences

The design requires 262,675 new training pairs and 525,350 H0/H1 rows per
representation, but no held-out–held-out distances. Cluster outcome evaluation
cannot begin until all predictions are immutable. Failure of the resource gate
triggers an equivalence-preserving optimization decision rather than scientific
truncation or outcome-guided narrowing.
