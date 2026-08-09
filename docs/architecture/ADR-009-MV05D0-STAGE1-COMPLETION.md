# ADR-009: MV5-D0 stage-1 completion and post-cache boundary

## Status

Accepted 2026-08-07 and superseded at the fold boundary by ADR-010. MV5-D0
stage 1 remains the accepted normalization input bundle.

## Context

ADR-008 stopped the monolithic raw-list route because deserialization required
about 45.74 GiB RSS. A subsequent local source search found exactly one existing
per-sample Seurat RDS for every frozen candidate. These files provide the same
RNA count matrices without loading the monolith.

## Decision

1. Use the existing per-sample Seurat RDS files as the source-migration route.
2. Canonicalize only RNA counts and sample-prefixed cell identifiers; do not use
   source metadata as an endpoint or selection variable.
3. Retain canonical sparse count hashes and runtime-complete v2 SCT identities.
4. Accept the 90 raw shards, 450 deterministic selections, and 450 matrix-only
   SCT caches as the complete MV5-D0 stage-1 input bundle.
5. Preserve exact all-active-level, separate-H0/H1 landscapes unchanged.
6. Stop after post-cache resource projection. Do not automatically launch the
   75 fold/seed jobs.
7. Keep gene topology, integration, full-matrix clustering, fusion, and all
   outcome evaluation behind their existing independent gates.

## Evidence

- 90/90 raw shards, zero failures, 53/53 exact comparisons against recovered
  monolithic samples, and 37 newly recovered existing-data samples;
- 450/450 deterministic selection summaries;
- 450/450 runtime-complete v2 SCT caches, zero failures;
- 2.562 normalization worker-hours, 1.81 GiB maximum RSS, 37.516-second maximum
  entry, and 2.992 GB total cache storage;
- independent validation of every record and file hash; and
- updated SCT cell-primary lower bound of 10.525 worker-hours.

## Consequences

The source-layout blocker is resolved without new data or a high-memory
exception. MV5-D0 stage 1 is complete, but G-MV5 remains open because no fold,
topology, distance, clustering, or biological endpoint was calculated. The next
decision is whether to authorize a narrowly specified, label-closed SCT
cell-primary fold stage under the measured resource envelope.
