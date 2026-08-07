# ADR-008: MV5-D0 source and normalization-identity gate

## Status

Superseded by ADR-009 on 2026-08-07 after a compliant existing per-sample
source route was found. This document remains the audit record for the rejected
monolithic-source route and the v2 identity requirement.

## Context

MV5-C2 authorized label-closed SCT cell-primary cache production under
30-minute/8-GiB per-entry and 40-GB storage guards. Actual inspection showed
that deserializing the monolithic raw-count RDS requires about 45.74 GiB RSS.
It also showed that rerunning SCT from identical inputs can differ from legacy
cache values while the legacy v1 key records only Seurat version.

## Decision

1. Reject monolithic source deserialization as a compliant production input.
2. Require one immutable raw shard per sample before normalization admission.
3. Allow an explicit, separately governed one-time migration exception only if
   no existing per-sample source can be located.
4. Replace v1 production identities with runtime-complete v2 identities.
5. Store the exact SCT data matrix instead of a full Seurat object.
6. Retain legacy caches for audit comparison, never as v2 cache hits.
7. Preserve the exact all-active-level, separate-H0/H1 landscape contract.
8. Keep labels closed and prohibit all downstream analysis until 450 v2 caches
   pass resource and completeness checks.

## Consequences

Per-entry SCT is feasible (six jobs, 3.361 GiB maximum RSS), and projected
matrix-cache storage falls to about 3.009 GiB. Production is nevertheless a
no-go because 37 of 90 sample shards were not recovered and source conversion
violates the current guard. This is an infrastructure/provenance limitation,
not evidence for or against the biological method and not yet a trigger for
new biological data.
