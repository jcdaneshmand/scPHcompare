# MV5-D0 normalization-cache specification v2

## Status

Accepted and completed on 2026-08-07. Production used the existing per-sample
Seurat sources, produced all 90 raw shards and 450 v2 SCT caches, and stopped
for the required post-cache feasibility projection.

## Frozen scientific input

- 90 existing-data samples across the 15 already-frozen studies;
- five seeds, `20260805` through `20260809`;
- exactly 384 deterministically selected cells per sample and seed;
- `Seurat::SCTransform(variable.features.n = 3000,
  return.only.var.genes = FALSE)`; and
- outcome labels closed through normalization, fold construction, topology,
  landscapes, distances, and label-free technical summaries.

## Source contract

Production workers must read one immutable sample shard, not deserialize the
monolithic historical list. Each shard records the sample ID, historical-source
SHA-256, counts, closed-label state, and false biological-outcome flag. Shards
must be atomically published and audited by file hash. A one-time conversion
that exceeds 8 GiB requires a separate owner-approved exception and may not be
misreported as satisfying the normal per-job cap.

## Cache identity v2

The key includes sample, seed, selected-cell hash, historical-source hash,
normalization settings, R/RNG identity, package versions, BLAS/LAPACK identity,
and thread settings. A legacy v1 cache is historical evidence only and is never
a v2 cache hit.

## Payload

The cache stores the complete sparse SCT `data` matrix and selected original
cell IDs. It does not retain the full Seurat object. Downstream feature ranking,
standardization, PCA, and mapping remain fit inside the training fold.

## Execution and failure policy

- at most two SCT workers;
- 1,800 elapsed seconds and 8 GiB process-tree RSS per entry;
- 40,000,000,000 bytes total SCT-cache storage;
- stop admitting new jobs after the first failure or observed cap breach;
- never overwrite stale, corrupt, partial, or identity-incompatible files;
- record every completed entry, failure, resource measurement, and hash; and
- stop after 450 caches for a new fold/runtime/storage projection.

## Prohibited work

Until the cache gate passes: no fold execution, PH, landscapes, distances,
clustering, integration, gene-view expansion, fusion, endpoints, label joins,
new-data acquisition, dependency installation, or default change.

The dissertation-aligned landscape definition remains separate H0/H1,
all active consecutive levels, and exact or error-controlled integration.
