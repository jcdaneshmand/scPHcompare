# MV5-D0 normalization-cache gate audit

| Field | Value |
|---|---|
| Date | 2026-08-07 |
| Contract | `mv05d0_normalization_cache_gate_v1` |
| Planned production scope | 90 samples × 5 seeds = 450 SCT caches |
| Outcome-label state | Closed |
| Downstream folds, PH, landscapes, clustering, endpoints | Not run |
| Decision | Stop before production; source staging and legacy identity gates failed |

## Outcome

MV5-D0 did not authorize the 450-entry production build. Two independent
pre-normalization gates failed before any fold, persistence diagram, landscape,
distance, clustering, or biological endpoint was computed.

First, the historical `expr_list_raw.Rds` is a monolithic compressed R object.
Deserializing it expanded to approximately 47,965,492 KiB RSS per process
(45.74 GiB), far above the frozen 8-GiB heavy-job guard. The initial
uncompressed sharding attempt reached its 1,804-second timeout and its R child
survived the client timeout. A concurrent process audit detected that stale
worker and the corrected worker, both above the guard; exactly those two PIDs
were terminated. No unrelated WSL process was changed.

Second, a six-sample matrix-only SCT pilot was resource-feasible and
reproducible in the current runtime, but did not exactly reproduce the legacy
MV5-C2 cache even though raw counts, selected cells, axes, and recorded Seurat
version were identical. The legacy v1 identity recorded too little of the
normalization runtime, allowing numerically different outputs to share the
same logical identity.

## Source-staging evidence

The rejected uncompressed attempt left 53 complete, ignored raw shards totaling
23,987,530,751 bytes before the worker was stopped. All six samples from the
earlier MV5-C pilot were among those shards. A bounded recovery audit reopened
only those six files, one at a time, verified 2.218 GiB of shard hashes, and
froze six label-closed seed-`20260805` selections. The remaining candidates
were not inferred, filled, or silently dropped.

The metadata's `File.Path` field is empty for all 90 frozen candidates, so an
existing per-sample upstream path is not currently available from the local
metadata. Completing source migration therefore requires either a separately
located per-sample source or an explicitly approved one-time high-memory
conversion outside the current 8-GiB contract. New biological data are not yet
triggered; this is a storage-layout problem in the existing data.

## Six-entry SCT matrix pilot

All six jobs completed with at most two workers and labels closed.

| Measurement | Result |
|---|---:|
| Entries | 6/6 |
| Worker-seconds | 175.43 s |
| Maximum process-tree RSS | 3.361 GiB |
| Matrix-cache storage | 41.08 MiB |
| Legacy full-object storage for the same entries | 351.81 MiB |
| Projected matrix-cache storage for 450 entries | 3.009 GiB |
| Frozen elapsed/RSS/storage caps | 1,800 s / 8 GiB / 40 GB |

The matrix payload retains the complete SCT `data` matrix used by the cached
fold worker and reduces these six files by 88.3%. This storage design remains
eligible, but its production identity must be v2.

## Legacy equivalence and identity correction

Across the six samples, axes matched but exact legacy equality failed. The
maximum absolute difference was `0.693147180559945`; for the diagnostic sample,
6,128 of 4,991,232 values (0.123%) differed. Differences were log-ratios such
as ±log(2). A clean single-worker rerun exactly equaled the two-worker output,
so queue concurrency was not the cause.

The implemented `mv05d0_sample_seed_sct_identity_v2` now includes R version,
all RNG kinds, Seurat/SeuratObject/sctransform/Matrix versions, BLAS/LAPACK
identifiers when reported, and OMP/OpenBLAS thread settings. Two independent
v2 builds produced identical cache keys, payload hashes, file hashes, and byte
sizes. Legacy v1 records remain readable for historical comparison but cannot
be reused as v2 production caches.

## Scientific boundary

The landscape contract is unchanged: H0 and H1 remain separate, every active
consecutive landscape level is retained, and integration remains exact or
error-controlled. No landscape was calculated because the normalization/source
gate is upstream of topology.

No tissue/approach outcome was opened, no method winner was selected, no new
data were added, and no full fold, PH, distance, clustering, fusion, or endpoint
job was run. PDFs, reviewer correspondence, raw shards, SCT caches, logs, and
`example_run.r` remain outside Git.

## Verification and process deviation

The complete test suite passes 500/500 and a clean staged source-package check
reports `Status: OK`. The first staged `R CMD build` copied the project's
`.Rprofile` and auto-bootstrapped `renv 1.1.5` into the WSL cache. This was an
unintended environment-side installation; it did not change the repository,
lockfile, or analysis package set. Subsequent package checking disabled user
and project profile activation and performed no restore. The deviation is
retained here rather than described as satisfying the no-install boundary.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Revise normalization provenance; landscape definition unchanged |
| Correctness demonstrated? | v2 cache reproducibility passes; legacy exact equivalence fails |
| Computation feasible? | Per-entry SCT yes; monolithic source staging no under current guard |
| Biological interpretation permitted? | Prohibited |
| Next action | MV5-D0a source-migration decision, then repeat the 450-cache gate |

## Next action

Do not proceed to the 75 fold/seed jobs. The next sprint must locate an existing
per-sample raw source or request an explicit, one-time high-memory conversion
exception with a separate cap, monitored single worker, compressed output, and
post-conversion hash audit. After all 90 source shards exist, regenerate all
450 selection identities and build only v2 matrix caches before reprojection.
