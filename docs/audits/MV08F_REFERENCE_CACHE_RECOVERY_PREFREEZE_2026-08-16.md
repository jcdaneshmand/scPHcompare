# MV8-F reference cache recovery prefreeze

Date: 2026-08-16

## Outcome

MV8-F authorizes reconstruction of the missing 450 accepted primary-reference
SCT caches and nothing downstream. The exact 90-sample by five-seed recovery
axis is frozen, all 90 retained individual sources are uniquely locatable, the
current normalization runtime equals the accepted runtime, the 170 retained
added-sample caches remain exact, and the full accepted 500-gene comparator
stack remains immutable and hash-valid.

The prospective contract is
[`MV08F_REFERENCE_CACHE_RECOVERY_PREFREEZE_V1.md`](../specifications/MV08F_REFERENCE_CACHE_RECOVERY_PREFREEZE_V1.md).
Machine-readable evidence is in
[`mv08f-recovery-prefreeze-evidence`](mv08f-recovery-prefreeze-evidence/).

## Frozen recovery axis

| Measurement | Result |
|---|---:|
| Primary samples | 90 |
| Seeds | 5 (`20260805`–`20260809`) |
| Selected cells per sample and seed | 384 |
| SCT caches to reconstruct | 450 |
| Minimum/maximum post-QC cells | 421 / 9,071 |
| Uniquely located retained sources | 90 / 90 |
| Outcome-label state | Closed |

The public recovery axis is the unmodified `primary90` subset of the accepted
MV7-FP manifest. Private execution inputs contain only sample IDs, expected
cell counts, the historical source digest, and closed-label flags; no source
path is published.

## Identity gates

The current WSL runtime exactly equals the runtime stored in an accepted
MV5-D0 v2 cache: R 4.4.1, Seurat 5.3.0, SeuratObject 5.0.2, sctransform 0.4.1,
Matrix 1.7-3, sequential future execution, and `1/1/1` OMP/OpenBLAS/MKL thread
identity. All 170 retained `added34` caches match their accepted file byte
sizes and SHA-256 values.

The immutable 500-gene comparator was rehashed in full:

| Artifact | Expected/live | Hash identity | Action |
|---|---:|---|---|
| Source bundles | 5 / 5 | Exact | Reuse |
| PH records | 1,240 / 1,240 | Exact | Reuse |
| Complete landscape groups | 20 / 20 | Exact | Reuse |

This means recovery is limited to upstream SCT inputs. The accepted 500-gene
source, PCA, typed views, persistence diagrams, and landscapes will not be
recomputed for the sensitivity comparison.

## Resources and stop boundary

Raw and SCT children retain the accepted 1,800-second and 8-GiB process-tree
caps with at most two workers. Aggregate private cache storage is capped at
40 GiB. Automatic retry is disabled; exact validated outputs may be resumed.

The independent validator passes 10/10 categories: recovery axis, source
coverage, live-cache state, runtime identity, immutable comparator, resource
caps, source freeze, decision boundary, artifact manifest, and privacy.

No 475-gene source or PCA fit, PH, landscape, distance, clustering, biological
label, HCA expression calculation, FASTQ download, or raw-read reprocessing is
authorized by this gate. Only a subsequent 450-of-450 exact comparison against
the accepted manifest may open the 475-gene calculation stage.

## Decision

Proceed with monitored, atomic reconstruction of 90 raw shards and 450 SCT
caches. Stop on any selected-cell, logical cache-key, payload, byte-size,
file-hash, resource, or label-boundary mismatch.
