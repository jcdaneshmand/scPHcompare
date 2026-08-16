# MV8-F reference cache recovery prefreeze v1

Date: 2026-08-16

## Purpose

MV8-F reconstructs only the 450 accepted primary-reference SCT caches that are
needed to calculate the prospective 475-gene harmonization sensitivity. It does
not refit the accepted 500-gene analysis, download HCA FASTQs, inspect outcome
labels, or compute a persistence diagram, landscape, clustering, or endpoint.

## Immutable inputs

The reconstruction axis is the `primary90` subset of the accepted MV7-FP cache
manifest: 90 samples, five seeds (`20260805` through `20260809`), 384 selected
cells per sample and seed, and 450 exact expected cache identities. Every
sample must have exactly one retained individual Seurat source and must remain
in the corrected primary-90 population.

The current WSL normalization runtime must exactly equal the runtime embedded
in an accepted MV5-D0 v2 cache, including R, Seurat, SeuratObject,
sctransform, Matrix, RNG, BLAS/LAPACK, future-plan, and thread identities. The
170 retained `added34` caches must match their accepted byte sizes and SHA-256
digests before recovery begins.

The accepted 500-gene comparator remains immutable. Its five source bundles,
1,240 PH records, and 20 complete landscape groups must match the public
MV7-H ledgers in file count, byte size, and SHA-256. A mismatch stops the stage;
the comparator is not regenerated implicitly.

## Authorized reconstruction

The accepted MV5-D0 v2 pipeline is reused without scientific alteration:

- exact per-sample RNA count extraction from the retained Seurat sources;
- deterministic 384-cell selections for the five accepted seeds;
- `SCTransform(variable.features.n = 3000, return.only.var.genes = FALSE)`;
- sequential normalization inside each child process;
- at most two monitored heavy children;
- atomic, resumable private raw and SCT cache writes; and
- outcome labels closed throughout.

Every completed cache must equal its accepted MV7-FP manifest row in selected
cell SHA-256, logical normalization cache key, payload contract, payload
SHA-256, file byte size, and file SHA-256. Approximate equality, an updated
runtime identity, or a newly generated replacement manifest is not acceptable.

## Resource and stop policy

Each raw or SCT child has an 1,800-second and 8-GiB process-tree cap. Aggregate
private raw plus SCT cache storage has a 40-GiB cap. No automatic retry is
authorized after a content, identity, resource, or label-boundary failure;
validated atomic outputs may be resumed after review.

The public recovery evidence contains identifiers, hashes, dimensions,
resources, and decisions only. It must not contain expression values, cell
barcodes, absolute private paths, tissue, approach, or biological outcomes.

## Continuation boundary

Only an independent 450-of-450 exact-identity validation may authorize the
prospective 475-gene source/PCA/typed-view stage. Cache recovery itself does not
authorize PH, landscapes, distance comparison, clustering, HCA expression
analysis, FASTQ download, or raw-read reprocessing.
