# MV5-AU corrected landscape matrix consumer prefreeze specification v1

Date frozen: 2026-08-12

Authorization: MV5-AT completion `9326726`

## Decision

The first corrected-matrix consumer shall be additive, default-off,
label-free average-linkage dendrogram construction, performed independently
for H0 and H1 within one explicitly bound topology view and iteration.

The implementation must first provide a verified read-only loader for a
completed `scph_corrected_landscape_artifact_v1` sidecar. The loader returns a
new versioned consumer bundle; it must never place corrected matrices in a
legacy field or accept a bare matrix without completion, provenance, cache,
sample-axis, dimension, iteration, and view validation.

## Why this consumer is first

Average-linkage `hclust(as.dist(matrix), method = "average")` needs no label,
outcome, cluster count, kernel scale, or stochastic seed. It retains the full
hierarchy and works directly with a validated dissimilarity matrix. H0 and H1
produce separate trees. The descriptive combined matrix is retained for
inspection but cannot drive the primary consumer.

PAM and spectral clustering remain closed because both require an explicit
cluster count; spectral additionally requires a kernel-scale rule. Ward.D2
and other linkage methods remain closed pending a comparison contract. Betti,
Euler, aggregated-landscape, and cross-iteration curve code consume sampled
functions rather than distance matrices and are schema-incompatible.

## Frozen interfaces

The later implementation may add two new public APIs:

1. `read_corrected_landscape_bundle(sidecar, iteration_id, view_id)`;
2. `corrected_landscape_average_trees(bundle)`.

`sidecar` must identify one completed artifact directory and its expected
input-set and matrix cache keys. `iteration_id` is one non-empty identifier.
`view_id` is exactly `cell_topology_v1` or `gene_topology_v1` and must be
explicitly supplied; it may not be guessed from filenames or sample counts.

The loaded bundle contains immutable H0, H1, and descriptive combined matrices,
sample IDs, input manifest identities, iteration/view IDs, scientific contract,
source cache keys, and a derived cache key. The tree result contains exactly
H0 and H1 average-linkage trees plus provenance. It contains no partition,
selected `k`, label, outcome, combined tree, or legacy alias.

## Validation and rollback

Implementation must test completion and shard hashes, class/schema rejection,
sidecar-key mismatch, view/iteration rejection, reordered/duplicated axes,
nonfinite/asymmetric/nonzero-diagonal matrices, separate H0/H1 identity,
permutation equivariance, deterministic tree serialization, unchanged source
artifacts, default-off equivalence, and absence of legacy/downstream writes.

Rollback is omission of the new explicit consumer control/API call. Existing
legacy clustering, visualization, Betti, modular, and cross-iteration paths
remain byte-for-byte unchanged.

## Abort boundary

Abort for scientific-contract drift, incomplete/corrupt artifacts, implicit
view inference, H0/H1 fusion, use of combined as primary, label/outcome access,
implicit `k`, legacy redirection, default change, artifact mutation, test/check
failure, or resource evidence outside the accepted MV5-AT envelope.

Passing this prefreeze authorizes only implementation and bounded label-free
smoke of the loader and two average-linkage trees. It does not authorize
partitions, evaluation, visualization, biological claims, new data, or
optimization/Rust.
