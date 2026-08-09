# MV5-F label-closed integration-induction resource gate specification v1

## Document control

| Field | Value |
|---|---|
| Contract ID | `mv05f_label_closed_integration_resource_gate_v1` |
| Date | 2026-08-08 |
| Status | Executed and independently validated |
| Existing-data axis | 90 samples, 15 studies, five seeds |
| Mapping axis | Reference-only SCT transfer projection, 500 genes, 384 cells, 30 PCs |
| Outcome state | Closed |
| Full integrated execution | Not authorized by this sprint |

## 1. Purpose and stop boundary

MV5-F determines whether the accepted MV5-B inductive mapping contract can be
reconstructed from the accepted MV5-D0/D1 sources and whether a future full
integrated cell-view execution is feasible. It does not test whether integration
improves tissue retrieval and does not use MV5-E results to alter the method.

This gate permits only source validation, fixed-panel reference/query
reconstruction, held-out mapping, ordered coordinate assembly, resource
measurement, independent validation, deterministic repeat, and projection.
It prohibits persistence homology, landscapes, distances, rankings, clustering,
gene topology, fusion, new data, label transfer, biological outcomes, and claim
promotion.

## 2. Frozen inherited axes

The future integrated route must preserve:

- the same 90 samples and 15 study-held-out folds;
- seeds `20260805:20260809`;
- the exact D1 training-selected 500-gene panel for each fold-seed;
- the exact D0-selected 384 cells per sample-seed;
- 30 reference PCs;
- separate H0 and H1 roles;
- the raw H0/H1 composite as descriptive only;
- the matched cell-energy and pseudobulk baseline roles; and
- the frozen MV5-E endpoint registry.

MV5-E tissue results cannot select folds, features, component weights, methods,
or exclusions. Tissue labels are not read during MV5-F.

## 3. Corrected real-data reconstruction

The matrix-only D0 caches intentionally omit Seurat SCT models, so MV5-F
reconstructs them from accepted raw shards while using D0 caches only to recover
and verify the selected-cell identities.

For one fold-seed group:

1. Validate the accepted D1 fold-cache hash and cache key.
2. Validate all 90 raw-shard hashes and all 90 seed-matched D0 SCT cache hashes.
3. Recover the exact selected cell IDs from each accepted D0 cache.
4. Subset raw counts to the D1 panel and selected cells.
5. Jointly SCTransform the training samples on the fixed panel.
6. Fit an exact 30-PC reference PCA on that fixed panel.
7. SCTransform each held-out sample independently.
8. Map each held-out sample with SCT transfer projection and integrated
   embeddings, without label transfer.
9. Assemble 90 ordered 384-by-30 coordinate matrices.

If a held-out sample lacks a frozen-panel feature, mapping uses the ordered
intersection of the frozen panel and query SCT assay. Missing genes are never
replaced and the full reference PCA remains fixed. Training samples must support
all 500 panel genes.

## 4. Immutability and identity

Every group identity includes the fold, fit scope, seed, training/query sample
axes, complete panel, dimensions, D1 payload identity, all raw-shard hashes, all
D0 SCT cache keys, mapping implementation identity, and runtime versions.

The full-panel reference identity includes cells, PCA embeddings, PCA loadings,
and SCT model count. It must be digest-identical before and after all query
mappings. Each query mapping has an immutable cache key, active feature list,
anchor count, and embedding hash. The private group payload contains coordinates
and mapping provenance only; timing remains outside the scientific payload so
clean rebuilds can be byte identical.

## 5. Outcome-blind pilot selection

Fold selection uses only the complete accepted D1 resource table across five
seeds. Seed `20260805` is executed for four structural roles:

| Role | Study | Training/query samples |
|---|---|---:|
| Minimum query / maximum reference | `SRA713577` | 89 / 1 |
| Maximum query / minimum reference | `SRA779509` | 65 / 25 |
| Maximum aggregate missing-feature burden | `SRA826293` | 82 / 8 |
| Median query size with nonzero missingness | `SRA749327` | 86 / 4 |

Ties are resolved by radix-sorted study ID. The pilot seed, roles, and selection
algorithm are frozen before execution.

## 6. Resource guards

The pilot uses one heavy worker with 1,800 seconds and 8 GiB process-tree RSS
per group, 7,200 seconds for the stage, and 10 GB total private result storage.

The full projection uses a 25% reserve, a 21.6-worker-hour planning cap, the
same 8-GiB RSS and 1,800-second group caps, and a 10-GB storage cap. Mapping is
projected from measured fixed-input, per-training-sample reference, per-query
SCT, per-query mapping, and assembly rates over the exact 75-group D1 axis.
Accepted D3/D4/D5 workloads are carried forward. D3/D4 receive the conservative
maximum of one and the MV5-C integrated-to-SCT finite-interval ratio.

## 7. Required validation

An independent validator must reconstruct without calling the production group
validator: payload and group identities; all ordered coordinate matrices; every
mapping and embedding identity; fixed-panel active subsets; reference
immutability; private file hashes; and absence of prohibited payloads.

Two clean admission-group rebuilds under the same implementation must be byte
identical. Focused and complete package tests and a source-package check must
pass without dependency or lockfile changes.

## 8. Decision rule

The decision is `go` only if the reserved complete projection passes all caps
and every pilot group passes identity, mapping, and scope checks. A `go` means
only that the owner may authorize a separate full label-closed integrated
execution goal. MV5-F records `full_execution_authorized=FALSE` and stops.
