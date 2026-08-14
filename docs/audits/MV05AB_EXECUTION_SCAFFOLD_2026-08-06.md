# MV5-A/MV5-B execution-scaffold audit

| Field | Value |
|---|---|
| Date | 2026-08-06 |
| Parent contract | `mv05_statistical_benchmark_contract_v1` |
| Split contract | `mv05_loso_manifest_v1` |
| Baseline contract | `mv05_matched_baseline_bundle_v1` |
| Mapping contract | `mv05_inductive_mapping_v1` |
| Outcome-label state | Closed |
| Biological outcomes | Not computed |
| Gate | MV5-A/MV5-B technical gate passes; G-MV5 remains open |

## Outcome

The frozen MV-05 design is now executable through the point immediately before
existing-data biological computation. The sprint produced immutable LOSO fold
manifests, three same-unit baseline implementations, analytical and full-shape
fixtures, and a deterministic reference-to-query mapping route.

## MV5-A: immutable folds

`mv05a-loso-execution-manifest-2026-08-06.csv` contains 2,232 rows: all 124
large-cohort samples repeated across 18 leave-one-study-out folds. Each fold has
one held-out study, a disjoint `training_reference`/`held_out_query` partition,
and a fold-specific training fit scope. Across folds, the held-out counts sum
to 124, so every sample is queried exactly once.

Study is present because it is the prespecified split variable. Tissue and
approach are absent, and every row states `outcome_label_state=closed`. The
manifest cannot be modified without invalidating its cache identity.

## MV5-A: matched baselines

Three baseline contracts now return the same typed, immutable sample-distance
bundle:

1. `cell_distribution_energy_shared_pca_v1` uses the square root of empirical
   V-statistic energy divergence over the same cell-by-PC point clouds used for
   cell topology.
2. `pseudobulk_shared_panel_euclidean_v1` uses Euclidean distance between gene
   means from the same training-scaled source matrices.
3. `gene_correlation_frobenius_v1` uses RMS Frobenius distance between matched
   gene-correlation matrices.

Each implementation validates common fit scope, representation, seed, and
axes, and records input cache keys, matrix SHA-256, formula ID, and bundle cache
key. Analytical fixtures passed. A separate synthetic scientific-shape fixture
used 500 genes, 384 cells per sample, 30 PCs, and three samples. All three
matrices were finite, symmetric, and exactly zero on the diagonal. The shared
PCA fit took 1.209 seconds; baseline computations took 0.066 seconds for cell
energy, 0.058 seconds for pseudobulk, and 0.245 seconds for gene correlation in
the measured run. These timings demonstrate bounded formula feasibility only,
not full-cohort cost.

## MV5-B: inductive mapping

The installed Seurat 5.3.0 APIs support the required architecture:
`FindTransferAnchors()` can project a reference PCA onto a query, and
`IntegrateEmbeddings()` can return corrected query embeddings. The project
wrapper deliberately omits `TransferData()` and UMAP projection because neither
is needed for sample geometry.

The deterministic fixture used:

- 80 reference cells and 40 disjoint held-out query cells;
- 200 features selected from the reference only;
- 10 reference-trained PCs;
- SCT residual recomputation under the reference model;
- exact RANN neighbor search; and
- no tissue, approach, class, or transferred labels.

Both repetitions found 111 anchors. They produced the same query-embedding
SHA-256 (`5a29aae7520ab9480be95db7df7cee23a1dda6f09e4ad8d3c83d0e7e22e8b627`)
and preserved reference identity SHA-256
(`117b6b7eead57d46a57abb77fbf30837f62711036ad92467139bace49d072e22`).

The SCT fixture emitted convergence warnings for the deliberately small
synthetic counts and used the installed native fallback because `glmGamPoi` is
not installed. The mapping nevertheless completed deterministically; no
dependency was installed or lockfile changed. These warnings must be reassessed
on the MV5-C existing-data feasibility run.

Focused MV5-A/MV5-B tests passed 33/33 expectations, including hand-calculated
energy and gene-correlation distances. The complete package suite passed
458/458 expectations.

The source package built in 489.4 seconds. `R CMD check --no-manual` completed
with `Status: OK` and zero errors, warnings, or notes.

Seurat's official reference documentation describes `pcaproject` as projecting
the reference PCA onto an scRNA-seq query and `IntegrateEmbeddings()` as adding
the corrected query embedding:

- <https://satijalab.org/seurat/reference/findtransferanchors>
- <https://satijalab.org/seurat/reference/integrateembeddings>

## Failure and activation boundary

`mv05_try_inductive_mapping_v1()` records any mapping exception as
`held_out_mapping_unavailable`; no all-sample integration substitute is allowed.
The successful synthetic fixture authorizes MV5-C feasibility only. It does not
authorize biological evaluation, a full 90-sample run, fusion, manuscript
claims, or changes to the existing public integration function.

## Next gate

MV5-C should execute one complete eligible existing-data tissue across all five
frozen subsample seeds and all prespecified single-view/baseline methods under
the existing caps. It must write and hash representations, mappings, diagrams,
distances, baseline matrices, and label-free predictions before opening outcome
labels, retain every failure, and produce a full-run resource projection.
