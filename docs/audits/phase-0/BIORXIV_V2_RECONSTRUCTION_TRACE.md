# P0-05 BioRxiv v2 Reconstruction Trace

## Objective

Map the preprint's reported analysis to exact repository inputs, configurations, intermediate objects, result tables, and figures. This is an evolving evidence ledger; unknowns remain explicit.

## Current disposition

`in_progress — exact run not yet reproducible`.

The 127-to-124 main-cohort and validation-cohort counts are now reconciled from the historical code, logs, manifests, and final PH artifacts. The checked-in matrix paths remain unavailable, the custom-input paths are placeholders, and no complete run command/configuration manifest has been found.

## Source-level expectations

The owner-held BioRxiv v2 manuscript and dissertation describe:

- a main multi-tissue analysis of 124 samples across eight tissues;
- a separate GSE120221 bone-marrow validation analysis of 25 samples;
- four principal configurations: Raw, SCT Individual, SCT Whole, and Integrated;
- H0/H1 persistent-homology summaries, including Betti/Euler/landscape analyses and clustering comparisons.

These expectations are used only to test repository correspondence; they do not prove which local files generated the paper.

## Candidate input metadata

`inst/extdata/inputs/metadata_MultiTissueAnalysis.csv` was added by commit `d4fb7e2d129cedbf0d2d8e00c48a97a060bfb2ca` (`analysis inputs`) on 2025-12-11.

| Property | Observed |
|---|---|
| Rows / unique sample names | 127 / 127 |
| Columns | `File Path`, `Sample Name`, `SRA`, `Tissue`, `Approach` |
| Tissues | 8 |
| SRA groups | 18 |
| Approaches | scRNA-seq and snRNA-seq |
| Matrix extension | 127 `.RData` paths |
| Paths resolving on this machine | 0 of 127 |

### Tissue counts

| Tissue | Rows |
|---|---:|
| bone marrow | 31 |
| colon | 13 |
| liver | 6 |
| pancreatic islets | 12 |
| pbmc | 12 |
| prostate | 16 |
| substantia nigra | 9 |
| testis | 28 |

### Count reconciliation

Historical `PH Pipeline V5.R` and current `R/PH_Calculation.R` both apply the same explicit rule after cell-level QC: discard a sample when its post-QC cell count is below `MIN_CELLS = 250`.

Historical `filtered_cells.csv` contains 124 samples. The three rows present in the 127-row public metadata but absent from that retained set are `SRS4386092`, `SRS4386091`, and `SRS4386107`. Dated pipeline logs identify them as `SRA850958_SRS4386092`, `SRA850958_SRS4386091`, and `SRA850958_SRS4386107`, with 169, 197, and 166 post-QC cells respectively. Each was explicitly skipped before PH. Final Raw, SCT Individual, SCT Whole, and Integrated PH lists each contain 124 named, non-null outputs.

The heterogeneous 127-row metadata contains 31 bone-marrow rows. Historical
`data/GSE120221/Sparse_RData` contains 25 separately processed validation
inputs; its pre- and post-filter tables each contain 25 samples, and all four
final validation PH lists contain 25 named, non-null outputs. MV8-B later
closed their public accessions and proved that these are exactly the 25
`SRA779509` libraries already among the 31 bone-marrow rows. The historical
processing flow is 25-to-25, but it is not an independent external cohort.

The code-first provenance implementation is recorded in [`../PROVENANCE_INSTRUMENTATION_2026-08-04.md`](../PROVENANCE_INSTRUMENTATION_2026-08-04.md). Exact historical evidence remains preserved in the local lineage repository rather than copying large artifacts here.

## Candidate configuration

The current `inst/extdata/iteration_config.csv` defines Raw, SCT Individual, SCT Whole, Seurat Integration, and Harmony Integration. The legacy result tables contain only Raw, SCT Individual, SCT Whole, and Integrated results, consistent with the four configurations reported in BioRxiv v2 and suggesting Harmony is a later software extension.

`inst/extdata/inputs/custom_iteration_inputs.R`, added in the same `analysis inputs` commit, is not an executable historical manifest: every object/matrix path is a `path/to/...` placeholder.

The untracked `example_run.r` is also not the paper command: it includes Harmony, points to a nonexistent metadata filename, and has no tracked provenance.

## Candidate legacy outputs

| File | Rows | Dataset labels | Dimensions/summaries | Integrity |
|---|---:|---|---|---|
| `MAIN_ALL_datasets_all_pairwise_statistics_meta_master.csv` | 3,600 | Integrated, Raw, SCT_Individual, SCT_Whole | dimension_0, dimension_1, euler, landscape | SHA-256 manifest recorded |
| `MAIN_Tissue all_cross_iteration_pairwise_stats.csv` | 320 | Cross-iteration tissue comparisons | Table-specific comparison fields | SHA-256 manifest recorded |
| `BONEMARROW_VALIDATION_ALL_datasets_all_pairwise_statistics_meta_master.csv` | 5,488 | Integrated/Raw/SCT Individual/SCT Whole bone marrow | dimension_0, dimension_1, euler, landscape | SHA-256 manifest recorded |
| `BONEMARROW_VALIDATION_Tissue all_cross_iteration_pairwise_stats.csv` | 40 | Cross-iteration validation comparisons | Table-specific comparison fields | SHA-256 manifest recorded |

The main all-datasets table includes approach, SRA, tissue, and one random-group-bootstrap comparison family. The validation table additionally includes sample comparisons.

### Statistical reconstruction warning

The main table contains 611 exact-zero permutation p-values; the validation table contains 5,328. These are legacy-output facts, not accepted inferential results. Phase 1 must reconstruct the resampling procedure and replace exact-zero Monte Carlo reporting with finite-sample estimates where appropriate.

## Missing links

| ID | Missing evidence | Why needed | Next search |
|---|---|---|---|
| B-01 | Resolved: exact 124-sample main manifest | Establish paper observation set | Historical `filtered_cells.csv`; three logged pre-PH exclusions |
| B-02 | Resolved: exact 25-sample GSE120221 manifest | Establish validation observation set | Historical 25-file input directory, pre/post-filter tables, and final PH lists |
| B-03 | Resolving expression-matrix locations or immutable data IDs | Re-run preprocessing/PH | Owner storage, historical machine paths, download/preparation scripts |
| B-04 | Package/session/hardware/seed configuration | Reproduce numerical results | Logs, saved session info, output timestamps, dissertation methods |
| B-05 | Seurat and PH intermediate objects | Distinguish recomputable from missing stages | Historical result directories and custom input objects |
| B-06 | Exact commands/parameters | Link software to outputs | Shell/R history, scripts, logs, README revisions |
| B-07 | Figure-to-table/output map | Regenerate paper figures | Figure source directories, plotting scripts, manuscript source |
| B-08 | Resolved: 127 -> 124 is post-QC exclusion; validation is independently 25 -> 25 | Reconcile counts | Code, dated logs, retained manifests, and final PH artifacts |

## Next P0-05 actions

1. Compare result-table values and labels with every BioRxiv v2 table/figure to determine which files are true paper outputs.
2. Recover file timestamps and generating-code revisions where possible.
3. Reconstruct exact commands and package/session/hardware/seed configuration where evidence permits.
4. Mark each paper figure/result as exact-match, plausible-match, conflicting, or missing.

## Acceptance condition

P0-05 is complete only when each reported BioRxiv v2 result maps to an exact input/configuration/output chain or is explicitly documented as unrecoverable.
