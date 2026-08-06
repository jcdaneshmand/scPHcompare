# scPHcompare

`scPHcompare` is an R package that wraps a collection of scripts for computing persistent homology (PH) on single-cell RNA-seq datasets and for comparing the resulting topological summaries across many samples. The code originated as several standalone R files; the package now exposes the main workflow through the function `run_unified_pipeline()`.

> **Scientific audit status:** the legacy/current PH route passes Seurat feature-by-cell matrices directly to a row-as-point PH API, while the dissertation intended cells as observations. Existing biological persistence diagrams and all distances, clustering, statistics, and figures derived from them are therefore historical reproduction artifacts, not validated corrected results. The compatibility pipeline remains available but is not an approved scientific default; see `docs/audits/LANDSCAPE_ORACLE_AND_DIAGRAM_ELIGIBILITY_2026-08-05.md`.

## Features

* Preprocessing and integration of single-cell datasets using **Seurat**
* Calculation of persistence diagrams for each sample
* Generation of distance matrices (Bottleneck, Spectral and Landscape) for PH objects
* Multiple clustering approaches including standard Seurat, k-means, hierarchical PH and spectral clustering
* Post-processing modules include:
  * Cluster comparison with statistical metrics
  * Betti, Euler, and Landscape curve analysis
  * Cross‑iteration comparisons (i.e, raw vs integrated)
* Helper utilities for plotting, caching and validating intermediate results
* `generate_toy_data()` helper for recreating a set of small synthetic datasets for demonstrations and tests

## Installation

After cloning the repository, restore the `renv` environment first so dependencies (including `devtools`) match the lockfile. Then install the package from the restored environment:

```r
# install.packages("renv")  # run once per machine
renv::restore()

library(devtools)            # installed from the lockfile
devtools::install("path/to/scPHcompare")
library(scPHcompare)
```

`renv::restore()` installs the specific Bioconductor and CRAN package versions recorded in `renv.lock`. R sessions started from the project root automatically use this library; alternatively call `renv::activate()` manually. If you add or upgrade packages, run `renv::snapshot()` to update the lockfile before committing.

## Input data

At present `scPHcompare` mainly supports expression matrices stored in `.RData` files. Datasets in this format can be obtained from resources such as [PanglaoDB](https://panglaodb.se/) for a wide range of tissues. Support for other standard matrix formats is planned for a future release.

## Quick start

Prepare a CSV file describing each dataset. The required columns are `File Path` and `Sample Name`. Additional metadata such as `SRA`, `Tissue` and `Approach` (scRNA‑seq or snRNA‑seq) will be preserved if provided. **At least one of `SRA`, `Tissue`, or `Approach` must be present to allow grouping and downstream analyses.** The expression matrices referenced by `File Path` must all have the same dimensionality so they can be aligned for merging and distance calculations (the pipeline pads missing features with zeros when needed but assumes each matrix represents the same gene space).

```r
results <- run_unified_pipeline(
  metadata_path = "metadata.csv",
  results_dir = "results",
  num_cores = 8,
  integration_methods = "seurat", # or include "harmony" to run the Harmony iteration
  run_cluster = TRUE,        # optional cluster metrics output
  run_betti = TRUE,          # optional Betti curve comparison
  run_cross_iteration = TRUE # cross-iteration comparisons
)
```

The wrapper executes preprocessing and PH calculation first. If the optional modules are enabled the function `run_postprocessing_pipeline()` performs clustering, visualisation and distance matrix generation. Results such as Seurat objects, plots and comparison tables are written to the directory specified by `results_dir`.

### Integration iterations and labels

`run_unified_pipeline()` can generate multiple iterations in one run. The available iteration labels and their associated assay names are listed in `inst/extdata/iteration_config.csv`:

* `Seurat Integration` → stored under the prefix `seurat_integration` with assay `integrated`
* `Harmony Integration` → stored under the prefix `harmony_integration` with assay `harmony`
* Additional rows describe the raw and SCT-derived iterations (`Raw`, `SCT_Individual`, `SCT_Whole`) that are produced alongside the integration outputs.

Select the integration strategy with the `integration_methods` argument (`"seurat"`, `"harmony"`, or both). The chosen method controls which integration iteration is generated, while the other baseline iterations (raw and SCT variants) remain available for comparison. The iteration labels from `iteration_config.csv` are reused when naming output files (for example, `seurat_integration_seurat_object.rds` or `harmony_integration_seurat_object.rds`).

Distance matrices can also be generated separately for each iteration by calling `process_iteration_calculate_matrices()` from `PH_PostProcessing_andAnalysis.R`.

### Providing custom Seurat objects or PH lists for an iteration

If you want to rerun the clustering and visualisation steps on your own Seurat
objects or PH lists, supply a named `custom_iteration_inputs` list **or a path**
to an R script that defines that list. A template containing every iteration is
available at `inst/extdata/custom_iteration_inputs_template.R`; copy it,
populate the file paths you wish to override (Seurat object, PD list, BDM, SDM,
landscape list or landscape L2 matrix), and either pass the resulting list
directly or point `custom_iteration_inputs` to the edited file. These helpers
are exported for convenience:

```r
# Load the packaged template (returns NULL if no paths are filled in)
load_custom_iteration_inputs_template()

# Import a file or list and receive a normalised override object
import_custom_iteration_inputs("./custom_iteration_inputs_template.R")
```

The overrides are applied before distance matrices are generated, ensuring all
clustering approaches are rerun on the substituted data while reusing any valid
precomputed matrices already on disk:

```r
custom_inputs <- list(
  "Seurat Integration" = list(
    seurat_object_path = "./overrides/seurat_integration_object.rds",
    ph_list_path = "./overrides/seurat_integration_ph_list.rds",
    bdm_matrix_path = "./overrides/seurat_integration_bdm.rds",
    sdm_matrix_path = "./overrides/seurat_integration_sdm.rds",
    landscape_list_path = "./overrides/seurat_integration_landscape_list.rds",
    landscape_l2_distance_matrix_path = "./overrides/seurat_integration_landscape_l2.rds"
  ),
  Raw = list(ph_list_path = "./overrides/raw_pd_list.rds")
)

results <- run_unified_pipeline(
  metadata_path = "metadata.csv",
  results_dir = "results",
  custom_iteration_inputs = custom_inputs,
  run_cluster = TRUE
)
```

Alternatively, point `custom_iteration_inputs` to the edited template file so
the pipeline can read the overrides directly from disk:

```r
results <- run_unified_pipeline(
  metadata_path = "metadata.csv",
  results_dir = "results",
  custom_iteration_inputs = "./overrides/custom_iteration_inputs_template.R",
  run_cluster = TRUE
)
```

## Deterministic analytical baseline

The maintained public smoke route runs real persistent homology on two tiny
point clouds with known topology and writes a stable scientific manifest plus
stage-level timing evidence:

```r
library(scPHcompare)
baseline <- run_toy_baseline(
  output_dir = "toy_baseline",
  seed = 20260805L
)
baseline$manifest
baseline$timings
```

The rotated square must produce three finite H0 features and one H1 feature;
the rotated line must produce three finite H0 features and no H1 feature. The
seed controls the rotation while preserving those distances, and the caller's
RNG state is restored after the run. `baseline_manifest.csv` contains the seed,
dependency versions, output digests, and expected feature counts. Timing and
wall-clock fields are kept in `stage_timings.csv` and `ph_attempt_log.csv` so
they do not make the scientific manifest nondeterministic. The manifest also
records the maximum error against the known H0/H1 birth/death values and the
accepted numerical tolerance (`1e-10`).

## Realistic production-route fixture

`run_realistic_fixture()` generates two redistributable 520-gene, 20-cell
synthetic samples and sends them through the actual sparse-RData loader, fixed
QC gates, Seurat construction, per-sample and merged SCTransform, Harmony
integration, and PH subprocess route. H0 remains the default; H1 is an explicit
option and is also covered by the locked reference contract:

```r
fixture <- run_realistic_fixture(
  output_dir = "realistic_fixture",
  seed = 20260805L
)
fixture$manifest
fixture$results$provenance$pipeline_metrics

h1_fixture <- run_realistic_fixture(
  output_dir = "realistic_fixture_h1",
  seed = 20260805L,
  max_dimension = 1L
)
```

The default seed is checked against the packaged reference contract in
`inst/extdata/realistic_fixture_reference.csv`. The contract covers sample and
feature counts, post-QC eligibility, output dimensions, finite values, H0/H1
interval counts, PH completion, and locked-environment hashes. The fixture selects only the
Harmony PH representation to bound smoke-test runtime; normal
`process_datasets_PH()` calls still compute every historical representation by
default. PH attempt provenance includes sampled monitor, child, descendant, and
process-tree peak RSS fields. These are poll-boundary sampled peaks, not
continuous operating-system maxima. The unified wrapper separately writes
stage-boundary process RSS to `pipeline_stage_metrics.csv`.

For a repeated H0/H1 comparison, run
`Rscript scripts/profile_realistic_ph_dimensions.R <empty-output-directory>`.
The profiler requires stable input/Harmony hashes, repeatable within-dimension
PH hashes, identical H0 output when H1 is enabled, and positive H1 counts.

## Legacy synthetic matrix generator

`generate_toy_data()` recreates 20 sparse 100-by-300 matrices spanning five
tissues, two sequencing approaches (`scRNA-seq` and `snRNA-seq`), and two SRA
identifiers. It is retained for matrix-loading and metadata demonstrations.
These matrices have only 100 features, below the full pipeline's default
500-feature-per-cell QC requirement, so they are **not** an end-to-end example
for `run_unified_pipeline()` and should not be represented as one.

```r
toy_files <- generate_toy_data()
loaded <- scPHcompare:::load_sparse_matrices(toy_files$matrices)
length(loaded$matrices)
```

All fixture families are synthetic and must not be used for biological
interpretation. The realistic fixture validates software routing and
reproducibility, not biological conclusions.

## Output overview

The pipeline creates several subfolders under `results_dir`. These are
initialised automatically when a step first writes output and are reused
in later runs:

* `seurat_objects/` – created by `run_postprocessing_pipeline()` when the
  processed Seurat object for an iteration is saved. Files are named
  `<iteration>_seurat_object.rds`, for example
  `results/seurat_objects/raw_seurat_object.rds`.
* `plots/` – generated by the visualisation routines. This directory
  contains UMAP plots, heatmaps and cluster comparison graphics. Within
  `withinIterationClusterComparison/` the cluster comparison module saves
  plots alongside raw and normalised metrics tables
  (`*_raw_comparison_metrics_with_pvals.csv` and
  `*_normalized_metrics_with_pvals.csv`). It also houses subfolders such
  as `betti_plots/`. Example files include
  `results/plots/UMAP_Plots_Raw_All_Clusters.pdf`,
  `results/plots/UMAP_Plots_sctInd_All_Clusters.pdf`,
  `results/plots/UMAP_Plots_sctWhole_All_Clusters.pdf`, and
  `results/plots/UMAP_Plots_Integrated_All_Clusters.pdf`.
* `BDM_progress_*` – temporary progress files written while computing the
  Bottleneck distance matrix in `process_iteration_calculate_matrices()`.
  They allow the computation to resume if interrupted, e.g.
  `results/BDM_progress_Raw.rds`.
* `plots/betti_plots/betti_cache/` – cache directory for Betti curve
  computations. Each dataset gets its own subfolder with hashed cache
  files such as
  `results/plots/betti_plots/betti_cache/Raw/cache_abcd1234.rds`.
* `plots/betti_plots/<dataset_name>/statistical_comparison_heatmaps/statistical_csvs/`
  – directory where pairwise-statistics tables (CSV) are saved for each
  comparison. Replace `<dataset_name>` with the name of the iteration.
* `cross_iteration_comparisons/` – output from the cross‑iteration
  module. Each grouping generates aggregated pairwise statistics, for
  example `Tissue_<tissue>_all_cross_iteration_pairwise_stats.csv`.

### Clustering metadata column naming

Clustering results added to the Seurat metadata follow consistent naming
patterns to ease downstream filtering. Columns created by
`apply_all_clustering_methods()` use lower‑case dataset names and a suffix
matching the grouping variable:

* K‑means clustering: `kmeans_cluster_<dataset>_<group>`
* Hierarchical clustering on PH distances:
  `hierarchical_cluster_<matrix>_ph_<dataset>_<group>_<linkage>` where
  `<matrix>` is `bdm`, `sdm`, or `landscape`
* Spectral clustering: `spectral_cluster_<matrix>_<dataset>_<group>`

Grouping suffixes correspond to `tissue`, `sra`, or `approach`. These
columns are consumed by `run_cluster_comparison()` when generating
comparison plots and metrics.

Pairwise statistics summarised across datasets for the manuscript are archived in the repository's `Result Tables/` directory. This folder contains the original analysis tables and should be preserved.

## Functions exported

The package exports the following user facing functions:

* `run_unified_pipeline()` – entry point that runs preprocessing and optional post‑processing modules
* `run_postprocessing_pipeline()` – standalone function for clustering and analyses on existing PH results
* `run_modular_analysis()` – runs cluster comparison, Betti curves and cross‑iteration analysis
* `process_datasets_PH()` – lower level function performing PH calculations on input datasets


Refer to the documentation in the `R/` directory for details on additional parameters.

## Citation and license

This project is released under the MIT license (see `DESCRIPTION`). If you use `scPHcompare` in your work please cite the package and the persistent homology methods accordingly.
