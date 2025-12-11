# scPHcompare

`scPHcompare` is an R package that wraps a collection of scripts for computing persistent homology (PH) on single-cell RNA-seq datasets and for comparing the resulting topological summaries across many samples. The code originated as several standalone R files; the package now exposes the main workflow through the function `run_unified_pipeline()`.

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

The package can be installed from a local checkout using `devtools`:

```r
# install.packages("devtools")
library(devtools)
install("path/to/scPHcompare")
library(scPHcompare)
```

`scPHcompare` depends on several Bioconductor and CRAN packages listed in the `DESCRIPTION` file. Make sure these are available in your R environment before running the pipeline.

## Dependency management with renv

This repository tracks its R package dependencies with [renv](https://rstudio.github.io/renv/). After cloning the repository, install `renv` if necessary and restore the lockfile to create a project-local library:

```r
# install.packages("renv")  # run once per machine
renv::restore()
```

`renv::restore()` installs the specific package versions recorded in `renv.lock`. R sessions started from the project root automatically use this library; alternatively call `renv::activate()` manually. If you add or upgrade packages, run `renv::snapshot()` to update the lockfile before committing.

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
directly or point `custom_iteration_inputs` to the edited file. The overrides
are applied before distance matrices are generated, ensuring all clustering
approaches are rerun on the substituted data while reusing any valid
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

When `custom_iteration_inputs` is left as `NULL`, the pipeline also checks the
packaged template at `inst/extdata/custom_iteration_inputs_template.R`. If you
save real file paths into that template (or a copy of it on your library path),
those overrides will be picked up automatically without needing to pass an
argument, letting you swap in custom Seurat objects or PH lists for specific
iterations with no changes to your pipeline call.

## Toy Example Data

To try the package without obtaining real datasets, a set of synthetic datasets can be generated directly in R:

```r
library(scPHcompare)
toy_files <- generate_toy_data()
```

`generate_toy_data()` recreates 20 sparse 100×300 matrices spanning five tissues, two sequencing approaches (`scRNA-seq` and `snRNA-seq`), and two SRA identifiers. The matrices and a corresponding `metadata.csv` file are written to `inst/extdata/toy/` and their paths are returned. The metadata can then be used to run the pipeline:

```r
results <- run_unified_pipeline(
  metadata_path = toy_files$metadata,
  results_dir = "toy_results",
  num_cores = 2
)
```

These toy datasets are randomly generated and extremely small. They are intended only for demonstrations and automated tests and should not be used for biological interpretation.

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
* `run_cross_iteration()` – helper to run only the cross‑iteration Betti curve comparisons on a set of prepared iterations
* `perform_integration()` – integrates a list of Seurat objects with checkpointing and anchor caching before PH calculation


Refer to the documentation in the `R/` directory for details on additional parameters.

## Citation and license

This project is released under the MIT license (see `DESCRIPTION`). If you use `scPHcompare` in your work please cite the package and the persistent homology methods accordingly.
