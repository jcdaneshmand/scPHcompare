# scPHcompare

`scPHcompare` is an R package that wraps a collection of scripts for computing persistent homology (PH) on single-cell RNA-seq datasets and for comparing the resulting topological summaries across many samples. The code originated as several standalone R files; the package now exposes the main workflow through the function `run_unified_pipeline()`.

## Features

* Preprocessing and integration of single-cell datasets using **Seurat**
  (integration with LIGER and MNN is planned for a future release)
* Calculation of persistence diagrams for each sample
* Generation of distance matrices (Bottleneck, Spectral and Landscape) for PH objects
* Multiple clustering approaches including standard Seurat, k-means, hierarchical PH and spectral clustering
* Optional post-processing modules for:
  * Cluster comparison with statistical metrics
  * Betti curve analysis
  * Cross‑iteration comparisons
* Helper utilities for plotting, caching and validating intermediate results

## Installation

The package can be installed from a local checkout using `devtools`:

```r
# install.packages("devtools")
library(devtools)
install("path/to/scPHcompare")
library(scPHcompare)
```

`scPHcompare` depends on several Bioconductor and CRAN packages listed in the `DESCRIPTION` file. Make sure these are available in your R environment before running the pipeline.

## Quick start

Prepare a CSV file describing each dataset. The required columns are `File Path` and `Sample Name`. Additional metadata such as `SRA`, `Tissue` and `Approach` (scRNA‑seq or snRNA‑seq) will be preserved if provided. Column names can be customised when calling the pipeline.

```r
results <- run_unified_pipeline(
  metadata_path = "metadata.csv",
  results_dir = "results",
  num_cores = 8,
  integration_method = "seurat",
  run_cluster = TRUE,        # optional cluster comparison
  run_cross_iteration = TRUE,# optional cross‑iteration analysis
  run_betti = TRUE           # optional Betti curve comparison
)
```

The wrapper executes preprocessing and PH calculation first. If the optional modules are enabled the function `run_postprocessing_pipeline()` performs clustering, visualisation and distance matrix generation. Results such as Seurat objects and plots are written to the directory specified by `results_dir`.

Distance matrices can also be generated separately for each iteration by calling `process_iteration_calculate_matrices()` from `PH_PostProcessing_andAnalysis.R`.

## Output overview

The pipeline creates a number of subfolders within `results_dir`:

* `seurat_objects/` – processed Seurat objects for each iteration
* `plots/` – UMAPs, heatmaps, cluster metrics and Betti curve graphics
* `BDM_progress_*` – progress logs for distance matrix calculations
* `plots/betti_plots/betti_cache/` – cached Betti curve data for reuse

Pairwise statistics summarised across datasets are provided in CSV files under the `Result Tables/` directory of this repository.

## Functions exported

The package exports the following user facing functions:

* `run_unified_pipeline()` – entry point that runs preprocessing and optional post‑processing modules
* `run_postprocessing_pipeline()` – standalone function for clustering and analyses on existing PH results
* `run_modular_analysis()` – helper to selectively run cluster comparison, Betti curves or cross‑iteration steps
* `process_datasets_PH()` – lower level function performing PH calculations on input datasets
* `run_cross_iteration()` – summarise metrics across different data iterations

Refer to the documentation in the `R/` directory for details on additional parameters.

## Citation and license

This project is released under the MIT license (see `DESCRIPTION`). If you use `scPHcompare` in your work please cite the package and the persistent homology methods accordingly.
