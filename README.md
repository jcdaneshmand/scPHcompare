# scPHcompare

This repository contains a set of R scripts for performing persistent homology
analysis and downstream comparison of single cell datasets.  The code base was
originally organised as several stand‐alone scripts.  A small wrapper
`unified_pipeline.R` is now provided to run the whole workflow from a single
entry point.

## Running the pipeline

1. Prepare a metadata CSV describing each sample.  The
   `process_datasets_PH` function in `run_ph_pipeline.R` expects columns
   `File Path` and `Sample Name` at minimum.
2. Call `run_unified_pipeline()` with the metadata file and an output
   directory.  Additional arguments are forwarded to `process_datasets_PH`.

```r
source("unified_pipeline.R")
results <- run_unified_pipeline(
  metadata_path = "metadata.csv",
  results_dir = "results",
  num_cores = 8,
  integration_method = "seurat",
  run_cluster = FALSE,        # optional
  run_cross_iteration = TRUE, # optional
  run_betti = FALSE           # optional
)
```

The wrapper runs the preprocessing and persistent homology analysis first.  If
`run_modular_analysis()` is available it can perform cluster comparison, Betti
curve tests, and cross-iteration analysis based on the flags forwarded from
`run_unified_pipeline()`.  Cross-iteration can still be executed separately when
`run_modular = FALSE` and `run_cross_iteration = TRUE`.  See the individual
scripts for the detailed options of each step.

Distance matrices (BDM, SDM and LDM) are generated during the post-processing
phase.  Use `process_iteration_calculate_matrices()` from
`PH_PostProcessing_andAnalysis.R` to compute these for each data iteration.  The
matrix routines were removed from `PH_Calculation.R` so that file only performs
persistence diagram generation.

