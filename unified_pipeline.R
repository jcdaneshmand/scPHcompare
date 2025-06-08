# unified_pipeline.R
# High level wrapper to run the entire scPHcompare workflow

# Load main pipeline steps
source("run_ph_pipeline.R")

# cross_iteration_analysis.R and Modular_Analysis_and_PostProcessing.R define
# many helper functions used for downstream analysis. They are sourced lazily
# inside the wrapper when those steps are requested.

#' Run the full scPHcompare pipeline
#'
#' @param metadata_path path to metadata CSV used by `process_datasets_PH`
#' @param results_dir output directory for all results
#' @param num_cores number of cores to use throughout
#' @param integration_method integration method for preprocessing (default 'seurat')
#' @param run_modular whether to run the modular post-processing stage
#' @param run_cross whether to run cross-iteration comparisons
#' @return list with results from each stage
run_unified_pipeline <- function(metadata_path,
                                 results_dir = "scph_results",
                                 num_cores = parallel::detectCores() %/% 2,
                                 integration_method = "seurat",
                                 run_modular = TRUE,
                                 run_cross = TRUE,
                                 ...) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. run preprocessing & PH analysis
  ph_dir <- file.path(results_dir, "preprocessing")
  ph_res <- process_datasets_PH(
    metadata_path = metadata_path,
    output_dir = ph_dir,
    num_cores = num_cores,
    integration_method = integration_method,
    ...
  )

  out <- list(preprocessing = ph_res)

  # 2. optional modular analysis step
  if (run_modular) {
    mod_env <- new.env()
    sys.source("Modular_Analysis_and_PostProcessing.R", envir = mod_env)
    if (exists("run_modular_analysis", envir = mod_env)) {
      out$modular <- mod_env$run_modular_analysis(results_dir)
    } else {
      message("run_modular_analysis() not found, skipping modular stage")
      out$modular <- NULL
    }
  }

  # 3. optional cross iteration step
  if (run_cross) {
    cross_env <- new.env()
    sys.source("cross_iteration_analysis.R", envir = cross_env)
    if (exists("run_cross_iteration", envir = cross_env)) {
      out$cross_iteration <- cross_env$run_cross_iteration(results_dir)
    } else {
      message("run_cross_iteration() not found, skipping cross iteration stage")
      out$cross_iteration <- NULL
    }
  }

  out
}

