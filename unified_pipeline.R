run_unified_pipeline <- function(metadata_path,
                                 results_dir = "results",
                                 num_cores = 8,
                                 integration_method = "seurat",
                                 run_modular = FALSE,
                                 run_cross_iteration = FALSE,
                                 ...) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required")
  }
  source("PH_Calculation.R")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)
  ph_results <- process_datasets_PH(metadata,
                                    integration_method = integration_method,
                                    num_cores = num_cores,
                                    ...)
  if (run_modular && exists("run_modular_analysis")) {
    try(run_modular_analysis(ph_results, results_dir = results_dir, ...), silent = TRUE)
  }
  if (run_cross_iteration && exists("run_cross_iteration")) {
    if (!is.null(ph_results$data_iterations)) {
      try(run_cross_iteration(ph_results$data_iterations,
                              results_folder = results_dir,
                              ...),
          silent = TRUE)
    } else {
      warning("Cross iteration requested but 'data_iterations' not found in results")
    }
  }
  ph_results
}
