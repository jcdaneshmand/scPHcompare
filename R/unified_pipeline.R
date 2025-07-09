run_unified_pipeline <- function(metadata_path,
                                 results_dir = "results",
                                 num_cores = 8,
                                 integration_method = "seurat",
                                 run_cluster = FALSE,
                                 run_modular = FALSE,
                                 run_cross_iteration = FALSE,
                                 run_betti = FALSE,
                                 ...) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required")
  }
  # Functions are available once the package is loaded
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)
  ph_results <- process_datasets_PH(metadata,
                                    integration_method = integration_method,
                                    num_cores = num_cores,
                                    ...)
  if (run_cluster || run_modular || run_cross_iteration || run_betti) {
    try(run_postprocessing_pipeline(ph_results,
                                    results_dir = results_dir,
                                    num_cores = num_cores,
                                    metadata_path = metadata_path,
                                    SRA_col = ph_results$SRA_col,
                                    Tissue_col = ph_results$Tissue_col,
                                    Approach_col = ph_results$Approach_col,
                                    ...),
        silent = TRUE)
  }
  ph_results
}
