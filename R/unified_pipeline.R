#' Run the complete PH pipeline
#'
#' This wrapper function first processes the datasets using
#' `process_datasets_PH()` and then optionally performs the
#' post-processing analyses provided by `run_postprocessing_pipeline()`.
#'
#' @param metadata_path Path to a CSV file containing dataset metadata.
#' @param results_dir Directory where results should be written.
#' @param num_cores Number of cores to use for parallel tasks.
#' @param integration_methods Character vector of integration methods to apply
#'   when processing datasets. Use "seurat" to produce the Seurat integration
#'   iteration, "harmony" to produce the Harmony iteration, or supply both to
#'   run them side by side (raw and SCT iterations are always generated
#'   alongside the chosen integration outputs).
#' @param run_cluster If `TRUE`, run clustering comparisons during
#'   post-processing.
#' @param run_cross_iteration If `TRUE`, perform cross-iteration analyses.
#' @param run_betti If `TRUE`, compute and compare Betti curves.
#' @param custom_iteration_inputs Optional named list mapping iteration labels or
#'   prefixes to lists with `seurat_object_path` and `ph_list_path` entries, or a
#'   path to an R script that defines such a list. When provided, matching
#'   iterations load the supplied Seurat objects and PH lists before recomputing
#'   distance matrices and clustering outputs. If `NULL`, the function looks for
#'   a populated template at `inst/extdata/custom_iteration_inputs_template.R`
#'   and will import it automatically when valid paths are detected.
#' @param ... Additional arguments passed to the underlying processing
#'   functions.
#'
#' @return The list produced by `process_datasets_PH()` that contains the
#'   processed iterations and metadata column names.
#'
#' @examples
#' \dontrun{
#' metadata <- read.csv("./data/metadata.csv")
#' run_unified_pipeline(metadata_path = "./data/metadata.csv")
#' }
#' @export
run_unified_pipeline <- function(metadata_path,
                                 results_dir = "results",
                                 num_cores = 8,
                                 integration_methods = "seurat",
                                 run_cluster = FALSE,
                                 run_betti = FALSE,
                                 run_cross_iteration = FALSE,
                                 custom_iteration_inputs = NULL,
                                 ...) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required")
  }
  custom_iteration_inputs <- import_custom_iteration_inputs(custom_iteration_inputs)
  # Functions are available once the package is loaded
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)
  ph_results <- process_datasets_PH(metadata,
                                    integration_methods = integration_methods,
                                    num_cores = num_cores,
                                    ...)
  if (run_cluster || run_betti || run_cross_iteration) {
    tryCatch(
      run_postprocessing_pipeline(ph_results,
                                  results_dir = results_dir,
                                  num_cores = num_cores,
                                  run_cluster = run_cluster,
                                  run_betti = run_betti,
                                  run_cross_iteration = run_cross_iteration,
                                  metadata_path = metadata_path,
                                  SRA_col = ph_results$SRA_col,
                                  Tissue_col = ph_results$Tissue_col,
                                  Approach_col = ph_results$Approach_col,
                                  custom_iteration_inputs = custom_iteration_inputs,
                                  ...),
      error = function(e) {
        stop(
          sprintf("Post-processing pipeline failed: %s", conditionMessage(e)),
          call. = FALSE
        )
      }
    )
  }
  ph_results
}
