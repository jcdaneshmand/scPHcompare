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
#' @param pipeline_metrics_file Optional path for the stage-level runtime and
#'   process-RSS CSV. Defaults to `pipeline_stage_metrics.csv` in `results_dir`.
#' @param corrected_landscape_control Optional explicit v1 control for additive,
#'   versioned corrected-landscape artifacts. `NULL` (default) performs no
#'   corrected-landscape work and preserves historical workflow behavior.
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
                                 pipeline_metrics_file = NULL,
                                 corrected_landscape_control = NULL,
                                 ...) {
  pipeline_started_at <- Sys.time()
  pipeline_rss_before <- current_process_rss()
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required")
  }
  custom_iteration_inputs <- import_custom_iteration_inputs(custom_iteration_inputs)
  # Functions are available once the package is loaded
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  if (is.null(pipeline_metrics_file)) {
    pipeline_metrics_file <- file.path(results_dir, "pipeline_stage_metrics.csv")
  }
  metadata_stage <- time_stage(
    "metadata_load",
    readr::read_csv(metadata_path, show_col_types = FALSE)
  )
  metadata <- metadata_stage$value
  timings <- metadata_stage$timing
  processing_stage <- time_stage(
    "ph_processing",
    process_datasets_PH(metadata,
                        integration_methods = integration_methods,
                        num_cores = num_cores,
                        ...)
  )
  ph_results <- processing_stage$value
  timings <- rbind(timings, processing_stage$timing)
  if (run_cluster || run_betti || run_cross_iteration ||
      !is.null(corrected_landscape_control)) {
    postprocessing_args <- c(list(
      ph_results = ph_results,
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
      corrected_landscape_control = corrected_landscape_control
    ), list(...))
    corrected_only <- !is.null(corrected_landscape_control) &&
      !run_cluster && !run_betti && !run_cross_iteration
    if (corrected_only) {
      postprocessing_args[c(
        "run_standard_seurat_clustering", "run_kmeans_clustering",
        "run_hierarchical_ph_clustering", "run_spectral_clustering",
        "run_visualizations", "run_sample_level_heatmap"
      )] <- rep(list(FALSE), 6L)
    }
    postprocessing_stage <- time_stage(
      "postprocessing",
      tryCatch(
        do.call(run_postprocessing_pipeline, postprocessing_args),
        error = function(e) {
          stop(
            sprintf("Post-processing pipeline failed: %s", conditionMessage(e)),
            call. = FALSE
          )
        }
      )
    )
    timings <- rbind(timings, postprocessing_stage$timing)
    if (!is.null(corrected_landscape_control) &&
        !is.null(postprocessing_stage$value)) {
      ph_results <- postprocessing_stage$value
    }
  }
  pipeline_finished_at <- Sys.time()
  timings <- rbind(
    timings,
    new_stage_timing(
      "pipeline_total",
      started_at = pipeline_started_at,
      finished_at = pipeline_finished_at,
      rss_before_bytes = pipeline_rss_before,
      rss_after_bytes = current_process_rss()
    )
  )
  dir.create(dirname(pipeline_metrics_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(timings, pipeline_metrics_file, row.names = FALSE, na = "")
  if (is.null(ph_results$provenance)) {
    ph_results$provenance <- list()
  }
  ph_results$provenance$pipeline_metrics <- timings
  ph_results$provenance$pipeline_metrics_file <- normalizePath(
    pipeline_metrics_file, winslash = "/", mustWork = TRUE
  )
  ph_results
}
