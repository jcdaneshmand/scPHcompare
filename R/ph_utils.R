# Helper utilities for the PH pipeline

# Load sparse matrices from a vector of file paths
load_sparse_matrices <- function(file_paths) {
  if (length(file_paths) == 0) {
    return(list(matrices = list(), sample_names = character(0)))
  }
  sparse_matrices <- lapply(file_paths, loadRData)
  sample_names <- gsub(".*/|\\.sparse.RData$", "", file_paths)
  names(sparse_matrices) <- sample_names
  list(matrices = sparse_matrices, sample_names = sample_names)
}

# Compute persistence diagrams for a list of expression matrices
compute_ph_batch <- function(expr_list, DIM, log_message, dataset_suffix, prefix,
                             max_cores = 6, memory_threshold = 0.25) {
  process_expression_list_with_monitoring(
    expr_list = expr_list,
    DIM = DIM,
    log_message = log_message,
    max_cores = max_cores,
    memory_threshold = memory_threshold,
    log_file = paste0("progress_log", prefix, dataset_suffix, ".csv"),
    results_file = paste0("intermediate_results", prefix, dataset_suffix, ".rds"),
    timeout_datasets = NULL
  )
}

# Save persistence diagram results to disk
save_ph_results <- function(result, expr_list, DIM, THRESHOLD,
                            dataset_suffix, prefix, log_message) {
  if (is.null(result)) {
    return(invisible(NULL))
  }
  PD_list <- result$PD_list
  thresholds <- result$thresholds
  names(PD_list) <- names(expr_list)
  names(thresholds) <- names(expr_list)
  tryCatch({
    saveRDS(PD_list,
            file = paste0("PD_list_dim", DIM, "_th", THRESHOLD,
                           prefix, dataset_suffix, ".Rds"))
    saveRDS(thresholds,
            file = paste0("thresholds_dim", DIM, prefix,
                           dataset_suffix, ".Rds"))
  }, error = function(e) {
    log_message(paste("Error saving persistence diagrams for", prefix, "data:",
                      e$message))
  })
}

# Assemble data iteration metadata for downstream processing
assemble_ph_results <- function(merged_unintegrated, integrated,
                                expr_list_raw, expr_list_sctInd,
                                expr_list_sctWhole, expr_list_integrated) {
  list(
    list(
      name = "Raw",
      seurat_obj = merged_unintegrated,
      assay = "RNA",
      bdm_matrix = "BDM_unintegrated_Raw.rds",
      sdm_matrix = "SDM_unintegrated_Raw.rds",
      pd_list = "PD_list_after_retries_unintegrated_Raw.rds",
      expr_list = expr_list_raw
    ),
    list(
      name = "SCT_Individual",
      seurat_obj = merged_unintegrated,
      assay = "SCT_Ind",
      bdm_matrix = "BDM_unintegrated_sctInd.rds",
      sdm_matrix = "SDM_unintegrated_sctInd.rds",
      pd_list = "PD_list_after_retries_unintegrated_sctInd.rds",
      expr_list = expr_list_sctInd
    ),
    list(
      name = "SCT_Whole",
      seurat_obj = merged_unintegrated,
      assay = "SCT",
      bdm_matrix = "BDM_unintegrated_sctWhole.rds",
      sdm_matrix = "SDM_unintegrated_sctWhole.rds",
      pd_list = "PD_list_after_retries_unintegrated_sctWhole.rds",
      expr_list = expr_list_sctWhole
    ),
    list(
      name = "Integrated",
      seurat_obj = integrated,
      assay = "integrated",
      bdm_matrix = "BDM_integrated.rds",
      sdm_matrix = "SDM_integrated.rds",
      pd_list = "PD_list_after_retries_integrated.rds",
      expr_list = expr_list_integrated
    )
  )
}

