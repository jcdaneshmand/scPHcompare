# Helper utilities for the PH pipeline

if (!exists("SEURAT_INTEGRATION_LABEL")) {
  SEURAT_INTEGRATION_LABEL <- "Seurat Integration"
}
if (!exists("SEURAT_INTEGRATION_PREFIX")) {
  SEURAT_INTEGRATION_PREFIX <- "seurat_integration"
}
if (!exists("HARMONY_INTEGRATION_LABEL")) {
  HARMONY_INTEGRATION_LABEL <- "Harmony Integration"
}
if (!exists("HARMONY_INTEGRATION_PREFIX")) {
  HARMONY_INTEGRATION_PREFIX <- "harmony_integration"
}
if (!exists("HARMONY_ASSAY_NAME")) {
  HARMONY_ASSAY_NAME <- "harmony"
}

# Load the default iteration configuration shipped with the package
# This keeps iteration labels, prefixes, and assay names centralized so that
# downstream modules do not hard-code integration identifiers.
get_iteration_config <- function(config_path = system.file("extdata", "iteration_config.csv", package = "scPHcompare")) {
  if (!nzchar(config_path) || !file.exists(config_path)) {
    stop("Iteration configuration file is missing: ", config_path)
  }

  cfg <- utils::read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("label", "prefix", "assay")
  missing_cols <- setdiff(required_cols, colnames(cfg))
  if (length(missing_cols) > 0) {
    stop("Iteration configuration is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  cfg
}

# Retrieve the set of integration iteration labels from configuration or defaults
get_integration_labels <- function(iteration_cfg = NULL) {
  cfg <- iteration_cfg
  if (is.null(cfg)) {
    cfg <- tryCatch(get_iteration_config(), error = function(e) NULL)
  }

  if (!is.null(cfg) && "label" %in% names(cfg)) {
    unique(cfg$label)
  } else {
    unique(c(SEURAT_INTEGRATION_LABEL, HARMONY_INTEGRATION_LABEL))
  }
}

# Prioritize a preferred integration iteration while keeping other iterations intact
prioritize_integration_iterations <- function(data_iterations, preferred_integration = SEURAT_INTEGRATION_LABEL) {
  if (length(data_iterations) == 0) {
    return(data_iterations)
  }

  available_names <- vapply(data_iterations, function(iter) iter$name, character(1))
  integration_labels <- get_integration_labels()
  preferred_matches <- intersect(preferred_integration, integration_labels)

  preferred_idx <- which(available_names %in% preferred_matches)
  if (length(preferred_idx) == 0 && SEURAT_INTEGRATION_LABEL %in% available_names) {
    preferred_idx <- which(available_names == SEURAT_INTEGRATION_LABEL)[1]
  }

  if (length(preferred_idx) == 0) {
    return(data_iterations)
  }

  ordered_idx <- c(preferred_idx, setdiff(seq_along(data_iterations), preferred_idx))
  data_iterations[ordered_idx]
}

# Load sparse matrices from a vector of file paths
load_sparse_matrices <- function(file_paths) {
  if (length(file_paths) == 0) {
    return(list(matrices = list(), sample_names = character(0)))
  }
  missing <- file_paths[!file.exists(file_paths)]
  if (length(missing) > 0) {
    stop(
      sprintf(
        "The following files do not exist: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
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
                                expr_list_sctWhole, expr_list_integrated,
                                harmony = NULL, expr_list_harmony = NULL,
                                harmony_assay = HARMONY_ASSAY_NAME) {
  iterations <- list(
    list(
      name = "Raw",
      seurat_obj = merged_unintegrated,
      assay = "RNA",
      bdm_matrix = "BDM_unintegrated_Raw.rds",
      sdm_matrix = "SDM_unintegrated_Raw.rds",
      pd_list = "PD_list_unintegrated_Raw.rds",
      expr_list = expr_list_raw
    ),
    list(
      name = "SCT_Individual",
      seurat_obj = merged_unintegrated,
      assay = "SCT_Ind",
      bdm_matrix = "BDM_unintegrated_sctInd.rds",
      sdm_matrix = "SDM_unintegrated_sctInd.rds",
      pd_list = "PD_list_unintegrated_sctInd.rds",
      expr_list = expr_list_sctInd
    ),
    list(
      name = "SCT_Whole",
      seurat_obj = merged_unintegrated,
      assay = "SCT",
      bdm_matrix = "BDM_unintegrated_sctWhole.rds",
      sdm_matrix = "SDM_unintegrated_sctWhole.rds",
      pd_list = "PD_list_unintegrated_sctWhole.rds",
      expr_list = expr_list_sctWhole
    ),
    list(
      name = SEURAT_INTEGRATION_LABEL,
      seurat_obj = integrated,
      assay = "integrated",
      bdm_matrix = paste0("BDM_", SEURAT_INTEGRATION_PREFIX, ".rds"),
      sdm_matrix = paste0("SDM_", SEURAT_INTEGRATION_PREFIX, ".rds"),
      pd_list = paste0("PD_list_", SEURAT_INTEGRATION_PREFIX, ".rds"),
      expr_list = expr_list_integrated
    )
  )

  if (!is.null(harmony) && !is.null(expr_list_harmony)) {
    iterations <- c(iterations, list(list(
      name = HARMONY_INTEGRATION_LABEL,
      seurat_obj = harmony,
      assay = harmony_assay,
      bdm_matrix = paste0("BDM_", HARMONY_INTEGRATION_PREFIX, ".rds"),
      sdm_matrix = paste0("SDM_", HARMONY_INTEGRATION_PREFIX, ".rds"),
      pd_list = paste0("PD_list_", HARMONY_INTEGRATION_PREFIX, ".rds"),
      expr_list = expr_list_harmony
    )))
  }

  iterations
}

