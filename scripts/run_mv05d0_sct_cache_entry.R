#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: run_mv05d0_sct_cache_entry.R RAW_SAMPLE SELECTION_CSV ",
    "CACHE_DIR AUDIT_CSV SAMPLE_ID SEED", call. = FALSE
  )
}
raw_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
selection_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
cache_dir <- args[[3L]]
audit_path <- args[[4L]]
sample_id <- args[[5L]]
seed <- as.integer(args[[6L]])
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")
if (!requireNamespace("future", quietly = TRUE)) {
  stop("future is required to freeze sequential SCT execution.", call. = FALSE)
}
future::plan(future::sequential)

selection <- utils::read.csv(
  selection_path, stringsAsFactors = FALSE, check.names = FALSE
)
row <- selection[selection$sample_id == sample_id & selection$seed == seed,
                 , drop = FALSE]
if (nrow(row) != 1L || row$outcome_label_state != "closed" ||
    as.logical(row$biological_outcomes_computed) ||
    row$selected_cells != 384L || !seed %in% 20260805:20260809) {
  stop("Sample-seed selection is outside the frozen MV5-D0 contract.",
       call. = FALSE)
}
raw <- readRDS(raw_path)
raw_valid <- if (identical(raw$contract_id, "mv05d0_raw_sample_cache_v2")) {
  tryCatch(mv05d0_validate_raw_sample_cache_v2(raw), error = identity)
} else if (identical(raw$contract_id, "mv05d0_raw_sample_cache_v1") &&
           identical(raw$outcome_label_state, "closed") &&
           identical(raw$biological_outcomes_computed, FALSE)) {
  raw
} else {
  simpleError("unsupported raw cache contract")
}
if (inherits(raw_valid, "error") || !identical(raw$sample_id, sample_id)) {
  stop("Raw sample cache identity is invalid.", call. = FALSE)
}
selected <- select_matched_cells(colnames(raw$counts), n = 384L, seed = seed)
if (!identical(attr(selected, "selected_cell_sha256"),
               row$selected_cell_sha256[[1L]])) {
  stop("Recomputed selected-cell identity differs from the frozen summary.",
       call. = FALSE)
}
runtime <- mv05d0_normalization_runtime_v1()
seurat_version <- runtime$seurat_version
identity <- mv05d0_normalization_cache_identity_v2(
  sample_id, seed, row$selected_cell_sha256[[1L]],
  raw$historical_source_sha256, runtime,
  variable_features_n = 3000L, return_only_var_genes = FALSE
)
safe_id <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
cache_file <- paste0(safe_id, "__", seed, "__sct.rds")
cache_path <- file.path(cache_dir, cache_file)
disposition <- mv05d0_cache_disposition_v2(cache_path, identity$cache_key)
started <- proc.time()[["elapsed"]]
if (identical(disposition, "build_missing")) {
  counts <- raw$counts[, selected, drop = FALSE]
  colnames(counts) <- paste(sample_id, selected, sep = "__")
  object <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE), project = sample_id,
    min.cells = 0L, min.features = 0L
  )
  object$sample_id <- sample_id
  object <- Seurat::SCTransform(
    object, variable.features.n = 3000L, return.only.var.genes = FALSE,
    verbose = FALSE, seed.use = seed
  )
  sct_data <- Seurat::GetAssayData(object, assay = "SCT", layer = "data")
  record <- mv05d0_new_normalization_cache_record_v2(
    identity,
    list(
      payload_contract_id = "mv05d0_sct_data_matrix_v1",
      sct_data = sct_data,
      selected_cell_ids = as.character(selected)
    )
  )
  partial <- tempfile(pattern = paste0(cache_file, "."), tmpdir = cache_dir)
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, cache_path)) {
    unlink(partial)
    stop("Failed to atomically publish normalization cache.", call. = FALSE)
  }
  disposition <- "built_atomic"
} else {
  record <- readRDS(cache_path)
  mv05d0_sct_matrix_from_cache_v1(record)
}
elapsed <- proc.time()[["elapsed"]] - started
audit <- data.frame(
  contract_id = "mv05d0_sample_seed_sct_cache_audit_v1",
  sample_id = sample_id, seed = seed,
  selected_cells = length(record$payload$selected_cell_ids),
  selected_cell_sha256 = identity$selected_cell_sha256,
  normalization_cache_key = identity$cache_key,
  payload_contract_id = if (is.null(record$payload$payload_contract_id))
    "mv05c2_seurat_object_legacy_v1" else record$payload$payload_contract_id,
  payload_sha256 = record$payload_sha256,
  private_cache_file = cache_file,
  private_cache_size_bytes = unname(file.info(cache_path)$size),
  private_cache_sha256 = digest::digest(
    file = cache_path, algo = "sha256", serialize = FALSE
  ),
  disposition = disposition, operation_seconds = elapsed,
  seurat_version = seurat_version,
  sctransform_version = runtime$sctransform_version,
  matrix_version = runtime$matrix_version,
  r_version = runtime$r_version,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(audit, audit_path)
message("Completed MV5-D0 SCT cache entry: ", sample_id, " / ", seed)
