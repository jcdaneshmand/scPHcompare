#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "usage: build_mv05c2_normalization_cache.R RAW_CACHE CELL_MANIFEST ",
    "SOURCE_PREFLIGHT CACHE_DIR AUDIT_CSV SEED", call. = FALSE
  )
}

raw_cache_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cell_manifest_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source_preflight_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
cache_dir <- args[[4L]]
audit_path <- args[[5L]]
seed <- as.integer(args[[6L]])
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

if (length(seed) != 1L || is.na(seed) || !seed %in% 20260805:20260809) {
  stop("SEED is outside the frozen MV5 seed set.", call. = FALSE)
}
cells <- utils::read.csv(
  cell_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
source_preflight <- utils::read.csv(
  source_preflight_path, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "sample_id", "seed", "cell_id", "selected_cell_sha256",
  "outcome_label_state", "biological_outcomes_computed"
)
if (!all(required %in% names(cells)) ||
    any(c("tissue", "approach") %in% names(cells)) ||
    any(cells$outcome_label_state != "closed") ||
    any(as.logical(cells$biological_outcomes_computed)) ||
    nrow(source_preflight) != 1L) {
  stop("Input manifests violate the frozen label boundary.", call. = FALSE)
}
cells <- cells[cells$seed == seed, , drop = FALSE]
if (!nrow(cells)) stop("No selected cells exist for SEED.", call. = FALSE)
cell_groups <- split(cells, cells$sample_id)
if (any(vapply(cell_groups, nrow, integer(1L)) != 384L)) {
  stop("Every sample-seed cache requires exactly 384 selected cells.",
       call. = FALSE)
}

source_sha256 <- source_preflight$private_cache_sha256[[1L]]
if (!identical(
  unname(file.info(raw_cache_path)$size),
  as.numeric(source_preflight$private_cache_size_bytes[[1L]])
)) {
  stop("Raw cache size does not match the frozen source preflight.",
       call. = FALSE)
}
raw_cache <- NULL

seurat_version <- as.character(utils::packageVersion("Seurat"))
rows <- list()
for (sample_id in sort(names(cell_groups), method = "radix")) {
  selected <- cell_groups[[sample_id]]
  selected_ids <- as.character(selected$cell_id)
  selected_sha <- unique(selected$selected_cell_sha256)
  if (length(selected_sha) != 1L) {
    stop("Selected-cell identity failed for ", sample_id, call. = FALSE)
  }
  identity <- mv05c2_normalization_cache_identity_v1(
    sample_id = sample_id, seed = seed,
    selected_cell_sha256 = selected_sha,
    source_cache_sha256 = source_sha256,
    seurat_version = seurat_version,
    variable_features_n = 3000L, return_only_var_genes = FALSE
  )
  safe_id <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
  cache_file <- paste0(safe_id, "__", seed, "__sct.rds")
  cache_path <- file.path(cache_dir, cache_file)
  disposition <- mv05c2_cache_disposition_v1(cache_path, identity$cache_key)
  started <- proc.time()[["elapsed"]]
  if (identical(disposition, "build_missing")) {
    if (is.null(raw_cache)) {
      raw_cache <- readRDS(raw_cache_path)
      if (!identical(raw_cache$contract_id,
                     "mv05c_existing_data_raw_cache_v1")) {
        stop("Raw cache contract is invalid.", call. = FALSE)
      }
    }
    if (!sample_id %in% names(raw_cache$samples) ||
        !all(selected_ids %in% colnames(raw_cache$samples[[sample_id]]))) {
      stop("Raw sample/cell identity failed for ", sample_id, call. = FALSE)
    }
    counts <- raw_cache$samples[[sample_id]][, selected_ids, drop = FALSE]
    colnames(counts) <- paste(sample_id, selected_ids, sep = "__")
    object <- Seurat::CreateSeuratObject(
      counts = Matrix::Matrix(counts, sparse = TRUE), project = sample_id,
      min.cells = 0L, min.features = 0L
    )
    object$sample_id <- sample_id
    object <- Seurat::SCTransform(
      object, variable.features.n = 3000L, return.only.var.genes = FALSE,
      verbose = FALSE, seed.use = seed
    )
    record <- mv05c2_new_normalization_cache_record_v1(
      identity,
      list(sct_object = object, selected_cell_ids = selected_ids)
    )
    partial <- tempfile(
      pattern = paste0(cache_file, "."), tmpdir = cache_dir
    )
    on.exit(unlink(partial), add = TRUE)
    saveRDS(record, partial, compress = FALSE, version = 3)
    if (!file.rename(partial, cache_path)) {
      stop("Failed to atomically publish normalization cache: ", sample_id)
    }
    disposition <- "built_atomic"
  } else {
    record <- readRDS(cache_path)
    mv05c2_validate_normalization_cache_record_v1(record)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05c2_sample_seed_sct_cache_audit_v1",
    sample_id = sample_id, seed = seed,
    selected_cells = length(record$payload$selected_cell_ids),
    selected_cell_sha256 = identity$selected_cell_sha256,
    normalization_cache_key = identity$cache_key,
    payload_sha256 = record$payload_sha256,
    private_cache_file = cache_file,
    private_cache_size_bytes = unname(file.info(cache_path)$size),
    private_cache_sha256 = digest::digest(
      file = cache_path, algo = "sha256", serialize = FALSE
    ),
    disposition = disposition, operation_seconds = elapsed,
    seurat_version = seurat_version,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  rm(record)
  invisible(gc())
}
audit <- do.call(rbind, rows)
write_provenance_csv(audit, audit_path)
message("MV5-C2 normalization cache stage completed for seed ", seed, ".")
