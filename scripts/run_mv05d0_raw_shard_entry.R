#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    "usage: run_mv05d0_raw_shard_entry.R SOURCE_RDS RAW_DIR AUDIT_CSV ",
    "SAMPLE_ID EXPECTED_CELLS HISTORICAL_SHA256 OVERLAP_DIR SOURCE_ROOT",
    call. = FALSE
  )
}
source_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
raw_dir <- args[[2L]]
audit_path <- args[[3L]]
sample_id <- args[[4L]]
expected_cells <- as.integer(args[[5L]])
historical_sha <- args[[6L]]
overlap_dir <- args[[7L]]
source_root <- normalizePath(args[[8L]], winslash = "/", mustWork = TRUE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

if (!startsWith(source_path, paste0(source_root, "/")) ||
    !identical(tools::file_path_sans_ext(basename(source_path)), sample_id) ||
    length(expected_cells) != 1L || is.na(expected_cells) ||
    expected_cells < 384L) {
  stop("Individual source path or expected shape is invalid.", call. = FALSE)
}
source_sha <- digest::digest(
  file = source_path, algo = "sha256", serialize = FALSE
)
object <- readRDS(source_path)
if (!inherits(object, "Seurat") || !"RNA" %in% names(object@assays)) {
  stop("Individual source is not a Seurat object with an RNA assay.",
       call. = FALSE)
}
counts <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
rm(object)
invisible(gc())
cell_ids <- colnames(counts)
prefix <- paste0(sample_id, "_")
prefixed <- startsWith(cell_ids, prefix)
if (all(prefixed)) {
  colnames(counts) <- substring(cell_ids, nchar(prefix) + 1L)
} else if (any(prefixed)) {
  stop("Individual source mixes prefixed and unprefixed cell IDs.",
       call. = FALSE)
}
counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
if (ncol(counts) != expected_cells || is.null(rownames(counts)) ||
    is.null(colnames(counts)) || anyDuplicated(rownames(counts)) ||
    anyDuplicated(colnames(counts))) {
  stop("Individual source counts fail the frozen sample shape.", call. = FALSE)
}
record <- mv05d0_new_raw_sample_cache_v2(
  sample_id, counts, historical_sha, source_sha
)
safe_id <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
cache_file <- paste0(safe_id, "__raw.rds")
cache_path <- file.path(raw_dir, cache_file)
disposition <- "built_atomic"
if (file.exists(cache_path)) {
  existing <- readRDS(cache_path)
  mv05d0_validate_raw_sample_cache_v2(existing)
  if (!identical(existing$sample_id, record$sample_id) ||
      !identical(existing$historical_source_sha256,
                 record$historical_source_sha256) ||
      !identical(existing$individual_source_sha256,
                 record$individual_source_sha256) ||
      !identical(existing$counts_sha256, record$counts_sha256)) {
    stop("Existing raw sample cache is stale; refusing overwrite.",
         call. = FALSE)
  }
  record <- existing
  disposition <- "reuse_validated"
} else {
  partial <- tempfile(pattern = paste0(cache_file, "."), tmpdir = raw_dir)
  saveRDS(record, partial, compress = "gzip", version = 3)
  if (!file.rename(partial, cache_path)) {
    unlink(partial)
    stop("Failed to atomically publish raw sample cache.", call. = FALSE)
  }
}
overlap_path <- file.path(overlap_dir, cache_file)
overlap_disposition <- "not_available"
if (file.exists(overlap_path)) {
  legacy <- readRDS(overlap_path)
  if (!identical(legacy$contract_id, "mv05d0_raw_sample_cache_v1") ||
      !identical(legacy$sample_id, sample_id) ||
      !identical(legacy$outcome_label_state, "closed") ||
      !identical(legacy$biological_outcomes_computed, FALSE) ||
      !identical(mv05d0_count_matrix_sha256_v1(legacy$counts),
                 record$counts_sha256)) {
    stop("Individual source differs from recovered monolithic shard.",
         call. = FALSE)
  }
  overlap_disposition <- "exact_content_identity"
}
audit <- data.frame(
  contract_id = "mv05d0_individual_raw_shard_audit_v1",
  sample_id = sample_id, genes = nrow(record$counts),
  cells = ncol(record$counts),
  historical_source_sha256 = historical_sha,
  individual_source_file = basename(source_path),
  individual_source_size_bytes = unname(file.info(source_path)$size),
  individual_source_sha256 = source_sha,
  counts_sha256 = record$counts_sha256,
  private_raw_cache_file = cache_file,
  private_raw_cache_size_bytes = unname(file.info(cache_path)$size),
  private_raw_cache_sha256 = digest::digest(
    file = cache_path, algo = "sha256", serialize = FALSE
  ),
  disposition = disposition,
  recovered_monolithic_comparison = overlap_disposition,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(audit, audit_path)
message("Completed individual raw shard: ", sample_id)
