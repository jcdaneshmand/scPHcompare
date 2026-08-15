#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: validate_mv07f_upstream_production.R PREFREEZE_DIR PRIVATE_ROOT ",
       "PRODUCTION_DIR OUTPUT", call. = FALSE)
}
prefreeze <- args[[1L]]; private_root <- args[[2L]]; production <- args[[3L]]
output <- args[[4L]]
source("R/provenance_utils.R"); source("R/toy_baseline.R")
source("R/dual_view_topology.R"); source("R/mv05_resource_safe_execution.R")
source("R/mv07f_validation_utils.R")
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(as.character(x))) == "true"
readp <- function(name) read.csv(file.path(production, name),
  stringsAsFactors = FALSE, check.names = FALSE)
queue <- read.csv(file.path(prefreeze, "mv07f-upstream-queue.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
raw <- readp("mv07f-raw-production.csv"); sct <- readp("mv07f-sct-production.csv")
selection <- readp("mv07f-selection-summary.csv"); manifest <- readp("mv07f-cache-manifest.csv")
resource <- readp("mv07f-resource-summary.csv"); summary <- readp("mv07f-production-summary.csv")
provenance <- readp("mv07f-execution-provenance.csv")
axis_ok <- nrow(raw) == 34L && nrow(sct) == 170L && nrow(selection) == 170L &&
  nrow(manifest) == 204L && !anyDuplicated(raw$sample_id) &&
  !anyDuplicated(paste(sct$sample_id, sct$seed)) &&
  setequal(raw$sample_id, unique(sct$sample_id)) &&
  identical(sort(unique(sct$seed)), 20260805:20260809)
raw_ok <- logical(nrow(raw)); raw_records <- vector("list", nrow(raw))
for (i in seq_len(nrow(raw))) {
  path <- file.path(private_root, "raw", raw$private_cache_file[[i]])
  record <- readRDS(path); raw_records[[i]] <- record
  valid <- tryCatch({mv05d0_validate_raw_sample_cache_v2(record); TRUE},
                    error = function(e) FALSE)
  raw_ok[[i]] <- valid && record$sample_id == raw$sample_id[[i]] &&
    ncol(record$counts) == raw$observed_cells[[i]] &&
    nrow(record$counts) == raw$observed_genes[[i]] &&
    record$counts_sha256 == raw$counts_sha256[[i]] &&
    sha(path) == raw$private_cache_sha256[[i]] &&
    as.numeric(file.info(path)$size) == raw$private_cache_bytes[[i]]
}
names(raw_records) <- raw$sample_id
sct_ok <- logical(nrow(sct))
for (i in seq_len(nrow(sct))) {
  path <- file.path(private_root, "sct", sct$private_cache_file[[i]])
  record <- readRDS(path)
  valid <- tryCatch({mv05d0_validate_normalization_cache_record_v2(record); TRUE},
                    error = function(e) FALSE)
  matrix <- if (valid) mv05d0_sct_matrix_from_cache_v1(record) else NULL
  selected <- select_matched_cells(colnames(raw_records[[sct$sample_id[[i]]]]$counts),
    n = 384L, seed = sct$seed[[i]])
  sct_ok[[i]] <- valid && record$identity$sample_id == sct$sample_id[[i]] &&
    record$identity$seed == sct$seed[[i]] && ncol(matrix) == 384L &&
    all(is.finite(matrix@x)) &&
    attr(selected, "selected_cell_sha256") == sct$selected_cell_sha256[[i]] &&
    record$payload_sha256 == sct$payload_sha256[[i]] &&
    sha(path) == sct$private_cache_sha256[[i]] &&
    as.numeric(file.info(path)$size) == sct$private_cache_bytes[[i]]
  rm(record, matrix, selected); invisible(gc())
}
selection_ok <- nrow(selection) == 170L &&
  identical(selection$sample_id, sct$sample_id) && identical(selection$seed, sct$seed) &&
  identical(selection$selected_cell_sha256, sct$selected_cell_sha256) &&
  all(selection$selected_cells == 384L) &&
  !any(truth(selection$biological_outcomes_computed))
paths <- c(file.path(private_root, "raw", raw$private_cache_file),
           file.path(private_root, "sct", sct$private_cache_file))
actual_manifest_sha <- vapply(paths, sha, character(1L))
actual_manifest_bytes <- as.numeric(file.info(paths)$size)
manifest_ok <- nrow(manifest) == 204L && all(file.exists(paths)) &&
  mv07f_manifest_matches_v1(
    manifest$private_cache_sha256, actual_manifest_sha,
    manifest$private_cache_bytes, actual_manifest_bytes)
partial <- c(setdiff(list.files(file.path(private_root, "raw")), raw$private_cache_file),
             setdiff(list.files(file.path(private_root, "sct")), sct$private_cache_file))
resource_ok <- nrow(resource) == 1L && truth(resource$resource_gate_passed) &&
  resource$total_worker_seconds <= resource$aggregate_elapsed_cap_seconds &&
  resource$peak_process_tree_rss_bytes <= resource$rss_cap_bytes &&
  resource$unique_cache_storage_bytes <= resource$storage_cap_bytes
summary_ok <- nrow(summary) == 1L && summary$raw_jobs == 34L &&
  summary$sct_jobs == 170L && summary$selections == 170L &&
  summary$partial_files == 0L && !truth(summary$label_access) &&
  all(summary[c("panel_fit_jobs", "pca_jobs", "ph_jobs", "landscape_jobs",
    "outcome_jobs")] == 0L)
provenance_ok <- nrow(provenance) == 1L && truth(provenance$receipt_before_source_parse) &&
  !truth(provenance$label_access) && !truth(provenance$panel_fit) &&
  !truth(provenance$pca) && !truth(provenance$ph) &&
  !truth(provenance$landscape) && !truth(provenance$outcomes)
queue_ok <- nrow(queue) == 204L && sum(queue$stage == "raw") == 34L &&
  sum(queue$stage == "sct") == 170L
checks <- data.frame(
  contract_id = "mv07f_upstream_independent_validation_v2",
  category = c("queue_axis", "raw_cache_content", "sct_cache_content",
    "selection_identity", "cache_manifest", "partial_state", "resource_gate",
    "production_summary", "execution_provenance"),
  passed = c(axis_ok && queue_ok, all(raw_ok), all(sct_ok), selection_ok,
    manifest_ok, !length(partial), resource_ok, summary_ok, provenance_ok),
  detail = c("204 fixed jobs", "34 of 34 raw records", "170 of 170 SCT records",
    "170 deterministic 384-cell hashes", "204 file hashes and byte sizes",
    "zero unexpected cache files", "child and aggregate caps pass",
    "zero panel PCA PH landscape outcome jobs", "receipt before source parse"),
  stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV7-F production validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
if (file.exists(output)) stop("Refusing overwrite: ", output)
write.table(checks, output, sep = ",", row.names = FALSE, col.names = TRUE,
            quote = TRUE, na = "")
message("MV7-F independent production validation: 9/9 categories pass")
