#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: validate_mv08f_cache_recovery.R PREFREEZE_DIR PRIVATE_ROOT ",
       "RAW_METRICS SELECTION SCT_METRICS OUTPUT_DIR")
}
prefreeze <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
private_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
raw_metrics_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
selection_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
sct_metrics_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output <- args[[6L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(output, all.files = TRUE, no.. = TRUE))) {
  stop("MV8-F recovery output directory must be empty.")
}
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv05_resource_safe_execution.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) if (is.logical(value)) !is.na(value) & value else
  tolower(trimws(as.character(value))) == "true"
axis <- read.csv(file.path(prefreeze, "mv08f-recovery-axis.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
raw <- read.csv(raw_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
selection <- read.csv(selection_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
sct <- read.csv(sct_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
axis <- axis[order(axis$seed, axis$sample_id, method = "radix"),, drop = FALSE]
selection <- selection[order(selection$seed, selection$sample_id,
                             method = "radix"),, drop = FALSE]
sct <- sct[order(sct$seed, sct$sample_id, method = "radix"),, drop = FALSE]
raw <- raw[order(raw$sample_id, method = "radix"),, drop = FALSE]
if (nrow(axis) != 450L || nrow(raw) != 90L || nrow(selection) != 450L ||
    nrow(sct) != 450L || anyDuplicated(raw$sample_id) ||
    anyDuplicated(paste(sct$sample_id, sct$seed, sep = "\r")) ||
    !identical(selection$sample_id, axis$sample_id) ||
    !identical(as.integer(selection$seed), as.integer(axis$seed)) ||
    !identical(selection$selected_cell_sha256, axis$selected_cell_sha256) ||
    !identical(sct$sample_id, axis$sample_id) ||
    !identical(as.integer(sct$seed), as.integer(axis$seed))) {
  stop("MV8-F recovery axes are incomplete or stale.")
}
raw_paths <- file.path(private_root, "raw", raw$private_raw_cache_file)
sct_paths <- file.path(private_root, "sct", sct$private_cache_file)
if (any(!file.exists(c(raw_paths, sct_paths)))) {
  stop("One or more MV8-F private caches are absent.")
}
raw_ok <- logical(nrow(raw))
for (index in seq_len(nrow(raw))) {
  record <- readRDS(raw_paths[[index]])
  valid <- tryCatch({
    mv05d0_validate_raw_sample_cache_v2(record)
    TRUE
  }, error = function(error) FALSE)
  raw_ok[[index]] <- valid && record$sample_id == raw$sample_id[[index]] &&
    record$counts_sha256 == raw$counts_sha256[[index]] &&
    nrow(record$counts) == raw$genes[[index]] &&
    ncol(record$counts) == raw$cells[[index]] &&
    sha(raw_paths[[index]]) == raw$private_raw_cache_sha256[[index]] &&
    as.numeric(file.info(raw_paths[[index]])$size) ==
      raw$private_raw_cache_size_bytes[[index]]
  rm(record)
  if (index %% 10L == 0L) invisible(gc(FALSE))
}
sct_rows <- vector("list", nrow(sct))
for (index in seq_len(nrow(sct))) {
  record <- readRDS(sct_paths[[index]])
  valid <- tryCatch({
    mv05d0_validate_normalization_cache_record_v2(record)
    TRUE
  }, error = function(error) FALSE)
  matrix <- if (valid) mv05d0_sct_matrix_from_cache_v1(record) else NULL
  row <- axis[index,, drop = FALSE]
  passed <- valid && record$identity$sample_id == row$sample_id &&
    record$identity$seed == row$seed &&
    record$identity$selected_cell_sha256 == row$selected_cell_sha256 &&
    record$cache_key == row$normalization_cache_key &&
    record$payload$payload_contract_id == row$payload_contract_id &&
    record$payload_sha256 == row$payload_sha256 &&
    ncol(matrix) == 384L && all(is.finite(matrix@x)) &&
    as.numeric(file.info(sct_paths[[index]])$size) == row$private_cache_bytes &&
    sha(sct_paths[[index]]) == row$private_cache_sha256 &&
    sct$disposition[[index]] %in% c("built_atomic", "reuse_validated") &&
    sct$exit_status[[index]] == 0L &&
    sct$selected_cell_sha256[[index]] == row$selected_cell_sha256 &&
    sct$normalization_cache_key[[index]] == row$normalization_cache_key &&
    sct$payload_sha256[[index]] == row$payload_sha256 &&
    sct$private_cache_sha256[[index]] == row$private_cache_sha256
  sct_rows[[index]] <- data.frame(
    contract_id = "mv08f_recovered_cache_exact_validation_v1",
    sample_id = row$sample_id, seed = as.integer(row$seed),
    selected_cell_sha256 = row$selected_cell_sha256,
    normalization_cache_key = row$normalization_cache_key,
    payload_contract_id = row$payload_contract_id,
    payload_sha256 = row$payload_sha256,
    private_cache_file = row$private_cache_file,
    private_cache_bytes = as.numeric(row$private_cache_bytes),
    private_cache_sha256 = row$private_cache_sha256,
    sct_genes = if (is.null(matrix)) NA_integer_ else nrow(matrix),
    sct_cells = if (is.null(matrix)) NA_integer_ else ncol(matrix),
    finite_payload = !is.null(matrix) && all(is.finite(matrix@x)),
    exact_identity_passed = passed, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
  rm(record, matrix)
  if (index %% 10L == 0L) invisible(gc(FALSE))
}
validated <- do.call(rbind, sct_rows)
if (!all(raw_ok) || !all(validated$exact_identity_passed)) {
  stop("MV8-F recovered cache exact-identity validation failed.")
}
cache_paths <- c(raw_paths, sct_paths)
resource <- data.frame(
  contract_id = "mv08f_cache_recovery_resource_summary_v1",
  raw_jobs = nrow(raw), sct_jobs = nrow(sct),
  raw_worker_seconds = sum(raw$elapsed_seconds),
  sct_worker_seconds = sum(sct$elapsed_seconds),
  total_worker_seconds = sum(raw$elapsed_seconds) + sum(sct$elapsed_seconds),
  peak_process_tree_rss_bytes = max(c(raw$peak_process_tree_rss_bytes,
    sct$peak_process_tree_rss_bytes)),
  private_cache_storage_bytes = sum(as.numeric(file.info(cache_paths)$size)),
  child_elapsed_cap_seconds = 1800, child_rss_cap_bytes = 8 * 1024^3,
  aggregate_storage_cap_bytes = 40 * 1024^3,
  all_resource_caps_passed = all(raw$elapsed_seconds <= 1800) &&
    all(sct$elapsed_seconds <= 1800) &&
    max(c(raw$peak_process_tree_rss_bytes,
          sct$peak_process_tree_rss_bytes)) <= 8 * 1024^3 &&
    sum(as.numeric(file.info(cache_paths)$size)) <= 40 * 1024^3,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08f_cache_recovery_decision_v1",
  decision = if (all(raw_ok) && all(validated$exact_identity_passed) &&
    resource$all_resource_caps_passed) {
      "recovery_exact_authorize_475_source_prefreeze"
    } else "stop_recovery_failure",
  raw_exact = sum(raw_ok), sct_exact = sum(validated$exact_identity_passed),
  accepted_axis_records = nrow(axis), unexpected_cache_files = length(c(
    setdiff(list.files(file.path(private_root, "raw")), basename(raw_paths)),
    setdiff(list.files(file.path(private_root, "sct")), basename(sct_paths)))),
  panel475_source_jobs_authorized = 0L, ph_jobs_authorized = 0L,
  landscape_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  hca_fastq_download_authorized = FALSE, label_access_authorized = FALSE,
  outcome_jobs_authorized = 0L,
  next_gate = "MV8-G_475_source_PH_landscape_prefreeze",
  stringsAsFactors = FALSE)
if (decision$unexpected_cache_files != 0L ||
    decision$decision == "stop_recovery_failure") stop("MV8-F final gate failed.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
outputs <- list(
  "mv08f-raw-recovery.csv" = raw,
  "mv08f-recovered-cache-validation.csv" = validated,
  "mv08f-recovery-resource-summary.csv" = resource,
  "mv08f-recovery-decision.csv" = decision)
paths <- vapply(names(outputs), function(name) {
  path <- file.path(output, name)
  write_provenance_csv(outputs[[name]], path)
  path
}, character(1L))
manifest <- data.frame(
  contract_id = "mv08f_cache_recovery_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)),
  contains_expression = FALSE, contains_cell_barcode = FALSE,
  contains_absolute_private_path = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(manifest, file.path(output, "mv08f-artifact-manifest.csv"))
message("MV8-F recovery exact validation passed: 90 raw + 450 SCT caches")
