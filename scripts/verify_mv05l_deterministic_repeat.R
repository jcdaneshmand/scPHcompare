#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: verify_mv05l_deterministic_repeat.R PUBLIC_DIR REPEAT_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
public_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")

files <- c(
  "mv05l-input-lock-2026-08-10.csv",
  "mv05l-method-map-2026-08-10.csv",
  "mv05l-endpoint-compatibility-2026-08-10.csv",
  "mv05l-paired-query-endpoints-2026-08-10.csv",
  "mv05l-pseudobulk-identity-control-2026-08-10.csv",
  "mv05l-sample-estimands-2026-08-10.csv",
  "mv05l-tissue-estimands-2026-08-10.csv",
  "mv05l-macro-estimands-2026-08-10.csv",
  "mv05l-estimand-intervals-2026-08-10.csv",
  "mv05l-primary-contrasts-2026-08-10.csv",
  "mv05l-bootstrap-audit-2026-08-10.csv",
  "mv05l-randomization-audit-2026-08-10.csv",
  "mv05l-production-summary-2026-08-10.csv"
)
public_paths <- file.path(public_dir, files)
repeat_paths <- file.path(repeat_dir, files)
if (any(!file.exists(public_paths)) || any(!file.exists(repeat_paths))) {
  stop("One or more MV5-L deterministic-repeat artifacts are absent.",
       call. = FALSE)
}
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
public_sha <- vapply(public_paths, file_sha, character(1L))
repeat_sha <- vapply(repeat_paths, file_sha, character(1L))
result <- data.frame(
  contract_id = "mv05l_deterministic_repeat_v1",
  artifact = files,
  public_sha256 = unname(public_sha),
  repeat_sha256 = unname(repeat_sha),
  byte_identical = unname(public_sha == repeat_sha),
  public_size_bytes = unname(file.info(public_paths)$size),
  repeat_size_bytes = unname(file.info(repeat_paths)$size),
  endpoint_recomputations = 0L, reranking_operations = 0L,
  method_tuning_operations = 0L, method_selection_operations = 0L,
  tissue_selection_operations = 0L, clustering_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, optimization_jobs_executed = 0L,
  source_evaluations_modified = 0L,
  stringsAsFactors = FALSE
)
if (!all(result$byte_identical) ||
    any(result$public_size_bytes != result$repeat_size_bytes)) {
  stop("MV5-L production outputs are not byte deterministic.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("MV5-L deterministic repeat validation passed.")
