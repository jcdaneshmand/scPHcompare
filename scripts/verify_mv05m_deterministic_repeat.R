#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_mv05m_deterministic_repeat.R PUBLIC_DIR REPEAT_DIR OUTPUT_CSV",
       call. = FALSE)
}
public_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")
files <- c(
  "mv05m-selection-criteria-2026-08-10.csv",
  "mv05m-candidate-scores-2026-08-10.csv",
  "mv05m-axis-readiness-2026-08-10.csv",
  "mv05m-selected-next-sprint-2026-08-10.csv",
  "mv05m-training-pair-scope-2026-08-10.csv",
  "mv05m-evidence-registry-2026-08-10.csv",
  "mv05m-production-summary-2026-08-10.csv"
)
sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
public_paths <- file.path(public_dir, files)
repeat_paths <- file.path(repeat_dir, files)
if (any(!file.exists(public_paths)) || any(!file.exists(repeat_paths))) {
  stop("MV5-M repeat artifact missing.", call. = FALSE)
}
public_sha <- vapply(public_paths, sha, character(1L))
repeat_sha <- vapply(repeat_paths, sha, character(1L))
result <- data.frame(
  contract_id = "mv05m_deterministic_repeat_v1",
  artifact = files,
  public_sha256 = unname(public_sha), repeat_sha256 = unname(repeat_sha),
  byte_identical = unname(public_sha == repeat_sha),
  biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, technical_mixing_jobs_executed = 0L,
  robustness_jobs_executed = 0L, gene_view_jobs_executed = 0L,
  fusion_jobs_executed = 0L, new_data_jobs_executed = 0L,
  optimization_jobs_executed = 0L,
  stringsAsFactors = FALSE
)
if (!all(result$byte_identical)) {
  stop("MV5-M audit artifacts are not byte deterministic.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("MV5-M deterministic repeat validation passed.")
