#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: verify_mv05j_public_repeat.R PUBLIC_DIR REPEAT_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
public_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")
files <- c(
  "mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz",
  "mv05j-method-completion-2026-08-09.csv",
  "mv05j-group-bundle-index-2026-08-09.csv",
  "mv05j-method-registry-2026-08-09.csv",
  "mv05j-component-scale-disposition-2026-08-09.csv",
  "mv05j-public-assembly-summary-2026-08-09.csv"
)
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
public_paths <- file.path(public_dir, files)
repeat_paths <- file.path(repeat_dir, files)
if (!all(file.exists(public_paths)) || !all(file.exists(repeat_paths))) {
  stop("One or more MV5-J repeat artifacts are absent.", call. = FALSE)
}
public_hash <- vapply(public_paths, file_sha, character(1L))
repeat_hash <- vapply(repeat_paths, file_sha, character(1L))
result <- data.frame(
  contract_id = "mv05j_public_repeat_validation_v1",
  artifact = files, public_sha256 = public_hash,
  repeat_sha256 = repeat_hash,
  byte_identical = public_hash == repeat_hash,
  public_size_bytes = unname(file.info(public_paths)$size),
  repeat_size_bytes = unname(file.info(repeat_paths)$size),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (!all(result$byte_identical) ||
    any(result$public_size_bytes != result$repeat_size_bytes)) {
  stop("MV5-J public artifact assembly is not byte-reproducible.",
       call. = FALSE)
}
write_provenance_csv(result, output_path)
message("MV5-J public artifact repeat validation passed.")
