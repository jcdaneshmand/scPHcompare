#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: validate_mv05d0_matrix_cache_equivalence.R NEW_CACHE_DIR ",
    "LEGACY_CACHE_DIR OUTPUT_CSV", call. = FALSE
  )
}
new_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
legacy_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")
source("R/mv05_resource_safe_execution.R")

files <- sort(list.files(new_dir, pattern = "__sct\\.rds$"), method = "radix")
if (!length(files)) stop("No matrix-only SCT caches were found.", call. = FALSE)
rows <- list()
for (file in files) {
  new <- readRDS(file.path(new_dir, file))
  legacy <- readRDS(file.path(legacy_dir, file))
  observed <- mv05d0_sct_matrix_from_cache_v1(new)
  expected <- mv05d0_sct_matrix_from_cache_v1(legacy)
  if (!identical(dim(observed), dim(expected)) ||
      !identical(dimnames(observed), dimnames(expected))) {
    stop("Matrix-only and legacy cache axes differ: ", file, call. = FALSE)
  }
  difference <- observed - expected
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05d0_matrix_payload_equivalence_v1",
    sample_id = new$identity$sample_id, seed = new$identity$seed,
    genes = nrow(observed), cells = ncol(observed),
    maximum_absolute_difference = max(abs(difference)),
    exact_numeric_identity = identical(observed, expected),
    matrix_only_cache_bytes = unname(file.info(file.path(new_dir, file))$size),
    legacy_object_cache_bytes = unname(file.info(file.path(legacy_dir, file))$size),
    observed_sha256 = digest::digest(observed, algo = "sha256", serialize = TRUE),
    expected_sha256 = digest::digest(expected, algo = "sha256", serialize = TRUE),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
result <- do.call(rbind, rows)
write_provenance_csv(result, output_path)
if (any(result$maximum_absolute_difference != 0) ||
    !all(result$exact_numeric_identity)) {
  stop("Matrix-only cache is not exactly equivalent to legacy cache.",
       call. = FALSE)
}
message("Validated exact matrix-only equivalence for ", nrow(result),
        " cache entries.")
