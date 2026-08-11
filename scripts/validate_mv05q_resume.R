#!/usr/bin/env Rscript

options(warn = 2)
source("R/provenance_utils.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv05q_resume.R BEFORE_CSV AFTER_CSV OUTPUT_CSV",
       call. = FALSE)
}
read_value <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                              check.names = FALSE)
before <- read_value(args[[1L]])
after <- read_value(args[[2L]])
before <- before[order(before$relative_path, method = "radix"), ]
after <- after[order(after$relative_path, method = "radix"), ]
same_axis <- identical(before$relative_path, after$relative_path)
same_bytes <- same_axis && identical(before$bytes, after$bytes)
same_hash <- same_axis && identical(before$sha256, after$sha256)
same_time <- same_axis && identical(before$modified_unix, after$modified_unix)
result <- data.frame(
  contract_id = "mv05q_immutable_resume_validation_v1",
  artifacts_before = nrow(before), artifacts_after = nrow(after),
  identical_path_axis = same_axis, identical_bytes = same_bytes,
  identical_sha256 = same_hash, identical_modified_timestamps = same_time,
  rebuilt_artifacts = if (same_axis) sum(before$sha256 != after$sha256 |
    before$modified_unix != after$modified_unix) else NA_integer_,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  labels_opened = FALSE, stringsAsFactors = FALSE)
if (!same_axis || !same_bytes || !same_hash || !same_time ||
    result$rebuilt_artifacts != 0L) {
  stop("MV5-Q immutable resume validation failed.", call. = FALSE)
}
write_provenance_csv(result, args[[3L]])
message("MV5-Q immutable resume passed: artifacts=", nrow(before),
        " rebuilt=0.")
