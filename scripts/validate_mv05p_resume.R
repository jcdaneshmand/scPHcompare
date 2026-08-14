#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv05p_resume.R BEFORE_SNAPSHOT AFTER_SNAPSHOT PUBLIC_OUTPUT",
       call. = FALSE)
}
source("R/provenance_utils.R")
read_one <- function(path) utils::read.csv(
  normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
  check.names = FALSE
)
before <- read_one(args[[1L]])
after <- read_one(args[[2L]])
public_output <- args[[3L]]
comparison_columns <- c(
  "unit_family", "production_group_id", "unit_id", "output_sha256",
  "output_size_bytes", "output_mtime_epoch_seconds", "status_sha256",
  "status_size_bytes", "status_mtime_epoch_seconds", "outcome_label_state",
  "biological_outcomes_computed", "clustering_jobs_executed"
)
before <- before[order(before$unit_id), comparison_columns, drop = FALSE]
after <- after[order(after$unit_id), comparison_columns, drop = FALSE]
unchanged <- nrow(before) == 4565L && nrow(after) == 4565L &&
  !anyDuplicated(before$unit_id) && !anyDuplicated(after$unit_id) &&
  isTRUE(all.equal(before, after, check.attributes = FALSE, tolerance = 0))
result <- data.frame(
  contract_id = "mv05p_immutable_resume_validation_v1",
  validation_id = "immutable_resume_all_units",
  required_units = 4565L, before_units = nrow(before), after_units = nrow(after),
  unchanged_units = if (unchanged) 4565L else NA_integer_,
  rebuilt_units = if (unchanged) 0L else NA_integer_,
  hash_size_timestamp_unchanged = unchanged,
  zero_rebuild_passed = unchanged,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)
if (!unchanged) {
  stop("MV5-P immutable all-unit resume validation failed.", call. = FALSE)
}
write_provenance_csv(result, public_output)
message("Validated immutable resume across all 4,565 MV5-P units.")
