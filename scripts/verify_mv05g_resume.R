#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_mv05g_resume.R BEFORE_CSV AFTER_CSV OUTPUT_CSV",
       call. = FALSE)
}
before <- utils::read.csv(
  normalizePath(args[[1L]], winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE, check.names = FALSE
)
after <- utils::read.csv(
  normalizePath(args[[2L]], winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "group_id", "cache_key", "payload_sha256", "coordinate_set_sha256",
  "result_size_bytes", "result_file_sha256", "outcome_label_state",
  "biological_outcomes_computed"
)
pass <- is.data.frame(before) && is.data.frame(after) &&
  nrow(before) == 75L && nrow(after) == 75L &&
  all(required %in% names(before)) && all(required %in% names(after)) &&
  identical(before, after) && all(before$outcome_label_state == "closed") &&
  !any(as.logical(before$biological_outcomes_computed))
result <- data.frame(
  contract_id = "mv05g_resume_validation_v1", groups = nrow(before),
  cache_keys_identical = identical(before$cache_key, after$cache_key),
  payload_hashes_identical = identical(
    before$payload_sha256, after$payload_sha256
  ),
  coordinate_hashes_identical = identical(
    before$coordinate_set_sha256, after$coordinate_set_sha256
  ),
  file_hashes_identical = identical(
    before$result_file_sha256, after$result_file_sha256
  ),
  complete_snapshot_byte_equivalent = pass,
  resumed_groups_rebuilt = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
if (!pass) stop("MV5-G resume changed an immutable result.", call. = FALSE)
dir.create(dirname(args[[3L]]), recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
write_provenance_csv(result, args[[3L]])
message("MV5-G resume preserved all 75 private results exactly.")
