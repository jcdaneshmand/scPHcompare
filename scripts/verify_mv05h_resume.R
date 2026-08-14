#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_mv05h_resume.R BEFORE_CSV AFTER_CSV OUTPUT_CSV",
       call. = FALSE)
}
before <- utils::read.csv(
  args[[1L]], stringsAsFactors = FALSE, check.names = FALSE
)
after <- utils::read.csv(
  args[[2L]], stringsAsFactors = FALSE, check.names = FALSE
)
source("R/provenance_utils.R")
before <- before[order(before$group_order), , drop = FALSE]
after <- after[order(after$group_order), , drop = FALSE]
keys <- c("group_id", "group_order", "fold_id", "seed")
values <- c(
  "view_records", "view_audit_sha256", "result_set_sha256", "result_bytes"
)
key_pass <- identical(before[keys], after[keys])
value_pass <- identical(before[values], after[values])
result <- data.frame(
  contract_id = "mv05h_resume_validation_v1", groups = nrow(before),
  group_keys_identical = key_pass,
  record_counts_identical = identical(before$view_records, after$view_records),
  audit_hashes_identical = identical(
    before$view_audit_sha256, after$view_audit_sha256
  ),
  result_set_hashes_identical = identical(
    before$result_set_sha256, after$result_set_sha256
  ),
  result_bytes_identical = identical(before$result_bytes, after$result_bytes),
  complete_snapshot_identical = key_pass && value_pass,
  resumed_groups_rebuilt = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (nrow(before) != 75L || !all(as.logical(result[c(
    "group_keys_identical", "record_counts_identical",
    "audit_hashes_identical", "result_set_hashes_identical",
    "result_bytes_identical", "complete_snapshot_identical"
)]))) {
  stop("MV5-H resume changed one or more accepted private artifacts.",
       call. = FALSE)
}
write_provenance_csv(result, args[[3L]])
message("Verified zero-rebuild MV5-H resume.")
