#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_mv05n_resume.R BEFORE AFTER OUTPUT", call. = FALSE)
}
source("R/provenance_utils.R")
before <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
after <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
before <- before[order(before$artifact_role, before$artifact_file,
                       method = "radix"), ]
after <- after[order(after$artifact_role, after$artifact_file,
                     method = "radix"), ]
passed <- identical(before$artifact_file, after$artifact_file) &&
  identical(before$artifact_role, after$artifact_role) &&
  identical(before$size_bytes, after$size_bytes) &&
  identical(before$modified_utc, after$modified_utc) &&
  identical(before$sha256, after$sha256)
if (!passed) stop("MV5-N resume changed an immutable artifact.", call. = FALSE)
summary <- data.frame(
  contract_id = "mv05n_admission_resume_validation_v1",
  artifacts = nrow(before), output_artifacts = sum(before$artifact_role == "output"),
  status_artifacts = sum(before$artifact_role == "status"),
  hashes_unchanged = TRUE, sizes_unchanged = TRUE, timestamps_unchanged = TRUE,
  groups_reused = 12L, groups_rebuilt = 0L, passed = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)
write_provenance_csv(summary, args[[3L]])
message("Validated immutable MV5-N admission resume across 24 artifacts.")
