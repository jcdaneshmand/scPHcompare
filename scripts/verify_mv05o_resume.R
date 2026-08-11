#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_mv05o_resume.R BEFORE_CSV AFTER_CSV OUTPUT_CSV",
       call. = FALSE)
}
source("R/provenance_utils.R")
read_one <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
before <- read_one(args[[1L]])
after <- read_one(args[[2L]])
columns <- c("artifact_file", "artifact_role", "size_bytes", "modified_utc", "sha256")
before <- before[order(before$artifact_role, before$artifact_file), columns]
after <- after[order(after$artifact_role, after$artifact_file), columns]
passed <- identical(before, after)
if (!passed) stop("MV5-O resume changed immutable artifacts.", call. = FALSE)
result <- data.frame(
  contract_id = "mv05o_immutable_resume_validation_v1",
  artifacts = nrow(before),
  output_artifacts = sum(before$artifact_role == "output"),
  status_artifacts = sum(before$artifact_role == "status"),
  hashes_unchanged = TRUE, sizes_unchanged = TRUE,
  timestamps_unchanged = TRUE, units_rebuilt = 0L, passed = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)
if (file.exists(args[[3L]])) stop("Refusing to overwrite resume validation.",
                                  call. = FALSE)
write_provenance_csv(result, args[[3L]])
message("Validated MV5-O immutable resume across ", nrow(before), " artifacts.")
