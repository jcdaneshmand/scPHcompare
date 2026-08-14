#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: validate_mv05d3_group_resume.R MANIFEST GROUP_ID ",
    "FOLD_CACHE_DIR RESULT_DIR VIEW_AUDIT_CSV OUTPUT_CSV EXPECTED_VIEWS",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
group_id <- args[[2L]]
fold_cache_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
result_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
audit_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output_path <- args[[6L]]
expected_views <- as.integer(args[[7L]])
source("R/provenance_utils.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
files <- sort(list.files(result_dir, pattern = "[.]rds$", full.names = TRUE),
              method = "radix")
if (length(files) != expected_views) {
  stop("MV5-D3 resume validation found an incomplete result set.",
       call. = FALSE)
}
before <- vapply(c(audit_path, files), file_sha, character(1L))
started <- proc.time()[["elapsed"]]
output <- system2(
  Sys.which("Rscript"),
  shQuote(c(
    "--vanilla", "scripts/run_mv05d3_ph_group.R", manifest_path, group_id,
    fold_cache_dir, result_dir, audit_path, "validate_resume"
  )),
  stdout = TRUE, stderr = TRUE
)
elapsed <- proc.time()[["elapsed"]] - started
status <- attr(output, "status")
if (!is.null(status) && status != 0L) {
  stop("MV5-D3 resume-validation subprocess failed: ",
       paste(output, collapse = "\n"), call. = FALSE)
}
after <- vapply(c(audit_path, files), file_sha, character(1L))
audit <- utils::read.csv(
  audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
result <- data.frame(
  contract_id = "mv05d3_group_resume_validation_v1",
  group_id = group_id, expected_views = expected_views,
  validated_views = nrow(audit), resume_elapsed_seconds = elapsed,
  files_compared = length(before),
  immutable_file_hashes = identical(before, after),
  checkpoint_sha256 = file_sha(audit_path),
  result_set_sha256 = digest::digest(
    after[-1L], algo = "sha256", serialize = TRUE
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  landscape_jobs_executed = 0L, distance_jobs_executed = 0L,
  clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, stringsAsFactors = FALSE
)
if (nrow(audit) != expected_views || !result$immutable_file_hashes) {
  stop("MV5-D3 group resume was incomplete or mutated files.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("Validated immutable MV5-D3 group resume: ", group_id)
