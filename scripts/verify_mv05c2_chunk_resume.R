#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: verify_mv05c2_chunk_resume.R BEFORE_CSV OUTPUT_DIR STATUS_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
before_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
status_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_path <- args[[4L]]
source("R/provenance_utils.R")
before <- utils::read.csv(before_path, stringsAsFactors = FALSE)
paths <- c(
  list.files(output_dir, pattern = "\\.csv$", full.names = TRUE),
  list.files(status_dir, pattern = "\\.csv$", full.names = TRUE)
)
paths <- sort(paths, method = "radix")
after <- data.frame(
  artifact_file = basename(paths),
  artifact_role = ifelse(dirname(paths) == output_dir, "output", "status"),
  size_bytes_after_resume = unname(file.info(paths)$size),
  sha256_after_resume = vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)),
  stringsAsFactors = FALSE
)
keys <- paste(before$artifact_role, before$artifact_file, sep = "\r")
matched <- match(
  paste(after$artifact_role, after$artifact_file, sep = "\r"), keys
)
if (anyNA(matched) || nrow(after) != nrow(before)) {
  stop("Resume changed the chunk artifact set.", call. = FALSE)
}
result <- data.frame(
  contract_id = "mv05c2_chunk_resume_validation_v1",
  artifact_file = after$artifact_file,
  artifact_role = after$artifact_role,
  size_bytes_before_resume = before$size_bytes[matched],
  size_bytes_after_resume = after$size_bytes_after_resume,
  sha256_before_resume = before$sha256[matched],
  sha256_after_resume = after$sha256_after_resume,
  immutable_resume_passed = before$size_bytes[matched] ==
    after$size_bytes_after_resume &
    before$sha256[matched] == after$sha256_after_resume,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (!all(result$immutable_resume_passed)) {
  stop("Resume modified one or more immutable chunk artifacts.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("Validated immutable resume for ", nrow(result), " chunk artifacts.")
