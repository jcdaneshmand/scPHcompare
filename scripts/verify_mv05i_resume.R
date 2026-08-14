#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: verify_mv05i_resume.R BEFORE OUTPUT_DIR STATUS_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
status_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
source("R/provenance_utils.R")
paths <- sort(c(
  list.files(output_dir, pattern = "[.]csv$", recursive = TRUE,
             full.names = TRUE),
  list.files(status_dir, pattern = "__status[.]csv$", recursive = TRUE,
             full.names = TRUE)
), method = "radix")
roots <- ifelse(startsWith(paths, paste0(output_dir, "/")),
                output_dir, status_dir)
relative <- substring(paths, nchar(roots) + 2L)
after <- data.frame(
  artifact_file = relative,
  artifact_role = ifelse(roots == output_dir,
                         "distance_output", "chunk_status"),
  size_bytes = unname(file.info(paths)$size),
  sha256 = vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)), stringsAsFactors = FALSE
)
if (!file.exists(args[[1L]])) {
  write_provenance_csv(after, args[[1L]])
  message("Captured immutable MV5-I resume snapshot for ", nrow(after),
          " files.")
  quit(status = 0L, save = "no")
}
before <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
names(after)[names(after) == "size_bytes"] <- "size_bytes_after_resume"
names(after)[names(after) == "sha256"] <- "sha256_after_resume"
matched <- match(
  paste(after$artifact_role, after$artifact_file, sep = "\r"),
  paste(before$artifact_role, before$artifact_file, sep = "\r")
)
if (anyNA(matched) || nrow(after) != nrow(before)) {
  stop("MV5-I resume changed the artifact set.", call. = FALSE)
}
result <- data.frame(
  contract_id = "mv05i_chunk_resume_validation_v1",
  artifact_file = after$artifact_file, artifact_role = after$artifact_role,
  size_bytes_before_resume = before$size_bytes[matched],
  size_bytes_after_resume = after$size_bytes_after_resume,
  sha256_before_resume = before$sha256[matched],
  sha256_after_resume = after$sha256_after_resume,
  immutable_resume_passed = before$size_bytes[matched] ==
    after$size_bytes_after_resume &
    before$sha256[matched] == after$sha256_after_resume,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (!all(result$immutable_resume_passed)) {
  stop("MV5-I resume modified an immutable artifact.", call. = FALSE)
}
write_provenance_csv(result, args[[4L]])
message("Validated immutable MV5-I resume for ", nrow(result), " files.")
