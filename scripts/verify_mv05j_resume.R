#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: verify_mv05j_resume.R BUNDLE_DIR AUDIT_DIR SNAPSHOT_CSV ",
    "SUMMARY_CSV", call. = FALSE
  )
}
bundle_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
audit_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
snapshot_path <- args[[3L]]
summary_path <- args[[4L]]
source("R/provenance_utils.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
bundle_paths <- sort(
  list.files(bundle_dir, pattern = "[.]rds$", full.names = TRUE),
  method = "radix"
)
audit_paths <- sort(
  list.files(audit_dir, pattern = "[.]csv$", full.names = TRUE),
  method = "radix"
)
paths <- c(bundle_paths, audit_paths)
if (length(bundle_paths) != 75L || length(audit_paths) != 75L) {
  stop("MV5-J resume verification requires 75 bundles and 75 audits.",
       call. = FALSE)
}
current <- data.frame(
  relative_file = c(
    paste0("bundles/", basename(bundle_paths)),
    paste0("audits/", basename(audit_paths))
  ),
  sha256 = vapply(paths, file_sha, character(1L)),
  size_bytes = unname(file.info(paths)$size), stringsAsFactors = FALSE
)
if (!file.exists(snapshot_path)) {
  write_provenance_csv(current, snapshot_path)
  message("Captured the pre-resume MV5-J file snapshot.")
} else {
  before <- utils::read.csv(
    snapshot_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  unchanged <- nrow(before) == nrow(current) &&
    identical(before$relative_file, current$relative_file) &&
    identical(before$sha256, current$sha256) &&
    identical(as.numeric(before$size_bytes), as.numeric(current$size_bytes))
  result <- data.frame(
    contract_id = "mv05j_resume_validation_v1",
    tracked_private_files = nrow(current), bundle_files = 75L,
    audit_files = 75L, unchanged_files = sum(
      before$relative_file == current$relative_file &
        before$sha256 == current$sha256 &
        before$size_bytes == current$size_bytes
    ),
    exact_snapshot_unchanged = unchanged,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  if (!unchanged) stop("MV5-J resume changed a completed artifact.",
                       call. = FALSE)
  write_provenance_csv(result, summary_path)
  message("MV5-J immutable resume validation passed.")
}
