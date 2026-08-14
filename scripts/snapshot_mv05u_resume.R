#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: snapshot_mv05u_resume.R RESULT_ROOT OUTPUT_CSV",
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-U resume snapshots.", call. = FALSE)
}
source("R/provenance_utils.R")
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
files <- sort(list.files(root, recursive = TRUE, full.names = TRUE),
              method = "radix")
info <- file.info(files)
files <- files[!info$isdir]
info <- file.info(files)
relative <- substring(normalizePath(files, winslash = "/"), nchar(root) + 2L)
if (length(files) != 240L || anyDuplicated(relative)) {
  stop("MV5-U resume snapshot requires exactly 240 private files.",
       call. = FALSE)
}
result <- data.frame(
  contract_id = "mv05u_resume_snapshot_v1",
  private_relative_path = relative,
  sha256 = vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L)),
  bytes = as.numeric(info$size),
  modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(result, args[[2L]])
message("Snapshotted 240 immutable MV5-U private files.")
