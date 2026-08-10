#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: snapshot_mv05n_resume.R OUTPUT_DIR STATUS_DIR SNAPSHOT_OUTPUT",
       call. = FALSE)
}
source("R/provenance_utils.R")
paths <- sort(c(
  list.files(normalizePath(args[[1L]], mustWork = TRUE), full.names = TRUE),
  list.files(normalizePath(args[[2L]], mustWork = TRUE), full.names = TRUE)
), method = "radix")
if (length(paths) != 24L) stop("MV5-N resume snapshot needs 24 files.",
                               call. = FALSE)
info <- file.info(paths)
snapshot <- data.frame(
  contract_id = "mv05n_admission_resume_snapshot_v1",
  artifact_file = basename(paths),
  artifact_role = ifelse(grepl("__status[.]csv$", paths), "status", "output"),
  size_bytes = unname(info$size),
  modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
  sha256 = vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L)),
  stringsAsFactors = FALSE
)
write_provenance_csv(snapshot, args[[3L]])
