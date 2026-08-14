#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: snapshot_mv05o_artifacts.R OUTPUT_DIR STATUS_DIR SNAPSHOT_CSV",
       call. = FALSE)
}
source("R/provenance_utils.R")
paths <- sort(c(
  list.files(normalizePath(args[[1L]], mustWork = TRUE), full.names = TRUE),
  list.files(normalizePath(args[[2L]], mustWork = TRUE), full.names = TRUE)
), method = "radix")
if (!length(paths) || any(!file.exists(paths))) {
  stop("MV5-O snapshot requires completed artifacts.", call. = FALSE)
}
info <- file.info(paths)
snapshot <- data.frame(
  contract_id = "mv05o_immutable_artifact_snapshot_v1",
  artifact_file = basename(paths),
  artifact_role = ifelse(grepl("__status[.]csv$", paths), "status", "output"),
  size_bytes = unname(info$size),
  modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
  sha256 = vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L)),
  stringsAsFactors = FALSE
)
if (file.exists(args[[3L]])) stop("Refusing to overwrite MV5-O snapshot.",
                                  call. = FALSE)
write_provenance_csv(snapshot, args[[3L]])
