#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: snapshot_mv05af_results.R RESULT_ROOT OUTPUT_CSV",
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-AF snapshots.", call. = FALSE)
}
source("R/provenance_utils.R")
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
if (file.exists(output)) stop("Refusing to overwrite: ", output, call. = FALSE)
groups <- list.dirs(root, recursive = FALSE, full.names = TRUE)
if (length(groups) != 150L || any(grepl("^\\.", basename(groups)))) {
  stop("MV5-AF snapshot requires exactly 150 groups.", call. = FALSE)
}
files <- sort(list.files(root, recursive = TRUE, full.names = TRUE,
                         all.files = TRUE, no.. = TRUE), method = "radix")
info <- file.info(files)
files <- files[!info$isdir]
info <- info[!info$isdir, , drop = FALSE]
if (length(files) != 1650L) {
  stop("MV5-AF snapshot requires exactly 1,650 files.", call. = FALSE)
}
result <- data.frame(
  contract_id = "mv05af_result_snapshot_v1",
  relative_path = substring(files, nchar(root) + 2L),
  bytes = as.numeric(info$size),
  sha256 = vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE), character(1L)),
  modification_time_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
  labels_opened = FALSE, rankings_computed = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(result, output)
message("Snapshotted 1,650 MV5-AF files with timestamps.")


