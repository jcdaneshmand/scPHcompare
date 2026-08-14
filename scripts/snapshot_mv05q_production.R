#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-Q snapshots.", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: snapshot_mv05q_production.R PRODUCTION_ROOT OUTPUT_CSV",
       call. = FALSE)
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
paths <- sort(list.files(root, recursive = TRUE, full.names = TRUE),
              method = "radix")
paths <- paths[file.info(paths)$isdir %in% FALSE]
info <- file.info(paths)
snapshot <- data.frame(
  relative_path = substring(paths, nchar(root) + 2L),
  bytes = as.numeric(info$size), modified_unix = as.numeric(info$mtime),
  sha256 = vapply(paths, digest::digest, character(1L), algo = "sha256",
                  serialize = FALSE), stringsAsFactors = FALSE)
utils::write.csv(snapshot, args[[2L]], row.names = FALSE, quote = TRUE,
                 na = "")
message("MV5-Q snapshot recorded ", nrow(snapshot), " artifacts.")
