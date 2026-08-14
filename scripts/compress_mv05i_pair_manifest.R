#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: compress_mv05i_pair_manifest.R INPUT_CSV OUTPUT_CSV_GZ",
       call. = FALSE)
}
input <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
if (file.exists(output)) {
  stop("Refusing to overwrite an existing compressed pair manifest.",
       call. = FALSE)
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
source_connection <- file(input, open = "rb")
target_connection <- gzfile(output, open = "wb", compression = 9)
on.exit({
  try(close(source_connection), silent = TRUE)
  try(close(target_connection), silent = TRUE)
}, add = TRUE)
repeat {
  block <- readBin(source_connection, what = "raw", n = 1024L * 1024L)
  if (!length(block)) break
  writeBin(block, target_connection)
}
close(source_connection)
close(target_connection)
message("Compressed MV5-I pair manifest without changing CSV content.")
