#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: snapshot_mv05c2_chunk_artifacts.R OUTPUT_DIR STATUS_DIR OUTPUT_CSV",
       call. = FALSE)
}
output_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
status_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]
source("R/provenance_utils.R")
paths <- c(
  list.files(output_dir, pattern = "\\.csv$", full.names = TRUE),
  list.files(status_dir, pattern = "\\.csv$", full.names = TRUE)
)
paths <- sort(paths, method = "radix")
if (!length(paths)) stop("No chunk artifacts exist to snapshot.")
result <- data.frame(
  artifact_file = basename(paths),
  artifact_role = ifelse(dirname(paths) == output_dir, "output", "status"),
  size_bytes = unname(file.info(paths)$size),
  sha256 = vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)),
  stringsAsFactors = FALSE
)
write_provenance_csv(result, output_path)
