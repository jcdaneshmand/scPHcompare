#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: validate_mv05i_compressed_manifest.R INPUT_CSV INPUT_CSV_GZ OUTPUT",
    call. = FALSE
  )
}
source("R/provenance_utils.R")
source_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
gzip_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
temporary <- tempfile(fileext = ".csv")
connection <- gzfile(gzip_path, open = "rb")
on.exit({
  try(close(connection), silent = TRUE)
  unlink(temporary)
}, add = TRUE)
output <- file(temporary, open = "wb")
repeat {
  block <- readBin(connection, what = "raw", n = 1024L * 1024L)
  if (!length(block)) break
  writeBin(block, output)
}
close(connection)
close(output)
source_rows <- nrow(utils::read.csv(
  source_path, stringsAsFactors = FALSE, check.names = FALSE
))
decompressed_rows <- nrow(utils::read.csv(
  temporary, stringsAsFactors = FALSE, check.names = FALSE
))
result <- data.frame(
  contract_id = "mv05i_pair_manifest_compression_validation_v1",
  source_rows = source_rows, decompressed_rows = decompressed_rows,
  source_size_bytes = unname(file.info(source_path)$size),
  gzip_size_bytes = unname(file.info(gzip_path)$size),
  decompressed_size_bytes = unname(file.info(temporary)$size),
  source_sha256 = sha(source_path),
  gzip_sha256 = sha(gzip_path),
  decompressed_sha256 = sha(temporary),
  byte_identity_passed = sha(source_path) == sha(temporary),
  row_identity_passed = source_rows == 70700L &&
    decompressed_rows == source_rows,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (!result$byte_identity_passed || !result$row_identity_passed) {
  stop("Compressed MV5-I pair manifest does not reproduce source bytes.",
       call. = FALSE)
}
write_provenance_csv(result, args[[3L]])
message("Validated byte-identical MV5-I compressed pair manifest.")
