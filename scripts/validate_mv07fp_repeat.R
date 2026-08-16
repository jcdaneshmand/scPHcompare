#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv07fp_repeat.R ORIGINAL REPEAT OUTPUT")
}
original <- args[[1]]
repeat_dir <- args[[2]]
output <- args[[3]]
files <- c(
  "mv07fp-source-summary.csv", "mv07fp-eligibility.csv",
  "mv07fp-panel.csv", "mv07fp-seed-stability.csv",
  "mv07fp-panel-comparison.csv", "mv07fp-decision.csv",
  "mv07fp-artifact-manifest.csv"
)
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
rows <- lapply(files, function(file) {
  first <- file.path(original, file)
  second <- file.path(repeat_dir, file)
  if (!file.exists(first) || !file.exists(second)) {
    stop("Missing repeat artifact: ", file)
  }
  first_sha <- sha(first)
  second_sha <- sha(second)
  first_bytes <- as.numeric(file.info(first)$size)
  second_bytes <- as.numeric(file.info(second)$size)
  data.frame(
    contract_id = "mv07fp_repeat_validation_v1",
    file = file,
    original_bytes = first_bytes,
    repeat_bytes = second_bytes,
    bytes_equal = first_bytes == second_bytes,
    original_sha256 = first_sha,
    repeat_sha256 = second_sha,
    sha256_equal = identical(first_sha, second_sha),
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
if (!all(result$bytes_equal) || !all(result$sha256_equal)) {
  stop("MV7-FP deterministic repeat validation failed.")
}
write.csv(result, output, row.names = FALSE, na = "")
message("MV7-FP deterministic repeat validation: 7/7 pass")
