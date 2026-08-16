#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv07h_prefreeze_repeat.R FIRST SECOND OUTPUT")
}
first <- args[[1L]]; second <- args[[2L]]; output <- args[[3L]]
files_first <- sort(list.files(first, recursive = TRUE, all.files = TRUE,
                               no.. = TRUE), method = "radix")
files_second <- sort(list.files(second, recursive = TRUE, all.files = TRUE,
                                no.. = TRUE), method = "radix")
if (!identical(files_first, files_second) || !length(files_first)) {
  stop("MV7-H prefreeze repeat file set differs.")
}
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
checks <- data.frame(
  contract_id = "mv07h_prefreeze_repeat_validation_v1", file = files_first,
  first_bytes = as.numeric(file.info(file.path(first, files_first))$size),
  second_bytes = as.numeric(file.info(file.path(second, files_second))$size),
  first_sha256 = vapply(file.path(first, files_first), sha, character(1L)),
  second_sha256 = vapply(file.path(second, files_second), sha, character(1L)),
  stringsAsFactors = FALSE)
checks$bytes_equal <- checks$first_bytes == checks$second_bytes
checks$sha256_equal <- checks$first_sha256 == checks$second_sha256
if (!all(checks$bytes_equal & checks$sha256_equal)) {
  stop("MV7-H prefreeze repeat differs.")
}
write.csv(checks, output, row.names = FALSE, na = "")
message("MV7-H prefreeze repeat: ", nrow(checks), "/", nrow(checks), " pass")
