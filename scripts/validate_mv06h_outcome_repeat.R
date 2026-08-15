#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(
  "usage: validate_mv06h_outcome_repeat.R FIRST_DIR SECOND_DIR OUTPUT",
  call. = FALSE)
first <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
second <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- args[[3L]]
if (file.exists(output)) stop("Refusing to overwrite MV6-H repeat evidence.",
                              call. = FALSE)
files <- c(
  "mv06h-bootstrap-audit.csv", "mv06h-decision.csv",
  "mv06h-label-open-receipt.csv", "mv06h-label-source-provenance.csv",
  "mv06h-method-intervals.csv", "mv06h-method-summaries.csv",
  "mv06h-prediction-postcheck.csv", "mv06h-primary-contrasts.csv",
  "mv06h-production-summary.csv", "mv06h-query-method-outcomes.csv",
  "mv06h-randomization-audit.csv", "mv06h-sample-method-summaries.csv",
  "mv06h-tissue-method-summaries.csv")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
a <- file.path(first, files); b <- file.path(second, files)
if (any(!file.exists(c(a, b)))) stop("MV6-H repeat corpus incomplete.",
                                     call. = FALSE)
result <- data.frame(
  contract_id = "mv06h_outcome_repeat_v1", artifact = files,
  first_sha256 = vapply(a, sha, character(1L)),
  second_sha256 = vapply(b, sha, character(1L)),
  first_bytes = as.numeric(file.info(a)$size),
  second_bytes = as.numeric(file.info(b)$size), stringsAsFactors = FALSE)
result$byte_identical <- result$first_sha256 == result$second_sha256 &
  result$first_bytes == result$second_bytes
if (!all(result$byte_identical)) stop("MV6-H outcome repeat drift.", call. = FALSE)
write.csv(result, output, row.names = FALSE, na = "")
message("MV6-H outcome repeat passed: 13/13 byte-identical.")
