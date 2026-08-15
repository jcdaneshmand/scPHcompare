#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(
  "usage: validate_mv06h_prediction_prefreeze_repeat.R FIRST_DIR SECOND_DIR OUTPUT",
  call. = FALSE)
first <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
second <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output <- args[[3L]]
if (file.exists(output)) stop("Refusing to overwrite repeat evidence.", call. = FALSE)
files <- c(
  "mv06h-prediction-group-manifest.csv", "mv06h-source-manifest.csv",
  "mv06h-implementation-manifest.csv", "mv06h-method-registry.csv",
  "mv06h-endpoint-registry.csv", "mv06h-contrast-registry.csv",
  "mv06h-inference-registry.csv", "mv06h-label-firewall.csv",
  "mv06h-prediction-lock.csv")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
paths1 <- file.path(first, files); paths2 <- file.path(second, files)
if (any(!file.exists(c(paths1, paths2)))) stop("Repeat corpus incomplete.",
                                               call. = FALSE)
result <- data.frame(
  contract_id = "mv06h_prediction_prefreeze_repeat_v1", artifact = files,
  first_sha256 = vapply(paths1, sha, character(1L)),
  second_sha256 = vapply(paths2, sha, character(1L)),
  first_bytes = as.numeric(file.info(paths1)$size),
  second_bytes = as.numeric(file.info(paths2)$size),
  labels_opened = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE)
result$byte_identical <- result$first_sha256 == result$second_sha256 &
  result$first_bytes == result$second_bytes
if (!all(result$byte_identical)) stop("MV6-H prefreeze repeat drift.", call. = FALSE)
write.csv(result, output, row.names = FALSE, na = "")
message("MV6-H prediction prefreeze repeat passed: 9/9 byte-identical.")
