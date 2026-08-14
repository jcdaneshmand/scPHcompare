#!/usr/bin/env Rscript

options(warn = 2)
source("R/provenance_utils.R")
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for MV5-Q repeat validation.", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: validate_mv05q_repeat.R PRODUCTION_ROOT REPEAT_ROOT OUTPUT_CSV",
       call. = FALSE)
}
production_root <- args[[1L]]
repeat_root <- args[[2L]]
output_path <- args[[3L]]
queue <- utils::read.csv("docs/audits/mv05q-analysis-queue-2026-08-10.csv",
                         stringsAsFactors = FALSE, check.names = FALSE)
maximum <- max(queue$training_samples)
repeat_fold <- min(queue$fold_id[queue$training_samples == maximum])
selected <- queue[queue$fold_id == repeat_fold, , drop = FALSE]
selected <- selected[order(selected$execution_order), , drop = FALSE]
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
files <- c("candidate-partitions.csv", "stability.csv",
           "selected-partitions.csv", "heldout-assignments.csv")
rows <- do.call(rbind, lapply(seq_len(nrow(selected)), function(index) {
  item <- selected[index, , drop = FALSE]
  do.call(rbind, lapply(files, function(name) {
    relative <- file.path("groups", safe_name(item$analysis_group_id), name)
    first <- file.path(production_root, relative)
    second <- file.path(repeat_root, relative)
    if (!file.exists(first) || !file.exists(second)) {
      stop("MV5-Q repeat artifact is missing.", call. = FALSE)
    }
    first_hash <- file_sha(first)
    second_hash <- file_sha(second)
    first_bytes <- file.info(first)$size
    second_bytes <- file.info(second)$size
    data.frame(
      contract_id = "mv05q_clean_repeat_validation_v1",
      fold_id = item$fold_id, representation = item$representation,
      distance_id = item$distance_id, analysis_group_id = item$analysis_group_id,
      artifact_name = name, production_sha256 = first_hash,
      repeat_sha256 = second_hash, production_bytes = first_bytes,
      repeat_bytes = second_bytes,
      byte_identical = identical(first_hash, second_hash) &&
        identical(first_bytes, second_bytes),
      training_samples = item$training_samples,
      selection_rule = "maximum_training_samples_then_canonical_fold_id_v1",
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      labels_opened = FALSE, stringsAsFactors = FALSE)
  }))
}))
if (nrow(selected) != 10L || nrow(rows) != 40L ||
    any(!rows$byte_identical) || unique(rows$training_samples) != maximum) {
  stop("MV5-Q clean repeat validation failed.", call. = FALSE)
}
write_provenance_csv(rows, output_path)
message("MV5-Q clean repeat passed: fold=", repeat_fold,
        " groups=10 artifacts=40 byte_identical=40.")
