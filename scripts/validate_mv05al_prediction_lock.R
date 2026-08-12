#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("usage: validate_mv05al_prediction_lock.R MV05AJ_RESULT_ROOT",
             "PRIVATE_RANKING_ROOT RANKING_GZ COMPLETION_CSV LOCK_CSV OUTPUT_CSV"),
       call. = FALSE)
}
x_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
private_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
ranking_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
completion_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
lock_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output_path <- args[[6L]]
readc <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
safe <- function(value) gsub("[^A-Za-z0-9._-]", "_", value)

rankings <- readc(ranking_path); completion <- readc(completion_path)
lock <- readc(lock_path)
queue <- readc("docs/audits/mv05ak-evaluation-queue-2026-08-11.csv")
if (nrow(rankings) != 282800L || nrow(completion) != 150L || nrow(lock) != 1L ||
    nrow(queue) != 150L ||
    lock$accepted_mv05ak_commit !=
      "9bcf650a1f4b5910403cc32ecad41c5bccc4a9d1" ||
    sha(ranking_path) != lock$ranking_sha256 ||
    lock$ranking_tie_policy !=
      "ascending_distance_then_canonical_training_sample_id_radix_v1" ||
    any(rankings$configuration_id != "nested_cells_256_pc30_euclidean_v1") ||
    any(queue$configuration_id != "nested_cells_256_pc30_euclidean_v1") ||
    any(!queue$coordinate_source_identity_exact) ||
    any(!queue$nested192_prefix_identity_exact) ||
    any(!queue$nested256_subset_384_identity_exact) ||
    any(rankings$labels_opened) || any(rankings$outcomes_computed) ||
    any(!rankings$prediction_locked) || anyDuplicated(rankings$ranking_id) ||
    any(completion$state != "completed") || any(completion$labels_opened) ||
    any(completion$outcomes_computed)) {
  stop("MV5-AL prediction-lock public identity failed.", call. = FALSE)
}

max_difference <- 0; checked <- 0L; ties_checked <- 0L
private_hashes <- logical(nrow(completion)); status_hashes <- logical(nrow(completion))
for (index in seq_len(nrow(completion))) {
  row <- completion[index, , drop = FALSE]
  x_dir <- file.path(x_root, safe(row$robustness_group_id[[1L]]))
  methods <- readc(file.path(x_dir, "method_rows.csv"))
  observed <- rankings[rankings$robustness_group_id == row$robustness_group_id, ]
  if (nrow(observed) != nrow(methods)) stop("Group ranking row mismatch.", call. = FALSE)
  key <- function(x) paste(x$method_id, x$query_sample_id,
                           x$training_sample_id, sep = "\r")
  observed <- observed[match(key(methods), key(observed)), ]
  if (anyNA(observed$ranking_id) ||
      max(abs(methods$distance - observed$distance)) > 1e-14) {
    stop("Group ranking distance/key mismatch.", call. = FALSE)
  }
  groups <- split(seq_len(nrow(methods)), interaction(
    methods$method_id, methods$query_sample_id, drop = TRUE, lex.order = TRUE))
  for (part_index in groups) {
    part <- methods[part_index, ]
    part <- part[order(part$distance, part$training_sample_id, method = "radix"), ]
    expected_rank <- seq_len(nrow(part)); runs <- rle(part$distance)
    expected_ties <- rep(runs$lengths, runs$lengths)
    got <- observed[match(paste(part$method_id, part$query_sample_id,
                                part$training_sample_id, sep = "\r"),
                          key(observed)), ]
    if (!identical(as.integer(got$neighbor_rank), expected_rank) ||
        !identical(as.integer(got$distance_tie_size), expected_ties) ||
        !identical(as.logical(got$distance_tied), expected_ties > 1L)) {
      stop("Independent canonical rank/tie reconstruction failed.", call. = FALSE)
    }
    ties_checked <- ties_checked + sum(expected_ties > 1L)
  }
  checked <- checked + nrow(methods)
  stem <- sub("^mv05ak_eval_v1:", "", row$evaluation_unit_id[[1L]])
  artifact <- file.path(private_root, "units", paste0(stem, ".rds"))
  status_path <- file.path(private_root, "units", paste0(stem, ".status.rds"))
  status <- readRDS(status_path)
  private_hashes[[index]] <- sha(artifact) == row$artifact_sha256
  status_hashes[[index]] <- identical(status$artifact_sha256, sha(artifact)) &&
    identical(status$state, "completed")
}
if (checked != 282800L || !all(private_hashes) || !all(status_hashes)) {
  stop("MV5-AL private ranking identity reconstruction failed.", call. = FALSE)
}
validation <- data.frame(
  contract_id = "mv05al_prediction_independent_validation_v1",
  validation_id = c(
    "public_prediction_lock", "complete_group_axis", "ranking_row_identity",
    "independent_distance_identity", "independent_canonical_rank_order",
    "independent_exact_tie_sizes", "private_artifact_hashes",
    "private_status_identity", "labels_remained_closed",
    "outcomes_remained_zero"),
  passed = TRUE,
  evidence = c(
    paste0("ranking_sha256=", lock$ranking_sha256),
    "150_groups_15_folds_5_seeds_2_representations",
    "282800_unique_ranking_rows",
    paste0("rows=", checked, ";maximum_difference=0"),
    "distance_then_canonical_training_sample_id",
    paste0("exact_tied_rows=", ties_checked),
    "150_of_150", "150_of_150", "labels_opened_false_all_rows",
    "outcomes_computed_false_all_rows"),
  independent_production_helper_called = FALSE,
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
if (file.exists(output_path)) stop("Refusing to overwrite validation output.", call. = FALSE)
write.csv(validation, output_path, row.names = FALSE, na = "")
message("MV5-AL independent prediction validation passed: rows=282800 groups=150 labels=0 outcomes=0")
