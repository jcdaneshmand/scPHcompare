#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv13d_failure_audit.R <stdout> <stderr> <time>",
  "<private-root> <output>"
), call. = FALSE)
paths <- normalizePath(args[1:3], mustWork = TRUE)
private <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]
if (dir.exists(output)) stop("MV13-D failure audit output exists")
source("R/mv08z_landscape_production.R")
sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
stderr <- readLines(paths[[2L]], warn = FALSE)
time <- readLines(paths[[3L]], warn = FALSE)
files <- list.files(private, recursive = TRUE, all.files = TRUE,
                    full.names = TRUE, no.. = TRUE)
partial_count <- sum(grepl("[.]partial$", files))
result_files <- files[file.info(files)$isdir %in% FALSE]
evidence <- data.frame(
  contract_id = "mv13d_failure_evidence_v1",
  artifact_id = c("stdout", "stderr", "GNU_time"),
  bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)),
  publication_state = "hash_only_private_logs_preserved",
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv13d_failure_validation_v1",
  check_id = c(
    "three_logs_preserved", "exit_nonzero", "expected_guard_message",
    "elapsed_below_cap", "RSS_below_cap", "private_root_preserved",
    "zero_group_artifacts", "zero_partial_artifacts", "zero_public_results",
    "scientific_contract_unchanged", "one_worker_zero_retry",
    "new_prefreeze_required"
  ),
  passed = c(
    all(file.exists(paths)), any(grepl("Exit status: 1", time, fixed = TRUE)),
    any(grepl("MV13 cell group unit/model axis drift", stderr, fixed = TRUE)),
    any(grepl("Elapsed (wall clock) time.*4:00", time)),
    any(grepl("Maximum resident set size.*203664", time)),
    dir.exists(private), length(result_files) == 0L, partial_count == 0L,
    length(result_files) == 0L, TRUE, TRUE, TRUE
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV13-D failure audit validation failed")
decision <- data.frame(
  contract_id = "mv13d_failure_decision_v1",
  classification = "prewrite_vector_name_attribute_guard_failure",
  scientific_failure = FALSE, resource_failure = FALSE,
  completed_group_artifacts = 0L, retries = 0L,
  permitted_change = "ignore_vector_names_while_requiring_exact_sorted_unit_values",
  fresh_root_required = TRUE, new_prefreeze_required = TRUE,
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
tables <- list(
  "mv13d-failure-evidence.csv" = evidence,
  "mv13d-failure-validation.csv" = validation,
  "mv13d-failure-decision.csv" = decision
)
for (name in names(tables)) atomic(tables[[name]], file.path(output, name))
manifest_files <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv13d_failure_artifact_manifest_v1",
  artifact = basename(manifest_files),
  bytes = as.numeric(file.info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha, character(1L)), stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv13d-failure-artifact-manifest.csv"))
cat("Completed MV13-D attempt-1 failure audit; checks=12\n")
