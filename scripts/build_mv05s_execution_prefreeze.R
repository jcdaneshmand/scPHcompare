#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-S execution prefreeze.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv05s_execution_prefreeze.R AUDIT_DIR EXPECTED_HEAD",
       call. = FALSE)
}
audit_dir <- args[[1L]]
expected_head <- args[[2L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-S execution prefreeze must run at exact prospective-engine HEAD.",
       call. = FALSE)
}
paths <- c(
  mv05r_source_freeze = "docs/audits/mv05r-source-freeze-2026-08-10.csv",
  mv05r_queue = "docs/audits/mv05r-evaluation-queue-2026-08-10.csv",
  mv05q_selected = "docs/audits/mv05q-selected-training-partitions-2026-08-10.csv.gz",
  mv05q_heldout = "docs/audits/mv05q-heldout-assignments-2026-08-10.csv.gz",
  mv05s_spec = "docs/specifications/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_SPECIFICATION_V1.md",
  mv05s_code = "R/mv05s_outcome_execution.R",
  mv05s_runner = "scripts/run_mv05s_outcome_execution.R",
  mv05s_validator = "scripts/validate_mv05s_outcome_execution.R",
  mv05s_prefreeze_builder = "scripts/build_mv05s_execution_prefreeze.R",
  mv05s_tests = "tests/testthat/test-mv05s-outcome-execution.R")
if (any(!file.exists(paths))) stop("MV5-S execution source is missing.",
                                   call. = FALSE)
sha <- vapply(paths, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1L))
freeze <- data.frame(
  contract_id = "mv05s_execution_source_freeze_v1",
  source_id = names(paths), artifact_locator = unname(paths), sha256 = sha,
  bytes = as.numeric(file.info(paths)$size), prospective_engine_head = head,
  external = FALSE, outcomes_computed = FALSE, evaluation_executed = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE)
freeze$execution_freeze_sha256 <- .mv05s_digest(
  paste(freeze$artifact_locator, freeze$sha256, sep = "\r"))

queue <- utils::read.csv(paths[["mv05r_queue"]], stringsAsFactors = FALSE,
                         check.names = FALSE)
mv05r_validate_evaluation_queue_v1(queue)
audit <- data.frame(
  contract_id = "mv05s_execution_prefreeze_audit_v1",
  prospective_engine_head = head,
  frozen_sources = nrow(freeze), evaluation_units = nrow(queue),
  analysis_groups = length(unique(queue$analysis_group_id)),
  algorithms = length(unique(queue$algorithm_id)),
  endpoints = length(unique(queue$endpoint_id)),
  training_units = sum(queue$evaluation_scope ==
                         "overlapping_training_partition_alignment"),
  heldout_units = sum(queue$evaluation_scope ==
                        "heldout_label_prediction_from_frozen_training_cluster"),
  nmi_variant = "max", bootstrap_replicates = 2000L,
  bootstrap_seed = 20260810L, worker_limit = 1L,
  outcomes_computed = FALSE, evaluation_executed = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE)
if (audit$evaluation_units != 2400L || audit$analysis_groups != 150L ||
    audit$algorithms != 2L || audit$endpoints != 8L ||
    audit$training_units != 1800L || audit$heldout_units != 600L) {
  stop("MV5-S prospective execution axes drifted.", call. = FALSE)
}
write_provenance_csv(freeze, file.path(
  audit_dir, "mv05s-execution-source-freeze-2026-08-10.csv"))
write_provenance_csv(audit, file.path(
  audit_dir, "mv05s-execution-prefreeze-audit-2026-08-10.csv"))
message("MV5-S execution prefreeze passed: sources=", nrow(freeze),
        " units=", nrow(queue), " training=", audit$training_units,
        " heldout=", audit$heldout_units,
        " outcomes=0 evaluation=0 selection=0")
