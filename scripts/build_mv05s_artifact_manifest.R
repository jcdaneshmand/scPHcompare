#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05s_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
decision <- c(
  "docs/PROJECT_EVIDENCE.md",
  "docs/audits/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_2026-08-10.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/specifications/MV05S_PREDICTION_LOCKED_CLUSTERING_OUTCOME_EXECUTION_SPECIFICATION_V1.md")
implementation <- c(
  "R/mv05s_outcome_execution.R",
  "scripts/build_mv05s_execution_prefreeze.R",
  "scripts/run_mv05s_outcome_execution.R",
  "scripts/validate_mv05s_outcome_execution.R",
  "scripts/build_mv05s_artifact_manifest.R",
  "tests/testthat/test-mv05s-outcome-execution.R")
evidence <- file.path("docs/audits", c(
  "mv05s-bootstrap-audit-2026-08-10.csv",
  "mv05s-execution-prefreeze-audit-2026-08-10.csv",
  "mv05s-execution-source-freeze-2026-08-10.csv",
  "mv05s-first-pass-resource-summary-2026-08-10.csv",
  "mv05s-first-pass-unit-completion-2026-08-10.csv",
  "mv05s-heldout-macro-outcomes-2026-08-10.csv",
  "mv05s-heldout-seed-outcomes-2026-08-10.csv",
  "mv05s-heldout-unit-summary-2026-08-10.csv",
  "mv05s-independent-validation-2026-08-10.csv",
  "mv05s-resource-summary-2026-08-10.csv",
  "mv05s-source-validation-2026-08-10.csv",
  "mv05s-training-context-summary-2026-08-10.csv",
  "mv05s-training-seed-outcomes-2026-08-10.csv",
  "mv05s-training-unit-summary-2026-08-10.csv",
  "mv05s-unit-completion-2026-08-10.csv"))
upstream <- file.path("docs/audits", c(
  "mv05r-source-freeze-2026-08-10.csv",
  "mv05r-evaluation-queue-2026-08-10.csv",
  "mv05q-selected-training-partitions-2026-08-10.csv.gz",
  "mv05q-heldout-assignments-2026-08-10.csv.gz"))
artifacts <- data.frame(
  artifact_path = c(decision, implementation, evidence, upstream),
  artifact_role = c(rep("decision_document", length(decision)),
                    rep("implementation_or_test", length(implementation)),
                    rep("outcome_evidence", length(evidence)),
                    rep("accepted_immutable_source", length(upstream))),
  stringsAsFactors = FALSE)
if (any(!file.exists(artifacts$artifact_path))) {
  stop("Missing MV5-S artifact(s): ", paste(
    artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
    collapse = ", "), call. = FALSE)
}
artifacts$contract_id <- "mv05s_public_artifact_manifest_v1"
artifacts$size_bytes <- unname(file.info(artifacts$artifact_path)$size)
artifacts$sha256 <- vapply(artifacts$artifact_path, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1L))
artifacts$outcome_state <- ifelse(
  artifacts$artifact_role == "accepted_immutable_source",
  "accepted_immutable_source",
  ifelse(artifacts$artifact_role == "implementation_or_test",
         "prospectively_frozen_or_postexecution_manifest_code",
         "complete_secondary_outcome_reporting"))
artifacts$method_selection_executed <- FALSE
artifacts$public_safe <- TRUE
artifacts <- artifacts[c(
  "contract_id", "artifact_path", "artifact_role", "size_bytes", "sha256",
  "outcome_state", "method_selection_executed", "public_safe")]
write_provenance_csv(artifacts, args[[1L]])
message("Wrote MV5-S public artifact manifest: ", args[[1L]])
