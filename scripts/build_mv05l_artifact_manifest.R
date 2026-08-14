#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05l_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")

decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-022-MV05L-LOCKED-REPRESENTATION-COMPARISON.md",
  "docs/audits/MV05L_LOCKED_REPRESENTATION_COMPARISON_2026-08-10.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/plans/phase-4-primary-benchmark.md",
  "docs/specifications/MV05L_LOCKED_REPRESENTATION_COMPARISON_SPECIFICATION_V1.md"
)
implementation <- c(
  "R/mv05l_representation_comparison.R",
  "scripts/run_mv05l_representation_comparison.R",
  "scripts/validate_mv05l_representation_comparison.R",
  "scripts/verify_mv05l_deterministic_repeat.R",
  "scripts/build_mv05l_artifact_manifest.R",
  "tests/testthat/test-mv05l-locked-representation-comparison.R"
)
evidence <- file.path("docs/audits", c(
  "mv05l-bootstrap-audit-2026-08-10.csv",
  "mv05l-deterministic-repeat-2026-08-10.csv",
  "mv05l-endpoint-compatibility-2026-08-10.csv",
  "mv05l-estimand-intervals-2026-08-10.csv",
  "mv05l-independent-validation-2026-08-10.csv",
  "mv05l-input-lock-2026-08-10.csv",
  "mv05l-macro-estimands-2026-08-10.csv",
  "mv05l-method-map-2026-08-10.csv",
  "mv05l-paired-query-endpoints-2026-08-10.csv",
  "mv05l-primary-contrasts-2026-08-10.csv",
  "mv05l-production-summary-2026-08-10.csv",
  "mv05l-pseudobulk-identity-control-2026-08-10.csv",
  "mv05l-randomization-audit-2026-08-10.csv",
  "mv05l-sample-estimands-2026-08-10.csv",
  "mv05l-tissue-estimands-2026-08-10.csv"
))
upstream <- file.path("docs/audits", c(
  "mv05e-query-endpoints-2026-08-08.csv",
  "mv05e-prediction-lock-2026-08-08.csv",
  "mv05e-independent-validation-2026-08-08.csv",
  "mv05e-deterministic-repeat-2026-08-08.csv",
  "mv05e-artifact-manifest-2026-08-08.csv",
  "mv05k-query-endpoints-2026-08-10.csv",
  "mv05k-prediction-lock-2026-08-10.csv",
  "mv05k-independent-validation-2026-08-10.csv",
  "mv05k-deterministic-repeat-2026-08-10.csv",
  "mv05k-artifact-manifest-2026-08-10.csv"
))
artifacts <- data.frame(
  artifact_path = c(decision, implementation, evidence, upstream),
  artifact_role = c(
    rep("decision_document", length(decision)),
    rep("implementation_or_test", length(implementation)),
    rep("comparison_evidence", length(evidence)),
    rep("locked_source_evaluation", length(upstream))
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop("Missing MV5-L artifact(s): ", paste(
    artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
    collapse = ", "), call. = FALSE)
}
artifacts$contract_id <- "mv05l_public_artifact_manifest_v1"
artifacts$size_bytes <- unname(file.info(artifacts$artifact_path)$size)
artifacts$sha256 <- vapply(artifacts$artifact_path, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1L))
artifacts$outcome_state <- ifelse(
  artifacts$artifact_role == "locked_source_evaluation",
  "accepted_immutable_source_evaluation",
  ifelse(artifacts$artifact_role == "implementation_or_test",
         "not_applicable_code",
         "locked_result_comparison_not_fully_outcome_blind")
)
artifacts$source_evaluations_modified <- FALSE
artifacts$public_safe <- TRUE
artifacts <- artifacts[, c(
  "contract_id", "artifact_path", "artifact_role", "size_bytes", "sha256",
  "outcome_state", "source_evaluations_modified", "public_safe"
)]
write_provenance_csv(artifacts, args[[1L]])
message("Wrote MV5-L public artifact manifest: ", args[[1L]])
