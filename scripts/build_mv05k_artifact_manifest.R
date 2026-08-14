#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05k_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")

decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-021-MV05K-PREDICTION-LOCKED-INTEGRATED-RETRIEVAL-EVALUATION.md",
  "docs/audits/MV05K_PREDICTION_LOCKED_INTEGRATED_RETRIEVAL_EVALUATION_2026-08-10.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/plans/phase-4-primary-benchmark.md",
  "docs/specifications/MV05K_PREDICTION_LOCKED_INTEGRATED_RETRIEVAL_EVALUATION_SPECIFICATION_V1.md"
)
implementation <- c(
  "R/mv05k_retrieval_evaluation.R",
  "scripts/run_mv05k_retrieval_evaluation.R",
  "scripts/validate_mv05k_retrieval_evaluation.R",
  "scripts/verify_mv05k_deterministic_repeat.R",
  "scripts/build_mv05k_artifact_manifest.R",
  "tests/testthat/test-mv05k-integrated-retrieval-evaluation.R"
)
evidence <- file.path("docs/audits", c(
  "mv05k-bootstrap-audit-2026-08-10.csv",
  "mv05k-deterministic-repeat-2026-08-10.csv",
  "mv05k-endpoint-dispositions-2026-08-10.csv",
  "mv05k-independent-validation-2026-08-10.csv",
  "mv05k-label-source-provenance-2026-08-10.csv",
  "mv05k-method-endpoint-summaries-2026-08-10.csv",
  "mv05k-method-intervals-2026-08-10.csv",
  "mv05k-paired-contrasts-2026-08-10.csv",
  "mv05k-prediction-lock-2026-08-10.csv",
  "mv05k-production-summary-2026-08-10.csv",
  "mv05k-query-endpoints-2026-08-10.csv",
  "mv05k-randomization-audit-2026-08-10.csv",
  "mv05k-sample-eligibility-2026-08-10.csv",
  "mv05k-sample-endpoint-summaries-2026-08-10.csv",
  "mv05k-seed-macro-endpoints-2026-08-10.csv",
  "mv05k-tissue-endpoint-summaries-2026-08-10.csv",
  "mv05k-tissue-seed-endpoints-2026-08-10.csv"
))
upstream <- file.path("docs/audits", c(
  "mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz",
  "mv05j-method-completion-2026-08-09.csv",
  "mv05j-group-bundle-index-2026-08-09.csv",
  "mv05j-method-registry-2026-08-09.csv",
  "mv05j-component-scale-disposition-2026-08-09.csv",
  "mv05j-public-assembly-summary-2026-08-09.csv",
  "mv05j-integrated-retrieval-evaluation-authorization-2026-08-09.csv"
))

artifacts <- data.frame(
  artifact_path = c(decision, implementation, evidence, upstream),
  artifact_role = c(
    rep("decision_document", length(decision)),
    rep("implementation_or_test", length(implementation)),
    rep("evaluation_evidence", length(evidence)),
    rep("frozen_prediction_input", length(upstream))
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop(
    "Missing MV5-K artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "), call. = FALSE
  )
}
artifacts$contract_id <- "mv05k_public_artifact_manifest_v1"
artifacts$size_bytes <- unname(file.info(artifacts$artifact_path)$size)
artifacts$sha256 <- vapply(artifacts$artifact_path, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1L))
artifacts$outcome_label_state <- ifelse(
  artifacts$artifact_role == "frozen_prediction_input", "closed",
  ifelse(artifacts$artifact_role == "implementation_or_test",
         "not_applicable_code", "opened_for_mv05k_integrated_retrieval_evaluation")
)
artifacts$biological_outcomes_computed <-
  artifacts$artifact_role %in% c("decision_document", "evaluation_evidence")
artifacts$public_safe <- TRUE
artifacts <- artifacts[, c(
  "contract_id", "artifact_path", "artifact_role", "size_bytes", "sha256",
  "outcome_label_state", "biological_outcomes_computed", "public_safe"
)]
write_provenance_csv(artifacts, args[[1L]])
message("Wrote MV5-K public artifact manifest: ", args[[1L]])
