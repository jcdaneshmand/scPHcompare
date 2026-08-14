#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05e_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")

decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-015-MV05E-PREDICTION-LOCKED-RETRIEVAL-EVALUATION.md",
  "docs/audits/MV05E_PREDICTION_LOCKED_RETRIEVAL_EVALUATION_2026-08-08.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/specifications/MV05E_PREDICTION_LOCKED_RETRIEVAL_EVALUATION_SPECIFICATION_V1.md"
)
implementation <- c(
  "R/mv05e_retrieval_evaluation.R",
  "scripts/run_mv05e_retrieval_evaluation.R",
  "scripts/validate_mv05e_retrieval_evaluation.R",
  "scripts/verify_mv05e_deterministic_repeat.R",
  "scripts/build_mv05e_artifact_manifest.R",
  "tests/testthat/test-mv05e-retrieval-evaluation.R"
)
evidence <- file.path("docs/audits", c(
  "mv05e-bootstrap-audit-2026-08-08.csv",
  "mv05e-deterministic-repeat-2026-08-08.csv",
  "mv05e-endpoint-dispositions-2026-08-08.csv",
  "mv05e-independent-validation-2026-08-08.csv",
  "mv05e-label-source-provenance-2026-08-08.csv",
  "mv05e-method-endpoint-summaries-2026-08-08.csv",
  "mv05e-method-intervals-2026-08-08.csv",
  "mv05e-paired-contrasts-2026-08-08.csv",
  "mv05e-prediction-lock-2026-08-08.csv",
  "mv05e-production-summary-2026-08-08.csv",
  "mv05e-query-endpoints-2026-08-08.csv",
  "mv05e-randomization-audit-2026-08-08.csv",
  "mv05e-sample-eligibility-2026-08-08.csv",
  "mv05e-sample-endpoint-summaries-2026-08-08.csv",
  "mv05e-seed-macro-endpoints-2026-08-08.csv",
  "mv05e-tissue-endpoint-summaries-2026-08-08.csv",
  "mv05e-tissue-seed-endpoints-2026-08-08.csv"
))
upstream <- file.path("docs/audits", c(
  "mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz",
  "mv05d5-method-completion-2026-08-08.csv",
  "mv05d5-group-bundle-index-2026-08-08.csv",
  "mv05d5-public-assembly-summary-2026-08-08.csv"
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
    "Missing MV5-E artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "), call. = FALSE
  )
}
artifacts$contract_id <- "mv05e_public_artifact_manifest_v1"
artifacts$size_bytes <- unname(file.info(artifacts$artifact_path)$size)
artifacts$sha256 <- vapply(artifacts$artifact_path, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1L))
artifacts$outcome_label_state <- ifelse(
  artifacts$artifact_role == "frozen_prediction_input", "closed",
  ifelse(artifacts$artifact_role == "implementation_or_test",
         "not_applicable_code", "opened_for_mv05e_retrieval_evaluation")
)
artifacts$biological_outcomes_computed <-
  artifacts$artifact_role %in% c("decision_document", "evaluation_evidence")
artifacts$public_safe <- TRUE
artifacts <- artifacts[, c(
  "contract_id", "artifact_path", "artifact_role", "size_bytes", "sha256",
  "outcome_label_state", "biological_outcomes_computed", "public_safe"
)]
write_provenance_csv(artifacts, args[[1L]])
message("Wrote MV5-E public artifact manifest: ", args[[1L]])
