#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05m_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-023-MV05M-BENCHMARK-GAP-GATE.md",
  "docs/audits/MV05M_POST_RETRIEVAL_BENCHMARK_GAP_GATE_2026-08-10.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/plans/phase-4-primary-benchmark.md",
  "docs/specifications/MV05M_POST_RETRIEVAL_BENCHMARK_GAP_GATE_SPECIFICATION_V1.md"
)
implementation <- c(
  "R/mv05m_benchmark_gap_gate.R",
  "scripts/build_mv05m_benchmark_gap_gate.R",
  "scripts/validate_mv05m_benchmark_gap_gate.R",
  "scripts/verify_mv05m_deterministic_repeat.R",
  "scripts/build_mv05m_artifact_manifest.R",
  "tests/testthat/test-mv05m-benchmark-gap-gate.R"
)
evidence <- file.path("docs/audits", c(
  "mv05m-axis-readiness-2026-08-10.csv",
  "mv05m-candidate-scores-2026-08-10.csv",
  "mv05m-deterministic-repeat-2026-08-10.csv",
  "mv05m-evidence-registry-2026-08-10.csv",
  "mv05m-independent-validation-2026-08-10.csv",
  "mv05m-production-summary-2026-08-10.csv",
  "mv05m-selected-next-sprint-2026-08-10.csv",
  "mv05m-selection-criteria-2026-08-10.csv",
  "mv05m-training-pair-scope-2026-08-10.csv"
))
upstream <- c(
  "docs/PROJECT_EVIDENCE.md",
  "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md",
  "docs/audits/MV05E_PREDICTION_LOCKED_RETRIEVAL_EVALUATION_2026-08-08.md",
  "docs/audits/MV05K_PREDICTION_LOCKED_INTEGRATED_RETRIEVAL_EVALUATION_2026-08-10.md",
  "docs/audits/MV05L_LOCKED_REPRESENTATION_COMPARISON_2026-08-10.md",
  "R/mv05d4_landscape_production.R",
  "R/mv05i_integrated_landscape_production.R",
  "R/mv05_benchmark_contract.R"
)
artifacts <- data.frame(
  artifact_path = c(decision, implementation, evidence, upstream),
  artifact_role = c(
    rep("decision_document", length(decision)),
    rep("implementation_or_test", length(implementation)),
    rep("gate_evidence", length(evidence)),
    rep("public_upstream_evidence", length(upstream))
  ), stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop("Missing MV5-M artifact(s): ", paste(
    artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
    collapse = ", "), call. = FALSE)
}
artifacts$contract_id <- "mv05m_public_artifact_manifest_v1"
artifacts$size_bytes <- unname(file.info(artifacts$artifact_path)$size)
artifacts$sha256 <- vapply(artifacts$artifact_path, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1L))
artifacts$biological_outcomes_computed <- FALSE
artifacts$confidential_source_included <- FALSE
artifacts$public_safe <- TRUE
artifacts <- artifacts[, c(
  "contract_id", "artifact_path", "artifact_role", "size_bytes", "sha256",
  "biological_outcomes_computed", "confidential_source_included", "public_safe"
)]
write_provenance_csv(artifacts, args[[1L]])
message("Wrote MV5-M public artifact manifest: ", args[[1L]])
