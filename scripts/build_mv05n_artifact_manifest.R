#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05n_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")

decision <- c(
  "PROJECT_PLAN.md",
  "docs/PROJECT_EVIDENCE.md",
  "docs/architecture/ADR-024-MV05N-LABEL-CLOSED-CLUSTERING-GATE.md",
  "docs/audits/MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_2026-08-10.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/specifications/MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_SPECIFICATION_V1.md"
)
implementation <- c(
  "R/mv05n_clustering_gate.R",
  "scripts/build_mv05n_training_pair_inventory.R",
  "scripts/stage_mv05n_admission_inputs.R",
  "scripts/mv05n_landscape_admission.py",
  "scripts/validate_mv05n_admission.R",
  "scripts/run_mv05n_baseline_admission.R",
  "scripts/validate_mv05n_baseline_admission.R",
  "scripts/snapshot_mv05n_resume.R",
  "scripts/verify_mv05n_resume.R",
  "scripts/build_mv05n_artifact_manifest.R",
  "tests/testthat/test-mv05n-clustering-gate.R"
)
evidence <- file.path("docs/audits", c(
  "mv05n-training-pair-identity-summary-2026-08-10.csv",
  "mv05n-training-pair-group-inventory-2026-08-10.csv",
  "mv05n-training-pair-chunk-inventory-2026-08-10.csv",
  "mv05n-admission-requests-2026-08-10.csv",
  "mv05n-admission-r-oracle-2026-08-10.csv",
  "mv05n-admission-completion-2026-08-10.csv",
  "mv05n-admission-resources-2026-08-10.csv",
  "mv05n-full-matrix-resource-projection-2026-08-10.csv",
  "mv05n-admission-exact-repeat-2026-08-10.csv",
  "mv05n-admission-resume-validation-2026-08-10.csv",
  "mv05n-baseline-admission-validation-2026-08-10.csv",
  "mv05n-baseline-admission-resources-2026-08-10.csv",
  "mv05n-combined-resource-projection-2026-08-10.csv",
  "mv05n-canonical-package-validation-2026-08-10.csv"
))
upstream <- c(
  "R/mv05_benchmark_contract.R",
  "R/mv05d4_landscape_production.R",
  "R/mv05i_integrated_landscape_production.R",
  "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md",
  "docs/audits/MV05M_POST_RETRIEVAL_BENCHMARK_GAP_GATE_2026-08-10.md"
)

artifacts <- data.frame(
  artifact_path = c(decision, implementation, evidence, upstream),
  artifact_role = c(
    rep("decision_document", length(decision)),
    rep("implementation_or_test", length(implementation)),
    rep("gate_evidence", length(evidence)),
    rep("public_upstream_evidence", length(upstream))
  ),
  stringsAsFactors = FALSE
)
missing <- !file.exists(artifacts$artifact_path)
if (any(missing)) {
  stop("Missing MV5-N artifact(s): ",
       paste(artifacts$artifact_path[missing], collapse = ", "),
       call. = FALSE)
}
artifacts$contract_id <- "mv05n_public_artifact_manifest_v1"
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
message("Wrote MV5-N public artifact manifest: ", args[[1L]])
