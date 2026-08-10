#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05h_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-018-MV05H-INTEGRATED-CELL-PH-PRODUCTION.md",
  "docs/audits/MV05H_INTEGRATED_CELL_PH_PRODUCTION_2026-08-09.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md"
)
evidence <- c(
  "docs/audits/mv05h-admission-independent-summary-2026-08-09.csv",
  "docs/audits/mv05h-admission-independent-validation-2026-08-09.csv",
  "docs/audits/mv05h-admission-resources-2026-08-09.csv",
  "docs/audits/mv05h-completion-summary-2026-08-09.csv",
  "docs/audits/mv05h-deterministic-repeat-detail-2026-08-09.csv",
  "docs/audits/mv05h-deterministic-repeat-summary-2026-08-09.csv",
  "docs/audits/mv05h-failures-2026-08-09.csv",
  "docs/audits/mv05h-independent-group-validation-2026-08-09.csv",
  "docs/audits/mv05h-independent-validation-summary-2026-08-09.csv",
  "docs/audits/mv05h-integrated-landscape-authorization-2026-08-09.csv",
  "docs/audits/mv05h-integrated-ph-manifest-2026-08-09.csv",
  "docs/audits/mv05h-integrated-ph-resources-2026-08-09.csv",
  "docs/audits/mv05h-post-ph-projection-components-2026-08-09.csv",
  "docs/audits/mv05h-post-ph-projection-summary-2026-08-09.csv",
  "docs/audits/mv05h-resume-after-2026-08-09.csv",
  "docs/audits/mv05h-resume-before-2026-08-09.csv",
  "docs/audits/mv05h-resume-validation-2026-08-09.csv",
  "docs/audits/mv05h-view-resources-2026-08-09.csv"
)
scripts <- c(
  "scripts/assemble_mv05h_completion.R",
  "scripts/build_mv05h_artifact_manifest.R",
  "scripts/build_mv05h_integrated_ph_manifest.R",
  "scripts/build_mv05h_post_ph_projection.R",
  "scripts/monitor_mv05h_ph_groups.R",
  "scripts/run_mv05h_ph_group.R",
  "scripts/run_mv05h_production_queue.sh",
  "scripts/snapshot_mv05h_private_results.R",
  "scripts/validate_mv05h_integrated_ph.R",
  "scripts/verify_mv05h_deterministic_repeat.R",
  "scripts/verify_mv05h_resume.R"
)
artifacts <- data.frame(
  artifact_path = c(
    decision, "docs/plans/WORK_LOG.md", "docs/PROJECT_EVIDENCE.md",
    "docs/specifications/MV05H_INTEGRATED_CELL_PH_PRODUCTION_SPECIFICATION_V1.md",
    "R/mv05h_integrated_ph_production.R", evidence, scripts,
    "tests/testthat/test-mv05h-integrated-ph-production.R"
  ),
  artifact_role = c(
    rep("decision_document", length(decision)), "audit_log", "evidence_ledger",
    "specification", "implementation", rep("evidence", length(evidence)),
    rep("execution_script", length(scripts)), "test"
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop("Missing MV5-H public artifact.", call. = FALSE)
}
artifacts$contract_id <- "mv05h_public_artifact_manifest_v1"
artifacts$size_bytes <- unname(file.info(artifacts$artifact_path)$size)
artifacts$sha256 <- vapply(artifacts$artifact_path, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1L))
artifacts$outcome_label_state <- "closed"
artifacts$biological_outcomes_computed <- FALSE
artifacts <- artifacts[, c(
  "contract_id", "artifact_path", "artifact_role", "size_bytes", "sha256",
  "outcome_label_state", "biological_outcomes_computed"
)]
write_provenance_csv(artifacts, args[[1L]])
message("Wrote MV5-H public artifact manifest.")
