#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05g_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-017-MV05G-LABEL-CLOSED-INTEGRATED-COORDINATE-PRODUCTION.md",
  "docs/audits/MV05G_LABEL_CLOSED_INTEGRATED_COORDINATE_PRODUCTION_2026-08-09.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md"
)
evidence <- c(
  "docs/audits/mv05g-deterministic-repeat-2026-08-09.csv",
  "docs/audits/mv05g-independent-group-validation-2026-08-09.csv",
  "docs/audits/mv05g-independent-validation-summary-2026-08-09.csv",
  "docs/audits/mv05g-integrated-coordinate-manifest-2026-08-08.csv",
  "docs/audits/mv05g-integrated-coordinate-resources-2026-08-08.csv",
  "docs/audits/mv05g-integrated-ph-authorization-2026-08-09.csv",
  "docs/audits/mv05g-post-coordinate-projection-components-2026-08-09.csv",
  "docs/audits/mv05g-post-coordinate-projection-summary-2026-08-09.csv",
  "docs/audits/mv05g-resume-after-2026-08-09.csv",
  "docs/audits/mv05g-resume-before-2026-08-09.csv",
  "docs/audits/mv05g-resume-validation-2026-08-09.csv"
)
scripts <- c(
  "scripts/build_mv05g_artifact_manifest.R",
  "scripts/build_mv05g_full_manifest.R",
  "scripts/build_mv05g_post_coordinate_projection.R",
  "scripts/run_mv05g_production_queue.sh",
  "scripts/snapshot_mv05g_private_results.R",
  "scripts/validate_mv05g_coordinate_production.R",
  "scripts/verify_mv05g_deterministic_repeat.R",
  "scripts/verify_mv05g_resume.R",
  "scripts/monitor_mv05f_mapping_pilot.R",
  "scripts/run_mv05f_mapping_group.R"
)
artifacts <- data.frame(
  artifact_path = c(
    decision, "docs/plans/WORK_LOG.md", "docs/PROJECT_EVIDENCE.md",
    "docs/specifications/MV05G_LABEL_CLOSED_INTEGRATED_COORDINATE_PRODUCTION_SPECIFICATION_V1.md",
    "R/mv05g_coordinate_production.R", evidence, scripts,
    "tests/testthat/test-mv05g-coordinate-production.R"
  ),
  artifact_role = c(
    rep("decision_document", length(decision)), "audit_log", "evidence_ledger",
    "specification", "implementation", rep("evidence", length(evidence)),
    rep("execution_script", length(scripts)), "test"
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop("Missing MV5-G public artifact.", call. = FALSE)
}
artifacts$contract_id <- "mv05g_public_artifact_manifest_v1"
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
message("Wrote MV5-G public artifact manifest.")
