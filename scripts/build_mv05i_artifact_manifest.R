#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05i_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
artifacts <- data.frame(
  artifact_path = c(
    "PROJECT_PLAN.md",
    "R/mv05i_integrated_landscape_production.R",
    "R/mv05i_post_projection.R",
    "docs/architecture/ADR-019-MV05I-INTEGRATED-CELL-LANDSCAPE-DISTANCES.md",
    "docs/audits/MV05I_INTEGRATED_CELL_LANDSCAPE_DISTANCES_2026-08-09.md",
    "docs/audits/mv05i-admission-min-completion-2026-08-09.csv",
    "docs/audits/mv05i-admission-min-r-oracle-2026-08-09.csv",
    "docs/audits/mv05i-admission-max-completion-2026-08-09.csv",
    "docs/audits/mv05i-admission-max-r-oracle-2026-08-09.csv",
    "docs/audits/mv05i-admission-resume-validation-2026-08-09.csv",
    "docs/audits/mv05i-integrated-cell-landscape-chunk-manifest-2026-08-09.csv",
    "docs/audits/mv05i-integrated-cell-landscape-pair-manifest-2026-08-09.csv.gz",
    "docs/audits/mv05i-pair-manifest-compression-validation-2026-08-09.csv",
    "docs/audits/mv05i-group-input-audit-2026-08-09.csv",
    "docs/audits/mv05i-production-group-resources-2026-08-09.csv",
    "docs/audits/mv05i-integrated-landscape-completion-2026-08-09.csv",
    "docs/audits/mv05i-independent-chunk-validation-2026-08-09.csv",
    "docs/audits/mv05i-independent-group-validation-2026-08-09.csv",
    "docs/audits/mv05i-exact-maximum-group-repeat-2026-08-09.csv",
    "docs/audits/mv05i-matrix-boundary-validation-2026-08-09.csv",
    "docs/audits/mv05i-production-resume-validation-2026-08-09.csv",
    "docs/audits/mv05i-integrated-cell-landscape-component-distances-2026-08-09.csv.gz",
    "docs/audits/mv05i-post-distance-projection-2026-08-09.csv",
    "docs/audits/mv05i-integrated-retrieval-authorization-2026-08-09.csv",
    "docs/audits/mv05i-verification-summary-2026-08-09.csv",
    "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
    "docs/plans/WORK_LOG.md",
    "docs/specifications/MV05I_INTEGRATED_CELL_LANDSCAPE_DISTANCE_SPECIFICATION_V1.md",
    "scripts/build_mv05i_artifact_manifest.R",
    "scripts/build_mv05i_pair_manifest.R",
    "scripts/build_mv05i_post_projection.R",
    "scripts/compress_mv05i_pair_manifest.R",
    "scripts/monitor_mv05i_landscape_groups.R",
    "scripts/mv05i_landscape_group.py",
    "scripts/stage_mv05i_group_inputs.R",
    "scripts/validate_mv05i_admission.R",
    "scripts/validate_mv05i_compressed_manifest.R",
    "scripts/validate_mv05i_production.R",
    "scripts/verify_mv05i_resume.R",
    "tests/testthat/test-mv05i-integrated-landscape-production.R"
  ),
  artifact_role = c(
    "decision_document", "implementation", "implementation",
    "decision_document", "decision_document",
    rep("evidence", 20L),
    "decision_document", "audit_log", "specification",
    rep("execution_script", 11L), "test"
  ),
  stringsAsFactors = FALSE
)
if (nrow(artifacts) != length(artifacts$artifact_role) ||
    any(!file.exists(artifacts$artifact_path))) {
  stop(
    "Missing MV5-I artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "), call. = FALSE
  )
}
artifacts$contract_id <- "mv05i_public_artifact_manifest_v1"
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
message("Wrote MV5-I public artifact manifest: ", args[[1L]])
