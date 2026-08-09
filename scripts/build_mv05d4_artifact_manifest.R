#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05d4_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
artifacts <- data.frame(
  artifact_path = c(
    "PROJECT_PLAN.md",
    "R/mv05d4_landscape_production.R",
    "R/mv05d4_post_projection.R",
    "docs/architecture/ADR-013-MV05D4-CELL-LANDSCAPE-DISTANCES.md",
    "docs/audits/MV05D4_CELL_LANDSCAPE_DISTANCES_2026-08-07.md",
    "docs/audits/mv05d4-admission-completion-2026-08-07.csv",
    "docs/audits/mv05d4-admission-max-completion-2026-08-07.csv",
    "docs/audits/mv05d4-admission-max-r-oracle-2026-08-07.csv",
    "docs/audits/mv05d4-admission-r-oracle-2026-08-07.csv",
    "docs/audits/mv05d4-admission-resume-validation-2026-08-07.csv",
    "docs/audits/mv05d4-cell-landscape-chunk-manifest-2026-08-07.csv",
    "docs/audits/mv05d4-cell-landscape-completion-2026-08-07.csv",
    "docs/audits/mv05d4-cell-landscape-component-distances-2026-08-07.csv.gz",
    "docs/audits/mv05d4-cell-landscape-pair-manifest-2026-08-07.csv.gz",
    "docs/audits/mv05d4-exact-group-repeat-2026-08-07.csv",
    "docs/audits/mv05d4-group-input-audit-2026-08-07.csv",
    "docs/audits/mv05d4-independent-chunk-validation-2026-08-07.csv",
    "docs/audits/mv05d4-independent-group-validation-2026-08-07.csv",
    "docs/audits/mv05d4-measured-primary-projection-2026-08-07.csv",
    "docs/audits/mv05d4-monitor-launch-failure-2026-08-07.csv",
    "docs/audits/mv05d4-production-group-resources-2026-08-07.csv",
    "docs/audits/mv05d4-production-v3-group-resources-2026-08-07.csv",
    "docs/audits/mv05d4-timing-separation-correction-2026-08-07.csv",
    "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
    "docs/plans/WORK_LOG.md",
    "docs/specifications/MV05D4_CELL_LANDSCAPE_DISTANCE_SPECIFICATION_V1.md",
    "scripts/build_mv05d4_artifact_manifest.R",
    "scripts/build_mv05d4_pair_manifest.R",
    "scripts/build_mv05d4_post_projection.R",
    "scripts/compress_mv05d4_pair_manifest.R",
    "scripts/monitor_mv05d4_landscape_groups.R",
    "scripts/mv05d4_landscape_group.py",
    "scripts/stage_mv05d4_group_inputs.R",
    "scripts/validate_mv05d4_admission.R",
    "scripts/validate_mv05d4_production.R",
    "scripts/verify_mv05d4_resume.R",
    "tests/testthat/test-mv05d4-landscape-production.R"
  ),
  artifact_role = c(
    "decision_document", "implementation", "implementation",
    "decision_document", "decision_document", rep("evidence", 18L),
    "decision_document", "audit_log", "specification",
    rep("execution_script", 10L), "test"
  ),
  stringsAsFactors = FALSE
)
if (nrow(artifacts) != length(artifacts$artifact_role) ||
    any(!file.exists(artifacts$artifact_path))) {
  stop(
    "Missing MV5-D4 artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "), call. = FALSE
  )
}
artifacts$contract_id <- "mv05d4_public_artifact_manifest_v1"
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
message("Wrote MV5-D4 public artifact manifest: ", args[[1L]])
