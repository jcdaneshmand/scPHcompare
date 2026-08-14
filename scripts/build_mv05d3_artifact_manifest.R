#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05d3_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
artifacts <- data.frame(
  artifact_path = c(
    "PROJECT_PLAN.md",
    "R/mv05d3_ph_production.R",
    "R/mv05d3_post_ph_projection.R",
    "docs/architecture/ADR-012-MV05D3-FULL-CELL-PH-PRODUCTION.md",
    "docs/audits/MV05D3_FULL_CELL_PH_PRODUCTION_2026-08-07.md",
    "docs/audits/mv05d3-cell-ph-completion-2026-08-07.csv",
    "docs/audits/mv05d3-exact-group-repeat-2026-08-07.csv",
    "docs/audits/mv05d3-full-cell-ph-manifest-2026-08-07.csv",
    "docs/audits/mv05d3-group-resources-2026-08-07.csv",
    "docs/audits/mv05d3-group-resume-validation-2026-08-07.csv",
    "docs/audits/mv05d3-independent-group-validation-2026-08-07.csv",
    "docs/audits/mv05d3-measured-primary-projection-2026-08-07.csv",
    "docs/audits/mv05d3-view-resources-2026-08-07.csv",
    "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
    "docs/plans/WORK_LOG.md",
    "docs/specifications/MV05D3_FULL_CELL_PH_PRODUCTION_SPECIFICATION_V1.md",
    "scripts/build_mv05d3_artifact_manifest.R",
    "scripts/build_mv05d3_full_ph_manifest.R",
    "scripts/monitor_mv05d3_ph_groups.R",
    "scripts/run_mv05d3_ph_group.R",
    "scripts/validate_mv05d3_group_resume.R",
    "scripts/validate_mv05d3_production.R",
    "tests/testthat/test-mv05d2-ph-profiling.R",
    "tests/testthat/test-mv05d3-ph-production.R"
  ),
  artifact_role = c(
    "decision_document", "implementation", "implementation",
    "decision_document", "decision_document", rep("evidence", 8L),
    "decision_document", "audit_log", "specification",
    rep("execution_script", 6L), rep("test", 2L)
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop(
    "Missing MV5-D3 artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "),
    call. = FALSE
  )
}
artifacts$contract_id <- "mv05d3_public_artifact_manifest_v1"
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
message("Wrote MV5-D3 public artifact manifest: ", args[[1L]])
