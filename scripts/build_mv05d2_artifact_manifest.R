#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05d2_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
artifacts <- data.frame(
  artifact_path = c(
    "PROJECT_PLAN.md",
    "R/mv05d2_ph_profiling.R",
    "R/mv05d2_post_ph_projection.R",
    "docs/architecture/ADR-011-MV05D2-BOUNDED-CELL-PH-PROFILING.md",
    "docs/audits/MV05D2_BOUNDED_CELL_PH_PROFILING_2026-08-07.md",
    "docs/audits/mv05d2-cell-ph-pilot-completion-2026-08-07.csv",
    "docs/audits/mv05d2-cell-ph-pilot-manifest-2026-08-07.csv",
    "docs/audits/mv05d2-cell-ph-pilot-resources-2026-08-07.csv",
    "docs/audits/mv05d2-cell-ph-profile-summary-2026-08-07.csv",
    "docs/audits/mv05d2-cell-ph-repeat-resources-2026-08-07.csv",
    "docs/audits/mv05d2-cross-engine-validation-2026-08-07.csv",
    "docs/audits/mv05d2-exact-repeat-validation-2026-08-07.csv",
    "docs/audits/mv05d2-full-cell-ph-projection-2026-08-07.csv",
    "docs/audits/mv05d2-post-ph-primary-projection-2026-08-07.csv",
    "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
    "docs/specifications/MV05D2_BOUNDED_CELL_PH_PROFILING_SPECIFICATION_V1.md",
    "scripts/build_mv05d2_artifact_manifest.R",
    "scripts/build_mv05d2_ph_pilot_manifest.R",
    "scripts/build_mv05d2_post_ph_projection.R",
    "scripts/monitor_mv05d2_ph_pilot.R",
    "scripts/run_mv05d2_ph_entry.R",
    "scripts/validate_mv05d2_ph_pilot.R",
    "tests/testthat/test-mv05d2-ph-profiling.R"
  ),
  artifact_role = c(
    "decision_document", "implementation", "implementation", "decision_document",
    "decision_document", rep("evidence", 9L), "decision_document",
    "specification", rep("execution_script", 6L), "test"
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop(
    "Missing MV5-D2 artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "),
    call. = FALSE
  )
}
artifacts$contract_id <- "mv05d2_public_artifact_manifest_v1"
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
message("Wrote MV5-D2 public artifact manifest: ", args[[1L]])
