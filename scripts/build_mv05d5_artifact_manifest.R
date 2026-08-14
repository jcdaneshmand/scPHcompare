#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05d5_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")

decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-014-MV05D5-CELL-RETRIEVAL-INPUTS.md",
  "docs/audits/MV05D5_CELL_RETRIEVAL_INPUTS_2026-08-08.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md"
)
evidence <- c(
  "docs/audits/mv05d5-admission-group-2026-08-08.csv",
  "docs/audits/mv05d5-admission-validation-2026-08-08.csv",
  "docs/audits/mv05d5-cell-retrieval-rankings-2026-08-08.csv.gz",
  "docs/audits/mv05d5-component-scale-disposition-2026-08-08.csv",
  "docs/audits/mv05d5-group-bundle-index-2026-08-08.csv",
  "docs/audits/mv05d5-independent-group-validation-2026-08-08.csv",
  "docs/audits/mv05d5-independent-validation-summary-2026-08-08.csv",
  "docs/audits/mv05d5-mean-profile-staging-2026-08-08.csv",
  "docs/audits/mv05d5-method-completion-2026-08-08.csv",
  "docs/audits/mv05d5-method-registry-2026-08-08.csv",
  "docs/audits/mv05d5-production-group-resources-2026-08-08.csv",
  "docs/audits/mv05d5-public-assembly-summary-2026-08-08.csv",
  "docs/audits/mv05d5-public-repeat-validation-2026-08-08.csv",
  "docs/audits/mv05d5-resume-validation-2026-08-08.csv"
)
scripts <- c(
  "scripts/assemble_mv05d5_public_artifacts.R",
  "scripts/build_mv05d5_artifact_manifest.R",
  "scripts/run_mv05d5_group.R",
  "scripts/run_mv05d5_queue.R",
  "scripts/stage_mv05d5_mean_profiles.R",
  "scripts/validate_mv05d5_admission.R",
  "scripts/validate_mv05d5_production.R",
  "scripts/verify_mv05d5_public_repeat.R",
  "scripts/verify_mv05d5_resume.R"
)
artifacts <- data.frame(
  artifact_path = c(
    decision, "docs/plans/WORK_LOG.md", "NAMESPACE",
    "docs/specifications/MV05D5_CELL_RETRIEVAL_INPUT_SPECIFICATION_V1.md",
    "R/mv05d5_retrieval_inputs.R", evidence, scripts,
    "tests/testthat/test-mv05d5-retrieval-inputs.R"
  ),
  artifact_role = c(
    rep("decision_document", length(decision)), "audit_log", "implementation",
    "specification",
    "implementation", rep("evidence", length(evidence)),
    rep("execution_script", length(scripts)), "test"
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop(
    "Missing MV5-D5 artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "), call. = FALSE
  )
}
artifacts$contract_id <- "mv05d5_public_artifact_manifest_v1"
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
message("Wrote MV5-D5 public artifact manifest: ", args[[1L]])
