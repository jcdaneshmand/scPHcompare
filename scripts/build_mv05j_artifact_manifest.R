#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05j_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")

decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-020-MV05J-INTEGRATED-CELL-RETRIEVAL-INPUTS.md",
  "docs/audits/MV05J_INTEGRATED_CELL_RETRIEVAL_INPUTS_2026-08-09.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/specifications/MV05J_INTEGRATED_CELL_RETRIEVAL_INPUT_SPECIFICATION_V1.md"
)
evidence <- c(
  "docs/audits/mv05j-admission-maximum-validation-2026-08-09.csv",
  "docs/audits/mv05j-admission-representative-validation-2026-08-09.csv",
  "docs/audits/mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz",
  "docs/audits/mv05j-component-scale-disposition-2026-08-09.csv",
  "docs/audits/mv05j-group-bundle-index-2026-08-09.csv",
  "docs/audits/mv05j-independent-group-validation-2026-08-09.csv",
  "docs/audits/mv05j-independent-validation-summary-2026-08-09.csv",
  "docs/audits/mv05j-integrated-retrieval-evaluation-authorization-2026-08-09.csv",
  "docs/audits/mv05j-method-completion-2026-08-09.csv",
  "docs/audits/mv05j-method-registry-2026-08-09.csv",
  "docs/audits/mv05j-production-group-resources-2026-08-09.csv",
  "docs/audits/mv05j-public-assembly-summary-2026-08-09.csv",
  "docs/audits/mv05j-public-repeat-validation-2026-08-09.csv",
  "docs/audits/mv05j-resume-validation-2026-08-09.csv"
)
scripts <- c(
  "scripts/assemble_mv05j_public_artifacts.R",
  "scripts/build_mv05j_artifact_manifest.R",
  "scripts/build_mv05j_resource_authorization.R",
  "scripts/run_mv05j_group.R",
  "scripts/run_mv05j_queue.R",
  "scripts/validate_mv05j_admission.R",
  "scripts/validate_mv05j_production.R",
  "scripts/verify_mv05j_public_repeat.R",
  "scripts/verify_mv05j_resume.R"
)
artifacts <- data.frame(
  artifact_path = c(
    decision, "R/mv05j_integrated_retrieval_inputs.R", evidence, scripts,
    "tests/testthat/test-mv05j-integrated-retrieval-inputs.R"
  ),
  artifact_role = c(
    rep("decision_document", length(decision)),
    "implementation", rep("evidence", length(evidence)),
    rep("execution_script", length(scripts)), "test"
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop(
    "Missing MV5-J artifact(s): ",
    paste(artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
          collapse = ", "), call. = FALSE
  )
}
artifacts$contract_id <- "mv05j_public_artifact_manifest_v1"
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
message("Wrote MV5-J public artifact manifest: ", args[[1L]])
