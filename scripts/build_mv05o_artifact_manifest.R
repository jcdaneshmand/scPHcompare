#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05o_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
decision <- c(
  "PROJECT_PLAN.md", "docs/PROJECT_EVIDENCE.md",
  "docs/architecture/ADR-025-MV05O-FULL-MATRIX-PRODUCTION-PREFREEZE.md",
  "docs/audits/MV05O_FULL_MATRIX_PRODUCTION_PREFREEZE_2026-08-10.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/WORK_LOG.md",
  "docs/specifications/MV05O_LABEL_CLOSED_FULL_MATRIX_PRODUCTION_PREFREEZE_SPECIFICATION_V1.md"
)
implementation <- c(
  "R/mv05o_production_prefreeze.R",
  "scripts/build_mv05o_prefreeze.R", "scripts/validate_mv05o_prefreeze.R",
  "scripts/verify_mv05o_prefreeze_repeat.R",
  "scripts/stage_mv05o_group_inputs.R", "scripts/mv05o_landscape_group.py",
  "scripts/run_mv05o_baseline_group.R", "scripts/snapshot_mv05o_artifacts.R",
  "scripts/verify_mv05o_resume.R",
  "scripts/validate_mv05o_landscape_runner_fixture.R",
  "scripts/build_mv05o_artifact_manifest.R",
  "tests/testthat/test-mv05o-production-prefreeze.R"
)
evidence <- file.path("docs/audits", c(
  "mv05o-source-freeze-2026-08-10.csv",
  "mv05o-production-group-queue-2026-08-10.csv",
  "mv05o-landscape-chunk-queue-2026-08-10.csv",
  "mv05o-baseline-group-queue-2026-08-10.csv",
  "mv05o-validation-plan-2026-08-10.csv",
  "mv05o-abort-rules-2026-08-10.csv",
  "mv05o-prefreeze-summary-2026-08-10.csv",
  "mv05o-independent-validation-2026-08-10.csv",
  "mv05o-deterministic-repeat-2026-08-10.csv",
  "mv05o-landscape-runner-fixture-validation-2026-08-10.csv",
  "mv05o-runner-resume-validation-2026-08-10.csv",
  "mv05o-canonical-package-validation-2026-08-10.csv"
))
upstream <- c(
  "docs/specifications/MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_SPECIFICATION_V1.md",
  "docs/audits/MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_2026-08-10.md",
  "docs/audits/mv05n-training-pair-identity-summary-2026-08-10.csv",
  "docs/audits/mv05n-combined-resource-projection-2026-08-10.csv"
)
artifacts <- data.frame(
  artifact_path = c(decision, implementation, evidence, upstream),
  artifact_role = c(rep("decision_document", length(decision)),
                    rep("implementation_or_test", length(implementation)),
                    rep("prefreeze_evidence", length(evidence)),
                    rep("accepted_upstream_evidence", length(upstream))),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop("Missing MV5-O artifact(s): ", paste(
    artifacts$artifact_path[!file.exists(artifacts$artifact_path)],
    collapse = ", "), call. = FALSE)
}
artifacts$contract_id <- "mv05o_public_artifact_manifest_v1"
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
message("Wrote MV5-O public artifact manifest: ", args[[1L]])
