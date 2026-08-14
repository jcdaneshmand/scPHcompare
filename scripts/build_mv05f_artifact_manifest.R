#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05f_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
source("R/provenance_utils.R")
decision <- c(
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-016-MV05F-LABEL-CLOSED-INTEGRATION-RESOURCE-GATE.md",
  "docs/audits/MV05F_LABEL_CLOSED_INTEGRATION_RESOURCE_GATE_2026-08-08.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md"
)
evidence <- c(
  "docs/audits/mv05f-deterministic-repeat-2026-08-08.csv",
  "docs/audits/mv05f-full-execution-decision-2026-08-08.csv",
  "docs/audits/mv05f-full-workload-components-2026-08-08.csv",
  "docs/audits/mv05f-full-workload-summary-2026-08-08.csv",
  "docs/audits/mv05f-independent-group-validation-2026-08-08.csv",
  "docs/audits/mv05f-independent-validation-summary-2026-08-08.csv",
  "docs/audits/mv05f-mapping-pilot-manifest-2026-08-08.csv",
  "docs/audits/mv05f-mapping-pilot-resources-2026-08-08.csv"
)
scripts <- c(
  "scripts/build_mv05f_artifact_manifest.R",
  "scripts/build_mv05f_mapping_pilot_manifest.R",
  "scripts/build_mv05f_resource_projection.R",
  "scripts/monitor_mv05f_mapping_pilot.R",
  "scripts/run_mv05f_mapping_group.R",
  "scripts/validate_mv05f_mapping_pilot.R",
  "scripts/verify_mv05f_deterministic_repeat.R"
)
artifacts <- data.frame(
  artifact_path = c(
    decision, "docs/plans/WORK_LOG.md",
    "docs/specifications/MV05F_LABEL_CLOSED_INTEGRATION_RESOURCE_GATE_SPECIFICATION_V1.md",
    "R/mv05f_integration_gate.R", evidence, scripts,
    "tests/testthat/test-mv05f-integration-gate.R"
  ),
  artifact_role = c(
    rep("decision_document", length(decision)), "audit_log", "specification",
    "implementation", rep("evidence", length(evidence)),
    rep("execution_script", length(scripts)), "test"
  ),
  stringsAsFactors = FALSE
)
if (any(!file.exists(artifacts$artifact_path))) {
  stop("Missing MV5-F public artifact.", call. = FALSE)
}
artifacts$contract_id <- "mv05f_public_artifact_manifest_v1"
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
message("Wrote MV5-F public artifact manifest.")
