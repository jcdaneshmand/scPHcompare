#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05ab_artifact_manifest.R AUDIT_DIR", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
audit_dir_arg <- gsub("\\\\", "/", args[[1L]])
audit_dir <- normalizePath(audit_dir_arg, winslash = "/", mustWork = TRUE)
paths <- c(
  "R/mv05_benchmark_execution.R",
  "R/mv05_inductive_mapping.R",
  "scripts/build_mv05ab_execution_scaffold.R",
  "scripts/run_mv05_inductive_fixture.R",
  "scripts/build_mv05ab_artifact_manifest.R",
  "tests/testthat/test-mv05-benchmark-execution.R",
  "docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md",
  "docs/architecture/ADR-005-MV05AB-EXECUTION-SCAFFOLD.md",
  "docs/audits/MV05AB_EXECUTION_SCAFFOLD_2026-08-06.md",
  file.path(audit_dir_arg, c(
    "mv05a-loso-execution-manifest-2026-08-06.csv",
    "mv05a-loso-execution-summary-2026-08-06.csv",
    "mv05a-baseline-analytical-validation-2026-08-06.csv",
    "mv05a-baseline-scientific-shape-feasibility-2026-08-06.csv",
    "mv05b-inductive-mapping-feasibility-2026-08-06.csv"
  ))
)
roles <- c(
  "execution_contract", "inductive_mapping_contract",
  "execution_evidence_builder", "mapping_fixture_runner",
  "artifact_manifest_builder", "contract_tests", "parent_specification",
  "decision_record", "audit_report", rep("generated_technical_evidence", 5L)
)
stopifnot(length(paths) == length(roles), all(file.exists(paths)))
manifest <- data.frame(
  artifact_path = gsub("\\\\", "/", paths),
  sha256 = vapply(paths, digest::digest, character(1L),
                  file = TRUE, algo = "sha256"),
  bytes = as.numeric(file.info(paths)$size),
  artifact_role = roles,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest,
  file.path(audit_dir, "mv05ab-artifact-manifest-2026-08-06.csv"),
  row.names = FALSE
)
message("Built MV5-A/MV5-B artifact manifest with ", nrow(manifest), " rows.")
