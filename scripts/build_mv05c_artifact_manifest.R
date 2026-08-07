#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05c_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
output <- args[[1L]]

fixed <- c(
  "NAMESPACE", "R/dual_view_topology.R", "R/mv05_benchmark_execution.R",
  "tests/testthat/test-mv05-benchmark-execution.R",
  "man/new_cell_projection_source.Rd",
  "man/construct_frozen_cell_topology_view.Rd",
  "PROJECT_PLAN.md",
  "docs/architecture/ADR-006-MV05C-EXISTING-DATA-FEASIBILITY.md",
  "docs/audits/MV05C_EXISTING_DATA_FEASIBILITY_2026-08-06.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/phase-4-primary-benchmark.md", "docs/plans/WORK_LOG.md",
  "scripts/build_mv05c_execution_manifest.R",
  "scripts/run_mv05c_existing_data_job.R",
  "scripts/run_mv05c_job_queue.sh", "scripts/assemble_mv05c_inputs.R",
  "scripts/mv05c_landscape_batch.py",
  "scripts/build_mv05c_label_free_artifacts.R",
  "scripts/build_mv05c_artifact_manifest.R"
)
audit_files <- list.files(
  "docs/audits", pattern = "^mv05c-", recursive = TRUE,
  full.names = TRUE, ignore.case = TRUE
)
paths <- sort(unique(c(fixed, audit_files)), method = "radix")
paths <- paths[file.exists(paths) & normalizePath(paths, winslash = "/") !=
                 normalizePath(output, winslash = "/", mustWork = FALSE)]
if (!length(paths) || any(grepl(
  "(^|/)(private|docs/private)(/|$)|\\.pdf$|pasted-text|example_run\\.r$",
  paths, ignore.case = TRUE
))) {
  stop("MV5-C public artifact scope is empty or contains forbidden material.")
}
role <- ifelse(
  grepl("^docs/audits/mv05c-|MV05C_EXISTING", paths), "audit_evidence",
  ifelse(grepl("^docs/architecture", paths), "architecture_decision",
  ifelse(grepl("^docs/plans|PROJECT_PLAN", paths), "auditable_plan",
  ifelse(grepl("^tests/", paths), "test",
  ifelse(grepl("^man/|NAMESPACE", paths), "public_api_documentation", "code"))))
)
manifest <- data.frame(
  contract_id = "mv05c_public_artifact_manifest_v1",
  path = gsub("\\\\", "/", paths), role = role,
  size_bytes = unname(file.info(paths)$size),
  sha256 = vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L)),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (anyDuplicated(manifest$path) || any(!is.finite(manifest$size_bytes)) ||
    any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
  stop("MV5-C artifact manifest failed identity validation.")
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, output, row.names = FALSE)
message("Hashed ", nrow(manifest), " public MV5-C artifacts.")
