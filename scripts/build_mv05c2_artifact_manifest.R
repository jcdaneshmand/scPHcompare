#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05c2_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
output_path <- args[[1L]]
source("R/provenance_utils.R")

fixed <- c(
  "R/mv05_resource_safe_execution.R",
  "tests/testthat/test-mv05-resource-safe-execution.R",
  "docs/specifications/MV05C2_RESOURCE_SAFE_EXECUTION_SPECIFICATION_V1.md",
  "docs/architecture/ADR-007-MV05C2-RESOURCE-SAFE-EXECUTION.md",
  "docs/audits/MV05C2_RESOURCE_SAFE_EXECUTION_2026-08-06.md",
  "PROJECT_PLAN.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/phase-4-primary-benchmark.md",
  "docs/plans/WORK_LOG.md"
)
scripts <- sort(list.files(
  "scripts", pattern = "mv05c2", full.names = TRUE
), method = "radix")
audits <- sort(list.files(
  "docs/audits", pattern = "^mv05c2.*\\.csv$", full.names = TRUE
), method = "radix")
paths <- unique(c(fixed, scripts, audits))
paths <- paths[normalizePath(
  paths, winslash = "/", mustWork = FALSE
) != normalizePath(output_path, winslash = "/", mustWork = FALSE)]
if (!length(paths) || any(!file.exists(paths)) ||
    any(grepl("(?i)(\\.pdf$|example_run|reviewer|pasted-text|^tmp/)",
              paths, perl = TRUE))) {
  stop("MV5-C2 public artifact scope is missing or contains private files.",
       call. = FALSE)
}
role <- ifelse(
  startsWith(paths, "R/") | startsWith(paths, "scripts/"), "implementation",
  ifelse(startsWith(paths, "tests/"), "verification",
         ifelse(grepl("\\.csv$", paths), "evidence", "decision_document"))
)
result <- data.frame(
  contract_id = "mv05c2_public_artifact_manifest_v1",
  artifact_path = gsub("\\\\", "/", paths),
  artifact_role = role,
  size_bytes = unname(file.info(paths)$size),
  sha256 = vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)),
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
result <- result[order(result$artifact_path, method = "radix"), , drop = FALSE]
write_provenance_csv(result, output_path)
message("Wrote MV5-C2 public artifact manifest for ", nrow(result), " files.")
