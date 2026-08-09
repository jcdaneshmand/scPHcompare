#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv05d1_artifact_manifest.R OUTPUT_CSV", call. = FALSE)
}
output_path <- args[[1L]]
source("R/provenance_utils.R")
fixed <- c(
  "R/mv05_resource_safe_execution.R",
  "R/mv05d1_post_fold_projection.R",
  "tests/testthat/test-mv05-resource-safe-execution.R",
  "docs/specifications/MV05D1_LABEL_CLOSED_CELL_FOLD_SPECIFICATION_V1.md",
  "docs/architecture/ADR-009-MV05D0-STAGE1-COMPLETION.md",
  "docs/architecture/ADR-010-MV05D1-LABEL-CLOSED-CELL-FOLDS.md",
  "docs/audits/MV05D1_LABEL_CLOSED_CELL_FOLD_COMPLETION_2026-08-07.md",
  "PROJECT_PLAN.md",
  "docs/plans/DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md",
  "docs/plans/phase-4-primary-benchmark.md",
  "docs/plans/WORK_LOG.md"
)
scripts <- sort(list.files(
  "scripts", pattern = "mv05d1", full.names = TRUE
), method = "radix")
audits <- sort(list.files(
  "docs/audits", pattern = "^mv05d1.*\\.csv$", full.names = TRUE
), method = "radix")
paths <- unique(c(fixed, scripts, audits))
paths <- paths[normalizePath(paths, winslash = "/", mustWork = FALSE) !=
                 normalizePath(output_path, winslash = "/", mustWork = FALSE)]
if (!length(paths) || any(!file.exists(paths)) ||
    any(grepl("(?i)(\\.pdf$|example_run|reviewer|pasted-text|^tmp/)",
              paths, perl = TRUE))) {
  stop("MV5-D1 public artifact scope is missing or contains private files.",
       call. = FALSE)
}
csvs <- paths[grepl("\\.csv$", paths)]
for (path in csvs) {
  value <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (any(tolower(names(value)) %in% c("tissue", "approach")) ||
      ("outcome_label_state" %in% names(value) &&
       any(value$outcome_label_state != "closed")) ||
      ("biological_outcomes_computed" %in% names(value) &&
       any(as.logical(value$biological_outcomes_computed)))) {
    stop("MV5-D1 public CSV violates the label firewall: ", path,
         call. = FALSE)
  }
}
role <- ifelse(
  startsWith(paths, "R/") | startsWith(paths, "scripts/"), "implementation",
  ifelse(startsWith(paths, "tests/"), "verification",
         ifelse(grepl("\\.csv$", paths), "evidence", "decision_document"))
)
result <- data.frame(
  contract_id = "mv05d1_public_artifact_manifest_v1",
  artifact_path = gsub("\\\\", "/", paths), artifact_role = role,
  size_bytes = unname(file.info(paths)$size),
  sha256 = vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L)),
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
result <- result[order(result$artifact_path, method = "radix"), , drop = FALSE]
write_provenance_csv(result, output_path)
message("Wrote MV5-D1 public artifact manifest for ", nrow(result), " files.")
