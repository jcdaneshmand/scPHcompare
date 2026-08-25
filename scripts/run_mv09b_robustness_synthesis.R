#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: run_mv09b_robustness_synthesis.R <prefreeze>",
  "<mv08zy-public-summary> <output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
source_path <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
execution_head <- tolower(trimws(args[[4L]]))
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV9-B output")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv09_robustness_synthesis.R")
.mv08z_verify_manifest(prefreeze, "mv09a-artifact-manifest.csv")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
contract <- readc(file.path(prefreeze, "mv09a-contract.csv"))
implementation <- readc(file.path(prefreeze, "mv09a-implementation-bindings.csv"))
decision <- readc(file.path(prefreeze, "mv09a-decision.csv"))
truth <- .mv08z_truth
if (nrow(contract) != 1L || execution_head != contract$execution_head ||
    sha(source_path) != contract$source_summary_sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256) ||
    !truth(decision$execution_authorized_after_commit)) {
  stop("MV9-B execution binding drift")
}
result <- mv09_build_robustness_synthesis_v1(readc(source_path))
artifacts <- list(
  "mv09b-aggregate-comparisons.csv" = result$aggregate,
  "mv09b-plot-data.csv" = result$plot_data,
  "mv09b-internal-seed-summary.csv" = result$internal_seed_summary,
  "mv09b-external-singleton.csv" = result$external_singleton,
  "mv09b-dimension-delta.csv" = result$dimension_delta,
  "mv09b-dimension-delta-summary.csv" = result$dimension_delta_summary
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
receipt <- data.frame(
  contract_id = "mv09b_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", aggregate_rows = nrow(result$aggregate),
  plot_rows = nrow(result$plot_data),
  internal_summary_rows = nrow(result$internal_seed_summary),
  external_singleton_rows = nrow(result$external_singleton),
  dimension_delta_rows = nrow(result$dimension_delta),
  dimension_summary_rows = nrow(result$dimension_delta_summary),
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  clustering_jobs = 0L, fusion_jobs = 0L, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv09b-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv09b-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv09b_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv09b-artifact-manifest.csv"))
message("Completed MV9-B aggregate robustness synthesis; rows=40")
