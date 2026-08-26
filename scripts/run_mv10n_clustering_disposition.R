#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: run_mv10n_clustering_disposition.R <prefreeze> <mv10h-synthesis>",
  "<output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
synthesis <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
execution_head <- tolower(trimws(args[[4L]]))
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-N output")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv10m_clustering_disposition.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10m-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10m-contract.csv"))
implementation <- readc(file.path(prefreeze, "mv10m-implementation-bindings.csv"))
decision <- readc(file.path(prefreeze, "mv10m-decision.csv"))
schemas <- readc(file.path(prefreeze, "mv10m-source-schema.csv"))
if (nrow(contract) != 1L || execution_head != contract$execution_head ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256) ||
    !.mv08z_truth(decision$execution_authorized_after_commit) ||
    !all(vapply(file.path(synthesis, schemas$file), sha, character(1L)) ==
           schemas$sha256)) {
  stop("MV10-N execution binding drift")
}
values <- lapply(schemas$file, function(name) readc(file.path(synthesis, name)))
names(values) <- schemas$source_id
result <- do.call(mv10m_build_methodological_disposition_v1, values)
artifacts <- list(
  "mv10n-selected-primary-seed-metrics.csv" =
    result$selected_primary_seed_metrics,
  "mv10n-primary-summary.csv" = result$primary_summary,
  "mv10n-representation-context.csv" = result$representation_context,
  "mv10n-method-sensitivity.csv" = result$method_sensitivity,
  "mv10n-disposition.csv" = result$disposition
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
receipt <- data.frame(
  contract_id = "mv10n_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", output_tables = 5L,
  selected_seed_rows = nrow(result$selected_primary_seed_metrics),
  primary_summary_rows = nrow(result$primary_summary),
  representation_context_rows = nrow(result$representation_context),
  method_sensitivity_rows = nrow(result$method_sensitivity),
  disposition_rows = nrow(result$disposition),
  labels_used = FALSE, outcomes_used = FALSE,
  biological_interpretation = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(output, "mv10n-terminal-receipt.csv"))
files <- sort(setdiff(list.files(output), "mv10n-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10n_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10n-artifact-manifest.csv"))
message("Completed MV10-N methodological disposition; tables=5")
