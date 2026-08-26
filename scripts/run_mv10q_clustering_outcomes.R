#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv10q_clustering_outcomes.R <prefreeze> <private-partitions>",
  "<mv07d> <mv07e> <public-output> <private-output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
partitions_path <- normalizePath(args[[2L]], mustWork = TRUE)
mv07d <- normalizePath(args[[3L]], mustWork = TRUE)
mv07e <- normalizePath(args[[4L]], mustWork = TRUE)
public_output <- args[[5L]]; private_output <- args[[6L]]
execution_head <- tolower(trimws(args[[7L]]))
for (path in c(public_output, private_output)) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)))
    stop("refusing to overwrite MV10-Q output")
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}
source("R/mv08z_landscape_production.R")
source("R/mv05s_outcome_execution.R")
source("R/mv10p_clustering_outcomes.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10p-artifact-manifest.csv")
decision <- readc(file.path(prefreeze, "mv10p-decision.csv"))
contract <- readc(file.path(prefreeze, "mv10p-contract.csv"))
implementation <- readc(file.path(prefreeze, "mv10p-implementation-bindings.csv"))
sources <- readc(file.path(prefreeze, "mv10p-source-freeze.csv"))
queue <- readc(file.path(prefreeze, "mv10p-queue.csv"))
endpoints <- readc(file.path(prefreeze, "mv10p-endpoints.csv"))
if (nrow(decision) != 1L || !.mv08z_truth(decision$execution_authorized_after_commit) ||
    nrow(contract) != 1L || execution_head != contract$execution_head ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) == implementation$sha256) ||
    !all(file.exists(sources$artifact)) ||
    !all(vapply(sources$artifact, sha, character(1L)) == sources$sha256) ||
    sha(partitions_path) != sources$sha256[[1L]]) {
  stop("MV10-Q execution binding drift")
}
selected <- mv10p_select_outcome_partitions_v1(readc(partitions_path))
reconciliation <- readc(file.path(mv07d, "mv07d-sample-reconciliation.csv"))
approach <- readc(file.path(mv07e, "mv07e-canonical-approach.csv"))
reconciliation <- reconciliation[
  as.logical(reconciliation$corrected_descriptive_124), , drop = FALSE
]
reconciliation <- reconciliation[order(reconciliation$sample_id), , drop = FALSE]
approach <- approach[order(approach$sample_id), , drop = FALSE]
metadata <- data.frame(
  contract_id = "mv10q_private_metadata_join_v1",
  sample_id = approach$sample_id, tissue = approach$tissue,
  study = approach$study, canonical_approach = approach$canonical_approach,
  corrected_primary_90 = as.logical(approach$corrected_primary_90),
  canonical_approach_source = approach$canonical_source,
  historical_heuristic_approach_used = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE
)
started <- proc.time()[["elapsed"]]
result <- mv10p_evaluate_outcomes_v1(selected, metadata, queue)
nonestimable <- endpoints[
  endpoints$execution_status == "structurally_not_estimable_single_class", ,
  drop = FALSE
]
structural <- data.frame(
  contract_id = "mv10q_outcome_structural_status_v1",
  endpoint_id = nonestimable$endpoint_id,
  population_id = nonestimable$population_id,
  label_axis = nonestimable$label_axis,
  status = nonestimable$execution_status, samples = 90L, label_classes = 1L,
  metric_rows_computed = 0L, p_value_computed = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE
)
elapsed <- proc.time()[["elapsed"]] - started
provenance <- data.frame(
  contract_id = "mv10q_outcome_provenance_v1", execution_head = execution_head,
  selected_partition_source_sha256 = sha(partitions_path),
  metadata_samples = nrow(metadata), evaluation_units = nrow(result$unit_summaries),
  seed_metric_rows = nrow(result$seed_metrics),
  contingency_rows = nrow(result$contingency),
  structural_nonestimable_endpoints = nrow(structural),
  elapsed_seconds = elapsed, workers = 1L, retries = 0L,
  p_values_computed = FALSE, method_selection_executed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
private_artifacts <- list("mv10q-private-metadata-join.csv" = metadata,
                          "mv10q-private-contingency.csv" = result$contingency)
for (name in names(private_artifacts)) atomic(
  private_artifacts[[name]], file.path(private_output, name)
)
private_manifest <- data.frame(
  contract_id = "mv10q_private_manifest_v1",
  artifact = names(private_artifacts),
  bytes = as.numeric(file.info(file.path(private_output,
                                         names(private_artifacts)))$size),
  sha256 = vapply(file.path(private_output, names(private_artifacts)),
                  sha, character(1L)), stringsAsFactors = FALSE
)
atomic(private_manifest, file.path(private_output, "mv10q-private-manifest.csv"))
public_artifacts <- list(
  "mv10q-seed-metrics.csv" = result$seed_metrics,
  "mv10q-unit-summaries.csv" = result$unit_summaries,
  "mv10q-structural-status.csv" = structural,
  "mv10q-provenance.csv" = provenance
)
for (name in names(public_artifacts)) atomic(
  public_artifacts[[name]], file.path(public_output, name)
)
receipt <- data.frame(
  contract_id = "mv10q_terminal_receipt_v1", execution_head = execution_head,
  completion_state = "complete", evaluation_units = nrow(result$unit_summaries),
  seed_metric_rows = nrow(result$seed_metrics),
  private_contingency_rows = nrow(result$contingency),
  public_bytes = sum(file.info(file.path(public_output,
                                         names(public_artifacts)))$size),
  private_bytes = sum(private_manifest$bytes), elapsed_seconds = elapsed,
  workers = 1L, retries = 0L, labels_opened_after_commit = TRUE,
  p_values_computed = FALSE, method_selection_executed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
atomic(receipt, file.path(public_output, "mv10q-terminal-receipt.csv"))
files <- sort(setdiff(list.files(public_output), "mv10q-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10q_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public_output, files))$size),
  sha256 = vapply(file.path(public_output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(public_output, "mv10q-artifact-manifest.csv"))
message("Completed MV10-Q clustering outcomes; units=300")
