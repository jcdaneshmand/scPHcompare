#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv07i_outcome_entry.R OUTCOME_PREFREEZE MV7D_PREFREEZE ",
       "MV7E_PREFREEZE SELECTED_PARTITIONS OUTPUT_DIR", call. = FALSE)
}
source("R/mv05s_outcome_execution.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
sha256 <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
prefreeze <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07d <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07e <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
selected_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[5L]]
decision <- read_csv(file.path(prefreeze, "mv07i-outcome-decision.csv"))
queue <- read_csv(file.path(prefreeze, "mv07i-outcome-queue.csv"))
partitions <- read_csv(file.path(prefreeze, "mv07i-outcome-partitions.csv"))
endpoints <- read_csv(file.path(prefreeze, "mv07i-outcome-endpoints.csv"))
input_manifest <- read_csv(file.path(
  prefreeze, "mv07i-outcome-input-manifest.csv"))
selected <- read_csv(selected_path)
if (nrow(decision) != 1L || decision$decision !=
      "authorize_MV7I_descriptive_outcome_execution_only" ||
    !as.logical(decision$labels_authorized_for_next_stage) ||
    !as.logical(decision$outcomes_authorized_for_next_stage) ||
    as.logical(decision$claims_authorized) || nrow(queue) != 120L ||
    nrow(partitions) != 12L || nrow(endpoints) != 6L || nrow(selected) != 7440L ||
    !identical(unique(partitions$selected_partition_sha256),
               sha256(selected_path)) ||
    nrow(input_manifest) != 8L || !all(file.exists(input_manifest$path)) ||
    !identical(tolower(unname(vapply(
      input_manifest$path, sha256, character(1L)))),
      tolower(input_manifest$sha256)) ||
    any(selected$outcome_label_state != "closed") ||
    any(as.logical(selected$biological_outcomes_computed))) {
  stop("MV7-I outcome execution admission is stale.", call. = FALSE)
}

named_files <- c(
  metadata_join = "metadata-join.csv", contingency = "contingency-long.csv",
  seed_metrics = "seed-metrics.csv", unit_summaries = "unit-summaries.csv",
  structural_status = "structural-status.csv", provenance = "provenance.csv")
if (dir.exists(output_dir)) {
  status_path <- file.path(output_dir, "status.csv")
  if (!file.exists(status_path)) stop("Existing MV7-I outcome output is incomplete.")
  status <- read_csv(status_path)
  paths <- stats::setNames(file.path(output_dir, named_files), names(named_files))
  if (!all(file.exists(paths))) stop("Existing MV7-I outcome artifact is missing.")
  expected <- unname(unlist(status[paste0(
    names(named_files), "_sha256")], use.names = FALSE))
  observed <- unname(vapply(paths, sha256, character(1L)))
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      !identical(tolower(observed), tolower(expected))) {
    stop("Existing MV7-I outcome output is stale.", call. = FALSE)
  }
  message("Reused complete MV7-I descriptive outcome artifacts.")
  quit(save = "no", status = 0L)
}
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = "mv07i_outcome__partial__",
                    tmpdir = dirname(output_dir))
dir.create(partial)
started <- proc.time()[["elapsed"]]

reconciliation <- read_csv(file.path(
  mv07d, "mv07d-sample-reconciliation.csv"))
approach <- read_csv(file.path(mv07e, "mv07e-canonical-approach.csv"))
reconciliation <- reconciliation[as.logical(
  reconciliation$corrected_descriptive_124), , drop = FALSE]
reconciliation <- reconciliation[order(
  reconciliation$sample_id, method = "radix"), ]
approach <- approach[order(approach$sample_id, method = "radix"), ]
if (nrow(reconciliation) != 124L || nrow(approach) != 124L ||
    !identical(reconciliation$sample_id, approach$sample_id) ||
    !identical(reconciliation$tissue, approach$tissue) ||
    !identical(reconciliation$study, approach$study) ||
    !setequal(selected$sample_id, reconciliation$sample_id)) {
  stop("MV7-I outcome metadata axis drifted.", call. = FALSE)
}
metadata <- data.frame(
  contract_id = "mv07i_outcome_private_metadata_join_v1",
  sample_id = approach$sample_id, tissue = approach$tissue,
  study = approach$study, canonical_approach = approach$canonical_approach,
  corrected_primary_90 = as.logical(approach$corrected_primary_90),
  canonical_approach_source = approach$canonical_source,
  historical_heuristic_approach_used = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE)

seed_rows <- list(); summary_rows <- list(); contingency_rows <- list()
seed_cursor <- 0L; summary_cursor <- 0L; contingency_cursor <- 0L
for (unit_index in seq_len(nrow(queue))) {
  unit <- queue[unit_index, , drop = FALSE]
  population_ids <- if (unit$population_id == "full124_descriptive") {
    metadata$sample_id
  } else {
    metadata$sample_id[metadata$corrected_primary_90]
  }
  metric_id <- if (unit$metric_id == "normalized_mutual_information_max")
    "normalized_mutual_information" else unit$metric_id
  unit_values <- numeric(5L)
  for (seed_index in seq_along(20260805:20260809)) {
    seed <- (20260805:20260809)[[seed_index]]
    part <- selected[
      selected$representation_id == unit$representation_id &
        selected$algorithm_id == unit$algorithm_id & selected$seed == seed,
      c("sample_id", "cluster"), drop = FALSE]
    part <- part[part$sample_id %in% population_ids, , drop = FALSE]
    part <- part[order(part$sample_id, method = "radix"), ]
    metadata_index <- match(part$sample_id, metadata$sample_id)
    labels <- metadata[[unit$label_axis]][metadata_index]
    if (nrow(part) != unit$expected_samples_per_seed || anyNA(metadata_index) ||
        anyNA(labels) || anyDuplicated(part$sample_id)) {
      stop("MV7-I outcome unit sample/label axis is malformed.", call. = FALSE)
    }
    metric <- mv05s_metric_v1(metric_id, part$cluster, labels)
    unit_values[[seed_index]] <- metric$estimate
    seed_cursor <- seed_cursor + 1L
    seed_rows[[seed_cursor]] <- data.frame(
      contract_id = "mv07i_outcome_seed_metric_v1",
      evaluation_unit_id = unit$evaluation_unit_id,
      execution_order = unit$execution_order, endpoint_id = unit$endpoint_id,
      population_id = unit$population_id, label_axis = unit$label_axis,
      representation_id = unit$representation_id,
      algorithm_id = unit$algorithm_id, algorithm_role = unit$algorithm_role,
      selected_k = unit$selected_k, metric_id = unit$metric_id, seed = seed,
      estimate = metric$estimate, samples = nrow(part),
      label_classes = length(unique(labels)),
      partition_clusters = length(unique(part$cluster)), status = metric$status,
      p_value_computed = FALSE, method_selection_executed = FALSE,
      stringsAsFactors = FALSE)
    if (unit$metric_id == "adjusted_rand_index") {
      cells <- as.data.frame(table(cluster = part$cluster,
                                   label_value = labels),
                            stringsAsFactors = FALSE)
      contingency_cursor <- contingency_cursor + 1L
      contingency_rows[[contingency_cursor]] <- data.frame(
        contract_id = "mv07i_outcome_private_contingency_v1",
        endpoint_id = unit$endpoint_id, population_id = unit$population_id,
        label_axis = unit$label_axis,
        representation_id = unit$representation_id,
        algorithm_id = unit$algorithm_id, selected_k = unit$selected_k,
        seed = seed, cluster = cells$cluster, label_value = cells$label_value,
        samples = cells$Freq, method_selection_executed = FALSE,
        stringsAsFactors = FALSE)
    }
  }
  summary_cursor <- summary_cursor + 1L
  summary_rows[[summary_cursor]] <- data.frame(
    contract_id = "mv07i_outcome_unit_summary_v1",
    evaluation_unit_id = unit$evaluation_unit_id,
    execution_order = unit$execution_order, endpoint_id = unit$endpoint_id,
    population_id = unit$population_id, label_axis = unit$label_axis,
    representation_id = unit$representation_id,
    algorithm_id = unit$algorithm_id, algorithm_role = unit$algorithm_role,
    selected_k = unit$selected_k, metric_id = unit$metric_id,
    seed_mean = mean(unit_values), seed_median = median(unit_values),
    seed_minimum = min(unit_values), seed_maximum = max(unit_values),
    seed_jackknife_se = mv05s_jackknife_se_v1(unit_values),
    completed_seeds = sum(is.finite(unit_values)), expected_seeds = 5L,
    status = if (all(is.finite(unit_values))) "completed" else
      "degenerate_label_or_partition_metric_not_identifiable",
    p_value_computed = FALSE, method_selection_executed = FALSE,
    stringsAsFactors = FALSE)
}
seed_metrics <- do.call(rbind, seed_rows); rownames(seed_metrics) <- NULL
unit_summaries <- do.call(rbind, summary_rows); rownames(unit_summaries) <- NULL
contingency <- do.call(rbind, contingency_rows); rownames(contingency) <- NULL
nonestimable <- endpoints[
  endpoints$execution_status == "structurally_not_estimable_single_class", ]
structural_status <- data.frame(
  contract_id = "mv07i_outcome_structural_status_v1",
  endpoint_id = nonestimable$endpoint_id,
  population_id = nonestimable$population_id,
  label_axis = nonestimable$label_axis,
  status = nonestimable$execution_status, samples = 90L, label_classes = 1L,
  metric_rows_computed = 0L, p_value_computed = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE)
if (nrow(seed_metrics) != 600L || nrow(unit_summaries) != 120L ||
    any(!is.finite(seed_metrics$estimate)) ||
    any(seed_metrics$status != "completed") ||
    any(unit_summaries$status != "completed") ||
    any(seed_metrics$p_value_computed) || any(unit_summaries$p_value_computed) ||
    any(seed_metrics$method_selection_executed) ||
    any(unit_summaries$method_selection_executed) || nrow(structural_status) != 1L) {
  stop("MV7-I descriptive outcome result contract failed.", call. = FALSE)
}
provenance <- data.frame(
  contract_id = "mv07i_outcome_provenance_v1",
  selected_partition_sha256 = sha256(selected_path),
  metadata_samples = nrow(metadata), evaluation_units = nrow(unit_summaries),
  seed_metric_rows = nrow(seed_metrics), contingency_rows = nrow(contingency),
  structural_nonestimable_endpoints = nrow(structural_status),
  p_values_computed = FALSE, method_selection_executed = FALSE,
  historical_heuristic_approach_used = FALSE,
  outcomes_computed = TRUE, stringsAsFactors = FALSE)

paths <- stats::setNames(file.path(partial, named_files), names(named_files))
write.csv(metadata, paths[["metadata_join"]], row.names = FALSE, na = "")
write.csv(contingency, paths[["contingency"]], row.names = FALSE, na = "")
write.csv(seed_metrics, paths[["seed_metrics"]], row.names = FALSE, na = "")
write.csv(unit_summaries, paths[["unit_summaries"]], row.names = FALSE, na = "")
write.csv(structural_status, paths[["structural_status"]],
          row.names = FALSE, na = "")
write.csv(provenance, paths[["provenance"]], row.names = FALSE, na = "")
elapsed <- proc.time()[["elapsed"]] - started
hashes <- vapply(paths, sha256, character(1L))
status <- data.frame(
  contract_id = "mv07i_outcome_status_v1", completion_state = "complete",
  elapsed_seconds = elapsed,
  metadata_join_sha256 = hashes[["metadata_join"]],
  contingency_sha256 = hashes[["contingency"]],
  seed_metrics_sha256 = hashes[["seed_metrics"]],
  unit_summaries_sha256 = hashes[["unit_summaries"]],
  structural_status_sha256 = hashes[["structural_status"]],
  provenance_sha256 = hashes[["provenance"]],
  p_values_computed = FALSE, method_selection_executed = FALSE,
  outcomes_computed = TRUE, stringsAsFactors = FALSE)
write.csv(status, file.path(partial, "status.csv"), row.names = FALSE, na = "")
if (!file.rename(partial, output_dir)) {
  unlink(partial, recursive = TRUE)
  stop("Failed to atomically publish MV7-I outcome artifacts.")
}
message("Completed MV7-I descriptive outcome artifacts.")
