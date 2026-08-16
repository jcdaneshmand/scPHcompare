#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: validate_mv08g_sources.R PREFREEZE PRIVATE_ROOT EXECUTION_EVIDENCE OUTPUT")
}
prefreeze <- args[[1L]]; private_root <- args[[2L]]
execution <- args[[3L]]; output <- args[[4L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G source validation output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
manifest <- read.csv(file.path(prefreeze, "mv08g-cache-manifest.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
panel <- read.csv(file.path(prefreeze, "mv08g-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
queue <- read.csv(file.path(prefreeze, "mv08g-source-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
metrics <- read.csv(file.path(execution, "mv08g-source-metrics.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
repeat_validation <- read.csv(file.path(execution,
  "mv08g-source-repeat-validation.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
execution_decision <- read.csv(file.path(execution, "mv08g-source-decision.csv"),
                               stringsAsFactors = FALSE, check.names = FALSE)
mv07h_validate_cache_manifest_v1(manifest)
mv08g_validate_common475_panel_v1(panel)
if (nrow(queue) != 5L || nrow(metrics) != 5L ||
    nrow(repeat_validation) != 1L || !truth(repeat_validation$byte_identical) ||
    execution_decision$decision != "source_complete_await_independent_validation") {
  stop("MV8-G source execution evidence is incomplete.")
}
rows <- list(); total_views <- c(cell_topology_v1 = 0L, gene_topology_v1 = 0L)
for (index in seq_len(nrow(queue))) {
  seed <- as.integer(queue$seed[[index]])
  path <- file.path(private_root, "source",
                    paste0("mv08g__", seed, "__source.rds"))
  metric <- metrics[metrics$seed == seed, , drop = FALSE]
  if (nrow(metric) != 1L || !file.exists(path) ||
      metric$disposition != "completed" || metric$output_sha256 != sha(path) ||
      metric$output_bytes != as.numeric(file.info(path)$size) ||
      metric$elapsed_seconds > metric$elapsed_cap_seconds ||
      metric$peak_process_tree_rss_bytes > metric$rss_cap_bytes) {
    stop("MV8-G source metric or artifact identity drift for seed ", seed)
  }
  record <- readRDS(path)
  mv08g_validate_source_record_v1(record)
  part <- manifest[manifest$seed == seed, , drop = FALSE]
  part <- part[order(part$sample_id, method = "radix"), , drop = FALSE]
  if (record$identity$seed != seed ||
      !identical(record$identity$sample_ids, part$sample_id) ||
      !identical(unname(record$identity$input_cache_keys),
                 unname(part$normalization_cache_key)) ||
      !identical(unname(record$identity$selected_cell_sha256),
                 unname(part$selected_cell_sha256)) ||
      !identical(record$panel$feature_id, panel$feature_id) ||
      any(!is.finite(record$center)) || any(!is.finite(record$scale)) ||
      any(record$scale <= sqrt(.Machine$double.eps))) {
    stop("MV8-G source transform identity drift for seed ", seed)
  }
  view_counts <- table(unlist(lapply(record$views, names), use.names = FALSE))
  if (!identical(as.integer(view_counts[.mv08g_views]), c(124L, 124L))) {
    stop("MV8-G source typed-view balance failed for seed ", seed)
  }
  total_views <- total_views + view_counts[.mv08g_views]
  rows[[index]] <- data.frame(
    contract_id = "mv08g_source_independent_validation_v1", seed = seed,
    source_sha256 = sha(path), source_bytes = as.numeric(file.info(path)$size),
    source_cache_key = record$cache_key, samples = length(record$views),
    panel_genes = nrow(record$panel), fit_cells = record$identity$fit_cells,
    pca_components = record$pca_model$n_components,
    cell_views = unname(view_counts[["cell_topology_v1"]]),
    gene_views = unname(view_counts[["gene_topology_v1"]]),
    cache_axis_exact = TRUE, selected_cell_axis_exact = TRUE,
    common475_axis_exact = TRUE, transform_finite = TRUE,
    resource_caps_passed = TRUE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
}
validated <- do.call(rbind, rows)
if (nrow(validated) != 5L || sum(total_views) != 1240L ||
    any(validated$samples != 124L) || any(validated$panel_genes != 475L) ||
    any(validated$pca_components != 30L) ||
    anyDuplicated(validated$source_cache_key)) {
  stop("MV8-G independent source validation is incomplete.")
}
summary <- data.frame(
  contract_id = "mv08g_source_validation_summary_v1",
  source_jobs = nrow(validated), source_repeat_jobs = 1L,
  samples_per_seed = 124L, typed_views = sum(total_views),
  cell_views = unname(total_views[["cell_topology_v1"]]),
  gene_views = unname(total_views[["gene_topology_v1"]]),
  maximum_process_tree_rss_bytes = max(metrics$peak_process_tree_rss_bytes),
  aggregate_elapsed_seconds = sum(metrics$elapsed_seconds),
  all_resource_caps_passed = all(validated$resource_caps_passed),
  exact_repeat_passed = truth(repeat_validation$byte_identical),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_source_validation_decision_v1",
  decision = "source_exact_authorize_PH_execution_prefreeze_only",
  source_jobs_exact = 5L, source_repeat_jobs_exact = 1L,
  ph_jobs_authorized = 0L, prospective_ph_jobs = 1240L,
  ripserr_primary_jobs_authorized = 0L,
  gene_gudhi_fallback_authorized_only_after_rss_cap = TRUE,
  landscape_jobs_authorized = 0L, matched_shift_jobs_authorized = 0L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "MV8-G_PH_execution_prefreeze",
  stringsAsFactors = FALSE)
outputs <- list(
  "mv08g-source-independent-validation.csv" = validated,
  "mv08g-source-validation-summary.csv" = summary,
  "mv08g-source-validation-decision.csv" = decision)
paths <- vapply(names(outputs), function(name) {
  path <- file.path(output, name); write_provenance_csv(outputs[[name]], path); path
}, character(1L))
artifact_manifest <- data.frame(
  contract_id = "mv08g_source_validation_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(artifact_manifest,
  file.path(output, "mv08g-source-validation-artifact-manifest.csv"))
message("MV8-G source validation passed: five exact bundles, 1,240 typed views")
