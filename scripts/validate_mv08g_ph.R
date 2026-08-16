#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: validate_mv08g_ph.R PREFREEZE SOURCE_ROOT PH_ROOT EXECUTION_EVIDENCE OUTPUT")
}
prefreeze <- args[[1L]]; source_root <- args[[2L]]; ph_root <- args[[3L]]
execution <- args[[4L]]; output <- args[[5L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G PH validation output must be empty.")
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
queue <- read.csv(file.path(prefreeze, "mv08g-ph-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
metrics <- read.csv(file.path(execution, "mv08g-ph-metrics.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
attempts <- read.csv(file.path(execution, "mv08g-ph-engine-attempts.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
repeats <- read.csv(file.path(execution, "mv08g-ph-repeat-validation.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
execution_decision <- read.csv(file.path(execution, "mv08g-ph-decision.csv"),
                               stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(queue) != 1240L || nrow(metrics) != 1240L || nrow(repeats) != 4L ||
    anyDuplicated(metrics$job_id) || !all(truth(repeats$byte_identical)) ||
    execution_decision$decision != "full_PH_complete_await_independent_validation" ||
    execution_decision$aggregate_elapsed_seconds >
      execution_decision$aggregate_elapsed_cap_seconds ||
    execution_decision$private_storage_bytes >
      execution_decision$aggregate_storage_cap_bytes) {
  stop("MV8-G PH execution evidence is incomplete.")
}
metrics <- metrics[match(queue$job_id, metrics$job_id), , drop = FALSE]
if (anyNA(metrics$job_id)) stop("MV8-G PH metric axis differs from its queue.")
source_cache <- new.env(parent = emptyenv())
get_source <- function(seed) {
  key <- as.character(seed)
  if (!exists(key, source_cache, inherits = FALSE)) {
    record <- readRDS(file.path(source_root, paste0("mv08g__", seed,
                                                   "__source.rds")))
    mv08g_validate_source_record_v1(record)
    assign(key, record, source_cache)
  }
  get(key, source_cache, inherits = FALSE)
}
rows <- vector("list", nrow(queue))
for (index in seq_len(nrow(queue))) {
  row <- queue[index, , drop = FALSE]; metric <- metrics[index, , drop = FALSE]
  path <- file.path(ph_root, paste0("mv08g__", row$seed, "__", row$sample_id,
    "__", row$view_id, "__ph.rds"))
  if (!file.exists(path) || metric$disposition != "completed" ||
      metric$output_sha256 != sha(path) ||
      metric$output_bytes != as.numeric(file.info(path)$size) ||
      metric$elapsed_seconds > metric$elapsed_cap_seconds ||
      metric$peak_process_tree_rss_bytes > metric$rss_cap_bytes) {
    stop("MV8-G PH resource or file identity drift at job ", row$job_id)
  }
  source_record <- get_source(row$seed)
  view <- source_record$views[[row$sample_id]][[row$view_id]]
  record <- readRDS(path)
  mv08g_validate_ph_record_v1(record, view)
  if (record$topology_result$provenance$ph_engine != metric$ph_engine) {
    stop("MV8-G selected engine provenance mismatch at job ", row$job_id)
  }
  diagram <- record$topology_result$diagram
  finite_h0 <- sum(diagram[, "dimension"] == 0 & is.finite(diagram[, "death"]) &
                     diagram[, "death"] > diagram[, "birth"])
  finite_h1 <- sum(diagram[, "dimension"] == 1 & is.finite(diagram[, "death"]) &
                     diagram[, "death"] > diagram[, "birth"])
  if (finite_h0 != length(view$point_ids) - 1L ||
      record$topology_result$provenance$essential_h0_count != 1L ||
      !isTRUE(record$h0_mst_oracle$passed) ||
      record$identity$common475_axis_sha256 !=
        "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba") {
    stop("MV8-G PH interval/MST contract failed at job ", row$job_id)
  }
  rows[[index]] <- data.frame(
    contract_id = "mv08g_ph_independent_validation_v1", job_id = row$job_id,
    seed = row$seed, sample_id = row$sample_id, view_id = row$view_id,
    selected_engine = metric$ph_engine, ph_sha256 = sha(path),
    ph_bytes = as.numeric(file.info(path)$size), point_count = length(view$point_ids),
    finite_h0_intervals = finite_h0, finite_h1_intervals = finite_h1,
    essential_h0_intervals = 1L,
    mst_maximum_absolute_error = record$h0_mst_oracle$maximum_absolute_error,
    mst_tolerance = record$h0_mst_oracle$tolerance, mst_passed = TRUE,
    resource_cap_passed = TRUE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
  rm(record, diagram)
  if (index %% 100L == 0L) invisible(gc(FALSE))
}
validated <- do.call(rbind, rows)
fallback_attempts <- attempts[attempts$ph_engine == "TDA_ripsDiag_GUDHI", ]
rss_primary <- attempts[attempts$ph_engine == "ripserr" &
  attempts$disposition == "rss_cap_exceeded", ]
fallback_ok <- nrow(fallback_attempts) == nrow(rss_primary) &&
  all(fallback_attempts$view_id == "gene_topology_v1") &&
  all(fallback_attempts$disposition == "completed") &&
  setequal(fallback_attempts$job_id, rss_primary$job_id)
summary <- data.frame(
  contract_id = "mv08g_ph_validation_summary_v1", ph_records = nrow(validated),
  cell_records = sum(validated$view_id == "cell_topology_v1"),
  gene_records = sum(validated$view_id == "gene_topology_v1"),
  finite_h0_intervals = sum(validated$finite_h0_intervals),
  finite_h1_intervals = sum(validated$finite_h1_intervals),
  mst_oracles_passed = sum(validated$mst_passed),
  ripserr_selected_records = sum(validated$selected_engine == "ripserr"),
  gudhi_selected_records = sum(validated$selected_engine == "TDA_ripsDiag_GUDHI"),
  fallback_policy_exact = fallback_ok, exact_repeats = sum(truth(repeats$byte_identical)),
  all_resource_caps_passed = all(validated$resource_cap_passed),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
if (summary$ph_records != 1240L || summary$cell_records != 620L ||
    summary$gene_records != 620L || summary$mst_oracles_passed != 1240L ||
    summary$exact_repeats != 4L || !summary$fallback_policy_exact) {
  stop("MV8-G PH independent closure failed.")
}
decision <- data.frame(
  contract_id = "mv08g_ph_validation_decision_v1",
  decision = "full_PH_exact_authorize_landscape_execution_prefreeze_only",
  ph_records_exact = 1240L, ph_repeats_exact = 4L,
  landscape_jobs_authorized = 0L, prospective_landscape_groups = 20L,
  matched_shift_jobs_authorized = 0L, prospective_matched_shift_groups = 20L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "MV8-G_landscape_execution_prefreeze", stringsAsFactors = FALSE)
outputs <- list(
  "mv08g-ph-independent-validation.csv" = validated,
  "mv08g-ph-validation-summary.csv" = summary,
  "mv08g-ph-validation-decision.csv" = decision)
paths <- vapply(names(outputs), function(name) {
  path <- file.path(output, name); write_provenance_csv(outputs[[name]], path); path
}, character(1L))
artifact_manifest <- data.frame(
  contract_id = "mv08g_ph_validation_artifact_manifest_v1",
  file = basename(paths), bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(artifact_manifest,
  file.path(output, "mv08g-ph-validation-artifact-manifest.csv"))
message("MV8-G PH validation passed: 1,240 records; landscape prefreeze only")
