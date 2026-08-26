#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv10r_clustering_outcome_closure.R <prefreeze>",
  "<private-partitions> <mv07d> <mv07e> <public-production>",
  "<private-production> <repeat-private> <output-dir> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
partitions_path <- normalizePath(args[[2L]], mustWork = TRUE)
mv07d <- normalizePath(args[[3L]], mustWork = TRUE)
mv07e <- normalizePath(args[[4L]], mustWork = TRUE)
production <- normalizePath(args[[5L]], mustWork = TRUE)
private_production <- normalizePath(args[[6L]], mustWork = TRUE)
repeat_private <- args[[7L]]; output <- args[[8L]]
execution_head <- tolower(trimws(args[[9L]]))
for (path in c(repeat_private, output)) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)))
    stop("refusing to overwrite MV10-R output")
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}
source("R/mv08z_landscape_production.R")
source("R/mv05s_outcome_execution.R")
source("R/mv10p_clustering_outcomes.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10p-artifact-manifest.csv")
.mv08z_verify_manifest(production, "mv10q-artifact-manifest.csv")
queue <- readc(file.path(prefreeze, "mv10p-queue.csv"))
endpoints <- readc(file.path(prefreeze, "mv10p-endpoints.csv"))
contract <- readc(file.path(prefreeze, "mv10p-contract.csv"))
receipt <- readc(file.path(production, "mv10q-terminal-receipt.csv"))
selected <- mv10p_select_outcome_partitions_v1(readc(partitions_path))
approach <- readc(file.path(mv07e, "mv07e-canonical-approach.csv"))
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
repeat_result <- mv10p_evaluate_outcomes_v1(selected, metadata, queue)
nonestimable <- endpoints[
  endpoints$execution_status == "structurally_not_estimable_single_class", ,
  drop = FALSE
]
repeat_structural <- data.frame(
  contract_id = "mv10q_outcome_structural_status_v1",
  endpoint_id = nonestimable$endpoint_id,
  population_id = nonestimable$population_id,
  label_axis = nonestimable$label_axis,
  status = nonestimable$execution_status, samples = 90L, label_classes = 1L,
  metric_rows_computed = 0L, p_value_computed = FALSE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE
)
repeat_private_artifacts <- list(
  "mv10q-private-metadata-join.csv" = metadata,
  "mv10q-private-contingency.csv" = repeat_result$contingency
)
for (name in names(repeat_private_artifacts)) atomic(
  repeat_private_artifacts[[name]], file.path(repeat_private, name)
)
public_map <- c(seed_metrics = "mv10q-seed-metrics.csv",
                unit_summaries = "mv10q-unit-summaries.csv",
                structural = "mv10q-structural-status.csv")
expected <- list(seed_metrics = repeat_result$seed_metrics,
                 unit_summaries = repeat_result$unit_summaries,
                 structural = repeat_structural)
rehash <- do.call(rbind, lapply(seq_along(public_map), function(index) {
  id <- names(public_map)[[index]]
  observed <- readc(file.path(production, public_map[[index]]))
  data.frame(
    contract_id = "mv10r_public_rehash_v1", output_order = index,
    artifact_id = id, rows = nrow(observed), expected_rows = nrow(expected[[id]]),
    exact_columns = identical(names(observed), names(expected[[id]])),
    numeric_repeat = isTRUE(all.equal(observed, expected[[id]], tolerance = 1e-14,
                                      check.attributes = FALSE)),
    stringsAsFactors = FALSE
  )
}))
private_files <- names(repeat_private_artifacts)
private_rehash <- data.frame(
  contract_id = "mv10r_private_rehash_v1", artifact = private_files,
  production_sha256 = vapply(file.path(private_production, private_files),
                             sha, character(1L)),
  repeat_sha256 = vapply(file.path(repeat_private, private_files),
                         sha, character(1L)), stringsAsFactors = FALSE
)
private_rehash$byte_identical <-
  private_rehash$production_sha256 == private_rehash$repeat_sha256
public_names <- names(readc(file.path(production, "mv10q-seed-metrics.csv")))
privacy <- !any(c("sample_id", "label_value", "true_label", "predicted_label") %in%
                  public_names)
validation <- data.frame(
  contract_id = "mv10r_validation_v1",
  check_id = c("prefreeze_manifest", "production_manifest", "terminal_complete",
               "execution_head", "three_hundred_units", "fifteen_hundred_rows",
               "one_structural_status", "three_public_recomputations",
               "all_public_numeric_repeat", "two_private_artifacts",
               "private_byte_repeat", "metadata_124", "selected_18600",
               "three_stacks", "two_dimensions", "five_methods", "five_seeds",
               "selected_K_2_3", "five_endpoints", "two_metrics",
               "all_finite", "all_completed", "no_p_values",
               "no_method_selection", "no_H0_H1_pooling", "public_privacy",
               "one_worker", "zero_retries", "storage_caps", "claim_firewall"),
  passed = c(TRUE, TRUE, receipt$completion_state == "complete",
             receipt$execution_head == execution_head &&
               contract$execution_head == execution_head,
             receipt$evaluation_units == 300L, receipt$seed_metric_rows == 1500L,
             nrow(repeat_structural) == 1L, nrow(rehash) == 3L,
             all(rehash$numeric_repeat & rehash$exact_columns),
             nrow(private_rehash) == 2L, all(private_rehash$byte_identical),
             nrow(metadata) == 124L, nrow(selected) == 18600L,
             length(unique(selected$stack_id)) == 3L,
             length(unique(selected$homology_dimension)) == 2L,
             length(unique(selected$method_id)) == 5L,
             length(unique(selected$seed)) == 5L,
             all(selected$k[selected$homology_dimension == "H0"] == 2L) &&
               all(selected$k[selected$homology_dimension == "H1"] == 3L),
             length(unique(queue$endpoint_id)) == 5L,
             length(unique(queue$metric_id)) == 2L,
             all(is.finite(repeat_result$seed_metrics$estimate)),
             all(repeat_result$seed_metrics$status == "completed"),
             !any(repeat_result$seed_metrics$p_value_computed),
             !any(repeat_result$seed_metrics$method_selection_executed),
             !any(grepl("H0_H1", repeat_result$seed_metrics$homology_dimension)),
             privacy, receipt$workers == 1L, receipt$retries == 0L,
             receipt$public_bytes <= 16 * 1024^2 &&
               receipt$private_bytes <= 16 * 1024^2,
             !receipt$biological_claims && !receipt$manuscript_claims),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-R validation failed")
decision <- data.frame(
  contract_id = "mv10r_decision_v1",
  decision = "close_complete_descriptive_clustering_outcomes",
  evaluation_units = 300L, seed_metric_rows = 1500L,
  next_stage = "separate_complete_outcome_review_prefreeze",
  method_selection_state = "closed", biological_claims_state = "closed",
  manuscript_claims_state = "closed", stringsAsFactors = FALSE
)
artifacts <- list("mv10r-public-rehash.csv" = rehash,
                  "mv10r-private-rehash.csv" = private_rehash,
                  "mv10r-validation.csv" = validation,
                  "mv10r-decision.csv" = decision)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-R clustering-outcome closure", "",
  "All 300 descriptive units and 1,500 seed metrics independently reproduce.",
  "Private metadata/contingency artifacts repeat byte-identically; public",
  "outputs remain aggregate-only. Interpretation and claims remain closed."
), file.path(output, "MV10R_CLUSTERING_OUTCOME_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10r-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10r_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10r-artifact-manifest.csv"))
message("Closed MV10-R clustering outcomes; checks=30")
