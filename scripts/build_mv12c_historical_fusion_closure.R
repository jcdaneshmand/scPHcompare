#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv12c_historical_fusion_closure.R <prefreeze>",
  "<matrix-bundle> <private> <public> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
bundle_path <- normalizePath(args[[2L]], mustWork = TRUE)
private <- normalizePath(args[[3L]], mustWork = TRUE)
public <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]; execution_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output)) stop("MV12-C output exists")
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
source("R/mv12_historical_fusion.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv; truth <- .mv08z_truth
.mv08z_verify_manifest(prefreeze, "mv12a-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv12b-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv12a-contract.csv"))
source_binding <- readc(file.path(prefreeze, "mv12a-source-binding.csv"))
receipt <- readc(file.path(public, "mv12b-terminal-receipt.csv"))
saved_paths <- c(
  assignments = file.path(private, "mv12b-sample-partitions.csv"),
  scales = file.path(public, "mv12b-distance-scales.csv"),
  catalog = file.path(public, "mv12b-matrix-catalog.csv"),
  quality = file.path(public, "mv12b-partition-quality.csv"),
  stability = file.path(public, "mv12b-seed-stability.csv"),
  consensus = file.path(public, "mv12b-consensus-diagnostics.csv"),
  detail = file.path(public, "mv12b-primary-decision-detail.csv"),
  disposition = file.path(public, "mv12b-disposition.csv")
)
bundle <- readRDS(bundle_path)
repeat_result <- mv12_fit_fusion_grid_v1(bundle)
repeat_stability <- mv10_seed_stability_v1(repeat_result$assignments)
repeat_consensus <- mv12_consensus_diagnostics_v1(repeat_result$assignments)
repeat_disposition <- mv12_fusion_decision_v1(repeat_stability,
                                               repeat_consensus)
values <- list(
  assignments = repeat_result$assignments, scales = repeat_result$scales,
  catalog = repeat_result$catalog, quality = repeat_result$quality,
  stability = repeat_stability, consensus = repeat_consensus,
  detail = repeat_disposition$detail,
  disposition = repeat_disposition$decision
)
temp <- tempfile("mv12c-repeat-"); dir.create(temp)
repeat_paths <- setNames(file.path(temp, paste0(names(values), ".csv")),
                         names(values))
for (name in names(values)) atomic(values[[name]], repeat_paths[[name]])
exact <- vapply(names(values), function(name) {
  sha(saved_paths[[name]]) == sha(repeat_paths[[name]])
}, logical(1L))
checks <- c(
  execution_head = execution_head == contract$execution_head,
  source_hash = sha(bundle_path) == source_binding$sha256,
  terminal_complete = identical(receipt$completion_state, "complete"),
  fifty_matrices = receipt$matrices == 50L,
  five_hundred_fits = receipt$partition_fits == 500L,
  sixty_two_thousand_assignments = receipt$private_assignment_rows == 62000L,
  twenty_scales = receipt$scale_rows == 20L,
  fifty_catalog = receipt$catalog_rows == 50L,
  five_hundred_quality = receipt$quality_rows == 500L,
  one_hundred_stability = receipt$stability_rows == 100L,
  three_hundred_consensus = receipt$consensus_rows == 300L,
  two_primary_details = receipt$primary_detail_rows == 2L,
  one_disposition = receipt$disposition_rows == 1L,
  all_scientific_artifacts_exact = all(exact),
  elapsed_cap = receipt$elapsed_seconds <= receipt$elapsed_cap_seconds,
  private_storage_cap = receipt$private_bytes <= receipt$private_storage_cap_bytes,
  public_storage_cap = receipt$public_bytes <= receipt$public_storage_cap_bytes,
  one_worker_zero_retry = receipt$workers == 1L && receipt$retries == 0L,
  H0_H1_separate = !truth(receipt$H0_H1_combined),
  labels_outcomes_closed = !truth(receipt$labels_used) &&
    !truth(receipt$outcomes_used),
  no_selection = !truth(receipt$method_or_weight_selected),
  claims_closed = !truth(receipt$biological_claims) &&
    !truth(receipt$manuscript_claims)
)
validation <- data.frame(
  contract_id = "mv12c_validation_v1", check_id = names(checks),
  passed = unname(checks), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV12-C closure failed")
repeat_audit <- data.frame(
  contract_id = "mv12c_artifact_repeat_v1", artifact_id = names(values),
  rows = vapply(values, nrow, integer(1L)),
  saved_sha256 = vapply(saved_paths, sha, character(1L)),
  repeat_sha256 = vapply(repeat_paths, sha, character(1L)),
  exact_repeat = unname(exact), stringsAsFactors = FALSE
)
decision <- repeat_disposition$decision
decision$historical_fusion_independently_closed <- TRUE
decision$next_action <- if (decision$option2_new_allqc_cell_topology_required) {
  "commit closure; begin prospectively governed newer all-QC cell topology"
} else {
  "commit clear negative closure; retain separate historical views"
}
dir.create(output, recursive = TRUE)
atomic(validation, file.path(output, "mv12c-validation.csv"))
atomic(repeat_audit, file.path(output, "mv12c-artifact-repeat.csv"))
atomic(repeat_disposition$detail,
       file.path(output, "mv12c-primary-decision-detail.csv"))
atomic(decision, file.path(output, "mv12c-decision.csv"))
writeLines(c(
  "# MV12-C historical fusion feasibility closure", "",
  "All eight scientific artifacts independently repeat byte-for-byte. The",
  "option-2 decision is the prospectively frozen stability-plus-consensus rule;",
  "labels, outcomes, method/weight selection, biology, and claims remain closed.",
  "", paste0("Validation: ", sum(validation$passed), "/",
             nrow(validation), " checks pass.")
), file.path(output, "MV12C_HISTORICAL_FUSION_CLOSURE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv12c-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv12c_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv12c-artifact-manifest.csv"))
cat("Completed MV12-C closure; checks=", nrow(validation), "\n", sep = "")
