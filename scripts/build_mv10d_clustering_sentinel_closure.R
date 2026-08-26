#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: build_mv10d_clustering_sentinel_closure.R <prefreeze>",
  "<mv07h-root> <mv08zu-private-root> <mv08zx-private-root>",
  "<sentinel-private> <sentinel-public> <output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
private <- normalizePath(args[[5L]], mustWork = TRUE)
public <- normalizePath(args[[6L]], mustWork = TRUE)
output <- args[[7L]]
if (dir.exists(output)) stop("MV10-D closure output already exists")

source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv10_clustering_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10b-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv10c-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10b-contract.csv"))
sentinel <- readc(file.path(prefreeze, "mv10b-sentinel-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv10b-implementation-bindings.csv"))
receipt <- readc(file.path(public, "mv10c-resource-receipt.csv"))
saved_quality <- readc(file.path(public, "mv10c-partition-quality.csv"))
job <- file.path(private, "job")
saved_partitions <- readc(file.path(job, "partitions.csv"))
worker_status <- readc(file.path(job, "status.csv"))
if (nrow(contract) != 1L || nrow(sentinel) != 1L ||
    receipt$execution_head != contract$execution_head ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) {
  stop("MV10-D closure binding drift", call. = FALSE)
}
loaded <- mv08zy_read_distance_stack_v1(
  sentinel, roots[[1L]], roots[[2L]], roots[[3L]]
)
matrix <- mv10_distance_matrix_v1(loaded$pairs)
fresh <- mv10_partition_grid_v1(matrix)
metadata <- data.frame(
  execution_head = contract$execution_head,
  catalog_order = sentinel$catalog_order,
  stack_id = sentinel$stack_id,
  representation_id = sentinel$representation_id,
  panel_id = sentinel$panel_id, seed = as.integer(sentinel$seed),
  homology_dimension = sentinel$homology_dimension,
  payload_set_sha256 = sentinel$payload_set_sha256,
  pair_axis_sha256 = sentinel$pair_axis_sha256,
  stringsAsFactors = FALSE
)
fresh_partitions <- cbind(
  metadata[rep(1L, nrow(fresh$partitions)), , drop = FALSE],
  fresh$partitions
)
fresh_quality <- cbind(
  metadata[rep(1L, nrow(fresh$quality)), , drop = FALSE], fresh$quality
)
rownames(fresh_partitions) <- NULL; rownames(fresh_quality) <- NULL
partition_identity <- isTRUE(all.equal(
  saved_partitions, fresh_partitions, tolerance = 0,
  check.attributes = FALSE
))
quality_identity <- isTRUE(all.equal(
  saved_quality, fresh_quality, tolerance = 1e-14,
  check.attributes = FALSE
))
private_files <- file.path(job, c("partitions.csv", "quality.csv", "status.csv"))
rehash <- data.frame(
  contract_id = "mv10d_private_rehash_v1",
  artifact = basename(private_files),
  bytes = as.numeric(file.info(private_files)$size),
  sha256 = vapply(private_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10d_validation_v1",
  check_id = c(
    "prefreeze_manifest", "sentinel_manifest", "execution_head",
    "implementation_bindings", "corrected_primary_H1_sentinel",
    "highest_frozen_seed", "source_payload", "source_pair_axis",
    "one_hundred_twenty_four_units", "seven_thousand_six_hundred_twenty_six_pairs",
    "forty_five_fits", "five_methods", "complete_k2_k10",
    "five_thousand_five_hundred_eighty_private_rows", "forty_five_quality_rows",
    "private_partition_exact_repeat", "quality_numeric_repeat",
    "worker_status_complete", "positive_resource_observation",
    "child_elapsed_cap", "child_rss_cap",
    "private_storage_cap", "one_worker_zero_retry", "empty_stderr",
    "label_outcome_firewall", "claim_firewall"
  ),
  passed = c(
    TRUE, TRUE, receipt$execution_head == contract$execution_head,
    all(vapply(implementation$file, sha, character(1L)) ==
          implementation$sha256),
    sentinel$stack_id == "allqc_residual_exact500" &&
      sentinel$homology_dimension == "H1",
    sentinel$seed == max(.mv10_required_seeds),
    loaded$payload_set_sha256 == sentinel$payload_set_sha256,
    loaded$pair_axis_sha256 == sentinel$pair_axis_sha256,
    nrow(matrix) == 124L, nrow(loaded$pairs) == 7626L,
    receipt$partition_fits == 45L,
    length(unique(saved_partitions$method_id)) == 5L,
    identical(sort(unique(saved_partitions$k)), .mv10_k_grid),
    nrow(saved_partitions) == 5580L, nrow(saved_quality) == 45L,
    partition_identity, quality_identity,
    worker_status$completion_state == "complete",
    receipt$elapsed_seconds > 0 && receipt$peak_process_tree_rss_bytes > 0,
    receipt$elapsed_seconds <= receipt$elapsed_cap_seconds,
    receipt$peak_process_tree_rss_bytes <= receipt$process_tree_rss_cap_bytes,
    receipt$private_bytes <= receipt$private_storage_cap_bytes,
    receipt$workers == 1L && receipt$retries == 0L,
    receipt$stderr_bytes == 0,
    !receipt$labels_used && !receipt$outcomes_used,
    !receipt$biological_claims && !receipt$manuscript_claims
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-D sentinel closure failed")
decision <- data.frame(
  contract_id = "mv10d_decision_v1",
  decision = "admit_full_label_closed_clustering_execution_after_commit",
  full_execution_authorized_after_commit = TRUE,
  maximum_workers = 1L, automatic_retries = 0L,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  biological_interpretation_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_stage = "mv10e_full_execution_then_mv10f_independent_closure",
  stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(rehash, file.path(output, "mv10d-private-rehash.csv"))
atomic(validation, file.path(output, "mv10d-validation.csv"))
atomic(decision, file.path(output, "mv10d-decision.csv"))
writeLines(c(
  "# MV10-D clustering sentinel closure", "",
  "The prospectively selected corrected-primary H1 matrix completed all five",
  "methods and K=2:10 under resource caps. Independent recomputation exactly",
  "reproduced private partitions and public quality. Full label-closed serial",
  "execution is authorized only after this closure is committed."
), file.path(output, "MV10D_CLUSTERING_SENTINEL_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10d-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10d_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10d-artifact-manifest.csv"))
cat("Closed MV10-D clustering sentinel; checks=26\n")
