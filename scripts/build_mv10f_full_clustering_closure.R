#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv10f_full_clustering_closure.R <prefreeze> <admission>",
  "<mv07h-root> <mv08zu-private-root> <mv08zx-private-root>",
  "<production-private> <production-public> <output> <execution-head>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
admission <- normalizePath(args[[2L]], mustWork = TRUE)
roots <- vapply(args[3:5], normalizePath, character(1L), mustWork = TRUE)
private <- normalizePath(args[[6L]], mustWork = TRUE)
public <- normalizePath(args[[7L]], mustWork = TRUE)
output <- args[[8L]]
execution_head <- tolower(trimws(args[[9L]]))
if (dir.exists(output)) stop("MV10-F closure output already exists")

source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv10_clustering_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(prefreeze, "mv10b-artifact-manifest.csv")
.mv08z_verify_manifest(admission, "mv10d-artifact-manifest.csv")
.mv08z_verify_manifest(public, "mv10e-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv10b-contract.csv"))
queue <- readc(file.path(prefreeze, "mv10b-workload-queue.csv"))
implementation <- readc(file.path(prefreeze, "mv10b-implementation-bindings.csv"))
admit <- readc(file.path(admission, "mv10d-decision.csv"))
receipt <- readc(file.path(public, "mv10e-terminal-receipt.csv"))
ledger <- readc(file.path(public, "mv10e-resource-ledger.csv"))
saved_assignments <- readc(file.path(private, "mv10e-sample-partitions.csv"))
saved_quality <- readc(file.path(public, "mv10e-partition-quality.csv"))
saved_stability <- readc(file.path(public, "mv10e-seed-stability.csv"))
saved_primary_k <- readc(file.path(public, "mv10e-primary-k-selection.csv"))
saved_agreement <- readc(file.path(public, "mv10e-method-agreement.csv"))
truth <- .mv08z_truth
if (execution_head != contract$execution_head ||
    receipt$execution_head != execution_head ||
    !truth(admit$full_execution_authorized_after_commit) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha, character(1L)) ==
           implementation$sha256)) stop("MV10-F closure binding drift")

fresh_partitions <- list(); fresh_quality <- list(); rehash <- list()
for (i in seq_len(nrow(queue))) {
  binding <- queue[i, , drop = FALSE]
  loaded <- mv08zy_read_distance_stack_v1(
    binding, roots[[1L]], roots[[2L]], roots[[3L]]
  )
  matrix <- mv10_distance_matrix_v1(loaded$pairs)
  grid <- mv10_partition_grid_v1(matrix)
  metadata <- data.frame(
    execution_head = execution_head,
    catalog_order = binding$catalog_order,
    stack_id = binding$stack_id,
    representation_id = binding$representation_id,
    panel_id = binding$panel_id, seed = as.integer(binding$seed),
    homology_dimension = binding$homology_dimension,
    payload_set_sha256 = binding$payload_set_sha256,
    pair_axis_sha256 = binding$pair_axis_sha256,
    stringsAsFactors = FALSE
  )
  expected_partitions <- cbind(
    metadata[rep(1L, nrow(grid$partitions)), , drop = FALSE], grid$partitions
  )
  expected_quality <- cbind(
    metadata[rep(1L, nrow(grid$quality)), , drop = FALSE], grid$quality
  )
  rownames(expected_partitions) <- NULL; rownames(expected_quality) <- NULL
  job <- file.path(private, "jobs", sprintf("job_%02d", i))
  job_partitions <- readc(file.path(job, "partitions.csv"))
  job_quality <- readc(file.path(job, "quality.csv"))
  job_status <- readc(file.path(job, "status.csv"))
  partition_equal <- isTRUE(all.equal(
    job_partitions, expected_partitions, tolerance = 0,
    check.attributes = FALSE
  ))
  quality_equal <- isTRUE(all.equal(
    job_quality, expected_quality, tolerance = 1e-14,
    check.attributes = FALSE
  ))
  if (!partition_equal || !quality_equal) {
    stop("MV10-F independent job reconstruction drift at order ", i)
  }
  fresh_partitions[[i]] <- expected_partitions
  fresh_quality[[i]] <- expected_quality
  rehash[[i]] <- data.frame(
    contract_id = "mv10f_private_job_rehash_v1", execution_order = i,
    catalog_order = binding$catalog_order,
    partitions_sha256 = sha(file.path(job, "partitions.csv")),
    quality_sha256 = sha(file.path(job, "quality.csv")),
    status_sha256 = sha(file.path(job, "status.csv")),
    partition_exact_repeat = partition_equal,
    quality_numeric_repeat = quality_equal,
    worker_complete = job_status$completion_state == "complete",
    stringsAsFactors = FALSE
  )
}
fresh_partitions <- do.call(rbind, fresh_partitions)
fresh_quality <- do.call(rbind, fresh_quality)
rownames(fresh_partitions) <- NULL; rownames(fresh_quality) <- NULL
fresh_stability <- mv10_seed_stability_v1(fresh_partitions)
fresh_primary_k <- mv10_select_primary_k_v1(fresh_partitions)
fresh_agreement <- mv10_method_agreement_v1(fresh_partitions)
rehash <- do.call(rbind, rehash)
assignment_equal <- isTRUE(all.equal(
  saved_assignments, fresh_partitions, tolerance = 0,
  check.attributes = FALSE
))
quality_equal <- isTRUE(all.equal(
  saved_quality, fresh_quality, tolerance = 1e-14,
  check.attributes = FALSE
))
stability_equal <- isTRUE(all.equal(
  saved_stability, fresh_stability, tolerance = 1e-14,
  check.attributes = FALSE
))
primary_k_equal <- isTRUE(all.equal(
  saved_primary_k, fresh_primary_k, tolerance = 1e-14,
  check.attributes = FALSE
))
agreement_equal <- isTRUE(all.equal(
  saved_agreement, fresh_agreement, tolerance = 1e-14,
  check.attributes = FALSE
))
ledger_hashes <- all(vapply(seq_len(nrow(ledger)), function(i) {
  job <- file.path(private, "jobs", sprintf("job_%02d", i))
  identical(sha(file.path(job, "partitions.csv")),
            ledger$partitions_sha256[[i]]) &&
    identical(sha(file.path(job, "quality.csv")),
              ledger$quality_sha256[[i]]) &&
    identical(sha(file.path(job, "status.csv")), ledger$status_sha256[[i]])
}, logical(1L)))
validation <- data.frame(
  contract_id = "mv10f_validation_v1",
  check_id = c(
    "prefreeze_manifest", "admission_manifest", "production_manifest",
    "execution_head", "implementation_bindings", "admission_authorized",
    "terminal_complete", "thirty_matrices", "one_thousand_three_hundred_fifty_fits",
    "one_hundred_sixty_seven_thousand_four_hundred_assignments",
    "one_thousand_three_hundred_fifty_quality_rows",
    "two_hundred_seventy_stability_rows", "two_primary_k_rows",
    "two_thousand_seven_hundred_agreement_rows", "three_representations",
    "five_seeds", "H0_H1_separate", "five_methods", "complete_k2_k10",
    "thirty_private_job_rehashes", "all_worker_complete",
    "job_partition_exact_repeat", "job_quality_numeric_repeat",
    "aggregate_assignment_exact_repeat", "aggregate_quality_numeric_repeat",
    "seed_stability_numeric_repeat", "primary_k_numeric_repeat",
    "method_agreement_numeric_repeat", "ledger_private_hashes",
    "positive_resource_observations", "all_children_elapsed_cap",
    "all_children_rss_cap", "aggregate_elapsed_cap",
    "private_storage_cap", "one_worker_zero_retry", "empty_stderr",
    "public_sample_identity_firewall", "label_outcome_firewall",
    "combination_firewall", "claim_firewall"
  ),
  passed = c(
    TRUE, TRUE, TRUE, execution_head == contract$execution_head,
    all(vapply(implementation$file, sha, character(1L)) ==
          implementation$sha256),
    truth(admit$full_execution_authorized_after_commit),
    receipt$completion_state == "complete", receipt$matrices == 30L,
    receipt$partition_fits == 1350L,
    nrow(saved_assignments) == 167400L, nrow(saved_quality) == 1350L,
    nrow(saved_stability) == 270L, nrow(saved_primary_k) == 2L,
    nrow(saved_agreement) == 2700L,
    length(unique(saved_assignments$stack_id)) == 3L,
    identical(sort(unique(saved_assignments$seed)), .mv10_required_seeds),
    setequal(saved_assignments$homology_dimension, c("H0", "H1")),
    length(unique(saved_assignments$method_id)) == 5L,
    identical(sort(unique(saved_assignments$k)), .mv10_k_grid),
    nrow(rehash) == 30L, all(rehash$worker_complete),
    all(rehash$partition_exact_repeat), all(rehash$quality_numeric_repeat),
    assignment_equal, quality_equal, stability_equal, primary_k_equal,
    agreement_equal, ledger_hashes,
    all(ledger$elapsed_seconds > 0) &&
      all(ledger$peak_process_tree_rss_bytes > 0),
    all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds),
    all(ledger$peak_process_tree_rss_bytes <=
          ledger$process_tree_rss_cap_bytes),
    receipt$elapsed_seconds <= receipt$aggregate_elapsed_cap_seconds,
    receipt$private_bytes <= receipt$private_storage_cap_bytes,
    receipt$workers == 1L && receipt$retries == 0L,
    receipt$stderr_bytes == 0,
    !"sample_id" %in% names(saved_quality) &&
      !"sample_id" %in% names(saved_stability) &&
      !"sample_id" %in% names(saved_agreement),
    !receipt$labels_used && !receipt$outcomes_used,
    !receipt$H0_H1_combined && !receipt$cell_gene_combined,
    !receipt$biological_claims && !receipt$manuscript_claims
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-F full closure failed")
decision <- data.frame(
  contract_id = "mv10f_decision_v1",
  decision = "close_full_label_closed_clustering_benchmark",
  next_stage = "separate_label_closed_descriptive_synthesis_prefreeze",
  result_interpretation_state = "closed_pending_separate_review",
  labels_outcomes_state = "closed", biological_interpretation_state = "closed",
  manuscript_claims_state = "closed", stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(rehash, file.path(output, "mv10f-private-job-rehash.csv"))
atomic(validation, file.path(output, "mv10f-validation.csv"))
atomic(decision, file.path(output, "mv10f-decision.csv"))
writeLines(c(
  "# MV10-F full clustering benchmark closure", "",
  "All 30 distance matrices and 1,350 partition fits were independently",
  "recomputed from immutable distance sources. Private assignments and all",
  "public aggregate outputs reproduce under the prospective contract. Results",
  "remain label-closed and uninterpreted pending a separate review prefreeze."
), file.path(output, "MV10F_FULL_CLUSTERING_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10f-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10f_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10f-artifact-manifest.csv"))
cat("Closed MV10-F full clustering benchmark; checks=40\n")
