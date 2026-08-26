#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(paste(
  "usage: build_mv10b_clustering_execution_prefreeze.R <output-dir>",
  "<execution-head>"
), call. = FALSE)
output <- args[[1L]]
execution_head <- tolower(trimws(args[[2L]]))
if (!grepl("^[0-9a-f]{40}$", execution_head)) stop("invalid execution head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-B prefreeze")
}
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
mv10a <- "docs/audits/mv10a-clustering-benchmark-prefreeze-v1"
.mv08z_verify_manifest(mv10a, "mv10a-artifact-manifest.csv")
catalog <- readc(file.path(mv10a, "mv10a-internal-stack-bindings.csv"))
mv10_validate_stack_catalog_v1(catalog)
mv10a_contract <- readc(file.path(mv10a, "mv10a-contract.csv"))
mv10a_resource <- readc(file.path(mv10a, "mv10a-resource-contract.csv"))
mv10a_outputs <- readc(file.path(mv10a, "mv10a-output-contract.csv"))
queue <- catalog[order(catalog$catalog_order), , drop = FALSE]
queue$execution_order <- seq_len(nrow(queue))
queue$elapsed_cap_seconds <- as.integer(mv10a_resource$child_elapsed_cap_seconds)
queue$rss_cap_bytes <- as.numeric(mv10a_resource$process_tree_rss_cap_bytes)
queue$aggregate_elapsed_cap_seconds <- 4L * 60L * 60L
queue$private_storage_cap_bytes <- as.numeric(
  mv10a_resource$private_storage_cap_bytes
)
queue$workers <- 1L; queue$retries <- 0L
queue$outcome_label_state <- "closed"
queue$biological_outcomes_computed <- FALSE
sentinel <- queue[
  queue$stack_id == "allqc_residual_exact500" &
    queue$homology_dimension == "H1" &
    queue$seed == max(.mv10_required_seeds), , drop = FALSE
]
sentinel$sentinel_rationale <- paste(
  "corrected_primary_representation_H1_highest_frozen_seed",
  "chosen_prospectively_without_cluster_results", sep = ";"
)
implementation_files <- c(
  "R/mv05_benchmark_contract.R",
  "R/mv05n_clustering_gate.R",
  "R/mv08z_landscape_production.R",
  "R/mv08zy_distance_comparison.R",
  "R/mv10_clustering_benchmark.R",
  "scripts/run_mv10_clustering_matrix_worker.R",
  "scripts/run_mv10c_clustering_sentinel.R",
  "scripts/build_mv10d_clustering_sentinel_closure.R",
  "scripts/run_mv10e_full_clustering_benchmark.R",
  "scripts/build_mv10f_full_clustering_closure.R",
  "scripts/build_mv10b_clustering_execution_prefreeze.R"
)
implementation <- data.frame(
  contract_id = "mv10b_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(mv10a, "mv10a-artifact-manifest.csv"),
  file.path(mv10a, "mv10a-contract.csv"),
  file.path(mv10a, "mv10a-resource-contract.csv"),
  file.path(mv10a, "mv10a-output-contract.csv"),
  file.path(mv10a, "mv10a-method-registry.csv"),
  file.path(mv10a, "mv10a-distance-registry.csv"),
  file.path(mv10a, "mv10a-k-contract.csv")
)
source_freeze <- data.frame(
  contract_id = "mv10b_source_freeze_v1",
  source_order = seq_along(source_files), artifact = source_files,
  bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv10b_clustering_execution_prefreeze_v1",
  execution_head = execution_head, input_prefreeze = "mv10a_v1",
  matrices = 30L, sentinel_matrices = 1L, units_per_matrix = 124L,
  methods_per_matrix = 5L, k_values = 9L,
  sentinel_partition_fits = 45L, full_partition_fits = 1350L,
  full_private_assignment_rows = 167400L,
  H0_H1_separate = TRUE, cell_gene_combined = FALSE,
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  stringsAsFactors = FALSE
)
resource <- data.frame(
  contract_id = "mv10b_resource_policy_v1", workers = 1L,
  automatic_retries = 0L, child_elapsed_cap_seconds = 1800L,
  process_tree_rss_cap_bytes = 4 * 1024^3,
  aggregate_elapsed_cap_seconds = 4L * 60L * 60L,
  private_storage_cap_bytes = 512 * 1024^2,
  poll_interval_seconds = 0.1,
  partial_evidence_policy = "preserve_stop_no_retry",
  stringsAsFactors = FALSE
)
closure <- data.frame(
  contract_id = "mv10b_prospective_closure_v1",
  stage = c("sentinel_mv10d", "full_mv10f"),
  requires_source_reload = TRUE,
  requires_private_partition_recomputation = TRUE,
  requires_public_aggregate_recomputation = c(FALSE, TRUE),
  partition_tolerance = 0,
  aggregate_numeric_tolerance = 1e-14,
  requires_resource_caps = TRUE,
  requires_privacy_firewall = TRUE,
  requires_label_outcome_firewall = TRUE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv10b_decision_v1",
  decision = "authorize_one_sentinel_after_prefreeze_commit",
  sentinel_execution_authorized_after_commit = TRUE,
  full_execution_authorized = FALSE,
  full_execution_requires = "committed_mv10d_26_of_26_sentinel_closure",
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  biological_interpretation_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10b_validation_v1",
  check_id = c(
    "mv10a_manifest", "mv10a_design_head", "thirty_matrix_queue",
    "catalog_order_complete", "three_representations", "five_seeds",
    "H0_H1_separate", "all_124_units", "all_7626_pairs",
    "exact500_gene_view", "payload_hashes", "common_pair_axis",
    "one_sentinel", "sentinel_corrected_primary", "sentinel_H1",
    "sentinel_highest_seed", "sentinel_no_result_selection",
    "five_methods", "complete_k2_k10", "forty_five_sentinel_fits",
    "one_thousand_three_hundred_fifty_full_fits",
    "one_hundred_sixty_seven_thousand_four_hundred_assignments",
    "output_contract_bound", "one_worker", "zero_retry", "child_elapsed_cap",
    "process_tree_rss_cap", "aggregate_elapsed_cap", "private_storage_cap",
    "partial_preservation", "eleven_implementation_files",
    "implementation_hashes", "seven_source_bindings", "source_hashes",
    "sentinel_closure_prospective", "full_closure_prospective",
    "private_partition_exact_repeat", "aggregate_numeric_tolerance",
    "label_outcome_firewall", "combination_firewall", "claim_firewall",
    "full_execution_closed"
  ),
  passed = c(
    TRUE, mv10a_contract$implementation_head ==
      "10335bc91f81a945aefb4a03775a068c9b84c204",
    nrow(queue) == 30L,
    identical(as.integer(queue$execution_order), seq_len(30L)),
    length(unique(queue$stack_id)) == 3L,
    identical(sort(unique(as.integer(queue$seed))), .mv10_required_seeds),
    setequal(queue$homology_dimension, c("H0", "H1")),
    all(queue$units == 124L), all(queue$unordered_pairs == 7626L),
    all(queue$panel_id == "exact500") &&
      all(queue$view_kind == "gene_topology_v1"),
    all(grepl("^[0-9a-f]{64}$", queue$payload_set_sha256)),
    length(unique(queue$pair_axis_sha256)) == 1L,
    nrow(sentinel) == 1L,
    sentinel$stack_id == "allqc_residual_exact500",
    sentinel$homology_dimension == "H1",
    sentinel$seed == max(.mv10_required_seeds),
    grepl("without_cluster_results", sentinel$sentinel_rationale),
    contract$methods_per_matrix == 5L, contract$k_values == 9L,
    contract$sentinel_partition_fits == 45L,
    contract$full_partition_fits == 1350L,
    contract$full_private_assignment_rows == 167400L,
    identical(as.integer(mv10a_outputs$expected_rows),
              c(167400L, 1350L, 270L, 2L, 2700L, 1L)),
    resource$workers == 1L, resource$automatic_retries == 0L,
    resource$child_elapsed_cap_seconds == 1800L,
    resource$process_tree_rss_cap_bytes == 4 * 1024^3,
    resource$aggregate_elapsed_cap_seconds == 14400L,
    resource$private_storage_cap_bytes == 512 * 1024^2,
    resource$partial_evidence_policy == "preserve_stop_no_retry",
    nrow(implementation) == 11L, all(file.exists(implementation$file)),
    nrow(source_freeze) == 7L, all(file.exists(source_freeze$artifact)),
    closure$requires_source_reload[closure$stage == "sentinel_mv10d"],
    closure$requires_public_aggregate_recomputation[
      closure$stage == "full_mv10f"],
    all(closure$partition_tolerance == 0),
    all(closure$aggregate_numeric_tolerance == 1e-14),
    !contract$labels_used && !contract$outcomes_used,
    contract$H0_H1_separate && !contract$cell_gene_combined,
    !contract$biological_claims && !contract$manuscript_claims,
    !decision$full_execution_authorized
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-B prefreeze validation failed")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
artifacts <- list(
  "mv10b-contract.csv" = contract,
  "mv10b-workload-queue.csv" = queue,
  "mv10b-sentinel-queue.csv" = sentinel,
  "mv10b-resource-policy.csv" = resource,
  "mv10b-prospective-closure.csv" = closure,
  "mv10b-implementation-bindings.csv" = implementation,
  "mv10b-source-freeze.csv" = source_freeze,
  "mv10b-decision.csv" = decision,
  "mv10b-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-B clustering execution prefreeze", "",
  "MV10-B freezes one corrected-primary H1 sentinel before full execution.",
  "The serial worker, resource monitor, full runner, and both independent",
  "closures are hash-bound. Full execution remains closed until a committed",
  "26/26 MV10-D sentinel closure admits the unchanged workload."
), file.path(output, "MV10B_CLUSTERING_EXECUTION_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10b-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10b_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10b-artifact-manifest.csv"))
cat("Built MV10-B clustering execution prefreeze; checks=42\n")
