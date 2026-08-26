#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv11b_cell_benchmark_prefreeze.R <matrix-bundle> <output>",
  "<execution-head> <prior-failure-audit>"
), call. = FALSE)
bundle_path <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
head <- tolower(trimws(args[[3L]]))
failure_audit <- normalizePath(args[[4L]], mustWork = TRUE)
if (dir.exists(output)) stop("MV11-B output already exists", call. = FALSE)

source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
if (length(head) != 1L || !grepl("^[0-9a-f]{40}$", head)) {
  stop("cannot bind MV11-B execution head", call. = FALSE)
}
expected_bundle_sha <-
  "beb58777197545ec7113898e6e1082cafb61f84b446de973fbdd5431c791774e"
expected_bundle_bytes <- 3599074
if (sha(bundle_path) != expected_bundle_sha ||
    as.numeric(file.info(bundle_path)$size) != expected_bundle_bytes) {
  stop("MV11-B source bundle drift", call. = FALSE)
}
bundle <- readRDS(bundle_path)
catalog <- mv11_cell_catalog_v1(bundle)
failure_files <- c(
  "MV11C_SENTINEL_ATTEMPT1_FAILURE_2026-08-25.md",
  "mv11c-failure-evidence.csv", "mv11c-failure-validation.csv"
)
failure_paths <- file.path(failure_audit, failure_files)
expected_failure_bytes <- c(1011, 486, 868)
expected_failure_sha <- c(
  "7685a91603612e5e9b72d6c79b995eb1ecf297a809d14df54d0d7d8993218979",
  "88bb1307a5295b933696cad84b77802fa1cead1ed7088469036e3d275699e9b6",
  "92012e6726602d796409ce47af342c202af7c5aeb74a482557ca38c1450f591a"
)
if (!all(file.exists(failure_paths)) ||
    !identical(as.numeric(file.info(failure_paths)$size),
               expected_failure_bytes) ||
    !identical(unname(vapply(failure_paths, sha, character(1L))),
               expected_failure_sha)) {
  stop("MV11-B prior failure evidence drift", call. = FALSE)
}
implementation_files <- c(
  "R/mv05_benchmark_contract.R", "R/mv05n_clustering_gate.R",
  "R/mv08z_landscape_production.R", "R/mv10_clustering_benchmark.R",
  "R/mv11_cell_benchmark.R", "scripts/run_mv11_cell_matrix_worker.R",
  "scripts/run_mv11c_cell_benchmark_sentinel.R",
  "scripts/build_mv11d_cell_benchmark_sentinel_closure.R",
  "scripts/run_mv11e_full_cell_benchmark.R",
  "scripts/build_mv11f_full_cell_benchmark_closure.R"
)
if (!all(file.exists(implementation_files))) {
  stop("MV11-B implementation stack is incomplete", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv11b_implementation_binding_v1",
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_binding <- data.frame(
  contract_id = "mv11b_source_binding_v1",
  source_id = "mv07i_historical_selectedfit_matrix_bundle",
  public_path = "private_historical_matrix_bundle_not_published",
  bytes = expected_bundle_bytes, sha256 = expected_bundle_sha,
  samples = 124L, seeds = 5L, admitted_components = "cell_H0;cell_H1",
  excluded_components = "cell_H0_H1_secondary;all_gene_components;all_medians",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
failure_binding <- data.frame(
  contract_id = "mv11b_prior_failure_binding_v2",
  artifact = failure_files, bytes = expected_failure_bytes,
  sha256 = expected_failure_sha, attempt = 1L,
  disposition = "fail_closed_before_clustering_output",
  stringsAsFactors = FALSE
)
queue <- catalog
queue$elapsed_cap_seconds <- 600
queue$rss_cap_bytes <- 4 * 1024^3
queue$aggregate_elapsed_cap_seconds <- 3600
queue$private_storage_cap_bytes <- 256 * 1024^2
queue$workers <- 1L
queue$automatic_retries <- 0L
queue$labels_allowed <- FALSE
queue$outcomes_allowed <- FALSE
queue$cross_view_comparison_allowed <- FALSE
sentinel_order <- queue$catalog_order[
  queue$homology_dimension == "H1" & queue$seed == max(.mv11_required_seeds)
]
contract <- data.frame(
  contract_id = "mv11b_cell_benchmark_prefreeze_v1",
  execution_head = head, execution_attempt = 2L,
  prior_attempt_disposition = "catalog_type_mismatch_before_output",
  recovery_change = "canonical_character_catalog_equality_only",
  matrices = 10L, samples_per_matrix = 124L,
  seeds = 5L, homology_dimensions = "H0;H1_separate",
  methods = 5L, k_grid = "2:10", partition_fits = 450L,
  private_assignment_rows = 55800L, public_quality_rows = 450L,
  public_stability_rows = 90L, public_primary_k_rows = 2L,
  public_method_agreement_rows = 900L,
  sentinel_catalog_order = sentinel_order,
  sentinel_partition_fits = 45L, sentinel_private_assignment_rows = 5580L,
  sentinel_execution_authorized_after_commit = TRUE,
  full_execution_authorized = FALSE, automatic_retries = 0L,
  H0_H1_combined = FALSE, cell_gene_combined = FALSE,
  labels_allowed = FALSE, outcomes_allowed = FALSE,
  inference_allowed = FALSE, biological_claims_allowed = FALSE,
  manuscript_claims_allowed = FALSE, stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv11b_decision_v1",
  sentinel_execution_authorized_after_commit = TRUE,
  full_execution_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  cross_view_comparison_authorized = FALSE,
  fusion_authorized = FALSE, biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_action = paste(
    "commit this prefreeze; run exactly one H1 seed sentinel; independently",
    "repeat it; authorize the ten-matrix execution only after exact closure"
  ), stringsAsFactors = FALSE
)
checks <- c(
  bundle_exact_hash = sha(bundle_path) == expected_bundle_sha,
  bundle_exact_bytes = as.numeric(file.info(bundle_path)$size) ==
    expected_bundle_bytes,
  bundle_outer_contract = identical(bundle$contract_id,
                                      "mv07i_matrix_bundle_v1"),
  prior_failure_three_files_bound = nrow(failure_binding) == 3L,
  prior_failure_hashes_exact = all(vapply(failure_paths, sha, character(1L)) ==
                                     failure_binding$sha256),
  execution_attempt_two = contract$execution_attempt == 2L,
  ten_matrices = nrow(queue) == 10L,
  five_seeds = identical(sort(unique(queue$seed)), .mv11_required_seeds),
  H0_H1_separate = setequal(queue$homology_dimension, c("H0", "H1")),
  cell_only = all(queue$view_kind == "cell_topology_v1"),
  exact500_only = all(queue$panel_id == "exact500"),
  complete_124_axes = all(queue$units == 124L),
  complete_pair_axes = all(queue$unordered_pairs == choose(124L, 2L)),
  five_methods = all(queue$methods == 5L),
  complete_k_grid = all(queue$k_grid == "2:10"),
  four_hundred_fifty_fits = contract$partition_fits == 450L,
  fifty_five_thousand_eight_hundred_assignments =
    contract$private_assignment_rows == 55800L,
  four_hundred_fifty_quality = contract$public_quality_rows == 450L,
  ninety_stability = contract$public_stability_rows == 90L,
  two_primary_k = contract$public_primary_k_rows == 2L,
  nine_hundred_agreement = contract$public_method_agreement_rows == 900L,
  one_worker = all(queue$workers == 1L),
  zero_retries = all(queue$automatic_retries == 0L),
  resource_caps_positive = all(queue$elapsed_cap_seconds > 0 &
                                 queue$rss_cap_bytes > 0),
  implementation_bound = all(implementation$bytes > 0 &
                               grepl("^[0-9a-f]{64}$", implementation$sha256)),
  sentinel_unique = length(sentinel_order) == 1L,
  full_execution_closed = !contract$full_execution_authorized,
  labels_outcomes_comparison_closed =
    all(!queue$labels_allowed & !queue$outcomes_allowed &
          !queue$cross_view_comparison_allowed)
)
validation <- data.frame(
  contract_id = "mv11b_validation_v1", check_id = names(checks),
  passed = unname(checks), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV11-B prefreeze validation failed: ",
       paste(validation$check_id[!validation$passed], collapse = ", "),
       call. = FALSE)
}
dir.create(output, recursive = TRUE)
atomic(contract, file.path(output, "mv11b-contract.csv"))
atomic(source_binding, file.path(output, "mv11b-source-binding.csv"))
atomic(failure_binding, file.path(output, "mv11b-prior-failure-binding.csv"))
atomic(queue, file.path(output, "mv11b-workload-queue.csv"))
atomic(implementation, file.path(output, "mv11b-implementation-bindings.csv"))
atomic(decision, file.path(output, "mv11b-decision.csv"))
atomic(validation, file.path(output, "mv11b-validation.csv"))
readme <- c(
  "# MV11-B matched historical cell benchmark prefreeze", "",
  paste0("Execution head: `", head, "`."), "",
  "This attempt-2 prospective contract binds the preserved fail-closed first",
  "attempt and admits exactly ten immutable historical cell",
  "H0/H1 matrices (124 samples, five seeds), the unchanged MV10 five-method",
  "K=2:10 grid, and one H1 sentinel after commit. It does not admit labels,",
  "outcomes, cross-view comparison, fusion, inference, biology, or claims.", "",
  paste0("Validation: ", sum(validation$passed), "/", nrow(validation),
         " checks pass.")
)
writeLines(readme, file.path(output, "MV11B_CELL_BENCHMARK_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv11b-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv11b_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv11b-artifact-manifest.csv"))
cat("Completed MV11-B prefreeze; checks=", nrow(validation), "\n", sep = "")
