#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv12a_historical_fusion_prefreeze.R <matrix-bundle>",
  "<output> <execution-head>"
), call. = FALSE)
bundle_path <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]; execution_head <- tolower(trimws(args[[3L]]))
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("invalid MV12-A output or execution head", call. = FALSE)
}
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv10_clustering_benchmark.R")
source("R/mv11_cell_benchmark.R")
source("R/mv12_historical_fusion.R")
sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
bundle_sha <-
  "beb58777197545ec7113898e6e1082cafb61f84b446de973fbdd5431c791774e"
bundle_bytes <- 3599074
if (sha(bundle_path) != bundle_sha ||
    as.numeric(file.info(bundle_path)$size) != bundle_bytes) {
  stop("MV12-A matrix bundle drift", call. = FALSE)
}
bundle <- readRDS(bundle_path)
mv11_validate_matrix_bundle_v1(bundle)
closure_files <- c(
  "docs/audits/mv11f-cell-benchmark-closure-v1/mv11f-artifact-manifest.csv",
  "docs/audits/mv11j-cross-view-closure-v1/mv11j-artifact-manifest.csv",
  "docs/audits/mv11l-cross-view-review-v1/mv11l-validation.csv"
)
if (!all(file.exists(closure_files))) stop("MV12-A closures absent")
closures <- data.frame(
  contract_id = "mv12a_closure_binding_v1",
  closure_id = c("cell_benchmark", "cross_view_agreement", "review"),
  file = closure_files,
  bytes = as.numeric(file.info(closure_files)$size),
  sha256 = vapply(closure_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv05_benchmark_contract.R", "R/mv05n_clustering_gate.R",
  "R/mv10_clustering_benchmark.R", "R/mv11_cell_benchmark.R",
  "R/mv12_historical_fusion.R",
  "scripts/run_mv12b_historical_fusion.R",
  "scripts/build_mv12c_historical_fusion_closure.R"
)
if (!all(file.exists(implementation_files))) stop("MV12-A implementation absent")
implementation <- data.frame(
  contract_id = "mv12a_implementation_binding_v1",
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_binding <- data.frame(
  contract_id = "mv12a_source_binding_v1",
  source_id = "mv07i_historical_dual_view_matrix_bundle",
  private_path_not_published =
    "tmp/mv07i-label-closed-v1/artifacts/matrix-bundle.rds",
  bytes = bundle_bytes, sha256 = bundle_sha, samples = 124L, seeds = 5L,
  admitted_components = "cell_H0;cell_H1;gene_H0;gene_H1",
  excluded_components = "H0_H1_composites;medians",
  labels_used = FALSE, outcomes_used = FALSE, stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv12a_historical_fusion_prefreeze_v1",
  execution_head = execution_head, matrices = 50L, samples = 124L,
  seeds = 5L, homology_dimensions = "H0;H1_separate",
  gene_weights = "0;0.25;0.5;0.75;1", primary_gene_weight = 0.5,
  methods = 5L, primary_method = .mv12_primary_method,
  common_k = "2;3", primary_dimension = "H1",
  partition_fits = 500L, private_assignment_rows = 62000L,
  public_scale_rows = 20L, public_catalog_rows = 50L,
  public_quality_rows = 500L, public_stability_rows = 100L,
  public_consensus_rows = 300L, public_primary_detail_rows = 2L,
  public_disposition_rows = 1L,
  elapsed_cap_seconds = 1800, private_storage_cap_bytes = 128 * 1024^2,
  public_storage_cap_bytes = 16 * 1024^2,
  workers = 1L, automatic_retries = 0L,
  execution_authorized_after_commit = TRUE,
  labels_allowed = FALSE, outcomes_allowed = FALSE,
  H0_H1_combination_allowed = FALSE, method_or_weight_selection_allowed = FALSE,
  biological_claims_allowed = FALSE, manuscript_claims_allowed = FALSE,
  stringsAsFactors = FALSE
)
trigger <- data.frame(
  contract_id = "mv12a_option2_trigger_v1",
  primary_axis = c("H1_PAM_K2", "H1_PAM_K3"), k = c(2L, 3L),
  stability_rule = "equal_weight_mean_seed_ARI_strictly_exceeds_both_components",
  consensus_rule =
    "equal_weight_min_component_ARI_minus_cell_gene_ARI_positive_in_at_least_3_of_5_seeds",
  primary_pass_requires_both_rules = TRUE,
  sensitivity_weights = "0.25;0.75_diagnostic_only",
  option2_required_unless =
    "both_primary_K_fail_and_no_sensitivity_weight_meets_both_rules",
  labels_used = FALSE, outcomes_used = FALSE, stringsAsFactors = FALSE
)
checks <- c(
  source_hash_exact = sha(bundle_path) == bundle_sha,
  source_bytes_exact = as.numeric(file.info(bundle_path)$size) == bundle_bytes,
  source_bundle_valid = identical(bundle$contract_id, "mv07i_matrix_bundle_v1"),
  three_closures_bound = nrow(closures) == 3L,
  closure_hashes_exact = all(vapply(closure_files, sha, character(1L)) ==
                               closures$sha256),
  implementation_bound = all(implementation$bytes > 0),
  fifty_matrices = contract$matrices == 50L,
  five_hundred_fits = contract$partition_fits == 500L,
  sixty_two_thousand_assignments = contract$private_assignment_rows == 62000L,
  public_cardinalities = all(unlist(contract[c(
    "public_scale_rows", "public_catalog_rows", "public_quality_rows",
    "public_stability_rows", "public_consensus_rows",
    "public_primary_detail_rows", "public_disposition_rows"
  )]) == c(20L, 50L, 500L, 100L, 300L, 2L, 1L)),
  H0_H1_separate = !contract$H0_H1_combination_allowed,
  equal_weight_primary = contract$primary_gene_weight == 0.5,
  PAM_primary = contract$primary_method == .mv12_primary_method,
  H1_primary = contract$primary_dimension == "H1",
  common_K_only = contract$common_k == "2;3",
  trigger_two_rows = nrow(trigger) == 2L,
  majority_consensus_rule = all(grepl("3_of_5", trigger$consensus_rule)),
  sensitivity_diagnostic_only = all(grepl("diagnostic_only",
                                           trigger$sensitivity_weights)),
  option2_rule_fixed = length(unique(trigger$option2_required_unless)) == 1L,
  one_worker_zero_retry = contract$workers == 1L &&
    contract$automatic_retries == 0L,
  positive_resource_caps = contract$elapsed_cap_seconds > 0 &&
    contract$private_storage_cap_bytes > 0 &&
    contract$public_storage_cap_bytes > 0,
  labels_outcomes_closed = !contract$labels_allowed && !contract$outcomes_allowed,
  selection_closed = !contract$method_or_weight_selection_allowed,
  claims_closed = !contract$biological_claims_allowed &&
    !contract$manuscript_claims_allowed,
  execution_after_commit = contract$execution_authorized_after_commit
)
validation <- data.frame(
  contract_id = "mv12a_validation_v1", check_id = names(checks),
  passed = unname(checks), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV12-A prefreeze validation failed")
decision <- data.frame(
  contract_id = "mv12a_decision_v1",
  historical_fusion_execution_authorized_after_commit = TRUE,
  option2_authorized_only_if_triggered = TRUE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  method_or_weight_selection_authorized = FALSE,
  biological_claims_authorized = FALSE, manuscript_claims_authorized = FALSE,
  next_action = paste(
    "commit prefreeze; execute 500 fixed fits once; independently repeat all",
    "scientific artifacts; apply frozen option-2 trigger without tuning"
  ), stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
atomic(contract, file.path(output, "mv12a-contract.csv"))
atomic(source_binding, file.path(output, "mv12a-source-binding.csv"))
atomic(closures, file.path(output, "mv12a-closure-bindings.csv"))
atomic(implementation, file.path(output, "mv12a-implementation-bindings.csv"))
atomic(trigger, file.path(output, "mv12a-option2-trigger.csv"))
atomic(decision, file.path(output, "mv12a-decision.csv"))
atomic(validation, file.path(output, "mv12a-validation.csv"))
writeLines(c(
  "# MV12-A historical fusion feasibility prefreeze", "",
  "This contract freezes separate H0/H1 median-normalized cell/gene distance",
  "fusion at five fixed weights. Equal-weight PAM H1 at K=2/3 is primary.",
  "Option 2 is triggered deterministically by stability and balanced consensus,",
  "never by labels or outcomes. No method/weight selection or claims are open.",
  "", paste0("Validation: ", sum(validation$passed), "/",
             nrow(validation), " checks pass.")
), file.path(output, "MV12A_HISTORICAL_FUSION_PREFREEZE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv12a-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv12a_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv12a-artifact-manifest.csv"))
cat("Completed MV12-A prefreeze; checks=", nrow(validation), "\n", sep = "")
