#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv13d_allqc_cell_full_prefreeze.R <mv13a> <mv13c>",
  "<private-locator> <sentinel-private> <output> <execution-head>"
), call. = FALSE)
mv13a <- normalizePath(args[[1L]], mustWork = TRUE)
mv13c <- normalizePath(args[[2L]], mustWork = TRUE)
locator <- normalizePath(args[[3L]], mustWork = TRUE)
sentinel <- normalizePath(args[[4L]], mustWork = TRUE)
output <- args[[5L]]; head <- tolower(trimws(args[[6L]]))
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", head)) {
  stop("invalid MV13-D output or execution head")
}
source("R/mv05_benchmark_contract.R"); source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R"); source("R/mv13_allqc_cell_topology.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
.mv08z_verify_manifest(mv13a, "mv13a-artifact-manifest.csv")
.mv08z_verify_manifest(mv13c, "mv13c-artifact-manifest.csv")
a_contract <- readc(file.path(mv13a, "mv13a-contract.csv"))
pca <- readc(file.path(mv13a, "mv13a-pca-queue.csv"))
views <- readc(file.path(mv13a, "mv13a-view-queue.csv"))
locator_binding <- readc(file.path(mv13a, "mv13a-private-locator-binding.csv"))
c_decision <- readc(file.path(mv13c, "mv13c-decision.csv"))
if (!c_decision$sentinel_independently_closed ||
    !c_decision$full_PCA_PH_prefreeze_eligible_next ||
    c_decision$full_execution_authorized_by_this_closure ||
    sha(locator) != locator_binding$sha256 || nrow(pca) != 7L ||
    nrow(views) != 636L) stop("MV13-D prerequisite drift")
saved <- readRDS(sentinel)
if (saved$contract_id != "mv13b_allqc_cell_sentinel_v1") {
  stop("MV13-D sentinel payload drift")
}
groups <- pca
groups$group_id <- paste(groups$dataset_scope, groups$panel_id, groups$seed,
                         sep = "__")
groups$unit_count <- ifelse(groups$dataset_scope == "internal124", 124L, 8L)
groups$adopt_closed_model <- groups$dataset_scope == "internal124" &
  groups$panel_id == "exact500" & groups$seed == 20260809L
groups$adopt_closed_view_count <- as.integer(groups$adopt_closed_model)
groups$new_model_count <- as.integer(!groups$adopt_closed_model)
groups$new_view_count <- groups$unit_count - groups$adopt_closed_view_count
groups$group_order <- seq_len(nrow(groups))
groups <- groups[, c(
  "group_order", "group_id", "dataset_scope", "panel_id", "seed",
  "representation_id", "pca_components", "unit_count",
  "adopt_closed_model", "adopt_closed_view_count", "new_model_count",
  "new_view_count"
)]
private_binding <- data.frame(
  contract_id = "mv13d_private_binding_v1",
  artifact_id = c("private_locator", "sentinel_payload"),
  bytes = as.numeric(file.info(c(locator, sentinel))$size),
  sha256 = vapply(c(locator, sentinel), sha, character(1L)),
  publication_state = "private_not_tracked", stringsAsFactors = FALSE
)
prior <- c(
  file.path(mv13a, "mv13a-artifact-manifest.csv"),
  file.path(mv13c, "mv13c-artifact-manifest.csv")
)
prior_binding <- data.frame(
  contract_id = "mv13d_prior_binding_v1", artifact_id = c("MV13_A", "MV13_C"),
  file = gsub("\\\\", "/", prior), bytes = as.numeric(file.info(prior)$size),
  sha256 = vapply(prior, sha, character(1L)), stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/dual_view_topology.R", "R/mv13_allqc_cell_topology.R",
  "scripts/build_mv13d_allqc_cell_full_prefreeze.R",
  "scripts/run_mv13d_allqc_cell_full.R",
  "scripts/build_mv13e_allqc_cell_full_closure.R"
)
if (!all(file.exists(implementation_files))) stop("MV13-D implementation absent")
implementation <- data.frame(
  contract_id = "mv13d_implementation_binding_v1", file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv13d_full_cell_topology_prefreeze_v1",
  execution_head = head, source_caches = 132L, pca_models = 7L,
  adopted_models = 1L, new_models = 6L, cell_views = 636L,
  adopted_views = 1L, new_views = 635L, dimension_records = 1272L,
  representation_id = .mv13_residual_id,
  cell_metric = "euclidean_shared_30PC_v1",
  ph_contract = "complete_VR_H0_H1_field2_threshold_minus1",
  elapsed_cap_seconds = 7200, rss_cap_bytes = 8 * 1024^3,
  private_storage_cap_bytes = 512 * 1024^2,
  public_storage_cap_bytes = 16 * 1024^2,
  workers = 1L, automatic_retries = 0L,
  full_execution_authorized = TRUE, independent_closure_required = TRUE,
  landscapes_authorized = FALSE, comparisons_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  biological_claims_authorized = FALSE, manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv13d_validation_v1",
  check_id = c(
    "MV13A_exact", "MV13C_exact", "sentinel_closed", "locator_exact",
    "sentinel_exact", "seven_groups", "six_new_models", "one_adopted_model",
    "six_hundred_thirty_five_new_views", "one_adopted_view",
    "six_hundred_thirty_six_total_views", "one_thousand_two_hundred_seventy_two_dimensions",
    "internal_five_models", "external_two_models", "selected384_shared30PC",
    "pearson_residual_only", "Ripserr_only", "one_worker_zero_retry",
    "positive_caps", "implementation_bound", "independent_closure_required",
    "landscapes_closed", "comparison_clustering_fusion_closed",
    "labels_outcomes_closed", "claims_closed"
  ),
  passed = c(
    nrow(a_contract) == 1L, nrow(c_decision) == 1L,
    c_decision$sentinel_independently_closed,
    sha(locator) == private_binding$sha256[1L],
    sha(sentinel) == private_binding$sha256[2L], nrow(groups) == 7L,
    sum(groups$new_model_count) == 6L, sum(groups$adopt_closed_model) == 1L,
    sum(groups$new_view_count) == 635L,
    sum(groups$adopt_closed_view_count) == 1L, sum(groups$unit_count) == 636L,
    contract$dimension_records == 1272L,
    sum(groups$dataset_scope == "internal124") == 5L,
    sum(groups$dataset_scope == "external8") == 2L,
    all(groups$pca_components == 30L),
    all(groups$representation_id == .mv13_residual_id),
    contract$ph_contract == "complete_VR_H0_H1_field2_threshold_minus1",
    contract$workers == 1L && contract$automatic_retries == 0L,
    all(unlist(contract[c("elapsed_cap_seconds", "rss_cap_bytes",
      "private_storage_cap_bytes", "public_storage_cap_bytes")]) > 0),
    nrow(implementation) == 5L, contract$independent_closure_required,
    !contract$landscapes_authorized,
    !any(c(contract$comparisons_authorized, contract$clustering_authorized,
           contract$fusion_authorized)),
    !contract$labels_authorized && !contract$outcomes_authorized,
    !contract$biological_claims_authorized &&
      !contract$manuscript_claims_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV13-D validation failed")
dir.create(output, recursive = TRUE)
tables <- list(
  "mv13d-contract.csv" = contract, "mv13d-group-queue.csv" = groups,
  "mv13d-private-bindings.csv" = private_binding,
  "mv13d-prior-bindings.csv" = prior_binding,
  "mv13d-implementation-bindings.csv" = implementation,
  "mv13d-validation.csv" = validation,
  "mv13d-decision.csv" = data.frame(
    contract_id = "mv13d_decision_v1", full_execution_authorized_after_commit = TRUE,
    independent_closure_required = TRUE, landscapes_authorized = FALSE,
    labels_authorized = FALSE, outcomes_authorized = FALSE,
    next_action = "commit prefreeze; run serial 7-group production once; independently refit all 7 groups",
    stringsAsFactors = FALSE)
)
for (name in names(tables)) atomic(tables[[name]], file.path(output, name))
files <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv13d_artifact_manifest_v1", artifact = basename(files),
  bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, sha, character(1L)), stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv13d-artifact-manifest.csv"))
cat("Completed MV13-D full-production prefreeze; checks=", nrow(validation), "\n",
    sep = "")
