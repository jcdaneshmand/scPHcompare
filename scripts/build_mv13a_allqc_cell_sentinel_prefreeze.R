#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "usage: build_mv13a_allqc_cell_sentinel_prefreeze.R <public-output>",
  "<private-output> <execution-head> <mv08p-private> <mv08pr-private>",
  "<mv08ps-private> <mv08o-internal-private> <mv08o-hca-private>"
), call. = FALSE)
public <- args[[1L]]; private <- args[[2L]]
head <- tolower(trimws(args[[3L]])); roots <- args[4:8]
if (dir.exists(public) || dir.exists(private) ||
    !grepl("^[0-9a-f]{40}$", head) || !all(dir.exists(roots))) {
  stop("invalid MV13-A output, head, or private roots", call. = FALSE)
}
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv13_allqc_cell_topology.R")
sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
binding_path <- paste0(
  "docs/audits/mv08r-full-topology-production-prefreeze-v1/",
  "mv08r-source-cache-bindings.csv"
)
panel_path <- paste0(
  "docs/audits/mv08e-reference-reconciliation-evidence/",
  "mv08e-common475-panel.csv"
)
bindings <- read.csv(binding_path, check.names = FALSE,
                     stringsAsFactors = FALSE)
panel <- read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(bindings) != 132L ||
    !identical(as.integer(bindings$source_order), 1:132) ||
    sum(bindings$dataset_scope == "internal124") != 124L ||
    sum(bindings$dataset_scope == "external8") != 8L ||
    nrow(panel) != 475L ||
    !identical(as.integer(panel$panel_order_475), 1:475) ||
    anyDuplicated(panel$feature_id)) stop("MV13-A source metadata drift")

candidates <- unlist(lapply(roots, function(root) list.files(
  file.path(root, "cache"), pattern = "[.]rds$", recursive = TRUE,
  full.names = TRUE
)), use.names = FALSE)
if (length(candidates) != 132L || anyDuplicated(normalizePath(candidates))) {
  stop("MV13-A requires exactly 132 primary cache candidates")
}
candidate_sha <- vapply(candidates, sha, character(1L))
locator <- lapply(seq_len(nrow(bindings)), function(i) {
  hit <- which(candidate_sha == tolower(bindings$cache_sha256[[i]]))
  if (length(hit) != 1L) stop("MV13-A cache lookup is not one-to-one at ", i)
  data.frame(
    source_order = bindings$source_order[[i]],
    dataset_scope = bindings$dataset_scope[[i]],
    unit_id = bindings$unit_id[[i]],
    fit_cells = bindings$fit_cells[[i]],
    private_cache_path = normalizePath(candidates[[hit]], mustWork = TRUE),
    cache_bytes = as.numeric(file.info(candidates[[hit]])$size),
    cache_sha256 = candidate_sha[[hit]], stringsAsFactors = FALSE
  )
})
locator <- do.call(rbind, locator)
if (!identical(as.numeric(locator$cache_bytes),
               as.numeric(bindings$cache_bytes))) {
  stop("MV13-A cache byte binding drift")
}
queue <- mv13_cell_topology_queue_v1(
  bindings$unit_id[bindings$dataset_scope == "internal124"],
  bindings$unit_id[bindings$dataset_scope == "external8"]
)
sentinel_unit <- bindings$unit_id[
  bindings$dataset_scope == "internal124" &
    bindings$fit_cells == max(bindings$fit_cells[bindings$dataset_scope ==
                                                  "internal124"])
]
if (length(sentinel_unit) != 1L) stop("MV13-A sentinel is ambiguous")
sentinel <- data.frame(
  contract_id = "mv13a_sentinel_v1", dataset_scope = "internal124",
  unit_id = sentinel_unit, seed = 20260809L, panel_id = "exact500",
  representation_id = .mv13_residual_id,
  role = "maximum_frozen_QC_cell_count_after_shared_PCA",
  pca_fit_units = 124L, pca_fit_cells = 124L * 384L,
  pca_components = 30L, ph_views = 1L, ph_max_dim = 1L,
  elapsed_cap_seconds = 1800, rss_cap_bytes = 8 * 1024^3,
  workers = 1L, retries = 0L, stringsAsFactors = FALSE
)
closures <- c(
  "docs/audits/mv08q-full-source-production-closure-v1/mv08q-artifact-manifest.csv",
  "docs/audits/mv08w-full-ph-production-closure-v1/mv08w-artifact-manifest.csv",
  "docs/audits/mv08zu-engine-v2-full-landscape-closure-v1/mv08zu-artifact-manifest.csv",
  "docs/audits/mv12c-historical-fusion-closure-v1/mv12c-artifact-manifest.csv"
)
if (!all(file.exists(closures))) stop("MV13-A prerequisite closure absent")
closure_binding <- data.frame(
  contract_id = "mv13a_closure_binding_v1", file = closures,
  bytes = as.numeric(file.info(closures)$size),
  sha256 = vapply(closures, sha, character(1L)), stringsAsFactors = FALSE
)
implementation <- c(
  "R/dual_view_topology.R", "R/mv13_allqc_cell_topology.R",
  "scripts/build_mv13a_allqc_cell_sentinel_prefreeze.R",
  "scripts/run_mv13b_allqc_cell_sentinel.R"
)
if (!all(file.exists(implementation))) stop("MV13-A implementation absent")
implementation_binding <- data.frame(
  contract_id = "mv13a_implementation_binding_v1", file = implementation,
  bytes = as.numeric(file.info(implementation)$size),
  sha256 = vapply(implementation, sha, character(1L)), stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv13a_allqc_cell_sentinel_prefreeze_v1",
  execution_head = head, source_caches = 132L, internal_units = 124L,
  external_units = 8L, pca_models = nrow(queue$pca),
  cell_views = nrow(queue$views), ph_records = nrow(queue$ph),
  internal_cell_views = sum(queue$views$dataset_scope == "internal124"),
  external_cell_views = sum(queue$views$dataset_scope == "external8"),
  representation_id = .mv13_residual_id,
  cell_metric = "euclidean_shared_30PC_v1",
  filtration = "complete_vietoris_rips_H0_H1_field2_threshold_minus1",
  landscape = paste(
    "finite_positive;essential_H0_excluded;all_consecutive_active_levels;",
    "exact_streamed_squared_L2;H0_H1_separate;no_grid;no_level_cap", sep = ""
  ),
  full_execution_authorized = FALSE, sentinel_authorized_after_commit = TRUE,
  landscapes_authorized = FALSE, comparisons_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE, stringsAsFactors = FALSE
)
public_binding <- bindings
public_binding$private_locator_state <- "private_hash_bound_not_published"
validation <- data.frame(
  contract_id = "mv13a_validation_v1",
  check_id = c(
    "one_hundred_thirty_two_sources", "internal124_external8",
    "cache_hashes_exact", "cache_bytes_exact", "private_locator_one_to_one",
    "common475_exact", "seven_PCA_models", "six_hundred_thirty_six_views",
    "one_thousand_two_hundred_seventy_two_PH", "internal_620_views",
    "external_16_views", "selected384", "shared_30PC",
    "pearson_residual_only", "H0_H1_separate", "landscape_definition",
    "maximum_fit_cell_sentinel", "one_worker_zero_retry", "resource_caps",
    "four_closures_bound", "implementation_bound", "full_execution_closed",
    "downstream_closed", "labels_outcomes_closed", "claims_closed"
  ),
  passed = c(
    nrow(bindings) == 132L, table(bindings$dataset_scope)[["internal124"]] == 124L &&
      table(bindings$dataset_scope)[["external8"]] == 8L,
    identical(tolower(locator$cache_sha256), tolower(bindings$cache_sha256)),
    identical(as.numeric(locator$cache_bytes), as.numeric(bindings$cache_bytes)),
    nrow(locator) == 132L && !anyDuplicated(locator$private_cache_path),
    nrow(panel) == 475L && all(panel$common475_axis_sha256 ==
      "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"),
    nrow(queue$pca) == 7L, nrow(queue$views) == 636L,
    nrow(queue$ph) == 1272L,
    sum(queue$views$dataset_scope == "internal124") == 620L,
    sum(queue$views$dataset_scope == "external8") == 16L,
    all(queue$views$selected_cells == 384L),
    all(queue$views$pca_components == 30L),
    all(queue$views$representation_id == .mv13_residual_id),
    identical(sort(unique(queue$ph$homology_dimension)), c("H0", "H1")),
    grepl("all_consecutive_active_levels", contract$landscape, fixed = TRUE),
    sentinel$unit_id == sentinel_unit && sentinel$seed == 20260809L,
    sentinel$workers == 1L && sentinel$retries == 0L,
    sentinel$elapsed_cap_seconds > 0 && sentinel$rss_cap_bytes == 8 * 1024^3,
    nrow(closure_binding) == 4L, nrow(implementation_binding) == 4L,
    !contract$full_execution_authorized,
    !any(c(contract$landscapes_authorized, contract$comparisons_authorized,
           contract$clustering_authorized, contract$fusion_authorized)),
    !contract$labels_authorized && !contract$outcomes_authorized,
    !contract$biological_claims_authorized &&
      !contract$manuscript_claims_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV13-A validation failed")
dir.create(private, recursive = TRUE); dir.create(public, recursive = TRUE)
atomic(locator, file.path(private, "mv13a-private-cache-locators.csv"))
locator_binding <- data.frame(
  contract_id = "mv13a_private_locator_binding_v1", rows = nrow(locator),
  bytes = as.numeric(file.info(file.path(private,
    "mv13a-private-cache-locators.csv"))$size),
  sha256 = sha(file.path(private, "mv13a-private-cache-locators.csv")),
  publication_state = "private_not_tracked", stringsAsFactors = FALSE
)
tables <- list(
  "mv13a-contract.csv" = contract,
  "mv13a-source-bindings.csv" = public_binding,
  "mv13a-private-locator-binding.csv" = locator_binding,
  "mv13a-pca-queue.csv" = queue$pca,
  "mv13a-view-queue.csv" = queue$views,
  "mv13a-ph-queue.csv" = queue$ph,
  "mv13a-sentinel.csv" = sentinel,
  "mv13a-common475-binding.csv" = data.frame(
    contract_id = "mv13a_common475_binding_v1", rows = nrow(panel),
    bytes = as.numeric(file.info(panel_path)$size), sha256 = sha(panel_path),
    axis_sha256 = unique(panel$common475_axis_sha256), stringsAsFactors = FALSE),
  "mv13a-closure-bindings.csv" = closure_binding,
  "mv13a-implementation-bindings.csv" = implementation_binding,
  "mv13a-validation.csv" = validation
)
for (name in names(tables)) atomic(tables[[name]], file.path(public, name))
manifest_files <- list.files(public, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv13a_artifact_manifest_v1",
  artifact = basename(manifest_files),
  bytes = as.numeric(file.info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha, character(1L)), stringsAsFactors = FALSE
)
atomic(manifest, file.path(public, "mv13a-artifact-manifest.csv"))
cat("Completed MV13-A prefreeze; checks=", nrow(validation), "\n", sep = "")
