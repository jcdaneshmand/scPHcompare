#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) stop(paste(
  "usage: build_mv17a_source_inventory_prefreeze.R",
  "<gene-bindings> <gene-sentinel-root> <gene-production-root>",
  "<gene-distance-root> <gene-distance-public-root>",
  "<cell-locator> <cell-group-root> <cell-axis-bindings>",
  "<cell-distance-root> <cell-distance-public-root> <output-dir>"
), call. = FALSE)

gene_bindings_path <- normalizePath(args[[1L]], mustWork = TRUE)
gene_roots <- c(
  mv08s_private_v3 = normalizePath(args[[2L]], mustWork = TRUE),
  mv08v_recovery_private_v2 = normalizePath(args[[3L]], mustWork = TRUE)
)
gene_distance_root <- normalizePath(args[[4L]], mustWork = TRUE)
gene_public_root <- normalizePath(args[[5L]], mustWork = TRUE)
cell_locator_path <- normalizePath(args[[6L]], mustWork = TRUE)
cell_group_root <- normalizePath(args[[7L]], mustWork = TRUE)
cell_axis_path <- normalizePath(args[[8L]], mustWork = TRUE)
cell_distance_root <- normalizePath(args[[9L]], mustWork = TRUE)
cell_public_root <- normalizePath(args[[10L]], mustWork = TRUE)
output <- args[[11L]]
if (dir.exists(output)) stop("MV17-A output root already exists", call. = FALSE)
dir.create(output, recursive = TRUE)

source("R/mv08z_landscape_production.R")
source("R/mv17_ph_calibration_localization_h2.R")
readc <- .mv08z_read_csv
atomic <- .mv08z_atomic_csv
sha <- .mv17a_sha256_file
execution_head <- tolower(trimws(Sys.getenv("MV17A_PREFREEZE_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV17A_PREFREEZE_HEAD must be the exact committed implementation head",
       call. = FALSE)
}

verify_manifest <- function(root, manifest_name) {
  manifest_path <- file.path(root, manifest_name)
  manifest <- readc(manifest_path)
  required <- c("artifact", "bytes", "sha256")
  if (!all(required %in% names(manifest))) stop("MV17-A manifest schema drift")
  paths <- file.path(root, manifest$artifact)
  if (any(!file.exists(paths)) ||
      !identical(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes)) ||
      !identical(tolower(vapply(paths, sha, character(1L))),
                 tolower(manifest$sha256))) {
    stop("MV17-A accepted closure manifest drift: ", root, call. = FALSE)
  }
  data.frame(
    source_closure = basename(root), artifacts = nrow(manifest),
    manifest_bytes = as.numeric(file.info(manifest_path)$size),
    manifest_sha256 = sha(manifest_path),
    artifact_set_sha256 = .mv17a_set_sha256(manifest$sha256),
    independently_rehashed = TRUE, stringsAsFactors = FALSE
  )
}

closure_specs <- list(
  c("docs/audits/mv08w-full-ph-production-closure-v1",
    "mv08w-artifact-manifest.csv"),
  c("docs/audits/mv08zu-engine-v2-full-landscape-closure-v1",
    "mv08zu-artifact-manifest.csv"),
  c("docs/audits/mv13e-allqc-cell-full-closure-v1",
    "mv13e-artifact-manifest.csv"),
  c("docs/audits/mv14-cell-landscape-closure-v1",
    "mv14-artifact-manifest.csv")
)
closure_binding <- do.call(rbind, lapply(closure_specs, function(spec) {
  verify_manifest(spec[[1L]], spec[[2L]])
}))
closure_binding$contract_id <- "mv17a_source_closure_binding_v1"
closure_binding <- closure_binding[c(
  "contract_id", "source_closure", "artifacts", "manifest_bytes",
  "manifest_sha256", "artifact_set_sha256", "independently_rehashed"
)]

gene_bindings <- readc(gene_bindings_path)
gene_ph <- mv17a_inventory_gene_ph_v1(gene_bindings, gene_roots)
gene_completion <- readc(file.path(gene_public_root,
                                   "mv08zt-chunk-completions.csv"))
gene_distance <- mv17a_inventory_delimited_chunks_v1(
  gene_distance_root, gene_completion, "gene_landscape_distance"
)
gene_expected <- readc(file.path(
  "docs/audits/mv08zu-engine-v2-full-landscape-closure-v1",
  "mv08zu-private-chunk-rehash.csv"
))
if (!identical(tolower(gene_completion$distances_sha256),
               tolower(gene_expected$distances_sha256)) ||
    gene_distance$artifacts != 628L || gene_distance$records != 152744L) {
  stop("MV17-A gene distance closure binding drift", call. = FALSE)
}

cell_locator <- readc(cell_locator_path)
locator_required <- c("dataset_scope", "private_cache_path", "cache_bytes",
                      "cache_sha256")
if (!all(locator_required %in% names(cell_locator)) || nrow(cell_locator) != 132L ||
    any(!file.exists(cell_locator$private_cache_path))) {
  stop("MV17-A cell source locator drift", call. = FALSE)
}
cell_cache_hash <- vapply(cell_locator$private_cache_path, sha, character(1L))
if (!identical(tolower(cell_cache_hash), tolower(cell_locator$cache_sha256)) ||
    !identical(as.numeric(file.info(cell_locator$private_cache_path)$size),
               as.numeric(cell_locator$cache_bytes))) {
  stop("MV17-A cell source cache rehash drift", call. = FALSE)
}
cell_group_paths <- sort(list.files(
  cell_group_root, pattern = "[.]rds$", full.names = TRUE
), method = "radix")
cell_groups <- mv17a_inventory_rds_groups_v1(cell_group_paths)
cell_axis <- readc(cell_axis_path)
axis_required <- c("artifact_file", "artifact_sha256", "homology_dimension",
                   "unit_id", "finite_intervals", "active_depth")
if (!all(axis_required %in% names(cell_axis)) || nrow(cell_axis) != 1272L ||
    !setequal(cell_axis$homology_dimension, c("H0", "H1")) ||
    any(cell_axis$active_depth < 0L)) {
  stop("MV17-A cell dimension-axis drift", call. = FALSE)
}
group_hash_map <- setNames(cell_groups$artifact_sha256,
                           basename(cell_group_paths))
if (!identical(tolower(unname(group_hash_map[cell_axis$artifact_file])),
               tolower(cell_axis$artifact_sha256))) {
  stop("MV17-A cell PH group binding drift", call. = FALSE)
}
cell_completion <- readc(file.path(cell_public_root,
                                   "mv14-chunk-completions.csv"))
cell_distance <- mv17a_inventory_delimited_chunks_v1(
  cell_distance_root, cell_completion, "cell_landscape_distance"
)
if (cell_distance$artifacts != 314L || cell_distance$records != 76372L) {
  stop("MV17-A cell distance corpus drift", call. = FALSE)
}

inventory <- rbind(
  data.frame(
    contract_id = "mv17a_source_inventory_v1", source_family = "gene_PH",
    artifact_kind = "PH_RDS", artifacts = gene_ph$artifacts,
    records = gene_ph$dimension_records, axes = gene_ph$selected_axis_identities,
    bytes = gene_ph$bytes, artifact_set_sha256 = gene_ph$artifact_set_sha256,
    schema_sha256 = gene_ph$schema_sha256, labels_opened = FALSE,
    outcomes_opened = FALSE, scientific_values_published = FALSE
  ),
  data.frame(
    contract_id = "mv17a_source_inventory_v1",
    source_family = "gene_landscape", artifact_kind = "exact_distance_CSV",
    artifacts = gene_distance$artifacts, records = gene_distance$records,
    axes = 28L, bytes = gene_distance$bytes,
    artifact_set_sha256 = gene_distance$artifact_set_sha256,
    schema_sha256 = gene_distance$schema_sha256, labels_opened = FALSE,
    outcomes_opened = FALSE, scientific_values_published = FALSE
  ),
  data.frame(
    contract_id = "mv17a_source_inventory_v1", source_family = "cell_source",
    artifact_kind = "selected384_source_cache_RDS", artifacts = nrow(cell_locator),
    records = nrow(cell_locator), axes = nrow(cell_locator),
    bytes = sum(as.numeric(cell_locator$cache_bytes)),
    artifact_set_sha256 = .mv17a_set_sha256(cell_cache_hash),
    schema_sha256 = digest::digest(paste(names(cell_locator), collapse = "|"),
                                   algo = "sha256", serialize = FALSE),
    labels_opened = FALSE, outcomes_opened = FALSE,
    scientific_values_published = FALSE
  ),
  data.frame(
    contract_id = "mv17a_source_inventory_v1", source_family = "cell_PH",
    artifact_kind = "PCA_view_PH_group_RDS", artifacts = nrow(cell_groups),
    records = sum(cell_groups$dimension_records), axes = sum(cell_groups$units),
    bytes = sum(cell_groups$artifact_bytes),
    artifact_set_sha256 = .mv17a_set_sha256(cell_groups$artifact_sha256),
    schema_sha256 = .mv17a_set_sha256(cell_groups$schema_sha256),
    labels_opened = FALSE, outcomes_opened = FALSE,
    scientific_values_published = FALSE
  ),
  data.frame(
    contract_id = "mv17a_source_inventory_v1",
    source_family = "cell_landscape", artifact_kind = "exact_distance_CSV",
    artifacts = cell_distance$artifacts, records = cell_distance$records,
    axes = 14L, bytes = cell_distance$bytes,
    artifact_set_sha256 = cell_distance$artifact_set_sha256,
    schema_sha256 = cell_distance$schema_sha256, labels_opened = FALSE,
    outcomes_opened = FALSE, scientific_values_published = FALSE
  )
)

axis_inventory <- data.frame(
  contract_id = "mv17a_axis_inventory_v1",
  source_family = c("gene", "cell"),
  diagram_artifacts = c(gene_ph$artifacts, sum(cell_groups$units)),
  dimension_records = c(gene_ph$dimension_records, nrow(cell_axis)),
  landscape_distance_pairs = c(gene_distance$records, cell_distance$records),
  selected_observations_per_view = 384L,
  feature_axis_min = c(gene_ph$point_count_min, min(cell_groups$feature_axis)),
  feature_axis_max = c(gene_ph$point_count_max, max(cell_groups$feature_axis)),
  PCA_models = c(0L, nrow(cell_groups)),
  PCA_components = c(0L, 30L),
  H0_H1_separate = TRUE,
  essential_H0_excluded = TRUE,
  private_identifiers_published = FALSE,
  stringsAsFactors = FALSE
)

estimands <- mv17a_estimand_registry_v1()
localization <- mv17a_localization_registry_v1()
fixtures <- mv17a_h2_fixture_registry_v1()
gates <- mv17a_stage_gate_registry_v1()
firewall <- mv17a_firewall_v1()
mv17a_validate_contract_v1(estimands, localization, fixtures, gates, firewall)

implementation_files <- c(
  "R/mv17_ph_calibration_localization_h2.R",
  "scripts/build_mv17a_source_inventory_prefreeze.R",
  "tests/testthat/test-mv17-ph-calibration-localization-h2.R",
  "tests/testthat/test-mv17a-source-inventory-prefreeze.R",
  "docs/plans/MV17_PH_CALIBRATION_LOCALIZATION_H2_SPRINT_PLAN.md"
)
implementation <- data.frame(
  contract_id = "mv17a_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv17a_source_inventory_estimand_prefreeze_v1",
  execution_head = execution_head,
  starting_head = "d9a74c012a5da992993c54dfaee4a852215b2ae0",
  gene_PH_artifacts = gene_ph$artifacts,
  gene_dimension_records = gene_ph$dimension_records,
  gene_landscape_pairs = gene_distance$records,
  cell_source_caches = nrow(cell_locator),
  cell_PCA_models = nrow(cell_groups),
  cell_views = sum(cell_groups$units),
  cell_dimension_records = nrow(cell_axis),
  cell_landscape_pairs = cell_distance$records,
  dimensions = "H0;H1_separate",
  essential_H0_policy = "exclude_infinite_interval",
  landscape_policy = "all_consecutive_active_levels_exact_no_grid_no_level_cap",
  values_published = FALSE, null_computation_authorized = FALSE,
  localization_computation_authorized = FALSE, H2_computation_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  biology_authorized = FALSE, manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)

validation <- data.frame(
  contract_id = "mv17a_validation_v1",
  check_id = c(
    "immutable_D250_start", "accepted_closure_manifests", "gene_PH_rehash",
    "gene_PH_axes", "gene_distance_rehash", "cell_source_rehash",
    "cell_PCA_view_PH_rehash", "cell_dimension_axis", "cell_distance_rehash",
    "H0_H1_separate", "essential_H0_policy", "estimands_complete",
    "localization_outputs_frozen", "H2_fixture_classes_frozen",
    "stage_gates_prospective", "implementation_bound", "public_aggregate_only",
    "zero_scientific_execution", "label_outcome_firewall",
    "downstream_firewall"
  ),
  passed = c(
    contract$starting_head == "d9a74c012a5da992993c54dfaee4a852215b2ae0",
    nrow(closure_binding) == 4L && all(closure_binding$independently_rehashed),
    gene_ph$artifacts == 1272L && gene_ph$dimension_records == 2544L,
    gene_ph$selected_axis_identities == 1272L &&
      gene_ph$point_count_min == 475L && gene_ph$point_count_max == 500L,
    gene_distance$artifacts == 628L && gene_distance$records == 152744L,
    nrow(cell_locator) == 132L,
    nrow(cell_groups) == 7L && sum(cell_groups$units) == 636L,
    nrow(cell_axis) == 1272L &&
      sum(cell_axis$homology_dimension == "H0") == 636L &&
      sum(cell_axis$homology_dimension == "H1") == 636L,
    cell_distance$artifacts == 314L && cell_distance$records == 76372L,
    all(axis_inventory$H0_H1_separate),
    all(axis_inventory$essential_H0_excluded),
    nrow(estimands) == 4L, nrow(localization) == 4L,
    nrow(fixtures) == 7L, nrow(gates) == 8L &&
      !any(gates$scientific_execution_authorized),
    nrow(implementation) == 5L && all(file.exists(implementation$file)),
    !any(inventory$scientific_values_published) &&
      !any(grepl("unit_id|sample_id|private_cache_path|output_file",
                 names(inventory))),
    !contract$null_computation_authorized &&
      !contract$localization_computation_authorized &&
      !contract$H2_computation_authorized,
    !contract$labels_authorized && !contract$outcomes_authorized,
    !contract$clustering_authorized && !contract$fusion_authorized &&
      !contract$biology_authorized && !contract$manuscript_claims_authorized
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV17-A source/estimand prefreeze failed")

decision <- data.frame(
  contract_id = "mv17a_decision_v1",
  MV17A_closed_after_commit = TRUE,
  authorized_next_implementation = "MV17-B;MV17-D;MV17-E",
  next_scientific_execution = "none_without_distinct_exact_head_prefreeze",
  real_null_calibration_authorized = FALSE,
  real_localization_authorized = FALSE,
  real_H2_authorized = FALSE,
  labels_outcomes_clustering_fusion_claims = "closed",
  stringsAsFactors = FALSE
)

artifacts <- list(
  "mv17a-contract.csv" = contract,
  "mv17a-source-closure-binding.csv" = closure_binding,
  "mv17a-source-inventory.csv" = inventory,
  "mv17a-axis-inventory.csv" = axis_inventory,
  "mv17a-calibration-estimands.csv" = estimands,
  "mv17a-localization-outputs.csv" = localization,
  "mv17a-h2-fixture-classes.csv" = fixtures,
  "mv17a-stage-gates.csv" = gates,
  "mv17a-firewall.csv" = firewall,
  "mv17a-implementation-binding.csv" = implementation,
  "mv17a-validation.csv" = validation,
  "mv17a-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV17-A source inventory and estimand prefreeze", "",
  "MV17-A independently rehashes the accepted MV8 gene and MV13/MV14 cell",
  "source families and publishes only aggregate identities, counts, hashes, and",
  "schemas. It defines separate cell/gene H0/H1 calibration estimands, H0/H1",
  "localization output classes, and H2 fixture classes.", "",
  "No null, localization, H2, label, outcome, clustering, fusion, biological,",
  "or manuscript-claim computation was run or authorized. After this audit is",
  "committed, only separately prefrozen MV17-B, MV17-D, and MV17-E qualification",
  "implementation is eligible."
), file.path(output, "MV17A_SOURCE_INVENTORY_ESTIMAND_PREFREEZE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv17a-artifact-manifest.csv"),
              method = "radix")
manifest <- data.frame(
  contract_id = "mv17a_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv17a-artifact-manifest.csv"))
message("Built MV17-A source/estimand prefreeze; checks=", nrow(validation))
