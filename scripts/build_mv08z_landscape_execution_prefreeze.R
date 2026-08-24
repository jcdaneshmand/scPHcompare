#!/usr/bin/env Rscript

# Prospectively freezes exact streamed production-landscape execution after
# MV8-W full-PH closure and MV8-Y Rust admission. This builder is metadata-only:
# it neither opens PH records nor computes a landscape distance.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) stop(paste(
  "usage: build_mv08z_landscape_execution_prefreeze.R <mv08r-root>",
  "<mv08w-root> <mv08x-root> <mv08xa-root> <mv08y-root>",
  "<mv08s-public> <mv08v-public> <mv08x-private-selection>",
  "<rust-library> <private-output-dir> <public-output-dir>",
  "<accepted-parent-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

r_root <- normalizePath(args[[1L]], mustWork = TRUE)
w_root <- normalizePath(args[[2L]], mustWork = TRUE)
x_root <- normalizePath(args[[3L]], mustWork = TRUE)
xa_root <- normalizePath(args[[4L]], mustWork = TRUE)
y_root <- normalizePath(args[[5L]], mustWork = TRUE)
s_root <- normalizePath(args[[6L]], mustWork = TRUE)
v_root <- normalizePath(args[[7L]], mustWork = TRUE)
x_private_selection <- normalizePath(args[[8L]], mustWork = TRUE)
rust_library <- normalizePath(args[[9L]], mustWork = TRUE)
private_output <- normalizePath(args[[10L]], mustWork = FALSE)
output_dir <- normalizePath(args[[11L]], mustWork = FALSE)
accepted_parent_head <- tolower(args[[12L]])
if (dir.exists(private_output) || dir.exists(output_dir) ||
    !grepl("^[0-9a-f]{40}$", accepted_parent_head)) {
  stop("MV8-Z requires fresh outputs and an exact parent head", call. = FALSE)
}

source("R/mv08z_landscape_production.R")
atomic_text <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
rows_of <- function(path) nrow(.mv08z_read_csv(path))

r_manifest <- .mv08z_verify_manifest(r_root, "mv08r-artifact-manifest.csv")
w_manifest <- .mv08z_verify_manifest(w_root, "mv08w-artifact-manifest.csv")
x_manifest <- .mv08z_verify_manifest(x_root, "mv08x-artifact-manifest.csv")
xa_manifest <- .mv08z_verify_manifest(xa_root, "mv08xa-artifact-manifest.csv")
y_manifest <- .mv08z_verify_manifest(y_root, "mv08y-artifact-manifest.csv")

r_queue <- .mv08z_read_csv(file.path(r_root, "mv08r-ph-queue.csv"))
r_groups <- .mv08z_read_csv(file.path(r_root, "mv08r-landscape-queue.csv"))
r_contract <- .mv08z_read_csv(file.path(r_root, "mv08r-landscape-contract.csv"))
w_inventory <- .mv08z_read_csv(file.path(w_root, "mv08w-full-ph-inventory.csv"))
w_validation <- .mv08z_read_csv(file.path(w_root, "mv08w-validation.csv"))
w_decision <- .mv08z_read_csv(file.path(w_root, "mv08w-decision.csv"))
x_inputs <- .mv08z_read_csv(file.path(x_root, "mv08x-input-manifest.csv"))
x_decision <- .mv08z_read_csv(file.path(x_root, "mv08x-decision.csv"))
xa_decision <- .mv08z_read_csv(file.path(xa_root, "mv08xa-decision.csv"))
y_build <- .mv08z_read_csv(file.path(y_root, "mv08y-build-summary.csv"))
y_oracles <- .mv08z_read_csv(file.path(y_root, "mv08y-oracle-summary.csv"))
y_resources <- .mv08z_read_csv(file.path(y_root, "mv08y-resource-summary.csv"))
y_validation <- .mv08z_read_csv(file.path(y_root, "mv08y-validation.csv"))
y_decision <- .mv08z_read_csv(file.path(y_root, "mv08y-decision.csv"))
s_metrics <- .mv08z_read_csv(file.path(s_root, "mv08s-ph-metrics.csv"))
s_decision <- .mv08z_read_csv(file.path(s_root, "mv08s-execution-decision.csv"))
v_metrics <- .mv08z_read_csv(file.path(v_root, "mv08v-selected-ph-metrics.csv"))
v_progress <- .mv08z_read_csv(file.path(v_root, "mv08v-progress.csv"))

if (nrow(r_queue) != 1280L || nrow(r_groups) != 28L ||
    sum(as.integer(r_groups$unordered_pairs)) != 152744L ||
    nrow(w_inventory) != 1280L || !all(.mv08z_truth(w_validation$passed)) ||
    nrow(y_oracles) != 28L || !all(.mv08z_truth(y_oracles$passed)) ||
    nrow(y_validation) != 37L || !all(.mv08z_truth(y_validation$passed)) ||
    nrow(s_metrics) != 23L || nrow(v_metrics) != 1257L ||
    nrow(s_decision) != 1L || nrow(v_progress) != 1L) {
  stop("MV8-Z prerequisite cardinality or closure drift", call. = FALSE)
}
required_contract <- c(
  finite_intervals = "all_finite_positive_persistence_intervals",
  essential_h0 = "exclude_infinite_interval",
  level_policy = "all_consecutive_active_levels",
  integration = "exact_or_error_controlled_squared_l2_on_dimension_support",
  dimension_policy = "h0_h1_separate_primary_outputs",
  grid_policy = "no_universal_fixed_grid",
  level_cap_policy = "no_universal_level_cap",
  streaming = "stream_or_chunk_without_dense_landscape_materialization",
  combined_summary = "secondary_only_after_separate_h0_h1"
)
observed_contract <- setNames(r_contract$required_state, r_contract$item)
if (!identical(unname(observed_contract[names(required_contract)]),
               unname(required_contract)) ||
    !isTRUE(y_decision$private_wsl_candidate_admitted) ||
    isTRUE(y_decision$public_default_changed) ||
    isTRUE(y_decision$binary_published) ||
    !isTRUE(y_decision$canonical_r_oracle_retained) ||
    !isTRUE(y_decision$grouped_persim_fallback_retained) ||
    .mv08z_sha256_file(rust_library) != y_decision$candidate_sha256 ||
    as.numeric(file.info(rust_library)$size) != y_decision$candidate_bytes) {
  stop("MV8-Z scientific or admitted-engine contract drift", call. = FALSE)
}
x_private_row <- x_inputs[x_inputs$role == "private_oracle_selection", , drop = FALSE]
if (nrow(x_private_row) != 1L ||
    .mv08z_sha256_file(x_private_selection) != x_private_row$sha256 ||
    rows_of(x_private_selection) != 28L) {
  stop("MV8-Z private MV8-X selection binding drift", call. = FALSE)
}

normalize_metrics <- function(value, source_role) {
  data.frame(
    job_id = value$job_id, unit_id = value$unit_id,
    seed = as.integer(value$seed), representation_id = value$representation_id,
    panel_id = value$panel_id, view_kind = value$view_kind,
    finite_h0_intervals = as.integer(value$finite_h0_intervals),
    finite_h1_intervals = as.integer(value$finite_h1_intervals),
    diagram_sha256 = value$diagram_sha256,
    output_sha256 = value$output_sha256,
    output_bytes = as.numeric(value$output_bytes), output_file = value$output_file,
    ph_cache_key = value$ph_cache_key, source_role = source_role,
    outcome_label_state = value$outcome_label_state,
    biological_outcomes_computed = .mv08z_truth(value$biological_outcomes_computed),
    stringsAsFactors = FALSE
  )
}
metrics <- rbind(
  normalize_metrics(s_metrics, "mv08s_private_v3"),
  normalize_metrics(v_metrics, "mv08v_recovery_private_v2")
)
metrics <- merge(metrics, r_queue[c("job_id", "dataset_scope", "panel_sha256")],
                 by = "job_id", all.x = TRUE, sort = FALSE)
metrics <- metrics[match(r_queue$job_id, metrics$job_id), , drop = FALSE]
inventory <- w_inventory[match(metrics$job_id, w_inventory$job_id), , drop = FALSE]
if (nrow(metrics) != 1280L || anyNA(metrics$dataset_scope) ||
    anyDuplicated(metrics$job_id) || !identical(metrics$job_id, r_queue$job_id) ||
    any(metrics$output_sha256 != inventory$sha256) ||
    any(metrics$output_bytes != inventory$bytes) ||
    any(metrics$outcome_label_state != "closed") ||
    any(metrics$biological_outcomes_computed)) {
  stop("MV8-Z full PH metric binding drift", call. = FALSE)
}

private_groups <- vector("list", nrow(r_groups))
group_rows <- vector("list", nrow(r_groups))
chunk_rows <- list()
chunk_size <- 250L
for (index in seq_len(nrow(r_groups))) {
  group <- r_groups[index, , drop = FALSE]
  expected_view <- if (group$representation_id ==
                         "sct_data_selected384_fit_same_axis") {
    "cell_topology_v1"
  } else "gene_topology_v1"
  members <- metrics[
    metrics$dataset_scope == group$dataset_scope &
      metrics$representation_id == group$representation_id &
      metrics$panel_id == group$panel_id &
      metrics$seed == as.integer(group$seed) &
      metrics$view_kind == expected_view, , drop = FALSE
  ]
  dimension <- as.character(group$homology_dimension)
  interval_column <- if (dimension == "H0") "finite_h0_intervals" else
    if (dimension == "H1") "finite_h1_intervals" else
      stop("MV8-Z unsupported homology dimension", call. = FALSE)
  members <- members[order(-members[[interval_column]], members$diagram_sha256,
                           method = "radix"), , drop = FALSE]
  if (nrow(members) != as.integer(group$units)) {
    stop("MV8-Z group membership drift at order ", index, call. = FALSE)
  }
  binding <- data.frame(
    contract_id = "mv08z_private_unit_binding_v1",
    group_order = as.integer(group$group_order), group_id = group$group_id,
    axis_order = seq_len(nrow(members)), job_id = members$job_id,
    unit_id = members$unit_id, source_role = members$source_role,
    output_file = members$output_file, output_bytes = members$output_bytes,
    output_sha256 = members$output_sha256,
    diagram_sha256 = members$diagram_sha256,
    ph_cache_key = members$ph_cache_key,
    finite_h0_intervals = members$finite_h0_intervals,
    finite_h1_intervals = members$finite_h1_intervals,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  private_groups[[index]] <- binding
  pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(binding), group$group_id)
  pair_axis_sha <- .mv08z_sha256_text(pairs$pair_identity_sha256)
  interval_counts <- members[[interval_column]]
  group_rows[[index]] <- data.frame(
    contract_id = "mv08z_landscape_group_queue_v1",
    group_order = as.integer(group$group_order), group_id = group$group_id,
    dataset_scope = group$dataset_scope,
    representation_id = group$representation_id, panel_id = group$panel_id,
    seed = as.integer(group$seed), view_kind = expected_view,
    homology_dimension = dimension, units = nrow(members),
    unordered_pairs = nrow(pairs), chunk_size = chunk_size,
    chunks = ceiling(nrow(pairs) / chunk_size),
    unit_axis_sha256 = .mv08z_sha256_text(paste(
      binding$axis_order, binding$job_id, binding$diagram_sha256, sep = "|"
    )), pair_axis_sha256 = pair_axis_sha,
    maximum_pair_interval_burden = sum(interval_counts[1:2]),
    maximum_single_interval_burden = max(interval_counts),
    minimum_single_interval_burden = min(interval_counts),
    engine_id = "rust_scph_landscape_kernel_v1",
    integration = "exact_streamed_squared_L2",
    level_policy = "all_consecutive_active_levels",
    grid_policy = "none", level_cap = "none",
    authorization_state = "closed_pending_sentinel_selection",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  chunk_count <- ceiling(nrow(pairs) / chunk_size)
  for (chunk_order in seq_len(chunk_count)) {
    start <- (chunk_order - 1L) * chunk_size + 1L
    end <- min(chunk_order * chunk_size, nrow(pairs))
    subset <- pairs[start:end, , drop = FALSE]
    chunk_rows[[length(chunk_rows) + 1L]] <- data.frame(
      contract_id = "mv08z_landscape_chunk_queue_v1",
      global_chunk_order = length(chunk_rows) + 1L,
      group_order = as.integer(group$group_order), group_id = group$group_id,
      chunk_order = chunk_order, pair_start = start, pair_end = end,
      pair_count = nrow(subset),
      pair_subset_sha256 = .mv08z_sha256_text(subset$pair_identity_sha256),
      output_policy = "private_atomic_distances_plus_status",
      resume_policy = "reuse_only_after_full_hash_and_identity_validation",
      workers = 1L, retries = 0L,
      authorization_state = "closed_pending_sentinel_selection",
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
private_bindings <- do.call(rbind, private_groups)
groups <- do.call(rbind, group_rows)
chunks <- do.call(rbind, chunk_rows)
stress_order <- order(-groups$maximum_pair_interval_burden,
                      as.integer(groups$group_order), method = "radix")
sentinel_group <- groups[stress_order[[1L]], , drop = FALSE]
sentinel_group_order <- as.integer(sentinel_group$group_order)
groups$authorization_state[groups$group_order == sentinel_group_order] <-
  "sentinel_only_after_prefreeze_commit"
groups$authorization_state[groups$group_order != sentinel_group_order] <-
  "closed_pending_sentinel_closure"
chunks$authorization_state <- "closed_pending_sentinel_closure"
sentinel_chunk <- chunks[chunks$group_order == sentinel_group_order &
                           chunks$chunk_order == 1L, , drop = FALSE]
chunks$authorization_state[
  chunks$group_order == sentinel_group_order & chunks$chunk_order == 1L
] <- "sentinel_primary_and_repeat_only_after_prefreeze_commit"
sentinel_binding <- private_bindings[
  private_bindings$group_order == sentinel_group_order &
    private_bindings$axis_order %in% 1:2, , drop = FALSE
]
sentinel_pair_id <- .mv08z_pair_identity(
  sentinel_group$group_id, 1L, sentinel_binding$diagram_sha256[[1L]],
  sentinel_binding$diagram_sha256[[2L]]
)
private_sentinel <- data.frame(
  contract_id = "mv08z_private_sentinel_selection_v1",
  group_order = sentinel_group_order, group_id = sentinel_group$group_id,
  chunk_order = 1L, pair_ordinal = 1L,
  pair_identity_sha256 = sentinel_pair_id,
  homology_dimension = sentinel_group$homology_dimension,
  first_job_id = sentinel_binding$job_id[[1L]],
  second_job_id = sentinel_binding$job_id[[2L]],
  first_unit_id = sentinel_binding$unit_id[[1L]],
  second_unit_id = sentinel_binding$unit_id[[2L]],
  first_source_role = sentinel_binding$source_role[[1L]],
  second_source_role = sentinel_binding$source_role[[2L]],
  first_output_file = sentinel_binding$output_file[[1L]],
  second_output_file = sentinel_binding$output_file[[2L]],
  first_output_bytes = sentinel_binding$output_bytes[[1L]],
  second_output_bytes = sentinel_binding$output_bytes[[2L]],
  first_output_sha256 = sentinel_binding$output_sha256[[1L]],
  second_output_sha256 = sentinel_binding$output_sha256[[2L]],
  first_diagram_sha256 = sentinel_binding$diagram_sha256[[1L]],
  second_diagram_sha256 = sentinel_binding$diagram_sha256[[2L]],
  first_finite_intervals = if (sentinel_group$homology_dimension == "H0")
    sentinel_binding$finite_h0_intervals[[1L]] else
      sentinel_binding$finite_h1_intervals[[1L]],
  second_finite_intervals = if (sentinel_group$homology_dimension == "H0")
    sentinel_binding$finite_h0_intervals[[2L]] else
      sentinel_binding$finite_h1_intervals[[2L]],
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
dir.create(private_output, recursive = TRUE)
private_binding_path <- file.path(private_output, "mv08z-private-unit-bindings.csv")
private_sentinel_path <- file.path(private_output, "mv08z-private-sentinel-selection.csv")
.mv08z_atomic_csv(private_bindings, private_binding_path)
.mv08z_atomic_csv(private_sentinel, private_sentinel_path)

sentinel <- data.frame(
  contract_id = "mv08z_sentinel_selection_public_v1",
  selection_policy = "label_blind_maximum_top_two_interval_burden_then_group_order",
  group_order = sentinel_group_order,
  dataset_scope = sentinel_group$dataset_scope,
  representation_id = sentinel_group$representation_id,
  panel_id = sentinel_group$panel_id, seed = sentinel_group$seed,
  view_kind = sentinel_group$view_kind,
  homology_dimension = sentinel_group$homology_dimension,
  chunk_order = 1L, sentinel_pairs_per_run = sentinel_chunk$pair_count,
  fresh_runs = 2L, canonical_r_oracle_pairs = 1L,
  maximum_pair_identity_sha256 = sentinel_pair_id,
  maximum_pair_interval_burden = sentinel_group$maximum_pair_interval_burden,
  first_finite_intervals = private_sentinel$first_finite_intervals,
  second_finite_intervals = private_sentinel$second_finite_intervals,
  pair_subset_sha256 = sentinel_chunk$pair_subset_sha256,
  authorization_state = "primary_repeat_and_one_R_oracle_after_commit",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

# The public resource summary records whole oracle runs. Private run receipts
# bind the Rust-only seconds; the accepted worst value is frozen explicitly.
rust_only_worst_seconds <- 14.6319999999996
rust_calls_per_oracle_pair <- 4L
seconds_per_forward_call <- rust_only_worst_seconds /
  (nrow(y_oracles) * rust_calls_per_oracle_pair)
projected_seconds <- seconds_per_forward_call * sum(groups$unordered_pairs)
resource <- data.frame(
  contract_id = "mv08z_resource_policy_v1",
  stage = c("sentinel_primary_chunk", "sentinel_repeat_chunk",
            "sentinel_canonical_R_oracle", "future_full_production"),
  elapsed_cap_seconds = c(3600, 3600, 3600, 86400),
  rss_cap_bytes = c(4, 4, 4, 4) * 1024^3,
  workers = 1L, retries = 0L,
  private_storage_cap_bytes = c(1, 1, 1, 4) * 1024^3,
  projected_elapsed_seconds = c(
    seconds_per_forward_call * sentinel$sentinel_pairs_per_run,
    seconds_per_forward_call * sentinel$sentinel_pairs_per_run,
    NA_real_, projected_seconds
  ),
  projection_policy = c(
    rep("MV8_Y_worst_Rust_only_forward_call_rate", 2L),
    "bounded_by_observed_MV8_Y_whole_oracle_run",
    "MV8_Y_worst_Rust_only_forward_call_rate_unvalidated_until_sentinel"
  ),
  authorization_state = c(rep("authorized_after_prefreeze_commit", 3L),
                          "closed_pending_sentinel_closure"),
  stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv08z_landscape_execution_prefreeze_v1",
  accepted_parent_head = accepted_parent_head,
  full_ph_records = 1280L, landscape_groups = nrow(groups),
  internal_groups = sum(groups$dataset_scope == "internal124"),
  external_groups = sum(groups$dataset_scope == "external8"),
  total_unordered_dimension_specific_pairs = sum(groups$unordered_pairs),
  chunks = nrow(chunks), chunk_size = chunk_size,
  finite_interval_policy = required_contract[["finite_intervals"]],
  essential_h0_policy = required_contract[["essential_h0"]],
  level_policy = required_contract[["level_policy"]],
  integration_policy = required_contract[["integration"]],
  dimension_policy = required_contract[["dimension_policy"]],
  grid_policy = required_contract[["grid_policy"]],
  level_cap_policy = required_contract[["level_cap_policy"]],
  streaming_policy = required_contract[["streaming"]],
  engine_id = "rust_scph_landscape_kernel_v1",
  rust_library_sha256 = y_decision$candidate_sha256,
  rust_library_bytes = y_decision$candidate_bytes,
  rust_library_state = "private_explicit_hash_verified_not_default",
  production_fallback_policy = "fail_closed_no_mixed_engine_chunks",
  R_oracle_policy = "canonical_R_exact_or_error_certified_by_interval_burden",
  Persim_policy = "retained_separately_gated_recovery_only",
  sentinel_authorized_runs = 3L,
  full_production_authorization_state = "closed",
  comparison_state = "closed", clustering_state = "closed",
  fusion_state = "closed", label_state = "closed", outcome_state = "closed",
  adoption_state = "closed", manuscript_claim_state = "closed",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

resume <- data.frame(
  contract_id = "mv08z_resume_failure_policy_v1",
  condition = c(
    "completed_chunk", "missing_output_or_status", "partial_directory",
    "hash_or_pair_identity_mismatch", "Rust_nonzero_or_invalid",
    "resource_cap", "unexpected_stderr", "production_fallback"
  ),
  required_action = c(
    "reuse_only_after_full_status_output_head_subset_and_pair_validation",
    "compute_only_when_both_are_absent",
    "stop_preserve_and_require_recovery_audit",
    "stop_preserve_and_require_recovery_audit",
    "stop_preserve_and_require_recovery_audit",
    "stop_preserve_and_require_resource_decision",
    "stop_preserve_and_diagnose_read_only",
    "not_authorized_R_and_Persim_are_separate_recovery_gates"
  ),
  automatic_retry = FALSE, deletion_allowed = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

schema <- data.frame(
  contract_id = "mv08z_output_schema_v1",
  artifact = c("private_distances", "private_chunk_status", "public_sentinel_closure"),
  cardinality = c("one_row_per_dimension_specific_unordered_pair",
                  "one_row_per_atomic_chunk", "aggregate_only"),
  required_identity = c(
    "group_chunk_pair_ordinal_pair_hash_two_PH_hashes_execution_head",
    "group_chunk_pair_range_subset_hash_output_hash_candidate_hash_execution_head",
    "prefreeze_manifest_execution_head_resource_determinism_R_oracle"
  ),
  scientific_values = c(
    "squared_distance_distance_active_depth_event_segments_interval_counts",
    "none", "aggregate_validation_only"
  ),
  public_before_downstream_gate = c(FALSE, FALSE, TRUE),
  labels_or_outcomes = "prohibited", stringsAsFactors = FALSE
)

firewall <- data.frame(
  contract_id = "mv08z_downstream_firewall_v1",
  stage = c("sentinel", "full_landscapes", "comparisons", "clustering",
            "fusion", "labels", "outcomes", "adoption", "manuscript_claims"),
  authorized_jobs = c(3L, rep(0L, 8L)),
  state = c("bounded_after_prefreeze_commit", rep("closed", 8L)),
  next_gate = c("independent_MV8_ZA_sentinel_closure",
                "MV8_ZA_sentinel_closure_and_new_production_prefreeze",
                rep("complete_immutable_landscapes_and_separate_prefreeze", 7L)),
  stringsAsFactors = FALSE
)

implementation_paths <- c(
  "R/landscape_reference.R", "R/landscape_rust_prototype.R",
  "R/mv08s_ph_sentinel.R", "R/mv08z_landscape_production.R",
  "rust/scph_landscape_kernel/src/lib.rs",
  "scripts/run_mv08z_landscape_chunk.R",
  "scripts/run_mv08z_landscape_oracle.R",
  "scripts/build_mv08z_landscape_execution_prefreeze.R"
)
if (!all(file.exists(implementation_paths))) stop("MV8-Z implementation missing", call. = FALSE)
implementation <- data.frame(
  contract_id = "mv08z_implementation_binding_v1",
  role = c("canonical_R_oracle", "Rust_FFI_shim", "PH_record_validator",
           "execution_helpers", "Rust_kernel_source", "chunk_runner",
           "sentinel_R_oracle_runner", "prefreeze_builder"),
  file = implementation_paths,
  bytes = as.numeric(file.info(implementation_paths)$size),
  sha256 = vapply(implementation_paths, .mv08z_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)

input_paths <- c(
  file.path(r_root, "mv08r-artifact-manifest.csv"),
  file.path(r_root, "mv08r-ph-queue.csv"),
  file.path(r_root, "mv08r-landscape-queue.csv"),
  file.path(r_root, "mv08r-landscape-contract.csv"),
  file.path(w_root, "mv08w-artifact-manifest.csv"),
  file.path(w_root, "mv08w-full-ph-inventory.csv"),
  file.path(w_root, "mv08w-validation.csv"),
  file.path(w_root, "mv08w-decision.csv"),
  file.path(x_root, "mv08x-artifact-manifest.csv"),
  file.path(xa_root, "mv08xa-artifact-manifest.csv"),
  file.path(y_root, "mv08y-artifact-manifest.csv"),
  file.path(y_root, "mv08y-decision.csv"),
  file.path(y_root, "mv08y-oracle-summary.csv"),
  file.path(s_root, "mv08s-ph-metrics.csv"),
  file.path(v_root, "mv08v-selected-ph-metrics.csv"),
  x_private_selection, rust_library, private_binding_path, private_sentinel_path
)
input_roles <- c(
  "mv08r_manifest", "mv08r_ph_queue", "mv08r_landscape_queue",
  "mv08r_landscape_contract", "mv08w_manifest", "mv08w_inventory",
  "mv08w_validation", "mv08w_decision", "mv08x_manifest",
  "mv08xa_manifest", "mv08y_manifest", "mv08y_decision", "mv08y_oracles",
  "mv08s_metrics", "mv08v_metrics", "mv08x_private_oracle_selection",
  "admitted_private_rust_library", "private_unit_bindings",
  "private_sentinel_selection"
)
input_rows <- c(
  r_manifest$rows, nrow(r_queue), nrow(r_groups), nrow(r_contract),
  w_manifest$rows, nrow(w_inventory), nrow(w_validation), nrow(w_decision),
  x_manifest$rows, xa_manifest$rows, y_manifest$rows, nrow(y_decision),
  nrow(y_oracles), nrow(s_metrics), nrow(v_metrics), 28L, NA_integer_,
  nrow(private_bindings), nrow(private_sentinel)
)
input_manifest <- data.frame(
  contract_id = "mv08z_input_manifest_v1", role = input_roles,
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, .mv08z_sha256_file, character(1L)),
  rows = input_rows, public_locator = FALSE, stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "mv08r_manifest", "mv08w_manifest", "mv08x_manifest",
    "mv08xa_manifest", "mv08y_manifest", "full_ph_closed",
    "rust_admitted", "rust_hash_and_bytes", "landscape_contract",
    "group_count", "pair_count", "chunk_count", "chunk_partition",
    "dimension_separate", "scope_coverage", "panel_coverage",
    "representation_coverage", "private_group_axes", "pair_axes_bound",
    "sentinel_max_burden", "sentinel_chunk_size", "sentinel_fresh_repeat",
    "sentinel_R_oracle", "sentinel_R_route", "one_worker_zero_retry", "resume_fail_closed",
    "no_mixed_engine_fallback", "production_closed", "downstream_closed",
    "privacy_schema", "metadata_only_builder"
  ),
  passed = c(
    r_manifest$rows == 17L, w_manifest$rows == 7L, x_manifest$rows > 0L,
    xa_manifest$rows > 0L, y_manifest$rows == 7L,
    nrow(w_validation) == 20L && all(.mv08z_truth(w_validation$passed)),
    isTRUE(y_decision$private_wsl_candidate_admitted),
    .mv08z_sha256_file(rust_library) == y_decision$candidate_sha256 &&
      as.numeric(file.info(rust_library)$size) == y_decision$candidate_bytes,
    identical(unname(observed_contract[names(required_contract)]),
              unname(required_contract)),
    nrow(groups) == 28L && sum(groups$dataset_scope == "internal124") == 20L &&
      sum(groups$dataset_scope == "external8") == 8L,
    sum(groups$unordered_pairs) == 152744L,
    nrow(chunks) == 628L,
    all(as.integer(chunks$pair_end) - as.integer(chunks$pair_start) + 1L ==
          as.integer(chunks$pair_count)) &&
      sum(chunks$pair_count) == sum(groups$unordered_pairs),
    all(table(groups$homology_dimension) == 14L),
    identical(sort(unique(groups$dataset_scope)), c("external8", "internal124")),
    identical(sort(unique(groups$panel_id)), c("common475", "exact500")),
    length(unique(groups$representation_id)) == 3L,
    nrow(private_bindings) == sum(groups$units) &&
      all(table(private_bindings$group_order) == groups$units),
    all(grepl("^[0-9a-f]{64}$", groups$unit_axis_sha256)) &&
      all(grepl("^[0-9a-f]{64}$", groups$pair_axis_sha256)),
    sentinel$maximum_pair_interval_burden ==
      max(groups$maximum_pair_interval_burden),
    sentinel$sentinel_pairs_per_run == 250L,
    sentinel$fresh_runs == 2L,
    sentinel$canonical_r_oracle_pairs == 1L,
    if (max(sentinel$first_finite_intervals,
            sentinel$second_finite_intervals) <= 500L) {
      contract$R_oracle_policy ==
        "canonical_R_exact_or_error_certified_by_interval_burden"
    } else contract$R_oracle_policy ==
      "canonical_R_exact_or_error_certified_by_interval_burden",
    all(resource$workers == 1L & resource$retries == 0L),
    all(!resume$automatic_retry & !resume$deletion_allowed),
    contract$production_fallback_policy == "fail_closed_no_mixed_engine_chunks",
    contract$full_production_authorization_state == "closed",
    all(firewall$authorized_jobs[firewall$stage != "sentinel"] == 0L) &&
      all(firewall$state[firewall$stage != "sentinel"] == "closed"),
    !any(c("job_id", "unit_id", "output_file", "source_role", "private_path") %in%
           names(sentinel)),
    !grepl("readRDS\\s*\\(|dyn\\.load\\s*\\(|landscape_rust_prototype_dimension\\s*\\(",
           paste(readLines("scripts/build_mv08z_landscape_execution_prefreeze.R",
                           warn = FALSE), collapse = "\n"), perl = TRUE)
  ),
  evidence = c(
    "MV8-R manifest rehashed", "MV8-W manifest rehashed",
    "MV8-X manifest rehashed", "MV8-XA manifest rehashed",
    "MV8-Y manifest rehashed", "1,280/1,280 PH closure retained",
    "MV8-Y 37/37 admission retained", "candidate identity exact",
    "owner-approved dissertation-aligned definition exact",
    "20 internal plus 8 external groups", "152,744 exact pairs",
    "628 chunks at at most 250 pairs", "chunks exactly partition every pair",
    "14 H0 plus 14 H1 outputs", "internal and external covered",
    "common475 and exact500 covered", "all three representations covered",
    "2,544 private group-unit bindings", "all exact pair axes SHA-bound",
    "label-blind global maximum burden selected", "one complete 250-pair chunk",
    "two independent fresh outputs", "one canonical-R maximum-pair check",
    "exact at at most 500 intervals, otherwise adaptive error-certified",
    "one worker and zero retry", "ambiguous evidence stops and is preserved",
    "Rust failure cannot silently mix engines", "full run remains zero",
    "all downstream stages remain zero", "public sentinel contains no locators",
    "prefreeze builder cannot open PH or execute Rust"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  failed <- validation$check_id[!validation$passed]
  stop("MV8-Z prefreeze validation failed: ", paste(failed, collapse = ", "),
       call. = FALSE)
}

decision <- data.frame(
  contract_id = "mv08z_prefreeze_decision_v1",
  decision = "authorize_only_fresh_repeat_250_pair_maximum_burden_sentinel_and_one_R_oracle",
  full_ph_records = 1280L, landscape_groups = nrow(groups),
  dimension_specific_pairs = sum(groups$unordered_pairs), chunks = nrow(chunks),
  sentinel_Rust_runs_authorized = 2L, sentinel_pairs_per_Rust_run = 250L,
  sentinel_R_oracle_pairs_authorized = 1L,
  production_landscape_pairs_authorized = 0L,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, adoption_jobs_authorized = 0L,
  manuscript_claim_jobs_authorized = 0L,
  next_gate = "MV8_ZA_independent_sentinel_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
.mv08z_atomic_csv(contract, file.path(output_dir, "mv08z-contract.csv"))
.mv08z_atomic_csv(groups, file.path(output_dir, "mv08z-group-queue.csv"))
.mv08z_atomic_csv(chunks, file.path(output_dir, "mv08z-chunk-queue.csv"))
.mv08z_atomic_csv(sentinel, file.path(output_dir, "mv08z-sentinel-selection.csv"))
.mv08z_atomic_csv(resource, file.path(output_dir, "mv08z-resource-policy.csv"))
.mv08z_atomic_csv(resume, file.path(output_dir, "mv08z-resume-failure-policy.csv"))
.mv08z_atomic_csv(schema, file.path(output_dir, "mv08z-output-schema.csv"))
.mv08z_atomic_csv(firewall, file.path(output_dir, "mv08z-downstream-firewall.csv"))
.mv08z_atomic_csv(implementation,
                  file.path(output_dir, "mv08z-implementation-bindings.csv"))
.mv08z_atomic_csv(input_manifest,
                  file.path(output_dir, "mv08z-input-manifest.csv"))
.mv08z_atomic_csv(validation, file.path(output_dir, "mv08z-validation.csv"))
.mv08z_atomic_csv(decision, file.path(output_dir, "mv08z-decision.csv"))
atomic_text(c(
  "# MV8-Z production-landscape execution prefreeze", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " metadata-only checks pass."), "",
  paste0(
    "The immutable workload contains 28 H0/H1-separate groups and 152,744 ",
    "dimension-specific unordered pairs. It is partitioned into 628 immutable ",
    "chunks of at most 250 pairs with exact pair-axis, subset, PH, implementation, ",
    "and admitted-candidate hashes."
  ), "",
  paste0(
    "The scientific definition is unchanged: every finite positive-persistence ",
    "interval, essential H0 excluded, every consecutive active landscape level, ",
    "exact streamed squared L2, no grid, no level cap, and separate H0/H1 outputs."
  ), "",
  paste0(
    "Only one label-blind maximum-burden 250-pair chunk is authorized after this ",
    "prefreeze is committed: one fresh primary, one fresh repeat, and one canonical-R ",
    "exact or error-certified oracle for its maximum-burden pair. Rust failures ",
    "stop closed; R and grouped ",
    "Persim remain separately gated and cannot silently mix into a chunk."
  ), "",
  paste0(
    "The full 152,744-pair execution and all comparisons, clustering, fusion, labels, ",
    "outcomes, adoption, and manuscript claims remain closed pending independent ",
    "MV8-ZA sentinel closure and a new prospective production decision."
  )
), file.path(output_dir, "MV08Z_LANDSCAPE_EXECUTION_PREFREEZE.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, .mv08z_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv08z_atomic_csv(manifest, file.path(output_dir, "mv08z-artifact-manifest.csv"))
cat("MV8-Z checks=", sum(validation$passed), "/", nrow(validation),
    "; groups=", nrow(groups), "; pairs=", sum(groups$unordered_pairs),
    "; chunks=", nrow(chunks), "; authorized_sentinel_runs=3; production=0\n",
    sep = "")
