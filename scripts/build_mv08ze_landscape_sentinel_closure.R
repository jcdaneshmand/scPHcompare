#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) stop(paste(
  "usage: build_mv08ze_landscape_sentinel_closure.R <mv08z-prefreeze>",
  "<mv08za-prefreeze> <mv08zb-prefreeze> <mv08zc-prefreeze>",
  "<mv08zd-prefreeze> <private-bindings> <private-sentinel>",
  "<mv08s-private> <mv08v-private> <rust-library>",
  "<execution-private> <execution-public> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

roots <- lapply(args[c(1L, 2L, 3L, 4L, 5L)], normalizePath, mustWork = TRUE)
names(roots) <- c("mv08z", "mv08za", "mv08zb", "mv08zc", "mv08zd")
binding_path <- normalizePath(args[[6L]], mustWork = TRUE)
sentinel_path <- normalizePath(args[[7L]], mustWork = TRUE)
s_root <- normalizePath(args[[8L]], mustWork = TRUE)
v_root <- normalizePath(args[[9L]], mustWork = TRUE)
rust_library <- normalizePath(args[[10L]], mustWork = TRUE)
private_root <- normalizePath(args[[11L]], mustWork = TRUE)
public_root <- normalizePath(args[[12L]], mustWork = TRUE)
output_dir <- normalizePath(args[[13L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZE output", call. = FALSE)

source("R/mv08z_landscape_production.R")
truth <- .mv08z_truth
read_csv <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
atomic_csv <- .mv08z_atomic_csv
atomic_text <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
verify_manifest <- function(root, name) {
  manifest <- read_csv(file.path(root, name))
  paths <- file.path(root, manifest$artifact)
  ok <- all(file.exists(paths)) &&
    all(as.numeric(file.info(paths)$size) == as.numeric(manifest$bytes)) &&
    all(vapply(paths, sha_file, character(1L)) == manifest$sha256)
  if (!ok) stop("MV8-ZE manifest drift: ", name, call. = FALSE)
  manifest
}
tree_files <- function(root) {
  values <- list.files(root, recursive = TRUE, full.names = TRUE,
                       all.files = TRUE, no.. = TRUE)
  values[!file.info(values)$isdir]
}
root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("MV8-ZE unknown PH source role", call. = FALSE)
}

manifest_names <- c(
  "mv08z-artifact-manifest.csv", "mv08za-artifact-manifest.csv",
  "mv08zb-artifact-manifest.csv", "mv08zc-artifact-manifest.csv",
  "mv08zd-artifact-manifest.csv"
)
manifests <- Map(verify_manifest, roots, manifest_names)
validations <- Map(function(root, prefix) {
  read_csv(file.path(root, paste0(prefix, "-validation.csv")))
}, roots, c("mv08z", "mv08za", "mv08zb", "mv08zc", "mv08zd"))

z_contract <- read_csv(file.path(roots$mv08z, "mv08z-contract.csv"))
z_decision <- read_csv(file.path(roots$mv08z, "mv08z-decision.csv"))
z_groups <- read_csv(file.path(roots$mv08z, "mv08z-group-queue.csv"))
z_chunks <- read_csv(file.path(roots$mv08z, "mv08z-chunk-queue.csv"))
z_inputs <- read_csv(file.path(roots$mv08z, "mv08z-input-manifest.csv"))
za_decision <- read_csv(file.path(roots$mv08za, "mv08za-decision.csv"))
recovery_decisions <- Map(function(root, prefix) {
  read_csv(file.path(root, paste0(prefix, "-decision.csv")))
}, roots[c("mv08zb", "mv08zc", "mv08zd")], c("mv08zb", "mv08zc", "mv08zd"))
recovery_failures <- Map(function(root, prefix) {
  read_csv(file.path(root, paste0(prefix, "-failure.csv")))
}, roots[c("mv08zb", "mv08zc", "mv08zd")], c("mv08zb", "mv08zc", "mv08zd"))

bindings_all <- read_csv(binding_path)
sentinel_private <- read_csv(sentinel_path)
sentinel_public <- read_csv(file.path(roots$mv08z, "mv08z-sentinel-selection.csv"))
ledger_path <- file.path(public_root, "mv08za-resource-ledger.csv")
progress_path <- file.path(public_root, "mv08za-progress.csv")
ledger <- read_csv(ledger_path)
progress <- read_csv(progress_path)

if (nrow(z_contract) != 1L || nrow(z_decision) != 1L ||
    nrow(za_decision) != 1L || nrow(sentinel_private) != 1L ||
    nrow(sentinel_public) != 1L || nrow(progress) != 1L || nrow(ledger) != 3L) {
  stop("MV8-ZE singleton/cardinality drift", call. = FALSE)
}
execution_head <- tolower(progress$execution_head)
group_order <- as.integer(sentinel_public$group_order)
chunk_order <- as.integer(sentinel_public$chunk_order)
group <- z_groups[as.integer(z_groups$group_order) == group_order, , drop = FALSE]
chunk <- z_chunks[as.integer(z_chunks$group_order) == group_order &
                    as.integer(z_chunks$chunk_order) == chunk_order, , drop = FALSE]
bindings <- bindings_all[as.integer(bindings_all$group_order) == group_order, , drop = FALSE]
if (nrow(group) != 1L || nrow(chunk) != 1L || nrow(bindings) != group$units) {
  stop("MV8-ZE sentinel group binding drift", call. = FALSE)
}

primary_dir <- file.path(private_root, "primary", .mv08z_safe_group(group_order),
                         .mv08z_safe_chunk(chunk_order))
repeat_dir <- file.path(private_root, "repeat", .mv08z_safe_group(group_order),
                        .mv08z_safe_chunk(chunk_order))
primary_path <- file.path(primary_dir, "distances.csv")
repeat_path <- file.path(repeat_dir, "distances.csv")
primary_status_path <- file.path(primary_dir, "status.csv")
repeat_status_path <- file.path(repeat_dir, "status.csv")
oracle_path <- file.path(private_root, "oracle.csv")
required_private <- c(primary_path, repeat_path, primary_status_path,
                      repeat_status_path, oracle_path)
if (!all(file.exists(required_private))) stop("MV8-ZE private execution incomplete", call. = FALSE)
primary <- read_csv(primary_path)
repeat_values <- read_csv(repeat_path)
primary_status <- read_csv(primary_status_path)
repeat_status <- read_csv(repeat_status_path)
oracle <- read_csv(oracle_path)

pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(bindings), group$group_id)
pairs <- pairs[pairs$pair_ordinal >= as.integer(chunk$pair_start) &
                 pairs$pair_ordinal <= as.integer(chunk$pair_end), , drop = FALSE]
pair_identity_ok <- nrow(pairs) == as.integer(chunk$pair_count) &&
  identical(primary$pair_ordinal, pairs$pair_ordinal) &&
  identical(repeat_values$pair_ordinal, pairs$pair_ordinal) &&
  identical(primary$pair_identity_sha256, pairs$pair_identity_sha256) &&
  identical(repeat_values$pair_identity_sha256, pairs$pair_identity_sha256) &&
  identical(primary$first_diagram_sha256, pairs$first_diagram_sha256) &&
  identical(primary$second_diagram_sha256, pairs$second_diagram_sha256) &&
  .mv08z_sha256_text(pairs$pair_identity_sha256) == chunk$pair_subset_sha256

scientific_columns <- setdiff(names(primary), "mode")
scientific_identity <- identical(primary[scientific_columns],
                                 repeat_values[scientific_columns])
strict_rows <- function(value, expected_mode) {
  nrow(value) == 250L && all(value$mode == expected_mode) &&
    all(value$execution_head == execution_head) &&
    all(value$group_order == group_order) && all(value$chunk_order == chunk_order) &&
    all(value$homology_dimension == group$homology_dimension) &&
    all(is.finite(value$squared_distance)) && all(value$squared_distance >= 0) &&
    all(abs(value$distance^2 - value$squared_distance) <=
          1e-12 * pmax(1, abs(value$squared_distance))) &&
    all(value$active_levels >= 1L) && all(value$event_segments >= 0L) &&
    all(truth(value$exact)) && all(truth(value$all_active_levels)) &&
    all(value$grid_points == 0L) && !any(truth(value$level_cap_applied)) &&
    all(value$rust_status == 0L) &&
    all(value$outcome_label_state == "closed") &&
    !any(truth(value$biological_outcomes_computed)) &&
    all(value$comparison_jobs == 0L) && all(value$clustering_jobs == 0L) &&
    all(value$fusion_jobs == 0L) && all(value$label_jobs == 0L) &&
    all(value$outcome_jobs == 0L)
}
status_ok <- function(value, path, mode) {
  nrow(value) == 1L && value$execution_head == execution_head &&
    value$mode == mode && value$completion_state == "complete" &&
    value$pair_count == 250L && value$pair_subset_sha256 == chunk$pair_subset_sha256 &&
    value$distances_sha256 == sha_file(path) &&
    as.numeric(value$distances_bytes) == as.numeric(file.info(path)$size) &&
    value$rust_library_sha256 == z_contract$rust_library_sha256 &&
    value$workers == 1L && value$retries == 0L && !truth(value$fallback_used) &&
    value$outcome_label_state == "closed" &&
    !truth(value$biological_outcomes_computed)
}

stage_paths <- c(primary_path, repeat_path, oracle_path)
ledger_ok <- identical(ledger$stage, c(
  "sentinel_primary_chunk", "sentinel_repeat_chunk", "sentinel_canonical_R_oracle"
)) && all(ledger$execution_head == execution_head) &&
  all(ledger$disposition == "completed") && all(ledger$exit_status == 0L) &&
  all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds) &&
  all(ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes) &&
  all(ledger$output_sha256 == vapply(stage_paths, sha_file, character(1L))) &&
  all(as.numeric(ledger$output_bytes) == as.numeric(file.info(stage_paths)$size)) &&
  all(ledger$stderr_class %in% c("empty", "expected_completion")) &&
  all(ledger$workers == 1L) && all(ledger$retries == 0L) &&
  all(ledger$outcome_label_state == "closed") &&
  !any(truth(ledger$biological_outcomes_computed))

logs <- file.path(private_root, "logs", paste0(ledger$stage, c(".stdout", ".stdout", ".stdout")))
errs <- file.path(private_root, "logs", paste0(ledger$stage, c(".stderr", ".stderr", ".stderr")))
log_sizes_ok <- all(file.exists(c(logs, errs))) &&
  all(as.numeric(file.info(logs)$size) == as.numeric(ledger$stdout_bytes)) &&
  all(as.numeric(file.info(errs)$size) == as.numeric(ledger$stderr_bytes))
stderr_text <- if (all(file.exists(errs))) vapply(errs, function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}, character(1L)) else rep("missing", 3L)
stderr_ok <- log_sizes_ok &&
  grepl("^Completed MV8-Z group_[0-9]+/chunk_[0-9]+; pairs=250$", trimws(stderr_text[[1L]])) &&
  grepl("^Completed MV8-Z group_[0-9]+/chunk_[0-9]+; pairs=250$", trimws(stderr_text[[2L]])) &&
  grepl("^MV8-Z sentinel oracle passed$", trimws(stderr_text[[3L]]))

oracle_pair <- primary[primary$pair_identity_sha256 == oracle$pair_identity_sha256, , drop = FALSE]
oracle_ok <- nrow(oracle) == 1L && oracle$execution_head == execution_head &&
  oracle$reference_route == "r_adaptive_certified" && truth(oracle$passed) &&
  truth(oracle$reference_error_estimate <= 1e-8) &&
  oracle$absolute_error <= oracle$acceptance_threshold &&
  oracle$rust_status == 0L && oracle$active_levels == oracle$expected_active_levels &&
  nrow(oracle_pair) == 1L &&
  abs(oracle_pair$squared_distance - oracle$candidate_squared_distance) <=
    oracle$acceptance_threshold &&
  oracle$outcome_label_state == "closed" &&
  !truth(oracle$biological_outcomes_computed)

source_paths <- vapply(seq_len(nrow(bindings)), function(index) {
  file.path(root_for(bindings$source_role[[index]]), bindings$output_file[[index]])
}, character(1L))
source_rehash_ok <- all(file.exists(source_paths)) &&
  all(as.numeric(file.info(source_paths)$size) == as.numeric(bindings$output_bytes)) &&
  all(vapply(source_paths, sha_file, character(1L)) == bindings$output_sha256)

input_binding_ok <- sha_file(binding_path) ==
  z_inputs$sha256[z_inputs$role == "private_unit_bindings"] &&
  sha_file(sentinel_path) ==
    z_inputs$sha256[z_inputs$role == "private_sentinel_selection"] &&
  sha_file(rust_library) == z_contract$rust_library_sha256 &&
  as.numeric(file.info(rust_library)$size) == as.numeric(z_contract$rust_library_bytes)
recovery_ok <- all(vapply(recovery_decisions, nrow, integer(1L)) == 1L) &&
  all(vapply(recovery_failures, nrow, integer(1L)) == 1L) &&
  all(vapply(recovery_decisions, function(value) !truth(value$scientific_contract_changed), logical(1L))) &&
  all(vapply(recovery_decisions, function(value) value$production_pairs_authorized == 0L, logical(1L))) &&
  all(vapply(recovery_failures, function(value) value$landscape_pair_outputs == 0L, logical(1L))) &&
  all(vapply(recovery_failures, function(value) !truth(value$later_children_started), logical(1L)))

private_files <- tree_files(private_root)
partials <- private_files[grepl("(__partial__|[.]partial$)", private_files)]
public_names <- sort(basename(tree_files(public_root)))
public_execution_ok <- identical(public_names, sort(c(
  "mv08za-progress.csv", "mv08za-resource-ledger.csv"
)))
progress_ok <- progress$state == "sentinel_execution_complete_closure_pending" &&
  progress$completed_child_processes == 3L && progress$Rust_chunks == 2L &&
  progress$pairs_per_Rust_chunk == 250L && progress$canonical_R_oracle_pairs == 1L &&
  progress$workers == 1L && progress$retries == 0L && progress$production_pairs == 0L &&
  progress$comparison_jobs == 0L && progress$clustering_jobs == 0L &&
  progress$fusion_jobs == 0L && progress$label_jobs == 0L &&
  progress$outcome_jobs == 0L && progress$outcome_label_state == "closed" &&
  !truth(progress$biological_outcomes_computed)

max_chunk_seconds <- max(ledger$elapsed_seconds[ledger$stage != "sentinel_canonical_R_oracle"])
seconds_per_pair <- max_chunk_seconds / 250
projected_seconds <- seconds_per_pair * z_contract$total_unordered_dimension_specific_pairs
planning_seconds <- 2 * projected_seconds
projection <- data.frame(
  contract_id = "mv08ze_landscape_projection_v1",
  observed_max_chunk_seconds = max_chunk_seconds,
  observed_pairs_per_chunk = 250L,
  observed_seconds_per_pair = seconds_per_pair,
  full_pair_count = z_contract$total_unordered_dimension_specific_pairs,
  measured_max_projection_seconds = projected_seconds,
  measured_max_projection_hours = projected_seconds / 3600,
  twofold_planning_seconds = planning_seconds,
  twofold_planning_hours = planning_seconds / 3600,
  observed_peak_process_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  private_sentinel_bytes = sum(as.numeric(file.info(private_files)$size)),
  projection_policy = "maximum_burden_chunk_linear_conservative",
  production_authorized = FALSE,
  stringsAsFactors = FALSE
)

checks <- data.frame(
  check_id = c(
    "mv08z_manifest", "mv08za_manifest", "mv08zb_manifest", "mv08zc_manifest",
    "mv08zd_manifest", "all_prefreeze_checks", "recovery_chain",
    "private_input_bindings", "execution_head", "sentinel_group_cardinality",
    "pair_identity_reconstruction", "primary_status", "repeat_status",
    "primary_exact_contract", "repeat_exact_contract", "scientific_repeat_identity",
    "public_resource_ledger", "private_log_receipts", "resource_caps",
    "one_worker_zero_retry", "execution_progress", "oracle_cardinality",
    "adaptive_R_certificate", "oracle_Rust_agreement", "oracle_pair_binding",
    "PH_source_rehash", "zero_private_partials", "public_execution_schema",
    "projection_finite", "full_production_closed", "downstream_firewall",
    "labels_outcomes_closed", "binary_private_explicit", "landscape_definition",
    "sentinel_scope_only"
  ),
  passed = c(
    nrow(manifests[[1L]]) == 13L, nrow(manifests[[2L]]) == 4L,
    nrow(manifests[[3L]]) == 5L, nrow(manifests[[4L]]) == 5L,
    nrow(manifests[[5L]]) == 5L,
    all(vapply(validations, function(value) all(truth(value$passed)), logical(1L))),
    recovery_ok, input_binding_ok,
    grepl("^[0-9a-f]{40}$", execution_head) && all(ledger$execution_head == execution_head),
    nrow(group) == 1L && group$units == 124L && nrow(pairs) == 250L,
    pair_identity_ok,
    status_ok(primary_status, primary_path, "sentinel_primary"),
    status_ok(repeat_status, repeat_path, "sentinel_repeat"),
    strict_rows(primary, "sentinel_primary"),
    strict_rows(repeat_values, "sentinel_repeat"), scientific_identity,
    ledger_ok, stderr_ok,
    all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds) &&
      all(ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes),
    all(ledger$workers == 1L) && all(ledger$retries == 0L), progress_ok,
    nrow(oracle) == 1L,
    oracle$reference_route == "r_adaptive_certified" &&
      truth(oracle$reference_error_estimate <= 1e-8),
    truth(oracle$passed) && oracle$absolute_error <= oracle$acceptance_threshold,
    oracle_ok, source_rehash_ok, length(partials) == 0L, public_execution_ok,
    is.finite(projected_seconds) && projected_seconds > 0 && planning_seconds > projected_seconds,
    z_contract$full_production_authorization_state == "closed" &&
      z_decision$production_landscape_pairs_authorized == 0L &&
      za_decision$production_pairs_authorized == 0L,
    all(progress[c("comparison_jobs", "clustering_jobs", "fusion_jobs",
                   "label_jobs", "outcome_jobs")] == 0L),
    progress$outcome_label_state == "closed" && !truth(progress$biological_outcomes_computed),
    z_contract$rust_library_state == "private_explicit_hash_verified_not_default" &&
      sha_file(rust_library) == z_contract$rust_library_sha256,
    z_contract$finite_interval_policy == "all_finite_positive_persistence_intervals" &&
      z_contract$essential_h0_policy == "exclude_infinite_interval" &&
      z_contract$level_policy == "all_consecutive_active_levels" &&
      z_contract$integration_policy ==
        "exact_or_error_controlled_squared_l2_on_dimension_support" &&
      z_contract$dimension_policy == "h0_h1_separate_primary_outputs" &&
      z_contract$grid_policy == "no_universal_fixed_grid" &&
      z_contract$level_cap_policy == "no_universal_level_cap" &&
      z_contract$streaming_policy ==
        "stream_or_chunk_without_dense_landscape_materialization",
    nrow(primary) == 250L && nrow(repeat_values) == 250L && nrow(oracle) == 1L &&
      progress$production_pairs == 0L
  ),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) {
  stop("MV8-ZE closure failed: ", paste(checks$check_id[!checks$passed], collapse = ", "),
       call. = FALSE)
}

execution_summary <- data.frame(
  contract_id = "mv08ze_execution_summary_v1", stage = ledger$stage,
  disposition = ledger$disposition, elapsed_seconds = ledger$elapsed_seconds,
  peak_process_tree_rss_bytes = ledger$peak_process_tree_rss_bytes,
  elapsed_cap_seconds = ledger$elapsed_cap_seconds,
  rss_cap_bytes = ledger$rss_cap_bytes, output_bytes = ledger$output_bytes,
  stderr_class = ledger$stderr_class, workers = ledger$workers,
  retries = ledger$retries, cap_passed = TRUE,
  stringsAsFactors = FALSE
)
scientific_summary <- data.frame(
  contract_id = "mv08ze_scientific_summary_v1",
  execution_head = execution_head, dataset_scope = group$dataset_scope,
  representation_id = group$representation_id, panel_id = group$panel_id,
  view_kind = group$view_kind, homology_dimension = group$homology_dimension,
  source_units_rehashed = nrow(bindings), pairs_per_chunk = nrow(primary),
  Rust_chunks = 2L, primary_repeat_scientifically_identical = scientific_identity,
  minimum_active_levels = min(primary$active_levels),
  maximum_active_levels = max(primary$active_levels),
  minimum_finite_interval_burden = min(primary$first_finite_intervals +
                                          primary$second_finite_intervals),
  maximum_finite_interval_burden = max(primary$first_finite_intervals +
                                          primary$second_finite_intervals),
  oracle_route = oracle$reference_route,
  oracle_reference_error_estimate = oracle$reference_error_estimate,
  oracle_absolute_engine_error = oracle$absolute_error,
  oracle_acceptance_threshold = oracle$acceptance_threshold,
  oracle_passed = truth(oracle$passed), exact = TRUE,
  all_active_consecutive_levels = TRUE, essential_H0_excluded = TRUE,
  H0_H1_separate = TRUE, uniform_grid_used = FALSE,
  universal_level_cap_used = FALSE, labels_outcomes_opened = FALSE,
  stringsAsFactors = FALSE
)
provenance <- data.frame(
  contract_id = "mv08ze_provenance_summary_v1",
  stage = names(roots), manifest_artifacts = vapply(manifests, nrow, integer(1L)),
  manifest_sha256 = vapply(seq_along(roots), function(index) {
    sha_file(file.path(roots[[index]], manifest_names[[index]]))
  }, character(1L)), stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv08ze_decision_v1",
  decision = "close_bounded_landscape_sentinel_and_require_new_full_production_prefreeze",
  sentinel_closed = TRUE, source_units_rehashed = nrow(bindings),
  Rust_pairs_independently_repeated = nrow(primary), canonical_R_oracle_pairs = 1L,
  full_production_authorized = FALSE, production_landscape_pairs_authorized = 0L,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, adoption_jobs_authorized = 0L,
  manuscript_claim_jobs_authorized = 0L,
  next_gate = "MV8_ZF_full_landscape_production_prefreeze",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(execution_summary, file.path(output_dir, "mv08ze-execution-summary.csv"))
atomic_csv(scientific_summary, file.path(output_dir, "mv08ze-scientific-summary.csv"))
atomic_csv(projection, file.path(output_dir, "mv08ze-projection.csv"))
atomic_csv(provenance, file.path(output_dir, "mv08ze-provenance-summary.csv"))
atomic_csv(checks, file.path(output_dir, "mv08ze-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08ze-decision.csv"))
atomic_text(c(
  "# MV8-ZE landscape sentinel closure", "",
  paste0("- Independent checks: ", nrow(checks), "/", nrow(checks), " pass."),
  "- Two fresh 250-pair maximum-burden Rust chunks are scientifically identical.",
  "- The maximum-burden pair passes the adaptive error-certified canonical-R oracle.",
  paste0("- Conservative measured-max full projection: ",
         format(round(projection$measured_max_projection_hours, 3), nsmall = 3),
         " worker-hours; twofold planning bound: ",
         format(round(projection$twofold_planning_hours, 3), nsmall = 3), " hours."),
  "- The dissertation-aligned all-active-level, exact streamed, H0/H1-separate definition is retained.",
  "- Full landscapes and every comparison/outcome stage remain closed pending MV8-ZF."
), file.path(output_dir, "MV08ZE_LANDSCAPE_SENTINEL_CLOSURE.md"))

artifacts <- sort(setdiff(basename(tree_files(output_dir)),
                          "mv08ze-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08ze-artifact-manifest.csv"))
cat("MV8-ZE landscape sentinel closure passed ", nrow(checks), "/",
    nrow(checks), "; full production remains closed\n", sep = "")
