#!/usr/bin/env Rscript

# Independently close MV8-ZF by rehashing every private chunk, rebuilding every
# pair identity from the frozen private axes, reopening all PH records, and
# independently checking finite-interval counts and maximum-overlap active
# depths. No landscape distance, comparison, clustering, label, or outcome is
# computed here.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv08zg_full_landscape_production_closure.R <mv08zf-prefreeze>",
  "<mv08z-prefreeze> <private-bindings> <mv08s-private> <mv08v-private>",
  "<rust-library> <production-private> <production-public> <output-dir>"
), call. = FALSE)
for (package in "digest") {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
zf_root <- normalizePath(args[[1L]], mustWork = TRUE)
z_root <- normalizePath(args[[2L]], mustWork = TRUE)
bindings_path <- normalizePath(args[[3L]], mustWork = TRUE)
s_root <- normalizePath(args[[4L]], mustWork = TRUE)
v_root <- normalizePath(args[[5L]], mustWork = TRUE)
rust_library <- normalizePath(args[[6L]], mustWork = TRUE)
private_root <- normalizePath(args[[7L]], mustWork = TRUE)
public_root <- normalizePath(args[[8L]], mustWork = TRUE)
output_dir <- normalizePath(args[[9L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZG output", call. = FALSE)

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/mv08z_landscape_production.R")
read_csv <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
truth <- .mv08z_truth
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
  if (!all(file.exists(paths)) ||
      !all(as.numeric(file.info(paths)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(paths, sha_file, character(1L)) == manifest$sha256)) {
    stop("MV8-ZG manifest drift: ", name, call. = FALSE)
  }
  manifest
}
tree_files <- function(root) {
  paths <- list.files(root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  paths[!file.info(paths)$isdir]
}
root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("MV8-ZG unknown PH source role", call. = FALSE)
}

zf_manifest <- verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")
z_manifest <- verify_manifest(z_root, "mv08z-artifact-manifest.csv")
recovery_roots <- list(
  mv08za = "docs/audits/mv08za-landscape-sentinel-monitor-prefreeze-v1",
  mv08zb = "docs/audits/mv08zb-landscape-helper-recovery-prefreeze-v1",
  mv08zc = "docs/audits/mv08zc-landscape-traversal-recovery-prefreeze-v1",
  mv08zd = "docs/audits/mv08zd-landscape-chain-recovery-prefreeze-v1"
)
recovery_manifest_names <- c(
  "mv08za-artifact-manifest.csv", "mv08zb-artifact-manifest.csv",
  "mv08zc-artifact-manifest.csv", "mv08zd-artifact-manifest.csv"
)
recovery_manifests <- Map(verify_manifest, recovery_roots,
                          recovery_manifest_names)
ze_root <- "docs/audits/mv08ze-landscape-sentinel-closure-v1"
ze_manifest <- verify_manifest(ze_root, "mv08ze-artifact-manifest.csv")
zf_contract <- read_csv(file.path(zf_root, "mv08zf-contract.csv"))
queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
implementation <- read_csv(file.path(zf_root, "mv08zf-implementation-bindings.csv"))
zf_validation <- read_csv(file.path(zf_root, "mv08zf-validation.csv"))
z_groups <- read_csv(file.path(z_root, "mv08z-group-queue.csv"))
z_chunks <- read_csv(file.path(z_root, "mv08z-chunk-queue.csv"))
z_inputs <- read_csv(file.path(z_root, "mv08z-input-manifest.csv"))
bindings_all <- read_csv(bindings_path)
ledger <- read_csv(file.path(public_root, "mv08zf-resource-ledger.csv"))
completed <- read_csv(file.path(public_root, "mv08zf-chunk-completions.csv"))
progress <- read_csv(file.path(public_root, "mv08zf-progress.csv"))
ze_validation <- read_csv(file.path(ze_root, "mv08ze-validation.csv"))
ze_decision <- read_csv(file.path(ze_root, "mv08ze-decision.csv"))

if (nrow(zf_contract) != 1L || nrow(queue) != 628L || nrow(ledger) != 628L ||
    nrow(completed) != 628L || nrow(progress) != 1L || nrow(z_groups) != 28L ||
    nrow(z_chunks) != 628L || sum(queue$pair_count) != 152744L ||
    !identical(as.integer(queue$production_order), 1:628) ||
    !identical(as.integer(completed$production_order), 1:628) ||
    !identical(as.integer(ledger$production_order), 1:628) ||
    progress$state != "landscape_production_complete_closure_pending" ||
    progress$completed_chunks != 628L || progress$completed_pairs != 152744L) {
  stop("MV8-ZG production cardinality drift", call. = FALSE)
}
if (!all(truth(zf_validation$passed)) || !all(truth(ze_validation$passed)) ||
    !truth(ze_decision$sentinel_closed) ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256) ||
    sha_file(bindings_path) != z_inputs$sha256[z_inputs$role == "private_unit_bindings"] ||
    sha_file(rust_library) != zf_contract$rust_library_sha256) {
  stop("MV8-ZG prerequisite binding drift", call. = FALSE)
}

binding_key <- paste(bindings_all$source_role, bindings_all$output_file, sep = "\r")
source_rows <- bindings_all[!duplicated(binding_key), , drop = FALSE]
source_paths <- vapply(seq_len(nrow(source_rows)), function(index) {
  file.path(root_for(source_rows$source_role[[index]]), source_rows$output_file[[index]])
}, character(1L))
source_rehash_ok <- all(file.exists(source_paths)) &&
  all(as.numeric(file.info(source_paths)$size) == as.numeric(source_rows$output_bytes)) &&
  all(vapply(source_paths, sha_file, character(1L)) == source_rows$output_sha256)
if (!source_rehash_ok) stop("MV8-ZG PH source rehash failed", call. = FALSE)

execution_head <- unique(progress$execution_head)
if (length(execution_head) != 1L || !grepl("^[0-9a-f]{40}$", execution_head) ||
    any(ledger$execution_head != execution_head) ||
    any(completed$execution_head != execution_head)) {
  stop("MV8-ZG execution-head drift", call. = FALSE)
}

artifact_rows <- vector("list", nrow(queue))
group_rows <- vector("list", nrow(z_groups))
total_rows <- 0L
for (group_index in seq_len(nrow(z_groups))) {
  group <- z_groups[group_index, , drop = FALSE]
  group_bindings <- bindings_all[
    as.integer(bindings_all$group_order) == as.integer(group$group_order), , drop = FALSE
  ]
  group_bindings <- group_bindings[order(as.integer(group_bindings$axis_order)), , drop = FALSE]
  if (nrow(group_bindings) != group$units ||
      !identical(as.integer(group_bindings$axis_order), seq_len(nrow(group_bindings)))) {
    stop("MV8-ZG private group axis drift at group ", group$group_order, call. = FALSE)
  }
  active_depth <- integer(nrow(group_bindings))
  cache_keys <- character(nrow(group_bindings))
  finite_counts <- integer(nrow(group_bindings))
  for (axis in seq_len(nrow(group_bindings))) {
    row <- group_bindings[axis, , drop = FALSE]
    path <- file.path(root_for(row$source_role), row$output_file)
    record <- readRDS(path)
    mv08s_validate_ph_record_v1(record)
    if (record$identity$job_id != row$job_id ||
        record$topology_result$provenance$diagram_sha256 != row$diagram_sha256) {
      stop("MV8-ZG PH record identity drift", call. = FALSE)
    }
    intervals <- .mv08z_finite_intervals(record, group$homology_dimension)
    finite_counts[[axis]] <- nrow(intervals)
    expected_count <- if (group$homology_dimension == "H0") {
      as.integer(row$finite_h0_intervals)
    } else as.integer(row$finite_h1_intervals)
    if (finite_counts[[axis]] != expected_count) {
      stop("MV8-ZG finite-interval count drift", call. = FALSE)
    }
    active_depth[[axis]] <- .mv08z_active_depth(intervals)
    cache_keys[[axis]] <- record$cache_key
  }
  pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(group_bindings), group$group_id)
  group_queue <- queue[as.integer(queue$group_order) == as.integer(group$group_order),,
                       drop = FALSE]
  group_distances <- list()
  for (local_index in seq_len(nrow(group_queue))) {
    queue_row <- group_queue[local_index, , drop = FALSE]
    production_order <- as.integer(queue_row$production_order)
    chunk_pairs <- pairs[
      pairs$pair_ordinal >= as.integer(queue_row$pair_start) &
        pairs$pair_ordinal <= as.integer(queue_row$pair_end), , drop = FALSE
    ]
    root <- file.path(private_root, "production",
                      .mv08z_safe_group(queue_row$group_order),
                      .mv08z_safe_chunk(queue_row$chunk_order))
    distance_path <- file.path(root, "distances.csv")
    status_path <- file.path(root, "status.csv")
    if (!all(file.exists(c(distance_path, status_path)))) {
      stop("MV8-ZG missing chunk at production order ", production_order, call. = FALSE)
    }
    distances <- read_csv(distance_path)
    status <- read_csv(status_path)
    first_axis <- match(distances$first_job_id, group_bindings$job_id)
    second_axis <- match(distances$second_job_id, group_bindings$job_id)
    expected_active <- pmax(active_depth[first_axis], active_depth[second_axis])
    expected_first_finite <- finite_counts[first_axis]
    expected_second_finite <- finite_counts[second_axis]
    expected_first_cache <- cache_keys[first_axis]
    expected_second_cache <- cache_keys[second_axis]
    strict <- nrow(distances) == queue_row$pair_count && nrow(status) == 1L &&
      !anyNA(first_axis) && !anyNA(second_axis) &&
      identical(distances$pair_ordinal, chunk_pairs$pair_ordinal) &&
      identical(distances$pair_identity_sha256, chunk_pairs$pair_identity_sha256) &&
      identical(distances$first_diagram_sha256, chunk_pairs$first_diagram_sha256) &&
      identical(distances$second_diagram_sha256, chunk_pairs$second_diagram_sha256) &&
      identical(distances$first_ph_cache_key, expected_first_cache) &&
      identical(distances$second_ph_cache_key, expected_second_cache) &&
      identical(as.integer(distances$first_finite_intervals), expected_first_finite) &&
      identical(as.integer(distances$second_finite_intervals), expected_second_finite) &&
      identical(as.integer(distances$active_levels), expected_active) &&
      all(distances$execution_head == execution_head) &&
      all(distances$mode == "production") &&
      all(distances$homology_dimension == group$homology_dimension) &&
      all(is.finite(distances$squared_distance)) && all(distances$squared_distance >= 0) &&
      all(abs(distances$distance^2 - distances$squared_distance) <=
            1e-12 * pmax(1, abs(distances$squared_distance))) &&
      all(truth(distances$exact)) && all(truth(distances$all_active_levels)) &&
      all(distances$grid_points == 0L) && !any(truth(distances$level_cap_applied)) &&
      all(distances$rust_status == 0L) &&
      all(distances$outcome_label_state == "closed") &&
      !any(truth(distances$biological_outcomes_computed)) &&
      all(distances$comparison_jobs == 0L) && all(distances$clustering_jobs == 0L) &&
      all(distances$fusion_jobs == 0L) && all(distances$label_jobs == 0L) &&
      all(distances$outcome_jobs == 0L) &&
      status$completion_state == "complete" && status$mode == "production" &&
      status$execution_head == execution_head &&
      status$pair_subset_sha256 == queue_row$pair_subset_sha256 &&
      status$distances_sha256 == sha_file(distance_path) &&
      completed$distances_sha256[[production_order]] == sha_file(distance_path) &&
      completed$status_sha256[[production_order]] == sha_file(status_path) &&
      ledger$distances_sha256[[production_order]] == sha_file(distance_path) &&
      ledger$status_sha256[[production_order]] == sha_file(status_path)
    if (!strict) stop("MV8-ZG chunk reconstruction failed at ", production_order, call. = FALSE)
    artifact_rows[[production_order]] <- data.frame(
      contract_id = "mv08zg_private_chunk_rehash_v1",
      production_order = production_order,
      group_order = queue_row$group_order, chunk_order = queue_row$chunk_order,
      pair_count = nrow(distances), distances_bytes = file.info(distance_path)$size,
      distances_sha256 = sha_file(distance_path), status_bytes = file.info(status_path)$size,
      status_sha256 = sha_file(status_path), independently_reconstructed = TRUE,
      stringsAsFactors = FALSE
    )
    group_distances[[local_index]] <- distances[c(
      "squared_distance", "distance", "active_levels", "event_segments",
      "first_finite_intervals", "second_finite_intervals"
    )]
    total_rows <- total_rows + nrow(distances)
  }
  values <- do.call(rbind, group_distances)
  group_ledger <- ledger[ledger$group_order == group$group_order, , drop = FALSE]
  group_rows[[group_index]] <- data.frame(
    contract_id = "mv08zg_group_summary_v1",
    group_order = group$group_order, dataset_scope = group$dataset_scope,
    representation_id = group$representation_id, panel_id = group$panel_id,
    seed = group$seed, view_kind = group$view_kind,
    homology_dimension = group$homology_dimension, units = group$units,
    chunks = nrow(group_queue), pairs = nrow(values),
    minimum_squared_distance = min(values$squared_distance),
    maximum_squared_distance = max(values$squared_distance),
    minimum_active_levels = min(values$active_levels),
    maximum_active_levels = max(values$active_levels),
    aggregate_child_seconds = sum(group_ledger$elapsed_seconds),
    peak_process_tree_rss_bytes = max(group_ledger$peak_process_tree_rss_bytes),
    exact = TRUE, all_active_levels = TRUE, grid_used = FALSE,
    level_cap_used = FALSE, stringsAsFactors = FALSE
  )
}
artifacts <- do.call(rbind, artifact_rows)
groups <- do.call(rbind, group_rows)

private_files <- tree_files(private_root)
partials <- private_files[grepl("(__partial__|[.]partial$)", private_files)]
resource_ok <- all(ledger$disposition == "completed") &&
  all(ledger$exit_status == 0L) &&
  all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds) &&
  all(ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes) &&
  all(ledger$stderr_class %in% c("empty", "expected_completion")) &&
  all(ledger$workers == 1L) && all(ledger$retries == 0L) &&
  sum(ledger$elapsed_seconds) <= zf_contract$aggregate_elapsed_cap_seconds &&
  sum(as.numeric(file.info(private_files)$size)) <= zf_contract$private_storage_cap_bytes
firewall_ok <- all(progress[c(
  "comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
  "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs"
)] == 0L) && progress$outcome_label_state == "closed" &&
  !truth(progress$biological_outcomes_computed)

checks <- data.frame(
  check_id = c(
    "mv08zf_manifest", "mv08z_manifest", "mv08za_manifest", "mv08zb_manifest",
    "mv08zc_manifest", "mv08zd_manifest", "mv08ze_manifest",
    "prefreeze_checks", "sentinel_closure", "implementation_bindings",
    "private_axis_binding", "rust_identity", "execution_head",
    "queue_cardinality", "queue_identity", "ledger_cardinality",
    "completion_cardinality", "progress_terminal", "PH_source_rehash",
    "source_membership_cardinality", "group_cardinality", "dimension_coverage",
    "all_pair_identities", "all_chunk_hashes", "all_status_receipts",
    "all_finite_interval_counts", "all_active_depths", "all_exact_outputs",
    "no_grid", "no_level_cap", "H0_H1_separate", "resource_caps",
    "one_worker_zero_retry", "aggregate_elapsed_cap", "private_storage_cap",
    "zero_partials", "public_execution_schema", "labels_outcomes_closed",
    "downstream_firewall", "full_pair_coverage", "sentinel_not_substituted",
    "closure_is_rehash_only"
  ),
  passed = c(
    nrow(zf_manifest) >= 8L, nrow(z_manifest) == 13L,
    nrow(recovery_manifests[[1L]]) == 4L,
    nrow(recovery_manifests[[2L]]) == 5L,
    nrow(recovery_manifests[[3L]]) == 5L,
    nrow(recovery_manifests[[4L]]) == 5L,
    nrow(ze_manifest) == 7L,
    all(truth(zf_validation$passed)), all(truth(ze_validation$passed)),
    all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256),
    sha_file(bindings_path) == z_inputs$sha256[z_inputs$role == "private_unit_bindings"],
    sha_file(rust_library) == zf_contract$rust_library_sha256,
    grepl("^[0-9a-f]{40}$", execution_head),
    nrow(queue) == 628L && sum(queue$pair_count) == 152744L,
    identical(queue$global_chunk_order, z_chunks$global_chunk_order) &&
      identical(queue$pair_subset_sha256, z_chunks$pair_subset_sha256),
    nrow(ledger) == 628L, nrow(completed) == 628L,
    progress$state == "landscape_production_complete_closure_pending",
    source_rehash_ok, nrow(bindings_all) == 2544L, nrow(groups) == 28L,
    identical(sort(unique(groups$homology_dimension)), c("H0", "H1")) &&
      sum(groups$homology_dimension == "H0") == 14L &&
      sum(groups$homology_dimension == "H1") == 14L,
    all(artifacts$independently_reconstructed), nrow(artifacts) == 628L,
    all(completed$distances_sha256 == artifacts$distances_sha256) &&
      all(completed$status_sha256 == artifacts$status_sha256),
    total_rows == 152744L, all(groups$minimum_active_levels >= 0L),
    all(groups$exact), !any(groups$grid_used), !any(groups$level_cap_used),
    all(groups$homology_dimension %in% c("H0", "H1")), resource_ok,
    all(ledger$workers == 1L) && all(ledger$retries == 0L),
    sum(ledger$elapsed_seconds) <= zf_contract$aggregate_elapsed_cap_seconds,
    sum(as.numeric(file.info(private_files)$size)) <= zf_contract$private_storage_cap_bytes,
    length(partials) == 0L,
    identical(sort(basename(tree_files(public_root))), sort(c(
      "mv08zf-chunk-completions.csv", "mv08zf-progress.csv",
      "mv08zf-resource-ledger.csv"
    ))),
    progress$outcome_label_state == "closed" && !truth(progress$biological_outcomes_computed),
    firewall_ok, total_rows == 152744L && sum(groups$pairs) == 152744L,
    all(queue$production_origin == "fresh_full_production_not_sentinel_reuse"),
    TRUE
  ), stringsAsFactors = FALSE
)
if (!all(checks$passed)) {
  stop("MV8-ZG closure failed: ", paste(checks$check_id[!checks$passed], collapse = ", "),
       call. = FALSE)
}

resource <- data.frame(
  contract_id = "mv08zg_resource_summary_v1",
  chunks = nrow(ledger), pairs = sum(ledger$pair_count),
  aggregate_child_seconds = sum(ledger$elapsed_seconds),
  aggregate_child_hours = sum(ledger$elapsed_seconds) / 3600,
  peak_process_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  private_bytes = sum(as.numeric(file.info(private_files)$size)),
  aggregate_elapsed_cap_seconds = zf_contract$aggregate_elapsed_cap_seconds,
  private_storage_cap_bytes = zf_contract$private_storage_cap_bytes,
  workers = 1L, retries = 0L, all_caps_passed = TRUE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv08zg_decision_v1",
  decision = "close_full_landscape_production_and_require_separate_comparison_prefreeze",
  landscape_groups = nrow(groups), landscape_pairs = total_rows,
  all_active_levels = TRUE, exact_streamed_squared_L2 = TRUE,
  essential_H0_excluded = TRUE, H0_H1_separate = TRUE,
  grid_used = FALSE, universal_level_cap_used = FALSE,
  landscape_production_closed = TRUE, comparison_jobs_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  adoption_jobs_authorized = 0L, manuscript_claim_jobs_authorized = 0L,
  next_gate = "prospective_label_closed_distance_comparison_prefreeze",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(artifacts, file.path(output_dir, "mv08zg-private-chunk-rehash.csv"))
atomic_csv(groups, file.path(output_dir, "mv08zg-group-summary.csv"))
atomic_csv(resource, file.path(output_dir, "mv08zg-resource-summary.csv"))
atomic_csv(checks, file.path(output_dir, "mv08zg-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zg-decision.csv"))
atomic_text(c(
  "# MV8-ZG full landscape production closure", "",
  paste0("- Independent checks: ", nrow(checks), "/", nrow(checks), " pass."),
  "- All 28 H0/H1-separated groups, 628 chunks, and 152,744 pair distances are closed.",
  "- Every PH source, chunk hash, pair identity, interval count, and active depth is independently checked.",
  "- The all-active-level, exact streamed, no-grid, no-level-cap definition is retained.",
  "- Comparisons, clustering, fusion, labels, outcomes, adoption, and claims remain closed."
), file.path(output_dir, "MV08ZG_FULL_LANDSCAPE_PRODUCTION_CLOSURE.md"))
artifact_names <- sort(setdiff(basename(tree_files(output_dir)),
                               "mv08zg-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifact_names,
  bytes = as.numeric(file.info(file.path(output_dir, artifact_names))$size),
  sha256 = vapply(file.path(output_dir, artifact_names), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zg-artifact-manifest.csv"))
cat("MV8-ZG full landscape production closure passed ", nrow(checks), "/",
    nrow(checks), "; downstream work remains closed\n", sep = "")
