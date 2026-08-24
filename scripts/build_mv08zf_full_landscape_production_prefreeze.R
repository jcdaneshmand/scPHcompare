#!/usr/bin/env Rscript

# Prospectively freeze complete streamed landscape production after MV8-ZE.
# This builder computes no landscape distance and opens no labels or outcomes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv08zf_full_landscape_production_prefreeze.R",
  "<private-bindings> <rust-library> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
bindings_path <- normalizePath(args[[1L]], mustWork = TRUE)
rust_library <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZF output", call. = FALSE)

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
    stop("MV8-ZF prerequisite manifest drift: ", name, call. = FALSE)
  }
  manifest
}

parent_head <- tolower(trimws(Sys.getenv("MV08ZF_PARENT_HEAD", unset = "")))
if (!nzchar(parent_head)) {
  parent_head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
}
if (!grepl("^[0-9a-f]{40}$", parent_head)) stop("cannot bind MV8-ZF parent HEAD", call. = FALSE)

roots <- list(
  mv08z = "docs/audits/mv08z-landscape-execution-prefreeze-v1",
  mv08za = "docs/audits/mv08za-landscape-sentinel-monitor-prefreeze-v1",
  mv08zb = "docs/audits/mv08zb-landscape-helper-recovery-prefreeze-v1",
  mv08zc = "docs/audits/mv08zc-landscape-traversal-recovery-prefreeze-v1",
  mv08zd = "docs/audits/mv08zd-landscape-chain-recovery-prefreeze-v1",
  mv08ze = "docs/audits/mv08ze-landscape-sentinel-closure-v1"
)
manifest_names <- c(
  "mv08z-artifact-manifest.csv", "mv08za-artifact-manifest.csv",
  "mv08zb-artifact-manifest.csv", "mv08zc-artifact-manifest.csv",
  "mv08zd-artifact-manifest.csv", "mv08ze-artifact-manifest.csv"
)
manifests <- Map(verify_manifest, roots, manifest_names)
z_contract <- read_csv(file.path(roots$mv08z, "mv08z-contract.csv"))
z_groups <- read_csv(file.path(roots$mv08z, "mv08z-group-queue.csv"))
z_chunks <- read_csv(file.path(roots$mv08z, "mv08z-chunk-queue.csv"))
z_inputs <- read_csv(file.path(roots$mv08z, "mv08z-input-manifest.csv"))
ze_validation <- read_csv(file.path(roots$mv08ze, "mv08ze-validation.csv"))
ze_decision <- read_csv(file.path(roots$mv08ze, "mv08ze-decision.csv"))
ze_projection <- read_csv(file.path(roots$mv08ze, "mv08ze-projection.csv"))
zd_bindings <- read_csv(file.path(roots$mv08zd, "mv08zd-implementation-bindings.csv"))
private_bindings <- read_csv(bindings_path)
if (nrow(z_contract) != 1L || nrow(z_groups) != 28L || nrow(z_chunks) != 628L ||
    nrow(ze_decision) != 1L || nrow(ze_projection) != 1L ||
    !all(truth(ze_validation$passed)) || !truth(ze_decision$sentinel_closed) ||
    sha_file(bindings_path) != z_inputs$sha256[z_inputs$role == "private_unit_bindings"] ||
    sha_file(rust_library) != z_contract$rust_library_sha256 ||
    nrow(private_bindings) != 2544L) {
  stop("MV8-ZF prerequisite closure drift", call. = FALSE)
}

queue <- z_chunks
queue$contract_id <- "mv08zf_full_landscape_queue_v1"
queue$production_order <- seq_len(nrow(queue))
queue$production_origin <- "fresh_full_production_not_sentinel_reuse"
queue$authorization_state <- "authorized_after_mv08zf_commit"
queue$workers <- 1L
queue$retries <- 0L
queue$resume_policy <- "strict_completed_prefix_only"
queue$comparison_jobs <- 0L
queue$clustering_jobs <- 0L
queue$fusion_jobs <- 0L
queue$label_jobs <- 0L
queue$outcome_jobs <- 0L
queue$adoption_jobs <- 0L
queue$manuscript_claim_jobs <- 0L
queue$outcome_label_state <- "closed"
queue$biological_outcomes_computed <- FALSE
queue <- queue[c(
  "contract_id", "production_order", "global_chunk_order", "group_order",
  "group_id", "chunk_order", "pair_start", "pair_end", "pair_count",
  "pair_subset_sha256", "output_policy", "production_origin", "workers",
  "retries", "resume_policy", "authorization_state", "comparison_jobs",
  "clustering_jobs", "fusion_jobs", "label_jobs", "outcome_jobs",
  "adoption_jobs", "manuscript_claim_jobs", "outcome_label_state",
  "biological_outcomes_computed"
)]

measured_projection <- as.numeric(ze_projection$measured_max_projection_seconds)
twofold_projection <- as.numeric(ze_projection$twofold_planning_seconds)
projected_private_bytes <- as.numeric(ze_projection$private_sentinel_bytes) /
  (2 * as.numeric(ze_projection$observed_pairs_per_chunk)) * sum(queue$pair_count)
twofold_private_bytes <- 2 * projected_private_bytes
aggregate_cap <- 40 * 3600
private_cap <- 1024^3
contract <- data.frame(
  contract_id = "mv08zf_full_landscape_production_v1",
  accepted_parent_head = parent_head, landscape_groups = nrow(z_groups),
  production_chunks = nrow(queue), production_pairs = sum(queue$pair_count),
  sentinel_pairs_reused_as_production = 0L,
  finite_interval_policy = z_contract$finite_interval_policy,
  essential_h0_policy = z_contract$essential_h0_policy,
  level_policy = z_contract$level_policy,
  integration_policy = z_contract$integration_policy,
  dimension_policy = z_contract$dimension_policy,
  grid_policy = z_contract$grid_policy,
  level_cap_policy = z_contract$level_cap_policy,
  streaming_policy = z_contract$streaming_policy,
  engine_id = z_contract$engine_id,
  rust_library_sha256 = z_contract$rust_library_sha256,
  rust_library_bytes = z_contract$rust_library_bytes,
  rust_library_state = z_contract$rust_library_state,
  child_elapsed_cap_seconds = 3600,
  child_rss_cap_bytes = 4 * 1024^3,
  aggregate_elapsed_cap_seconds = aggregate_cap,
  private_storage_cap_bytes = private_cap,
  minimum_free_disk_bytes = 10 * 1024^3,
  minimum_available_memory_bytes = 8 * 1024^3,
  workers = 1L, automatic_retries = 0L,
  fallback_policy = "none_stop_and_preserve",
  atomic_write = TRUE, resume_policy = "strict_completed_prefix_only",
  public_receipt_policy = "aggregate_chunk_and_resource_only",
  comparison_state = "closed", clustering_state = "closed",
  fusion_state = "closed", label_state = "closed", outcome_state = "closed",
  adoption_state = "closed", manuscript_claim_state = "closed",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
resource <- data.frame(
  contract_id = "mv08zf_resource_policy_v1",
  stage = c("landscape_chunk", "full_production"),
  elapsed_cap_seconds = c(contract$child_elapsed_cap_seconds,
                          contract$aggregate_elapsed_cap_seconds),
  rss_cap_bytes = c(contract$child_rss_cap_bytes, contract$child_rss_cap_bytes),
  workers = 1L, retries = 0L,
  private_storage_cap_bytes = contract$private_storage_cap_bytes,
  measured_projection_seconds = c(
    ze_projection$observed_max_chunk_seconds, measured_projection
  ),
  twofold_planning_seconds = c(
    2 * ze_projection$observed_max_chunk_seconds, twofold_projection
  ),
  authorization_state = c(
    "authorized_after_mv08zf_commit", "authorized_after_launch_headroom_gate"
  ), stringsAsFactors = FALSE
)
storage <- data.frame(
  contract_id = "mv08zf_storage_projection_v1",
  sentinel_private_bytes = ze_projection$private_sentinel_bytes,
  sentinel_Rust_chunks = 2L,
  projected_full_private_bytes = projected_private_bytes,
  twofold_projected_private_bytes = twofold_private_bytes,
  private_storage_cap_bytes = private_cap,
  minimum_free_disk_bytes = contract$minimum_free_disk_bytes,
  projection_policy = "sentinel_total_bytes_per_Rust_pair_linear_conservative",
  stringsAsFactors = FALSE
)
resume <- data.frame(
  contract_id = "mv08zf_resume_failure_policy_v1",
  condition = c(
    "fresh_launch", "completed_strict_prefix", "interrupt_with_unowned_logs",
    "child_timeout", "child_RSS_cap", "ordinary_child_failure",
    "aggregate_elapsed_cap", "private_storage_cap", "hash_or_identity_drift"
  ),
  required_action = c(
    "require_fresh_private_and_public_roots",
    "rehash_every_completed_chunk_then_resume_at_next_order",
    "stop_preserve_and_require_recovery_audit",
    "stop_preserve_and_require_new_decision",
    "stop_preserve_and_require_new_decision",
    "stop_preserve_and_require_recovery_audit",
    "stop_before_next_chunk_and_preserve",
    "stop_before_next_chunk_and_preserve",
    "stop_preserve_and_do_not_resume"
  ),
  automatic_retry = FALSE, deletion_allowed = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
closure <- data.frame(
  contract_id = "mv08zg_prospective_closure_requirement_v1",
  requirement = c(
    "rehash_all_prefreeze_manifests", "rehash_all_unique_PH_sources",
    "rehash_all_628_distance_and_status_files", "reconstruct_all_pair_identities",
    "reopen_all_PH_records", "reconstruct_all_finite_interval_counts",
    "reconstruct_all_maximum_overlap_active_depths", "validate_exact_distance_schema",
    "validate_resources_and_strict_prefix", "validate_privacy_and_firewalls"
  ),
  expected_count = c(7L, length(unique(paste(
    private_bindings$source_role, private_bindings$output_file, sep = "\r"
  ))), 1256L, 152744L, 2544L, 2544L, 152744L, 152744L, 628L, 1L),
  recompute_landscape_distance = FALSE,
  labels_or_outcomes = FALSE, stringsAsFactors = FALSE
)

implementation_paths <- c(
  "R/mv08z_landscape_production.R",
  "scripts/run_mv08z_landscape_chunk.R",
  "scripts/run_mv08zf_full_landscape_production.R",
  "scripts/build_mv08zg_full_landscape_production_closure.R",
  "scripts/build_mv08zf_full_landscape_production_prefreeze.R"
)
implementation_roles <- c(
  "landscape_helper", "chunk_worker", "serial_parent_runner",
  "independent_closure_builder", "prefreeze_builder"
)
if (!all(file.exists(implementation_paths))) stop("MV8-ZF implementation absent", call. = FALSE)
prior_worker <- zd_bindings[
  zd_bindings$file == "scripts/run_mv08z_landscape_chunk.R", , drop = FALSE
]
if (nrow(prior_worker) != 1L) stop("MV8-ZF prior worker binding absent", call. = FALSE)
prior_helper <- read_csv(file.path(
  roots$mv08z, "mv08z-implementation-bindings.csv"
))
prior_helper <- prior_helper[
  prior_helper$file == "R/mv08z_landscape_production.R", , drop = FALSE
]
if (nrow(prior_helper) != 1L ||
    prior_helper$sha256 != sha_file("R/mv08z_landscape_production.R")) {
  stop("MV8-ZF unchanged helper binding absent", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv08zf_implementation_binding_v1",
  role = implementation_roles, file = implementation_paths,
  old_sha256 = c(prior_helper$sha256, prior_worker$sha256,
                 rep(NA_character_, 3L)),
  bytes = as.numeric(file.info(implementation_paths)$size),
  sha256 = vapply(implementation_paths, sha_file, character(1L)),
  scientific_change = FALSE,
  safety_change = implementation_roles == "chunk_worker",
  stringsAsFactors = FALSE
)
input_roles <- c(names(roots), "private_unit_bindings", "private_rust_library")
input_paths <- c(
  vapply(seq_along(roots), function(index) {
    file.path(roots[[index]], manifest_names[[index]])
  }, character(1L)), bindings_path, rust_library
)
input_manifest <- data.frame(
  contract_id = "mv08zf_input_manifest_v1", role = input_roles,
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha_file, character(1L)),
  public_locator = c(
    paste0("public_audit_", names(roots)),
    "private_runtime_binding", "private_hash_verified_candidate"
  ), stringsAsFactors = FALSE
)

worker_text <- paste(readLines("scripts/run_mv08z_landscape_chunk.R", warn = FALSE),
                     collapse = "\n")
runner_text <- paste(readLines("scripts/run_mv08zf_full_landscape_production.R",
                               warn = FALSE), collapse = "\n")
closure_text <- paste(readLines("scripts/build_mv08zg_full_landscape_production_closure.R",
                                warn = FALSE), collapse = "\n")
validation <- data.frame(
  check_id = c(
    "mv08z_manifest", "mv08za_manifest", "mv08zb_manifest", "mv08zc_manifest",
    "mv08zd_manifest", "mv08ze_manifest", "sentinel_closure",
    "private_axis_binding", "rust_identity", "group_cardinality",
    "dimension_coverage", "queue_cardinality", "pair_cardinality",
    "queue_identity", "fresh_production_scope", "all_active_levels",
    "essential_H0_exclusion", "exact_or_error_controlled_integration",
    "H0_H1_separate", "no_grid", "no_level_cap", "streamed_chunks",
    "worker_production_guard", "worker_amendment_transition",
    "runner_exact_head", "runner_one_worker", "runner_zero_retry",
    "runner_strict_prefix", "runner_atomic_receipts", "child_resource_caps",
    "aggregate_elapsed_cap", "storage_projection", "launch_headroom",
    "closure_implementation", "closure_pair_reconstruction",
    "implementation_hashes", "public_schema", "downstream_firewall",
    "labels_outcomes_closed", "zero_execution"
  ),
  passed = c(
    nrow(manifests[[1L]]) == 13L, nrow(manifests[[2L]]) == 4L,
    nrow(manifests[[3L]]) == 5L, nrow(manifests[[4L]]) == 5L,
    nrow(manifests[[5L]]) == 5L, nrow(manifests[[6L]]) == 7L,
    all(truth(ze_validation$passed)) && truth(ze_decision$sentinel_closed),
    sha_file(bindings_path) == z_inputs$sha256[z_inputs$role == "private_unit_bindings"],
    sha_file(rust_library) == z_contract$rust_library_sha256,
    nrow(z_groups) == 28L, sum(z_groups$homology_dimension == "H0") == 14L &&
      sum(z_groups$homology_dimension == "H1") == 14L,
    nrow(queue) == 628L && identical(queue$production_order, 1:628),
    sum(queue$pair_count) == 152744L,
    identical(queue$global_chunk_order, z_chunks$global_chunk_order) &&
      identical(queue$pair_subset_sha256, z_chunks$pair_subset_sha256),
    all(queue$production_origin == "fresh_full_production_not_sentinel_reuse"),
    z_contract$level_policy == "all_consecutive_active_levels",
    z_contract$essential_h0_policy == "exclude_infinite_interval",
    z_contract$integration_policy ==
      "exact_or_error_controlled_squared_l2_on_dimension_support",
    z_contract$dimension_policy == "h0_h1_separate_primary_outputs",
    z_contract$grid_policy == "no_universal_fixed_grid",
    z_contract$level_cap_policy == "no_universal_level_cap",
    z_contract$streaming_policy ==
      "stream_or_chunk_without_dense_landscape_materialization",
    grepl("full production remains closed without MV8-ZF authorization",
          worker_text, fixed = TRUE) &&
      grepl("MV08ZF_PREFREEZE", worker_text, fixed = TRUE),
    implementation$old_sha256[implementation$role == "chunk_worker"] == prior_worker$sha256,
    grepl("MV08ZF_GIT_HEAD", runner_text, fixed = TRUE),
    grepl("workers = 1L", runner_text, fixed = TRUE),
    grepl("retries = 0L", runner_text, fixed = TRUE),
    grepl("completed strict prefix", runner_text, fixed = TRUE),
    grepl("atomic_csv", runner_text, fixed = TRUE),
    contract$child_elapsed_cap_seconds == 3600 &&
      contract$child_rss_cap_bytes == 4 * 1024^3,
    aggregate_cap > twofold_projection && aggregate_cap == 40 * 3600,
    twofold_private_bytes < private_cap && private_cap == 1024^3,
    contract$minimum_free_disk_bytes == 10 * 1024^3 &&
      contract$minimum_available_memory_bytes == 8 * 1024^3,
    grepl("independently checking finite-interval counts", closure_text, fixed = TRUE),
    grepl("pair_identity_sha256", closure_text, fixed = TRUE) &&
      grepl("expected_active", closure_text, fixed = TRUE),
    all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256),
    contract$public_receipt_policy == "aggregate_chunk_and_resource_only",
    all(queue[c("comparison_jobs", "clustering_jobs", "fusion_jobs",
                "label_jobs", "outcome_jobs", "adoption_jobs",
                "manuscript_claim_jobs")] == 0L),
    all(queue$outcome_label_state == "closed") &&
      !any(truth(queue$biological_outcomes_computed)),
    TRUE
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV8-ZF prefreeze failed: ",
       paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE)
}
decision <- data.frame(
  contract_id = "mv08zf_decision_v1",
  decision = "authorize_fresh_serial_full_landscape_production_after_commit_and_headroom_gate",
  full_production_authorized = TRUE,
  production_landscape_pairs_authorized = sum(queue$pair_count),
  production_chunks_authorized = nrow(queue), workers = 1L,
  automatic_retries = 0L, scientific_contract_changed = FALSE,
  production_safety_guard_corrected = TRUE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, adoption_jobs_authorized = 0L,
  manuscript_claim_jobs_authorized = 0L,
  next_gate = "MV8_ZG_independent_full_landscape_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(contract, file.path(output_dir, "mv08zf-contract.csv"))
atomic_csv(queue, file.path(output_dir, "mv08zf-production-queue.csv"))
atomic_csv(resource, file.path(output_dir, "mv08zf-resource-policy.csv"))
atomic_csv(storage, file.path(output_dir, "mv08zf-storage-projection.csv"))
atomic_csv(resume, file.path(output_dir, "mv08zf-resume-failure-policy.csv"))
atomic_csv(closure, file.path(output_dir, "mv08zg-prospective-closure.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08zf-implementation-bindings.csv"))
atomic_csv(input_manifest, file.path(output_dir, "mv08zf-input-manifest.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zf-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zf-decision.csv"))
atomic_text(c(
  "# MV8-ZF full landscape production prefreeze", "",
  paste0("- Validation: ", nrow(validation), "/", nrow(validation), " pass."),
  "- Workload: 28 H0/H1-separated groups, 628 fresh chunks, 152,744 pairs.",
  "- Definition: all finite positive intervals, all consecutive active levels, exact streamed squared L2, no grid or level cap.",
  paste0("- Conservative measured/twofold projections: ",
         round(measured_projection / 3600, 3), "/",
         round(twofold_projection / 3600, 3), " worker-hours; aggregate stop 40 hours."),
  paste0("- Twofold private-output projection: ", round(twofold_private_bytes / 1024^2, 1),
         " MiB inside a 1-GiB cap."),
  "- One worker, zero automatic retries, strict-prefix resume, and manifest-verified production authorization are required.",
  "- Comparisons, clustering, fusion, labels, outcomes, adoption, and claims remain closed."
), file.path(output_dir, "MV08ZF_FULL_LANDSCAPE_PRODUCTION_PREFREEZE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zf-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zf-artifact-manifest.csv"))
cat("MV8-ZF full landscape production prefreeze passed ", nrow(validation), "/",
    nrow(validation), "; zero production executed\n", sep = "")
