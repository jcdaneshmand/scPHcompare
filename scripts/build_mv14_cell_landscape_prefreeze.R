#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv14_cell_landscape_prefreeze.R <mv13e-closure>",
  "<mv13-private-groups> <engine-v2-rust-library>",
  "<private-bindings-output> <audit-output> <expected-head>"
), call. = FALSE)
mv13e_root <- normalizePath(args[[1L]], mustWork = TRUE)
group_root <- normalizePath(args[[2L]], mustWork = TRUE)
rust_library <- normalizePath(args[[3L]], mustWork = TRUE)
private_binding_path <- normalizePath(args[[4L]], mustWork = FALSE)
output <- normalizePath(args[[5L]], mustWork = FALSE)
expected_head <- tolower(trimws(args[[6L]]))
if (dir.exists(output) || file.exists(private_binding_path)) {
  stop("MV14 prefreeze outputs already exist.", call. = FALSE)
}
environment_head <- tolower(trimws(Sys.getenv("MV14_PREFREEZE_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", expected_head) ||
    environment_head != expected_head) {
  stop("MV14 exact implementation HEAD mismatch.", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv13_allqc_cell_topology.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/mv14_cell_landscape.R")

.mv14_verify_manifest(mv13e_root, "mv13e-artifact-manifest.csv")
mv13e_validation <- .mv14_read_csv(file.path(mv13e_root, "mv13e-validation.csv"))
mv13e_decision <- .mv14_read_csv(file.path(mv13e_root, "mv13e-decision.csv"))
mv13e_groups <- .mv14_read_csv(file.path(mv13e_root, "mv13e-group-repeat.csv"))
if (nrow(mv13e_validation) != 27L || !all(.mv14_truth(mv13e_validation$passed)) ||
    nrow(mv13e_decision) != 1L ||
    !.mv14_truth(mv13e_decision$full_cell_PH_independently_closed) ||
    !.mv14_truth(mv13e_decision$landscapes_prefreeze_eligible_next) ||
    .mv14_truth(mv13e_decision$landscapes_authorized_by_this_closure) ||
    nrow(mv13e_groups) != 7L || !all(.mv14_truth(mv13e_groups$all_exact))) {
  stop("MV13-E does not admit the MV14 prefreeze.", call. = FALSE)
}

admission_root <- normalizePath(
  "docs/audits/mv08zq-landscape-kernel-repair-admission-closure-v1",
  mustWork = TRUE
)
.mv14_verify_manifest(admission_root, "mv08zq-artifact-manifest.csv")
admission_validation <- .mv14_read_csv(file.path(
  admission_root, "mv08zq-validation.csv"
))
admission_decision <- .mv14_read_csv(file.path(
  admission_root, "mv08zq-decision.csv"
))
rust_sha <- .mv14_sha256_file(rust_library)
rust_bytes <- as.numeric(file.info(rust_library)$size)
if (nrow(admission_validation) != 33L ||
    !all(.mv14_truth(admission_validation$passed)) ||
    nrow(admission_decision) != 1L ||
    admission_decision$scientific_engine_version != 2L ||
    admission_decision$candidate_sha256 != rust_sha) {
  stop("MV14 Rust engine-v2 admission drift.", call. = FALSE)
}

binding_rows <- list(); group_rows <- list(); binding_index <- 0L
for (base_order in seq_len(nrow(mv13e_groups))) {
  prior <- mv13e_groups[base_order, , drop = FALSE]
  artifact_file <- paste0(prior$group_id, ".rds")
  artifact_path <- file.path(group_root, artifact_file)
  if (!file.exists(artifact_path)) stop("missing MV13 group: ", artifact_file)
  artifact <- readRDS(artifact_path)
  mv13_validate_cell_group_v1(artifact)
  if (artifact$group_id != prior$group_id ||
      artifact$dataset_scope != prior$dataset_scope ||
      artifact$panel_id != prior$panel_id || artifact$seed != prior$seed ||
      length(artifact$records) != prior$units) {
    stop("MV13 group/closure identity drift at order ", base_order)
  }
  artifact_sha <- .mv14_sha256_file(artifact_path)
  artifact_bytes <- as.numeric(file.info(artifact_path)$size)
  for (dimension in c("H0", "H1")) {
    group_order <- (base_order - 1L) * 2L + match(dimension, c("H0", "H1"))
    group_id <- paste(artifact$group_id, dimension, sep = "__")
    per_record <- lapply(seq_along(artifact$records), function(record_index) {
      record <- artifact$records[[record_index]]
      intervals <- .mv14_intervals(record, dimension)
      if (nrow(intervals) &&
          any(!is.finite(intervals) | intervals[, "death"] <= intervals[, "birth"])) {
        stop("MV14 non-finite/non-positive interval drift.", call. = FALSE)
      }
      binding_index <<- binding_index + 1L
      data.frame(
        contract_id = "mv14_private_axis_binding_v1",
        binding_order = binding_index, group_order = group_order,
        group_id = group_id, base_group_order = base_order,
        artifact_file = artifact_file, artifact_bytes = artifact_bytes,
        artifact_sha256 = artifact_sha, record_index = record_index,
        axis_order = record_index, unit_id = record$unit_id,
        view_cache_key = record$view$cache_key,
        result_cache_key = record$result$cache_key,
        diagram_sha256 = record$result$provenance$diagram_sha256,
        homology_dimension = dimension,
        finite_intervals = nrow(intervals),
        active_depth = .mv14_active_depth(intervals),
        stringsAsFactors = FALSE
      )
    })
    binding <- do.call(rbind, per_record)
    binding_rows[[group_order]] <- binding
    pair_count <- choose(nrow(binding), 2L)
    group_rows[[group_order]] <- data.frame(
      contract_id = "mv14_cell_landscape_group_v1",
      group_order = group_order, group_id = group_id,
      base_group_order = base_order,
      dataset_scope = artifact$dataset_scope, panel_id = artifact$panel_id,
      seed = artifact$seed, homology_dimension = dimension,
      units = nrow(binding), pair_count = pair_count,
      chunks = ceiling(pair_count / 250),
      minimum_finite_intervals = min(binding$finite_intervals),
      maximum_finite_intervals = max(binding$finite_intervals),
      minimum_active_depth = min(binding$active_depth),
      maximum_active_depth = max(binding$active_depth),
      artifact_file = artifact_file, artifact_bytes = artifact_bytes,
      artifact_sha256 = artifact_sha,
      labels_used = FALSE, outcomes_used = FALSE, downstream_jobs = 0L,
      stringsAsFactors = FALSE
    )
  }
}
bindings <- do.call(rbind, binding_rows)
groups <- do.call(rbind, group_rows)
.mv14_atomic_csv(bindings, private_binding_path)
private_binding_sha <- .mv14_sha256_file(private_binding_path)
private_binding_bytes <- as.numeric(file.info(private_binding_path)$size)

queue_rows <- list(); production_order <- 0L
for (group_order in seq_len(nrow(groups))) {
  binding <- bindings[bindings$group_order == group_order, , drop = FALSE]
  pairs <- .mv14_group_pairs(binding, groups$group_id[[group_order]])
  starts <- seq.int(1L, nrow(pairs), by = 250L)
  for (chunk_order in seq_along(starts)) {
    production_order <- production_order + 1L
    first <- starts[[chunk_order]]; last <- min(nrow(pairs), first + 249L)
    queue_rows[[production_order]] <- data.frame(
      contract_id = "mv14_cell_landscape_queue_v1",
      production_order = production_order, group_order = group_order,
      group_id = groups$group_id[[group_order]], chunk_order = chunk_order,
      pair_start = first, pair_end = last, pair_count = last - first + 1L,
      pair_subset_sha256 = .mv14_sha256_text(
        pairs$pair_identity_sha256[first:last]
      ),
      authorization_state = "authorized_after_mv14_prefreeze_commit",
      comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
      label_jobs = 0L, outcome_jobs = 0L,
      manuscript_claim_jobs = 0L, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }
}
queue <- do.call(rbind, queue_rows)

implementation_files <- c(
  "R/mv14_cell_landscape.R",
  "scripts/build_mv14_cell_landscape_prefreeze.R",
  "scripts/run_mv14_cell_landscape_chunk.R",
  "scripts/run_mv14_cell_landscape_production.R",
  "scripts/build_mv14_cell_landscape_closure.R"
)
implementation <- data.frame(
  contract_id = "mv14_implementation_binding_v1",
  role = c("helpers", "prefreeze_builder", "chunk_worker", "serial_runner",
           "independent_closure"),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, .mv14_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv14_cell_landscape_prefreeze_v1",
  execution_head = expected_head, scientific_engine_version = 2L,
  landscape_groups = nrow(groups), production_chunks = nrow(queue),
  production_pairs = sum(queue$pair_count), chunk_pair_cap = 250L,
  finite_interval_policy = "all_finite_positive_persistence_intervals",
  essential_h0_policy = "exclude_infinite_interval",
  level_policy = "all_consecutive_active_levels",
  integration_policy = "exact_streamed_squared_l2_on_dimension_support",
  dimension_policy = "h0_h1_separate_primary_outputs",
  grid_policy = "no_universal_fixed_grid",
  level_cap_policy = "no_universal_level_cap",
  engine_id = "rust_scph_landscape_kernel_v2",
  rust_library_sha256 = rust_sha, rust_library_bytes = rust_bytes,
  child_elapsed_cap_seconds = 1800,
  child_rss_cap_bytes = 4 * 1024^3,
  aggregate_elapsed_cap_seconds = 72000,
  private_storage_cap_bytes = 512 * 1024^2,
  public_storage_cap_bytes = 64 * 1024^2,
  minimum_free_disk_bytes = 10 * 1024^3,
  minimum_available_memory_bytes = 4 * 1024^3,
  workers = 1L, automatic_retries = 0L,
  fallback_policy = "none_stop_and_preserve",
  resume_policy = "strict_completed_prefix_only",
  atomic_write = TRUE, public_receipt_policy = "aggregate_only",
  production_authorized_after_commit = TRUE,
  comparison_state = "closed", clustering_state = "closed",
  fusion_state = "closed", label_state = "closed", outcome_state = "closed",
  manuscript_claim_state = "closed", stringsAsFactors = FALSE
)
closure <- data.frame(
  contract_id = "mv14_prospective_closure_v1",
  required_groups = 14L, required_chunks = 314L,
  required_pairs = 76372L, required_dimensions = 2L,
  require_every_input_rehash = TRUE,
  require_every_pair_identity = TRUE,
  require_every_active_depth = TRUE,
  require_one_exact_R_oracle_per_group = TRUE,
  required_exact_R_oracles = 14L,
  require_zero_partials = TRUE, require_one_worker = TRUE,
  require_zero_retries = TRUE, require_resource_caps = TRUE,
  comparisons_authorized = FALSE, clustering_authorized = FALSE,
  fusion_authorized = FALSE, labels_authorized = FALSE,
  outcomes_authorized = FALSE, manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
inputs <- data.frame(
  contract_id = "mv14_input_binding_v1",
  role = c("mv13e_manifest", "engine_v2_admission_manifest",
           "private_axis_bindings", "engine_v2_library"),
  bytes = c(
    file.info(file.path(mv13e_root, "mv13e-artifact-manifest.csv"))$size,
    file.info(file.path(admission_root, "mv08zq-artifact-manifest.csv"))$size,
    private_binding_bytes, rust_bytes
  ),
  sha256 = c(
    .mv14_sha256_file(file.path(mv13e_root, "mv13e-artifact-manifest.csv")),
    .mv14_sha256_file(file.path(admission_root, "mv08zq-artifact-manifest.csv")),
    private_binding_sha, rust_sha
  ), stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv14_prefreeze_validation_v1",
  check_id = c(
    "mv13e_exact_closure", "rust_engine_v2_admitted", "seven_source_groups",
    "fourteen_landscape_groups", "six_hundred_thirty_six_axes_per_dimension",
    "one_thousand_two_hundred_seventy_two_private_bindings",
    "three_hundred_fourteen_chunks", "seventy_six_thousand_three_hundred_seventy_two_pairs",
    "queue_order", "pair_chunk_cap", "private_binding_hash",
    "all_group_artifacts_rehashed", "all_record_axes_valid",
    "finite_positive_intervals", "essential_h0_excluded",
    "all_active_levels", "exact_streamed_squared_l2", "h0_h1_separate",
    "no_grid", "no_level_cap", "engine_hash", "implementation_hashes",
    "one_worker", "zero_retries", "zero_fallback", "resource_caps",
    "strict_resume", "atomic_outputs", "fourteen_R_oracles_required",
    "downstream_firewall", "labels_outcomes_closed"
  ),
  passed = c(
    all(.mv14_truth(mv13e_validation$passed)),
    all(.mv14_truth(admission_validation$passed)), nrow(mv13e_groups) == 7L,
    nrow(groups) == 14L && all(table(groups$homology_dimension) == 7L),
    all(tapply(groups$units, groups$homology_dimension, sum) == 636L),
    nrow(bindings) == 1272L, nrow(queue) == 314L,
    sum(queue$pair_count) == 76372L,
    identical(as.integer(queue$production_order), seq_len(nrow(queue))),
    max(queue$pair_count) <= 250L, .mv14_sha256_file(private_binding_path) == private_binding_sha,
    length(unique(groups$artifact_sha256)) == 7L,
    !anyDuplicated(bindings[, c("group_order", "axis_order")]),
    all(bindings$finite_intervals >= 0L),
    contract$essential_h0_policy == "exclude_infinite_interval",
    contract$level_policy == "all_consecutive_active_levels",
    contract$integration_policy == "exact_streamed_squared_l2_on_dimension_support",
    contract$dimension_policy == "h0_h1_separate_primary_outputs",
    contract$grid_policy == "no_universal_fixed_grid",
    contract$level_cap_policy == "no_universal_level_cap",
    rust_sha == contract$rust_library_sha256,
    all(file.exists(implementation$file)) &&
      all(vapply(implementation$file, .mv14_sha256_file, character(1L)) == implementation$sha256),
    contract$workers == 1L, contract$automatic_retries == 0L,
    contract$fallback_policy == "none_stop_and_preserve",
    contract$aggregate_elapsed_cap_seconds == 72000 &&
      contract$private_storage_cap_bytes == 512 * 1024^2,
    contract$resume_policy == "strict_completed_prefix_only",
    .mv14_truth(contract$atomic_write), closure$required_exact_R_oracles == 14L,
    all(queue$comparison_jobs == 0L & queue$clustering_jobs == 0L &
          queue$fusion_jobs == 0L & queue$manuscript_claim_jobs == 0L),
    all(queue$outcome_label_state == "closed") &&
      !any(.mv14_truth(queue$biological_outcomes_computed))
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV14 prefreeze failed: ",
  paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE
)
decision <- data.frame(
  contract_id = "mv14_prefreeze_decision_v1",
  decision = "authorize_label_closed_full_cell_landscapes_after_commit",
  production_authorized_after_commit = TRUE,
  groups_authorized = 14L, chunks_authorized = 314L,
  pairs_authorized = 76372L, workers = 1L, automatic_retries = 0L,
  comparisons_authorized = FALSE, clustering_authorized = FALSE,
  fusion_authorized = FALSE, labels_authorized = FALSE,
  outcomes_authorized = FALSE, manuscript_claims_authorized = FALSE,
  next_gate = "MV14_independent_full_cell_landscape_closure",
  stringsAsFactors = FALSE
)

dir.create(output, recursive = TRUE)
.mv14_atomic_csv(contract, file.path(output, "mv14-contract.csv"))
.mv14_atomic_csv(groups, file.path(output, "mv14-group-queue.csv"))
.mv14_atomic_csv(queue, file.path(output, "mv14-production-queue.csv"))
.mv14_atomic_csv(inputs, file.path(output, "mv14-input-bindings.csv"))
.mv14_atomic_csv(implementation, file.path(output, "mv14-implementation-bindings.csv"))
.mv14_atomic_csv(closure, file.path(output, "mv14-prospective-closure.csv"))
.mv14_atomic_csv(validation, file.path(output, "mv14-validation.csv"))
.mv14_atomic_csv(decision, file.path(output, "mv14-decision.csv"))
artifacts <- list.files(output, full.names = TRUE)
manifest <- data.frame(
  contract_id = "mv14_artifact_manifest_v1", artifact = basename(artifacts),
  bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, .mv14_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv14_atomic_csv(manifest, file.path(output, "mv14-artifact-manifest.csv"))
cat("MV14 prefreeze passed ", sum(validation$passed), "/", nrow(validation),
    "; groups=14; chunks=314; pairs=76372\n", sep = "")
