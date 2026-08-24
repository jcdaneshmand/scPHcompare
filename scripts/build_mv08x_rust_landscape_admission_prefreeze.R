#!/usr/bin/env Rscript

# Prospectively freeze the post-MV8-W Rust landscape-kernel rebuild and
# canonical-R admission corpus. This script selects metadata-only stress pairs;
# it does not open PH artifacts or compute any landscape distance.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv08x_rust_landscape_admission_prefreeze.R <mv08r-root>",
  "<mv08w-root> <mv08s-public> <mv08v-public>",
  "<private-selection.csv> <output-dir>"
), call. = FALSE)

r_root <- normalizePath(args[[1L]], mustWork = TRUE)
w_root <- normalizePath(args[[2L]], mustWork = TRUE)
s_root <- normalizePath(args[[3L]], mustWork = TRUE)
v_root <- normalizePath(args[[4L]], mustWork = TRUE)
private_selection_path <- normalizePath(args[[5L]], mustWork = FALSE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
if (file.exists(private_selection_path) || dir.exists(output_dir)) {
  stop("refusing to overwrite MV8-X prefreeze output", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest required", call. = FALSE)
}

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
truth <- function(value) {
  if (is.logical(value)) return(!is.na(value) & value)
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}
verify_manifest <- function(root, name) {
  manifest_path <- file.path(root, name)
  manifest <- read_csv(manifest_path)
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !all(vapply(paths, sha_file, character(1L)) == manifest$sha256) ||
      !all(as.numeric(file.info(paths)$size) == as.numeric(manifest$bytes))) {
    stop("bound manifest drift: ", name, call. = FALSE)
  }
  list(rows = nrow(manifest), sha256 = sha_file(manifest_path))
}

r_manifest <- verify_manifest(r_root, "mv08r-artifact-manifest.csv")
w_manifest <- verify_manifest(w_root, "mv08w-artifact-manifest.csv")
r_queue <- read_csv(file.path(r_root, "mv08r-ph-queue.csv"))
groups <- read_csv(file.path(r_root, "mv08r-landscape-queue.csv"))
landscape_contract <- read_csv(file.path(r_root, "mv08r-landscape-contract.csv"))
w_inventory <- read_csv(file.path(w_root, "mv08w-full-ph-inventory.csv"))
w_validation <- read_csv(file.path(w_root, "mv08w-validation.csv"))
w_decision <- read_csv(file.path(w_root, "mv08w-decision.csv"))
s_metrics <- read_csv(file.path(s_root, "mv08s-ph-metrics.csv"))
s_decision <- read_csv(file.path(s_root, "mv08s-execution-decision.csv"))
v_metrics <- read_csv(file.path(v_root, "mv08v-selected-ph-metrics.csv"))
v_progress <- read_csv(file.path(v_root, "mv08v-progress.csv"))

if (nrow(r_queue) != 1280L || nrow(groups) != 28L ||
    nrow(w_inventory) != 1280L || !all(truth(w_validation$passed)) ||
    nrow(w_decision) != 1L || w_decision$full_ph_records != 1280L ||
    nrow(s_metrics) != 23L || nrow(v_metrics) != 1257L ||
    nrow(s_decision) != 1L || nrow(v_progress) != 1L) {
  stop("MV8-X prerequisite cardinality drift", call. = FALSE)
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
  combined_summary = "secondary_only_after_separate_h0_h1",
  engine_gate = "rust_rebuild_hash_rebind_and_R_oracle_before_production"
)
observed_contract <- setNames(
  landscape_contract$required_state, landscape_contract$item
)
if (!identical(unname(observed_contract[names(required_contract)]),
               unname(required_contract)) ||
    !all(truth(landscape_contract$owner_approved)) ||
    !all(landscape_contract$outcome_label_state == "closed") ||
    any(truth(landscape_contract$biological_outcomes_computed))) {
  stop("MV8-X landscape contract drift", call. = FALSE)
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
    output_bytes = as.numeric(value$output_bytes),
    output_file = value$output_file, source_role = source_role,
    outcome_label_state = value$outcome_label_state,
    biological_outcomes_computed = truth(value$biological_outcomes_computed),
    stringsAsFactors = FALSE
  )
}
metrics <- rbind(
  normalize_metrics(s_metrics, "mv08s_private_v3"),
  normalize_metrics(v_metrics, "mv08v_recovery_private_v2")
)
metrics <- merge(
  metrics,
  r_queue[c("job_id", "dataset_scope", "panel_sha256")],
  by = "job_id", all.x = TRUE, all.y = FALSE, sort = FALSE
)
metrics <- metrics[match(r_queue$job_id, metrics$job_id), , drop = FALSE]
inventory <- w_inventory[match(metrics$job_id, w_inventory$job_id), , drop = FALSE]
if (nrow(metrics) != 1280L || anyNA(metrics$dataset_scope) ||
    anyDuplicated(metrics$job_id) || !identical(metrics$job_id, r_queue$job_id) ||
    any(metrics$output_sha256 != inventory$sha256) ||
    any(metrics$output_bytes != inventory$bytes) ||
    any(metrics$outcome_label_state != "closed") ||
    any(metrics$biological_outcomes_computed) ||
    !all(truth(inventory$independently_validated)) ||
    !all(truth(inventory$h0_mst_passed))) {
  stop("MV8-X full PH metric binding drift", call. = FALSE)
}

select_group <- function(group) {
  expected_view_kind <- if (
    group$representation_id == "sct_data_selected384_fit_same_axis"
  ) "cell_topology_v1" else "gene_topology_v1"
  candidates <- metrics[
    metrics$dataset_scope == group$dataset_scope &
      metrics$representation_id == group$representation_id &
      metrics$panel_id == group$panel_id &
      metrics$seed == as.integer(group$seed) &
      metrics$view_kind == expected_view_kind,
    , drop = FALSE
  ]
  if (nrow(candidates) != as.integer(group$units)) {
    stop("MV8-X landscape-group membership drift", call. = FALSE)
  }
  dimension <- as.character(group$homology_dimension)
  interval_column <- if (dimension == "H0") {
    "finite_h0_intervals"
  } else if (dimension == "H1") {
    "finite_h1_intervals"
  } else stop("unsupported MV8-X homology dimension", call. = FALSE)
  counts <- candidates[[interval_column]]
  ordering <- order(counts, candidates$diagram_sha256, method = "radix")
  first <- candidates[ordering[[1L]], , drop = FALSE]
  second <- candidates[ordering[[length(ordering)]], , drop = FALSE]
  if (first$job_id == second$job_id) {
    stop("MV8-X stress pair requires two distinct PH records", call. = FALSE)
  }
  pair_payload <- paste(
    "mv08x_oracle_pair_v1", group$group_order, group$dataset_scope,
    group$representation_id, group$panel_id, group$seed, dimension,
    first$diagram_sha256, second$diagram_sha256, sep = "|"
  )
  pair_sha <- digest::digest(pair_payload, algo = "sha256", serialize = FALSE)
  route <- if (max(counts[ordering[c(1L, length(ordering))]]) <= 500L) {
    "r_exact_breakpoint"
  } else "r_adaptive_certified"
  data.frame(
    contract_id = "mv08x_oracle_selection_private_v1",
    oracle_order = as.integer(group$group_order),
    oracle_id = sprintf("mv08x_oracle_%02d_%s", as.integer(group$group_order),
                        substr(pair_sha, 1L, 12L)),
    pair_identity_sha256 = pair_sha,
    group_id = group$group_id, dataset_scope = group$dataset_scope,
    representation_id = group$representation_id, panel_id = group$panel_id,
    panel_sha256 = first$panel_sha256, view_kind = expected_view_kind,
    seed = as.integer(group$seed),
    homology_dimension = dimension, group_units = as.integer(group$units),
    group_unordered_pairs = as.integer(group$unordered_pairs),
    reference_route = route, exact_max_intervals = 500L,
    absolute_tolerance = 1e-8, relative_tolerance = 1e-8,
    subdivisions = 200L,
    first_finite_intervals = first[[interval_column]],
    second_finite_intervals = second[[interval_column]],
    first_job_id = first$job_id, second_job_id = second$job_id,
    first_unit_id = first$unit_id, second_unit_id = second$unit_id,
    first_source_role = first$source_role,
    second_source_role = second$source_role,
    first_output_file = first$output_file,
    second_output_file = second$output_file,
    first_output_bytes = first$output_bytes,
    second_output_bytes = second$output_bytes,
    first_output_sha256 = first$output_sha256,
    second_output_sha256 = second$output_sha256,
    first_diagram_sha256 = first$diagram_sha256,
    second_diagram_sha256 = second$diagram_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
private_selection <- do.call(rbind, lapply(
  seq_len(nrow(groups)), function(index) select_group(groups[index, , drop = FALSE])
))
private_selection <- private_selection[
  order(private_selection$oracle_order), , drop = FALSE
]
atomic_csv(private_selection, private_selection_path)

public_selection <- private_selection[c(
  "contract_id", "oracle_order", "oracle_id", "pair_identity_sha256",
  "dataset_scope", "representation_id", "panel_id", "view_kind", "seed",
  "homology_dimension", "group_units", "group_unordered_pairs",
  "reference_route", "exact_max_intervals", "absolute_tolerance",
  "relative_tolerance", "subdivisions", "first_finite_intervals",
  "second_finite_intervals", "outcome_label_state",
  "biological_outcomes_computed"
)]
public_selection$contract_id <- "mv08x_oracle_selection_public_v1"

implementation_paths <- c(
  "rust/scph_landscape_kernel/Cargo.toml",
  "rust/scph_landscape_kernel/Cargo.lock",
  "rust/scph_landscape_kernel/src/lib.rs",
  "rust/scph_landscape_kernel/include/scph_landscape_kernel.h",
  "R/landscape_contract.R", "R/landscape_reference.R",
  "R/landscape_rust_prototype.R", "R/mv08s_ph_sentinel.R",
  "scripts/build_mv08x_rust_landscape_admission_prefreeze.R",
  "scripts/acquire_mv08x_rust_toolchain.sh",
  "scripts/run_mv08x_rust_rebuild.sh",
  "scripts/run_mv08x_rust_landscape_oracles.R",
  "scripts/build_mv08y_rust_landscape_admission_closure.R"
)
if (!all(file.exists(implementation_paths))) {
  stop("MV8-X implementation set incomplete", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv08x_implementation_binding_v1",
  role = c(
    "cargo_manifest", "cargo_lock", "rust_kernel", "c_abi_header",
    "r_landscape_contract", "r_canonical_oracle", "r_rust_shim",
    "ph_record_validator", "prefreeze_builder", "toolchain_acquirer",
    "rust_rebuilder", "oracle_runner", "closure_builder"
  ),
  file = implementation_paths,
  bytes = as.numeric(file.info(implementation_paths)$size),
  sha256 = vapply(implementation_paths, sha_file, character(1L)),
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
  file.path(s_root, "mv08s-ph-metrics.csv"),
  file.path(s_root, "mv08s-execution-decision.csv"),
  file.path(v_root, "mv08v-selected-ph-metrics.csv"),
  file.path(v_root, "mv08v-progress.csv"),
  private_selection_path
)
input_manifest <- data.frame(
  contract_id = "mv08x_input_manifest_v1",
  role = c(
    "mv08r_manifest", "mv08r_ph_queue", "mv08r_landscape_queue",
    "mv08r_landscape_contract", "mv08w_manifest", "mv08w_inventory",
    "mv08w_validation", "mv08w_decision", "mv08s_ph_metrics",
    "mv08s_decision", "mv08v_ph_metrics", "mv08v_progress",
    "private_oracle_selection"
  ),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha_file, character(1L)),
  rows = c(
    r_manifest$rows, nrow(r_queue), nrow(groups), nrow(landscape_contract),
    w_manifest$rows, nrow(w_inventory), nrow(w_validation), nrow(w_decision),
    nrow(s_metrics), nrow(s_decision), nrow(v_metrics), nrow(v_progress),
    nrow(private_selection)
  ), stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv08x_rust_landscape_admission_contract_v1",
  engine_id = "rust_scph_landscape_kernel_v1",
  scientific_object = "dimension_specific_exact_all_active_squared_l2",
  finite_interval_policy = required_contract[["finite_intervals"]],
  essential_h0_policy = required_contract[["essential_h0"]],
  level_policy = required_contract[["level_policy"]],
  integration_policy = required_contract[["integration"]],
  dimension_policy = required_contract[["dimension_policy"]],
  grid_policy = required_contract[["grid_policy"]],
  level_cap_policy = required_contract[["level_cap_policy"]],
  streaming_policy = required_contract[["streaming"]],
  r_oracle = "landscape_reference_v4",
  rust_toolchain = "1.97.1-x86_64-unknown-linux-gnu",
  rustup_installer_version = "1.29.0",
  rustup_installer_url = paste0(
    "https://static.rust-lang.org/rustup/archive/1.29.0/",
    "x86_64-unknown-linux-gnu/rustup-init"
  ),
  rustup_installer_sha256 =
    "4acc9acc76d5079515b46346a485974457b5a79893cfb01112423c89aeb5aa10",
  external_rust_crates = 0L, build_jobs = 1L, oracle_workers = 1L,
  automatic_retries = 0L, production_landscape_jobs = 0L,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

resource <- data.frame(
  contract_id = "mv08x_resource_policy_v1",
  stage = c("clean_build_a", "clean_build_b", "oracle_run_a", "oracle_run_b"),
  elapsed_cap_seconds = c(3600, 3600, 3600, 3600),
  rss_cap_bytes = c(8, 8, 12, 12) * 1024^3,
  workers = 1L, retries = 0L,
  private_storage_cap_bytes = c(4, 4, 4, 4) * 1024^3,
  stringsAsFactors = FALSE
)

acceptance <- data.frame(
  contract_id = "mv08x_acceptance_v1",
  gate = c(
    "toolchain_identity", "installer_hash", "source_lock_hashes",
    "format_unit_clippy", "two_clean_builds", "native_c_abi",
    "dependency_inventory", "analytical_fixtures", "failure_fallback",
    "all_28_groups", "r_oracle_equivalence", "reverse_invariant",
    "self_zero_invariant", "two_run_scientific_identity",
    "resource_caps", "aggregate_only_publication", "downstream_firewalls"
  ),
  required_state = c(
    "rustc_1.97.1_host_x86_64_unknown_linux_gnu",
    "official_rustup_1.29.0_sha256_exact",
    "all_prefrozen_implementation_hashes_exact", "all_pass",
    "byte_identical_candidate_library", "success_and_error_paths_pass",
    "recorded_and_no_unexpected_project_dependency", "all_pass",
    "missing_and_corrupt_library_use_canonical_R", "28_of_28",
    "exact_or_within_R_error_certificate", "bit_identical_and_counts_swap",
    "exact_zero", "byte_identical_normalized_results", "all_pass",
    "no_unit_accession_donor_or_private_path", "all_zero_and_closed"
  ), stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "mv08r_manifest", "mv08w_manifest", "full_ph_closed",
    "full_ph_identity", "landscape_contract", "group_cardinality",
    "dimension_coverage", "scope_coverage", "panel_coverage",
    "representation_coverage", "selection_cardinality", "selection_distinct",
    "selection_routes", "private_selection_bound", "implementation_bound",
    "toolchain_frozen", "resource_policy", "production_closed",
    "downstream_closed", "privacy_schema"
  ),
  passed = c(
    r_manifest$rows > 0L, w_manifest$rows == 7L,
    nrow(w_validation) == 20L && all(truth(w_validation$passed)),
    identical(sort(metrics$job_id), sort(w_inventory$job_id)),
    identical(unname(observed_contract[names(required_contract)]),
              unname(required_contract)),
    nrow(groups) == 28L && all(as.integer(groups$units) >= 2L),
    identical(sort(unique(groups$homology_dimension)), c("H0", "H1")) &&
      all(table(groups$homology_dimension) == 14L),
    identical(sort(unique(groups$dataset_scope)), c("external8", "internal124")),
    identical(sort(unique(groups$panel_id)), c("common475", "exact500")),
    length(unique(groups$representation_id)) == 3L,
    nrow(private_selection) == 28L && nrow(public_selection) == 28L,
    all(private_selection$first_job_id != private_selection$second_job_id),
    all(private_selection$reference_route %in%
          c("r_exact_breakpoint", "r_adaptive_certified")),
    sha_file(private_selection_path) == tail(input_manifest$sha256, 1L),
    nrow(implementation) == length(implementation_paths),
    contract$rust_toolchain == "1.97.1-x86_64-unknown-linux-gnu" &&
      nchar(contract$rustup_installer_sha256) == 64L,
    all(resource$workers == 1L) && all(resource$retries == 0L),
    contract$production_landscape_jobs == 0L,
    all(unlist(contract[c(
      "comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
      "outcome_jobs"
    )], use.names = FALSE) == 0) &&
      contract$outcome_label_state == "closed" &&
      !contract$biological_outcomes_computed,
    !any(c("job_id", "unit_id", "output_file", "source_role", "accession",
           "private_path") %in% names(public_selection))
  ),
  evidence = c(
    "MV8-R manifest independently rehashes", "MV8-W manifest independently rehashes",
    "MV8-W 20/20 closure remains accepted", "1,280 PH identities align",
    "owner-approved landscape definition is unchanged", "28 frozen groups",
    "14 H0 plus 14 H1 groups", "internal and external scopes covered",
    "common475 and exact500 covered", "all three representations covered",
    "one stress pair per group", "every stress pair uses distinct records",
    "reference route is fixed from finite-interval depth",
    "private locator selection is SHA-bound", "all execution sources are hash-bound",
    "official isolated Rust toolchain identity is frozen",
    "one worker, zero retry, bounded resources", "production landscapes remain zero",
    "all downstream and outcome jobs remain zero", "public selection is aggregate-only"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-X prefreeze validation failed", call. = FALSE)

decision <- data.frame(
  contract_id = "mv08x_prefreeze_decision_v1",
  decision = "authorize_isolated_toolchain_rebuild_and_28_group_R_oracle_admission_only",
  full_ph_records = 1280L, landscape_groups = 28L,
  oracle_pairs = nrow(public_selection), production_landscape_jobs = 0L,
  rust_binary_state = "absent_rebuild_required",
  next_gate = "MV8Y_independent_Rust_admission_closure",
  landscape_execution_authorized = FALSE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, validations_passed = sum(validation$passed),
  validations_total = nrow(validation), stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(contract, file.path(output_dir, "mv08x-contract.csv"))
atomic_csv(public_selection, file.path(output_dir, "mv08x-oracle-selection.csv"))
atomic_csv(resource, file.path(output_dir, "mv08x-resource-policy.csv"))
atomic_csv(acceptance, file.path(output_dir, "mv08x-acceptance.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08x-implementation-bindings.csv"))
atomic_csv(input_manifest, file.path(output_dir, "mv08x-input-manifest.csv"))
atomic_csv(validation, file.path(output_dir, "mv08x-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08x-decision.csv"))
atomic_text(c(
  "# MV8-X Rust landscape admission prefreeze", "",
  "**Result:** 20/20 prospective checks pass; rebuild/oracle admission only.", "",
  paste0(
    "MV8-X freezes one deterministic finite-interval stress pair for each of ",
    "the 28 post-MV8-W landscape groups. The corpus covers 14 H0 and 14 H1 ",
    "groups across both dataset scopes, both panels, and all three representations."
  ), "",
  paste0(
    "The isolated official Rust 1.97.1 toolchain may be reacquired and the ",
    "dependency-free kernel rebuilt twice. Admission requires canonical R ",
    "equivalence, reverse/self invariants, deterministic builds and runs, native ",
    "ABI/fallback checks, and all resource/privacy gates."
  ), "",
  paste0(
    "No production landscape, comparison, clustering, fusion, label, outcome, ",
    "adoption, or manuscript-claim job is authorized."
  )
), file.path(output_dir, "MV08X_RUST_LANDSCAPE_ADMISSION_PREFREEZE.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08x-artifact-manifest.csv"))
cat("MV8-X prefreeze checks=", sum(validation$passed), "/", nrow(validation),
    "; oracle_pairs=", nrow(public_selection), "; production_landscapes=0\n", sep = "")
