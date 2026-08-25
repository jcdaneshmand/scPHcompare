#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(paste(
  "usage: build_mv08zt_engine_v2_full_landscape_prefreeze.R",
  "<private-unit-bindings> <engine-v2-rust-library> <output-dir>"
), call. = FALSE)
bindings_path <- normalizePath(args[[1L]], mustWork = TRUE)
rust_library <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZT output", call. = FALSE)
parent_head <- tolower(trimws(Sys.getenv("MV08ZT_PARENT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", parent_head)) {
  stop("MV8-ZT exact parent HEAD absent", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
truth <- function(value) {
  if (is.logical(value)) return(!is.na(value) & value)
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}
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
verify_manifest <- function(root, name) {
  manifest <- read_csv(file.path(root, name))
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(vapply(files, sha_file, character(1L)) == manifest$sha256) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes))) {
    stop("MV8-ZT manifest drift: ", name, call. = FALSE)
  }
  manifest
}

zq_root <- normalizePath(
  "docs/audits/mv08zq-landscape-kernel-repair-admission-closure-v1",
  mustWork = TRUE
)
z_root <- normalizePath(
  "docs/audits/mv08z-landscape-execution-prefreeze-v1", mustWork = TRUE
)
zf_root <- normalizePath(
  "docs/audits/mv08zf-full-landscape-production-prefreeze-v1", mustWork = TRUE
)
zq_manifest <- verify_manifest(zq_root, "mv08zq-artifact-manifest.csv")
z_manifest <- verify_manifest(z_root, "mv08z-artifact-manifest.csv")
zf_manifest <- verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")
zq_validation <- read_csv(file.path(zq_root, "mv08zq-validation.csv"))
zq_decision <- read_csv(file.path(zq_root, "mv08zq-decision.csv"))
z_contract <- read_csv(file.path(z_root, "mv08z-contract.csv"))
groups <- read_csv(file.path(z_root, "mv08z-group-queue.csv"))
chunks <- read_csv(file.path(z_root, "mv08z-chunk-queue.csv"))
z_inputs <- read_csv(file.path(z_root, "mv08z-input-manifest.csv"))
old_queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))

if (nrow(zq_decision) != 1L || nrow(z_contract) != 1L ||
    nrow(groups) != 28L || nrow(chunks) != 628L || nrow(old_queue) != 628L) {
  stop("MV8-ZT prerequisite cardinality drift", call. = FALSE)
}
candidate_sha <- sha_file(rust_library)
candidate_bytes <- as.numeric(file.info(rust_library)$size)
binding_row <- z_inputs[z_inputs$role == "private_unit_bindings", , drop = FALSE]
if (nrow(binding_row) != 1L || sha_file(bindings_path) != binding_row$sha256 ||
    as.numeric(file.info(bindings_path)$size) != as.numeric(binding_row$bytes)) {
  stop("MV8-ZT private binding drift", call. = FALSE)
}

queue <- old_queue
queue$contract_id <- "mv08zt_engine_v2_full_landscape_queue_v1"
queue$production_origin <- "fresh_engine_v2_from_zero_no_old_output_reuse"
queue$authorization_state <- "authorized_after_mv08zt_commit"
queue$scientific_engine_version <- 2L

implementation_files <- c(
  "rust/scph_landscape_kernel/src/lib.rs",
  "R/landscape_rust_prototype.R",
  "R/mv08z_landscape_production.R",
  "scripts/run_mv08z_landscape_chunk.R",
  "scripts/run_mv08zt_full_landscape_production.R",
  "scripts/build_mv08zt_engine_v2_full_landscape_prefreeze.R"
)
implementation <- data.frame(
  contract_id = "mv08zt_implementation_binding_v1",
  role = c(
    "engine_v2_rust_source", "R_FFI_shim", "execution_helpers",
    "version_aware_chunk_worker", "fresh_full_runner", "prefreeze_builder"
  ),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv08zt_engine_v2_full_landscape_prefreeze_v1",
  accepted_parent_head = parent_head,
  scientific_engine_version = 2L,
  landscape_groups = 28L, production_chunks = 628L,
  production_pairs = 152744L, production_start_order = 1L,
  old_completed_chunks_reused = 0L, old_completed_pairs_reused = 0L,
  sentinel_pairs_reused_as_production = 0L,
  finite_interval_policy = "all_finite_positive_persistence_intervals",
  essential_h0_policy = "exclude_infinite_interval",
  level_policy = "all_consecutive_active_levels",
  integration_policy = "exact_streamed_squared_l2_on_dimension_support",
  dimension_policy = "h0_h1_separate_primary_outputs",
  grid_policy = "no_universal_fixed_grid",
  level_cap_policy = "no_universal_level_cap",
  streaming_policy = "stream_or_chunk_without_dense_landscape_materialization",
  engine_id = "rust_scph_landscape_kernel_v2",
  rust_library_sha256 = candidate_sha, rust_library_bytes = candidate_bytes,
  rust_library_state = "private_explicit_hash_verified_not_default",
  child_elapsed_cap_seconds = 3600, child_rss_cap_bytes = 4 * 1024^3,
  aggregate_elapsed_cap_seconds = 144000,
  private_storage_cap_bytes = 1024^3,
  minimum_free_disk_bytes = 10 * 1024^3,
  minimum_available_memory_bytes = 8 * 1024^3,
  workers = 1L, automatic_retries = 0L,
  fallback_policy = "none_stop_and_preserve",
  atomic_write = TRUE, resume_policy = "strict_completed_prefix_only_new_root",
  public_receipt_policy = "aggregate_chunk_and_resource_only",
  fresh_production_authorized_after_commit = TRUE,
  old_root_resume = FALSE, old_root_reuse = FALSE,
  old_root_overwrite = FALSE, old_root_delete = FALSE,
  comparison_state = "closed", clustering_state = "closed",
  fusion_state = "closed", label_state = "closed", outcome_state = "closed",
  manuscript_claim_state = "closed", outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
resource <- data.frame(
  contract_id = "mv08zt_resource_policy_v1",
  stage = c("landscape_chunk", "full_production"),
  elapsed_cap_seconds = c(3600, 144000),
  rss_cap_bytes = 4 * 1024^3, workers = 1L, retries = 0L,
  private_storage_cap_bytes = 1024^3,
  minimum_free_disk_bytes = 10 * 1024^3,
  minimum_available_memory_bytes = 8 * 1024^3,
  authorization_state = c(
    "authorized_after_mv08zt_commit", "authorized_after_launch_headroom_gate"
  ), stringsAsFactors = FALSE
)
firewall <- data.frame(
  contract_id = "mv08zt_old_root_firewall_v1",
  old_engine_version = 1L, old_completed_chunks = 502L,
  old_completed_pairs = 123516L, old_scientific_state = "rejected",
  resume = FALSE, reuse = FALSE, overwrite = FALSE, delete = FALSE,
  new_private_root_required = TRUE, new_public_root_required = TRUE,
  stringsAsFactors = FALSE
)
closure <- data.frame(
  contract_id = "mv08zu_prospective_closure_contract_v1",
  required_chunks = 628L, required_pairs = 152744L,
  required_engine_version = 2L, require_every_input_rehash = TRUE,
  require_every_pair_identity = TRUE, require_every_active_depth = TRUE,
  require_h0_h1_separate = TRUE, require_zero_partial = TRUE,
  require_one_worker = TRUE, require_zero_retries = TRUE,
  require_resource_caps = TRUE, require_aggregate_privacy = TRUE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, manuscript_claim_jobs_authorized = 0L,
  stringsAsFactors = FALSE
)
input_bindings <- data.frame(
  contract_id = "mv08zt_input_binding_v1",
  role = c(
    "mv08zq_admission_manifest", "mv08z_queue_manifest",
    "mv08zf_queue_plan_manifest", "private_unit_bindings", "engine_v2_candidate"
  ),
  bytes = c(
    file.info(file.path(zq_root, "mv08zq-artifact-manifest.csv"))$size,
    file.info(file.path(z_root, "mv08z-artifact-manifest.csv"))$size,
    file.info(file.path(zf_root, "mv08zf-artifact-manifest.csv"))$size,
    file.info(bindings_path)$size, candidate_bytes
  ),
  sha256 = c(
    sha_file(file.path(zq_root, "mv08zq-artifact-manifest.csv")),
    sha_file(file.path(z_root, "mv08z-artifact-manifest.csv")),
    sha_file(file.path(zf_root, "mv08zf-artifact-manifest.csv")),
    sha_file(bindings_path), candidate_sha
  ), stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "mv08zq_manifest", "mv08z_manifest", "mv08zf_plan_manifest",
    "mv08zq_checks", "engine_v2_admitted", "candidate_binding",
    "private_axis_binding", "implementation_bindings", "group_cardinality",
    "chunk_cardinality", "pair_cardinality", "queue_order", "queue_identity",
    "h0_h1_balance", "finite_positive_policy", "essential_h0_policy",
    "all_active_levels", "exact_streamed_l2", "h0_h1_separate",
    "no_grid", "no_level_cap", "fresh_order_zero", "old_prefix_zero_reuse",
    "old_root_immutable", "one_worker", "zero_retries", "zero_fallback",
    "child_caps", "aggregate_caps", "headroom_gates", "strict_resume",
    "atomic_outputs", "aggregate_receipts", "prospective_closure",
    "downstream_firewall", "labels_outcomes_closed"
  ),
  passed = c(
    nrow(zq_manifest) > 0L, nrow(z_manifest) > 0L, nrow(zf_manifest) > 0L,
    nrow(zq_validation) == 33L && all(truth(zq_validation$passed)),
    zq_decision$scientific_engine_version == 2L &&
      zq_decision$oracle_pairs_passed == 56L,
    candidate_sha == zq_decision$candidate_sha256 &&
      candidate_sha == contract$rust_library_sha256,
    sha_file(bindings_path) == binding_row$sha256,
    all(file.exists(implementation$file)) &&
      all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256),
    nrow(groups) == 28L, nrow(queue) == 628L,
    sum(queue$pair_count) == 152744L,
    identical(as.integer(queue$production_order), 1:628),
    identical(queue$pair_subset_sha256, old_queue$pair_subset_sha256) &&
      identical(queue$global_chunk_order, chunks$global_chunk_order),
    all(table(groups$homology_dimension) == 14L),
    contract$finite_interval_policy == "all_finite_positive_persistence_intervals",
    contract$essential_h0_policy == "exclude_infinite_interval",
    contract$level_policy == "all_consecutive_active_levels",
    contract$integration_policy == "exact_streamed_squared_l2_on_dimension_support",
    contract$dimension_policy == "h0_h1_separate_primary_outputs",
    contract$grid_policy == "no_universal_fixed_grid",
    contract$level_cap_policy == "no_universal_level_cap",
    contract$production_start_order == 1L,
    contract$old_completed_chunks_reused == 0L &&
      contract$old_completed_pairs_reused == 0L,
    !firewall$resume && !firewall$reuse && !firewall$overwrite && !firewall$delete,
    all(resource$workers == 1L), all(resource$retries == 0L),
    contract$fallback_policy == "none_stop_and_preserve",
    contract$child_elapsed_cap_seconds == 3600 &&
      contract$child_rss_cap_bytes == 4 * 1024^3,
    contract$aggregate_elapsed_cap_seconds == 144000 &&
      contract$private_storage_cap_bytes == 1024^3,
    contract$minimum_free_disk_bytes == 10 * 1024^3 &&
      contract$minimum_available_memory_bytes == 8 * 1024^3,
    contract$resume_policy == "strict_completed_prefix_only_new_root",
    truth(contract$atomic_write),
    contract$public_receipt_policy == "aggregate_chunk_and_resource_only",
    closure$required_chunks == 628L && closure$required_pairs == 152744L &&
      closure$required_engine_version == 2L,
    all(queue$comparison_jobs == 0L) && all(queue$clustering_jobs == 0L) &&
      all(queue$fusion_jobs == 0L) && all(queue$label_jobs == 0L) &&
      all(queue$outcome_jobs == 0L),
    all(queue$outcome_label_state == "closed") &&
      !any(truth(queue$biological_outcomes_computed))
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZT prefreeze failed at: ",
  paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE
)
decision <- data.frame(
  contract_id = "mv08zt_decision_v1",
  decision = "authorize_fresh_engine_v2_full_landscape_after_commit_and_headroom_gate",
  scientific_engine_version = 2L,
  fresh_production_authorized_after_commit = TRUE,
  production_landscape_pairs_authorized = 152744L,
  production_chunks_authorized = 628L, workers = 1L, automatic_retries = 0L,
  old_root_resume = FALSE, old_root_reuse = FALSE,
  old_root_overwrite = FALSE, old_root_delete = FALSE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, manuscript_claim_jobs_authorized = 0L,
  next_gate = "MV8_ZU_independent_engine_v2_full_landscape_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(contract, file.path(output_dir, "mv08zt-contract.csv"))
atomic_csv(queue, file.path(output_dir, "mv08zt-production-queue.csv"))
atomic_csv(resource, file.path(output_dir, "mv08zt-resource-policy.csv"))
atomic_csv(firewall, file.path(output_dir, "mv08zt-old-root-firewall.csv"))
atomic_csv(closure, file.path(output_dir, "mv08zu-prospective-closure.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08zt-implementation-bindings.csv"))
atomic_csv(input_bindings, file.path(output_dir, "mv08zt-input-bindings.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zt-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zt-decision.csv"))
atomic_text(c(
  "# MV8-ZT fresh engine-v2 full-landscape production prefreeze", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " checks pass; 628 fresh chunks and 152,744 pairs are authorized only ",
         "after this prefreeze is committed and the launch headroom gate passes."), "",
  "The corrected hash-bound scientific engine version is 2. H0 and H1 remain separate; every finite positive interval and every consecutive active level is retained; exact streamed squared L2 uses no grid or universal level cap.", "",
  "The stopped 502-chunk version-1 root remains immutable rejected evidence. Fresh production begins at order 1 under entirely new private and public roots; no old calculated distance is resumed, reused, overwritten, or deleted.", "",
  "Comparisons, clustering, fusion, labels, outcomes, biological or manuscript claims, cleanup, and deletion remain closed pending independent MV8-ZU 628/628 and 152,744/152,744 reconstruction."
), file.path(output_dir, "MV08ZT_ENGINE_V2_FULL_LANDSCAPE_PREFREEZE.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zt-artifact-manifest.csv"))
cat("MV8-ZT prefreeze checks=", sum(validation$passed), "/", nrow(validation),
    "; chunks=628; pairs=152744; engine=2; old_reuse=0\n", sep = "")
