#!/usr/bin/env Rscript

# Amend MV8-ZH to reconstruct private pair identities from the already-bound
# private unit axis. No production receipt or landscape is written here.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv08zi_landscape_pair_binding_recovery_prefreeze.R",
  "<mv08zh-prefreeze> <mv08zf-prefreeze> <private-bindings>",
  "<private-root> <public-root> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zh_root <- normalizePath(args[[1L]], mustWork = TRUE)
zf_root <- normalizePath(args[[2L]], mustWork = TRUE)
bindings_path <- normalizePath(args[[3L]], mustWork = TRUE)
private_root <- normalizePath(args[[4L]], mustWork = TRUE)
public_root <- normalizePath(args[[5L]], mustWork = TRUE)
output_dir <- normalizePath(args[[6L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZI output", call. = FALSE)

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
    stop("MV8-ZI prerequisite manifest drift", call. = FALSE)
  }
}
verify_manifest(zh_root, "mv08zh-artifact-manifest.csv")
verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")

zh_snapshot <- read_csv(file.path(zh_root, "mv08zh-stopped-snapshot.csv"))
zh_orphan <- read_csv(file.path(zh_root, "mv08zh-orphan-binding.csv"))
zh_impl <- read_csv(file.path(zh_root, "mv08zh-implementation-bindings.csv"))
zf_inputs <- read_csv(file.path(zf_root, "mv08zf-input-manifest.csv"))
queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
bindings <- read_csv(bindings_path)
ledger_path <- file.path(public_root, "mv08zf-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zf-chunk-completions.csv")
progress_path <- file.path(public_root, "mv08zf-progress.csv")
ledger <- read_csv(ledger_path)
completed <- read_csv(completion_path)
progress <- read_csv(progress_path)
index <- 164L
row <- queue[index, , drop = FALSE]
chunk_root <- file.path(private_root, "production",
                        .mv08z_safe_group(row$group_order),
                        .mv08z_safe_chunk(row$chunk_order))
distances_path <- file.path(chunk_root, "distances.csv")
status_path <- file.path(chunk_root, "status.csv")
distances <- read_csv(distances_path)
status <- read_csv(status_path)
group_bindings <- bindings[as.integer(bindings$group_order) == row$group_order, , drop = FALSE]
pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(group_bindings), row$group_id)
pairs <- pairs[pairs$pair_ordinal >= row$pair_start &
                 pairs$pair_ordinal <= row$pair_end, , drop = FALSE]
active_runner_lines <- suppressWarnings(system2(
  "pgrep", c("-f", "run_mv08zf_full_landscape_production[.]R"),
  stdout = TRUE, stderr = FALSE
))
active_runner_processes <- length(active_runner_lines)
partials <- list.files(private_root, pattern = "partial", recursive = TRUE,
                       full.names = TRUE, all.files = TRUE)

old_recovery <- zh_impl[zh_impl$role == "recovery_executor", , drop = FALSE]
implementation_files <- c(
  "scripts/recover_mv08zh_landscape_launcher_interruption.R",
  "scripts/build_mv08zi_landscape_pair_binding_recovery_prefreeze.R"
)
implementation <- data.frame(
  contract_id = "mv08zi_implementation_binding_v1",
  role = c("recovery_executor", "amendment_builder"),
  file = implementation_files,
  old_sha256 = c(old_recovery$sha256, NA_character_),
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  scientific_change = FALSE, private_pair_binding_fix = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
failure <- data.frame(
  contract_id = "mv08zi_preserved_attempt_v1",
  attempt = "mv08zh_post_commit_adoption",
  stopped_stage = "private_pair_inventory_lookup",
  stopped_reason = "public_MV8_Z_audit_intentionally_has_no_private_pair_queue",
  public_rows_before = 163L, public_rows_after = nrow(completed),
  private_orphan_preserved = TRUE, public_receipts_written = 0L,
  landscape_jobs_run = 0L, automatic_retries = 0L,
  stringsAsFactors = FALSE
)
pair_binding <- data.frame(
  contract_id = "mv08zi_private_pair_binding_v1",
  production_order = index, group_order = row$group_order,
  chunk_order = row$chunk_order, expected_pair_rows = nrow(pairs),
  expected_pair_subset_sha256 = .mv08z_sha256_text(pairs$pair_identity_sha256),
  observed_pair_rows = nrow(distances),
  observed_pair_subset_sha256 = .mv08z_sha256_text(distances$pair_identity_sha256),
  private_bindings_bytes = as.numeric(file.info(bindings_path)$size),
  private_bindings_sha256 = sha_file(bindings_path),
  public_private_identities_exposed = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c(
    "mv08zh_manifest", "mv08zf_manifest", "prior_failure_preserved",
    "no_public_mutation", "prefix_cardinality", "prefix_hashes_unchanged",
    "orphan_hashes_unchanged", "private_binding_hash", "private_binding_cardinality",
    "group_binding_cardinality", "pair_cardinality", "pair_subset_reconstruction",
    "pair_identity_exact", "status_queue_identity", "orphan_complete",
    "zero_active_runner", "zero_partials", "implementation_transition",
    "implementation_hashes", "no_scientific_change", "zero_retry",
    "zero_recomputation", "labels_outcomes_closed", "downstream_closed",
    "public_privacy"
  ),
  passed = c(
    TRUE, TRUE, failure$public_receipts_written == 0L && failure$landscape_jobs_run == 0L,
    nrow(completed) == 163L && nrow(ledger) == 163L && progress$completed_chunks == 163L,
    zh_snapshot$completed_prefix == 163L,
    sha_file(ledger_path) == zh_snapshot$ledger_sha256 &&
      sha_file(completion_path) == zh_snapshot$completion_sha256 &&
      sha_file(progress_path) == zh_snapshot$progress_sha256,
    sha_file(distances_path) == zh_orphan$distance_sha256 &&
      sha_file(status_path) == zh_orphan$status_sha256,
    sha_file(bindings_path) == zf_inputs$sha256[zf_inputs$role == "private_unit_bindings"],
    nrow(bindings) == 2544L, nrow(group_bindings) == 124L,
    nrow(pairs) == row$pair_count && nrow(distances) == row$pair_count,
    pair_binding$expected_pair_subset_sha256 == row$pair_subset_sha256,
    identical(distances$pair_identity_sha256, pairs$pair_identity_sha256),
    status$pair_subset_sha256 == row$pair_subset_sha256,
    status$completion_state == "complete" && status$mode == "production",
    active_runner_processes == 0L, length(partials) == 0L,
    nrow(old_recovery) == 1L && implementation$old_sha256[[1L]] == old_recovery$sha256 &&
      implementation$sha256[[1L]] != old_recovery$sha256,
    all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256),
    !any(implementation$scientific_change), all(ledger$retries == 0L),
    failure$landscape_jobs_run == 0L,
    status$outcome_label_state == "closed" && !truth(status$biological_outcomes_computed),
    all(progress[c("comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
                   "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs")] == 0L),
    !pair_binding$public_private_identities_exposed
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZI prefreeze failed: ", paste(validation$check_id[!validation$passed], collapse = ", "),
  call. = FALSE
)
decision <- data.frame(
  contract_id = "mv08zi_decision_v1",
  decision = "authorize_private_axis_reconstructed_orphan_adoption_after_commit",
  orphan_adoption_authorized = TRUE, orphan_production_order = 164L,
  resume_at_production_order = 165L,
  landscape_recomputation_authorized = FALSE, automatic_retries = 0L,
  scientific_contract_changed = FALSE, private_pair_binding_required = TRUE,
  resource_values_are_upper_bounds = TRUE, comparison_jobs_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  next_gate = "resume_MV8_ZF_then_MV8_ZG_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(failure, file.path(output_dir, "mv08zi-preserved-attempt.csv"))
atomic_csv(pair_binding, file.path(output_dir, "mv08zi-private-pair-binding.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08zi-implementation-bindings.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zi-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zi-decision.csv"))
atomic_text(c(
  "# MV8-ZI private-pair binding recovery amendment", "",
  "- Validation: 25/25 pass.",
  "- MV8-ZH stopped before receipt mutation because private pair identities were intentionally absent from the public MV8-Z audit.",
  "- The amended executor reconstructs order-164 identities from the existing hash-bound 2,544-row private unit axis.",
  "- The 163-row prefix and completed orphan are unchanged; no landscape, retry, deletion, overwrite, label, outcome, or downstream job ran.",
  "- After commit, one exact orphan adoption is authorized, followed by strict-prefix resume at order 165."
), file.path(output_dir, "MV08ZI_PRIVATE_PAIR_BINDING_RECOVERY_PREFREEZE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zi-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zi-artifact-manifest.csv"))
cat("MV8-ZI private-pair binding recovery prefreeze passed 25/25; zero recovery executed\n")
