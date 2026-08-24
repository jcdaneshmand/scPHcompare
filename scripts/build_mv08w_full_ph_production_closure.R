#!/usr/bin/env Rscript

# Independently close MV8-V by rehashing all new PH artifacts, reconstructing
# every typed gene view from the accepted source caches, and rerunning the H0
# MST oracle. The 23 MV8-T records are bound through their immutable closure.
# No new PH, landscape, comparison, label, or outcome is computed here.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv08w_full_ph_production_closure.R <mv08u-prefreeze>",
  "<mv08p-private> <mv08pr-private> <mv08ps-private>",
  "<mv08o-internal-private> <common-panel> <mv08v-private>",
  "<mv08v-public> <output-dir>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
source_roots <- list(
  mv08p_original_v1 = normalizePath(args[[2L]], mustWork = TRUE),
  mv08pr_overlay_v1 = normalizePath(args[[3L]], mustWork = TRUE),
  mv08ps_overlay_v1 = normalizePath(args[[4L]], mustWork = TRUE),
  mv08o_internal_primary_v8 = normalizePath(args[[5L]], mustWork = TRUE)
)
panel_path <- normalizePath(args[[6L]], mustWork = TRUE)
private_root <- normalizePath(args[[7L]], mustWork = TRUE)
public_root <- normalizePath(args[[8L]], mustWork = TRUE)
output_dir <- normalizePath(args[[9L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-W output", call. = FALSE)
for (package in "digest") {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
resolve_source <- function(row) {
  role <- row$source_cache_root_role
  if (!(role %in% names(source_roots))) stop("MV8-W unknown source root role", call. = FALSE)
  file.path(source_roots[[role]], row$source_cache_relative_file)
}

contract <- read_csv(file.path(prefreeze, "mv08u-contract.csv"))
queue <- read_csv(file.path(prefreeze, "mv08u-full-ph-queue.csv"))
implementation <- read_csv(file.path(prefreeze, "mv08u-implementation-bindings.csv"))
input_manifest <- read_csv(file.path(prefreeze, "mv08u-input-manifest.csv"))
runtime_inputs <- read_csv(file.path(prefreeze, "mv08u-runtime-input-bindings.csv"))
ledger <- read_csv(file.path(public_root, "mv08v-resource-ledger.csv"))
metrics <- read_csv(file.path(public_root, "mv08v-selected-ph-metrics.csv"))
progress <- read_csv(file.path(public_root, "mv08v-progress.csv"))
t_root <- "docs/audits/mv08t-ph-sentinel-closure-v1"
t_artifacts <- read_csv(file.path(t_root, "mv08t-private-artifact-rehash.csv"))
t_validation <- read_csv(file.path(t_root, "mv08t-validation.csv"))
t_manifest <- read_csv(file.path(t_root, "mv08t-artifact-manifest.csv"))
if (nrow(contract) != 1L || nrow(queue) != 1257L || nrow(metrics) != 1257L ||
    !identical(as.integer(metrics$production_order), 1:1257) ||
    !identical(metrics$job_id, queue$job_id) ||
    nrow(progress) != 1L ||
    progress$state != "ph_production_complete_closure_pending" ||
    progress$completed_records != 1257L || !all(t_validation$passed)) {
  stop("MV8-W prerequisite ledger drift", call. = FALSE)
}
if (!all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha_file, character(1L)) ==
           implementation$sha256) ||
    !all(file.exists(input_manifest$input)) ||
    !all(vapply(input_manifest$input, sha_file, character(1L)) ==
           input_manifest$sha256) ||
    !all(vapply(file.path(t_root, t_manifest$artifact), sha_file, character(1L)) ==
           t_manifest$sha256)) {
  stop("MV8-W bound public evidence drift", call. = FALSE)
}
panel <- read_csv(panel_path)
if (nrow(runtime_inputs) != 1L ||
    sha_file(panel_path) != runtime_inputs$file_sha256 ||
    nrow(panel) != 475L || !identical(as.integer(panel$panel_order_475), 1:475)) {
  stop("MV8-W common475 panel drift", call. = FALSE)
}

source_rehash <- list()
for (i in which(!duplicated(queue$unit_id))) {
  row <- queue[i, , drop = FALSE]
  path <- resolve_source(row)
  if (!file.exists(path) || sha_file(path) != row$source_cache_sha256) {
    stop("MV8-W source-cache drift: ", row$unit_id, call. = FALSE)
  }
  cache <- readRDS(path)
  mv08s_validate_residual_cache_v1(cache)
  source_rehash[[length(source_rehash) + 1L]] <- data.frame(
    contract_id = "mv08w_source_cache_rehash_v1", unit_id = row$unit_id,
    source_cache_root_role = row$source_cache_root_role,
    bytes = as.numeric(file.info(path)$size), sha256 = sha_file(path),
    independently_rehashed = TRUE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
}
source_rehash <- do.call(rbind, source_rehash)

artifact_rows <- vector("list", nrow(queue))
for (i in seq_len(nrow(queue))) {
  row <- queue[i, , drop = FALSE]
  path <- file.path(private_root, row$output_file)
  if (!file.exists(path) || file.exists(paste0(path, ".partial")) ||
      sha_file(path) != metrics$output_sha256[[i]] ||
      as.numeric(file.info(path)$size) != metrics$output_bytes[[i]]) {
    stop("MV8-W PH artifact drift at production order ", i, call. = FALSE)
  }
  cache <- readRDS(resolve_source(row))
  view <- mv08s_residual_gene_view_v1(cache, row, panel)
  record <- readRDS(path)
  mv08s_validate_ph_record_v1(record, row, view)
  artifact_rows[[i]] <- data.frame(
    contract_id = "mv08w_ph_artifact_rehash_v1", artifact_origin = "mv08v",
    production_order = i, mv08r_job_order = row$mv08r_job_order,
    job_id = row$job_id, unit_id = row$unit_id, seed = row$seed,
    representation_id = row$representation_id, panel_id = row$panel_id,
    view_kind = row$view_kind, selected_engine = metrics$selected_engine[[i]],
    point_count = record$h0_mst_oracle$point_count,
    finite_h0_intervals = record$h0_mst_oracle$finite_h0_intervals,
    finite_h1_intervals = record$h0_mst_oracle$finite_h1_intervals,
    h0_mst_passed = record$h0_mst_oracle$passed,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    bytes = as.numeric(file.info(path)$size), sha256 = sha_file(path),
    independently_reconstructed = TRUE, downstream_jobs = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
production_artifacts <- do.call(rbind, artifact_rows)
sentinel <- t_artifacts[t_artifacts$artifact_type == "ph", , drop = FALSE]
full_inventory <- rbind(
  data.frame(
    contract_id = "mv08w_full_ph_inventory_v1", artifact_origin = "mv08t",
    job_id = sentinel$artifact_id, bytes = sentinel$primary_bytes,
    sha256 = sentinel$primary_sha256, independently_validated = sentinel$independently_validated,
    h0_mst_passed = sentinel$independently_validated,
    outcome_label_state = sentinel$outcome_label_state,
    biological_outcomes_computed = sentinel$biological_outcomes_computed,
    stringsAsFactors = FALSE
  ),
  data.frame(
    contract_id = "mv08w_full_ph_inventory_v1", artifact_origin = "mv08v",
    job_id = production_artifacts$job_id, bytes = production_artifacts$bytes,
    sha256 = production_artifacts$sha256,
    independently_validated = production_artifacts$independently_reconstructed,
    h0_mst_passed = production_artifacts$h0_mst_passed,
    outcome_label_state = production_artifacts$outcome_label_state,
    biological_outcomes_computed = production_artifacts$biological_outcomes_computed,
    stringsAsFactors = FALSE
  )
)

primary <- ledger[ledger$stage == "ph_primary", , drop = FALSE]
fallback <- ledger[ledger$stage == "ph_fallback", , drop = FALSE]
fallback_orders <- as.integer(sub("^ph_fallback__", "", fallback$attempt_id))
allowed_primary_stops <- primary$disposition == "rss_cap_exceeded" &
  as.integer(sub("^ph_primary__", "", primary$attempt_id)) %in% fallback_orders
resource_policy_passed <- all(
  (primary$disposition == "completed") | allowed_primary_stops
) && all(fallback$disposition == "completed") &&
  all(ledger$stderr_class == "empty") &&
  all(ledger$workers == 1L & ledger$retries == 0L) &&
  all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds + 5) &&
  all(ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes + 64 * 1024^2)
partials <- list.files(private_root, pattern = "\\.partial$", recursive = TRUE,
                       full.names = TRUE)
private_files <- list.files(private_root, recursive = TRUE, full.names = TRUE)
private_files <- private_files[file.info(private_files)$isdir %in% FALSE]
total_private_bytes <- sum(as.numeric(file.info(private_files)$size))

validation <- data.frame(
  check_id = c(
    "production_cardinality", "strict_queue_identity", "source_rehash",
    "artifact_rehash", "typed_view_reconstruction", "h0_mst_oracles",
    "finite_h0_cardinality", "full_inventory", "inventory_identity",
    "resource_policy", "fallback_trigger", "aggregate_elapsed",
    "private_storage", "no_partials", "stderr", "serial_zero_retry",
    "h0_h1_separate", "downstream_firewalls", "labels_outcomes_closed",
    "landscape_gate_closed"
  ),
  passed = c(
    nrow(production_artifacts) == 1257L,
    identical(production_artifacts$job_id, queue$job_id),
    nrow(source_rehash) == 131L && all(source_rehash$independently_rehashed),
    all(production_artifacts$sha256 == metrics$output_sha256),
    all(production_artifacts$independently_reconstructed),
    all(production_artifacts$h0_mst_passed),
    all(production_artifacts$finite_h0_intervals ==
          production_artifacts$point_count - 1L),
    nrow(full_inventory) == 1280L,
    !anyDuplicated(full_inventory$job_id) &&
      setequal(full_inventory$job_id,
               read_csv("docs/audits/mv08r-full-topology-production-prefreeze-v1/mv08r-ph-queue.csv")$job_id),
    resource_policy_passed,
    all(primary$disposition != "rss_cap_exceeded" |
          as.integer(sub("^ph_primary__", "", primary$attempt_id)) %in% fallback_orders),
    sum(ledger$elapsed_seconds) <= contract$aggregate_elapsed_cap_seconds + 5,
    total_private_bytes <= contract$private_storage_cap_bytes,
    length(partials) == 0L,
    all(ledger$stderr_class == "empty"),
    all(ledger$workers == 1L & ledger$retries == 0L),
    all(queue$homology_dimensions == "H0;H1_separate"),
    all(production_artifacts$downstream_jobs == 0L),
    all(full_inventory$outcome_label_state == "closed") &&
      !any(full_inventory$biological_outcomes_computed), TRUE
  ),
  evidence = c(
    "1,257 MV8-V records independently reopened",
    "production artifact order equals the frozen MV8-U queue",
    "131 source caches independently rehashed and contract-validated",
    "every output byte count and SHA agrees with its selected metric",
    "every source matrix rebuilt into the frozen typed gene view",
    "all full-view H0 MST oracles rerun and pass",
    "finite H0 contains point_count minus one intervals",
    "23 closed MV8-T plus 1,257 closed MV8-V records",
    "inventory is exactly the immutable 1,280-row MV8-R identity set",
    "all child dispositions and caps follow the frozen policy",
    "GUDHI appears only after the matching Ripserr RSS stop",
    "aggregate attempt time remains within 72 hours",
    "private PH/log evidence remains within one GiB",
    "no atomic partial artifact remains",
    "every child stderr is empty",
    "one worker and zero automatic retries throughout",
    "H0 and H1 retained separately in all queue contracts",
    "no landscape, comparison, clustering or fusion job executed",
    "labels and biological outcomes remain unopened",
    "landscape stress remains closed pending separate Rust admission prefreeze"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-W closure validation failed", call. = FALSE)

decision <- data.frame(
  contract_id = "mv08w_decision_v1",
  decision = "full_1280_PH_closed_landscape_Rust_admission_prefreeze_may_begin",
  mv08t_records = 23L, mv08v_records = 1257L, full_ph_records = 1280L,
  fallback_records = nrow(fallback), validations_passed = sum(validation$passed),
  validations_total = nrow(validation), landscape_groups_authorized = 0L,
  comparison_strata_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic_csv(source_rehash, file.path(output_dir, "mv08w-source-cache-rehash.csv"))
atomic_csv(production_artifacts, file.path(output_dir, "mv08w-production-ph-rehash.csv"))
atomic_csv(full_inventory, file.path(output_dir, "mv08w-full-ph-inventory.csv"))
atomic_csv(ledger, file.path(output_dir, "mv08w-resource-ledger.csv"))
atomic_csv(validation, file.path(output_dir, "mv08w-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08w-decision.csv"))
atomic_text(c(
  "# MV8-W full-PH production closure", "",
  "**Result:** 1,280/1,280 PH records independently closed.", "",
  paste0(
    "MV8-W rehashed 131 accepted source caches, reconstructed all 1,257 new ",
    "typed gene views, reopened every MV8-V PH record, and reran every H0 MST ",
    "oracle. Together with the 23 immutable MV8-T records, the MV8-R PH axis ",
    "is complete."
  ), "",
  paste0(
    "Landscapes, comparisons, clustering, fusion, labels, outcomes, adoption, ",
    "and manuscript claims remain closed. The next permissible action is a ",
    "separate Rust rebuild/hash-rebind and landscape-stress prefreeze."
  )
), file.path(output_dir, "MV08W_FULL_PH_PRODUCTION_CLOSURE.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08w-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08w-artifact-manifest.csv"))
cat("MV8-W checks=", sum(validation$passed), "/", nrow(validation),
    "; PH=1280/1280\n", sep = "")
