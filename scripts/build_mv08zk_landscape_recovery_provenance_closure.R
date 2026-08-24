#!/usr/bin/env Rscript

# Bind the immutable MV8-ZG 152,744-pair scientific closure to MV8-ZH/MV8-ZI
# recovery provenance and publish truthful resource semantics for the adopted
# order-164 child. This companion never recomputes landscape distances.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: build_mv08zk_landscape_recovery_provenance_closure.R",
  "<mv08zg-closure> <mv08zf-prefreeze> <mv08zh-prefreeze>",
  "<mv08zi-prefreeze> <mv08zj-prefreeze> <production-public> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zg_root <- normalizePath(args[[1L]], mustWork = TRUE)
zf_root <- normalizePath(args[[2L]], mustWork = TRUE)
zh_root <- normalizePath(args[[3L]], mustWork = TRUE)
zi_root <- normalizePath(args[[4L]], mustWork = TRUE)
zj_root <- normalizePath(args[[5L]], mustWork = TRUE)
public_root <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZK output", call. = FALSE)

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
truth <- function(value) {
  if (is.logical(value)) return(value)
  tolower(as.character(value)) %in% c("true", "t", "1", "yes")
}
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
verify_manifest <- function(root, name) {
  path <- file.path(root, name)
  manifest <- read_csv(path)
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(files, sha_file, character(1L)) == manifest$sha256)) {
    stop("MV8-ZK prerequisite manifest drift: ", name, call. = FALSE)
  }
  data.frame(
    stage = sub("-artifact-manifest[.]csv$", "", name),
    artifacts = nrow(manifest), public_bytes = sum(as.numeric(manifest$bytes)),
    manifest_sha256 = sha_file(path), stringsAsFactors = FALSE
  )
}

audit_chain <- rbind(
  verify_manifest(zg_root, "mv08zg-artifact-manifest.csv"),
  verify_manifest(zf_root, "mv08zf-artifact-manifest.csv"),
  verify_manifest(zh_root, "mv08zh-artifact-manifest.csv"),
  verify_manifest(zi_root, "mv08zi-artifact-manifest.csv"),
  verify_manifest(zj_root, "mv08zj-artifact-manifest.csv")
)
zg_validation <- read_csv(file.path(zg_root, "mv08zg-validation.csv"))
zg_resource <- read_csv(file.path(zg_root, "mv08zg-resource-summary.csv"))
zg_decision <- read_csv(file.path(zg_root, "mv08zg-decision.csv"))
zf_contract <- read_csv(file.path(zf_root, "mv08zf-contract.csv"))
zf_queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
zh_orphan <- read_csv(file.path(zh_root, "mv08zh-orphan-binding.csv"))
zh_policy <- read_csv(file.path(zh_root, "mv08zh-recovery-policy.csv"))
zh_decision <- read_csv(file.path(zh_root, "mv08zh-decision.csv"))
zi_binding <- read_csv(file.path(zi_root, "mv08zi-private-pair-binding.csv"))
zi_decision <- read_csv(file.path(zi_root, "mv08zi-decision.csv"))
zj_contract <- read_csv(file.path(zj_root, "mv08zj-contract.csv"))
zj_validation <- read_csv(file.path(zj_root, "mv08zj-validation.csv"))
zj_implementation <- read_csv(file.path(zj_root, "mv08zj-implementation-bindings.csv"))
ledger <- read_csv(file.path(public_root, "mv08zf-resource-ledger.csv"))
completed <- read_csv(file.path(public_root, "mv08zf-chunk-completions.csv"))
progress <- read_csv(file.path(public_root, "mv08zf-progress.csv"))

adopted_order <- 164L
adopted_ledger <- ledger[ledger$production_order == adopted_order, , drop = FALSE]
adopted_completion <- completed[
  completed$production_order == adopted_order, , drop = FALSE
]
measured_ledger <- ledger[ledger$production_order != adopted_order, , drop = FALSE]

group_orders <- sort(unique(as.integer(ledger$group_order)))
group_resource <- do.call(rbind, lapply(group_orders, function(group_order) {
  rows <- ledger[ledger$group_order == group_order, , drop = FALSE]
  measured <- rows[rows$production_order != adopted_order, , drop = FALSE]
  upper <- rows[rows$production_order == adopted_order, , drop = FALSE]
  data.frame(
    contract_id = "mv08zk_group_resource_v1",
    group_order = group_order,
    chunks = nrow(rows), pairs = sum(as.integer(rows$pair_count)),
    measured_chunks = nrow(measured), upper_bound_chunks = nrow(upper),
    measured_child_seconds = sum(as.numeric(measured$elapsed_seconds)),
    adopted_elapsed_upper_bound_seconds = sum(as.numeric(upper$elapsed_seconds)),
    aggregate_child_seconds_upper_bound = sum(as.numeric(rows$elapsed_seconds)),
    measured_peak_rss_bytes = if (nrow(measured)) {
      max(as.numeric(measured$peak_process_tree_rss_bytes))
    } else 0,
    adopted_rss_upper_bound_bytes = if (nrow(upper)) {
      max(as.numeric(upper$peak_process_tree_rss_bytes))
    } else 0,
    overall_peak_rss_upper_bound_bytes = max(as.numeric(rows$peak_process_tree_rss_bytes)),
    telemetry_semantics = if (nrow(upper)) {
      "measured_except_order_164_conservative_upper_bound"
    } else "measured",
    stringsAsFactors = FALSE
  )
}))

resource <- data.frame(
  contract_id = "mv08zk_resource_interpretation_v1",
  chunks = nrow(ledger), pairs = sum(as.integer(ledger$pair_count)),
  measured_chunks = nrow(measured_ledger), upper_bound_chunks = nrow(adopted_ledger),
  measured_child_seconds_excluding_order_164 =
    sum(as.numeric(measured_ledger$elapsed_seconds)),
  order_164_elapsed_upper_bound_seconds =
    as.numeric(adopted_ledger$elapsed_seconds),
  aggregate_child_seconds_upper_bound = sum(as.numeric(ledger$elapsed_seconds)),
  measured_peak_rss_bytes_excluding_order_164 =
    max(as.numeric(measured_ledger$peak_process_tree_rss_bytes)),
  order_164_rss_upper_bound_bytes =
    as.numeric(adopted_ledger$peak_process_tree_rss_bytes),
  overall_peak_rss_upper_bound_bytes =
    max(as.numeric(ledger$peak_process_tree_rss_bytes)),
  private_bytes = zg_resource$private_bytes,
  aggregate_elapsed_cap_seconds = zf_contract$aggregate_elapsed_cap_seconds,
  private_storage_cap_bytes = zf_contract$private_storage_cap_bytes,
  upper_bounds_are_measurements = FALSE,
  all_caps_passed_conservatively =
    sum(as.numeric(ledger$elapsed_seconds)) <= zf_contract$aggregate_elapsed_cap_seconds &&
    max(as.numeric(ledger$peak_process_tree_rss_bytes)) <= zf_contract$child_rss_cap_bytes &&
    zg_resource$private_bytes <= zf_contract$private_storage_cap_bytes,
  workers = 1L, retries = 0L, stringsAsFactors = FALSE
)

downstream_columns <- c(
  "comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
  "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs"
)
current_hashes <- vapply(zj_implementation$file, sha_file, character(1L))
validation <- data.frame(
  check_id = c(
    "audit_chain", "MV8_ZG_validation", "MV8_ZG_pair_closure",
    "MV8_ZJ_validation", "implementation_bindings", "production_cardinality",
    "terminal_progress", "strict_order", "adopted_order_unique",
    "launcher_recovery_binding", "pair_recovery_binding",
    "ledger_distance_binding", "ledger_status_binding",
    "completion_distance_binding", "completion_status_binding",
    "upper_bound_semantics", "measured_receipt_count", "conservative_elapsed_cap",
    "conservative_rss_cap", "private_storage_cap", "one_worker_zero_retry",
    "zero_recomputation", "execution_head", "downstream_firewall",
    "aggregate_only_schema"
  ),
  passed = c(
    nrow(audit_chain) == 5L, all(truth(zg_validation$passed)),
    zg_decision$landscape_pairs == 152744L &&
      truth(zg_decision$landscape_production_closed),
    all(truth(zj_validation$passed)),
    identical(unname(current_hashes), zj_implementation$sha256),
    nrow(ledger) == 628L && nrow(completed) == 628L &&
      nrow(zf_queue) == 628L && sum(as.integer(ledger$pair_count)) == 152744L,
    progress$state == "landscape_production_complete_closure_pending" &&
      progress$completed_chunks == 628L && progress$completed_pairs == 152744L,
    identical(as.integer(ledger$production_order), 1:628) &&
      identical(as.integer(completed$production_order), 1:628),
    nrow(adopted_ledger) == 1L && nrow(adopted_completion) == 1L &&
      adopted_ledger$production_order == 164L,
    zh_orphan$production_order == 164L && zh_decision$orphan_production_order == 164L &&
      zh_policy$orphan_order == 164L,
    zi_binding$production_order == 164L && zi_decision$orphan_production_order == 164L &&
      zi_binding$observed_pair_subset_sha256 == zh_orphan$pair_subset_sha256,
    adopted_ledger$distances_sha256 == zh_orphan$distance_sha256,
    adopted_ledger$status_sha256 == zh_orphan$status_sha256,
    adopted_completion$distances_sha256 == zh_orphan$distance_sha256,
    adopted_completion$status_sha256 == zh_orphan$status_sha256,
    truth(zh_orphan$upper_bound_not_measurement) &&
      !truth(zj_contract$upper_bounds_are_measurements) &&
      !truth(resource$upper_bounds_are_measurements) &&
      adopted_ledger$elapsed_seconds == zh_orphan$parent_elapsed_upper_bound_seconds &&
      adopted_ledger$peak_process_tree_rss_bytes == zh_orphan$parent_rss_upper_bound_bytes,
    nrow(measured_ledger) == 627L && nrow(adopted_ledger) == 1L,
    resource$aggregate_child_seconds_upper_bound <=
      resource$aggregate_elapsed_cap_seconds,
    resource$overall_peak_rss_upper_bound_bytes <= zf_contract$child_rss_cap_bytes,
    resource$private_bytes <= resource$private_storage_cap_bytes,
    all(ledger$workers == 1L) && all(ledger$retries == 0L) &&
      progress$workers == 1L && progress$retries == 0L,
    !truth(zh_policy$landscape_recomputation) &&
      !truth(zh_decision$landscape_recomputation_authorized) &&
      !truth(zi_decision$landscape_recomputation_authorized),
    length(unique(ledger$execution_head)) == 1L &&
      unique(ledger$execution_head) == progress$execution_head,
    all(as.numeric(unlist(progress[downstream_columns], use.names = FALSE)) == 0) &&
      progress$outcome_label_state == "closed" &&
      !truth(progress$biological_outcomes_computed),
    !any(c("unit_id", "sample_id", "accession", "donor_id", "source_path",
           "private_path", "private_file") %in%
         c(names(audit_chain), names(group_resource), names(resource)))
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZK recovery-provenance closure failed: ",
  paste(validation$check_id[!validation$passed], collapse = ", "),
  call. = FALSE
)

decision <- data.frame(
  contract_id = "mv08zk_decision_v1",
  decision = "bind_recovery_provenance_to_full_MV8_ZG_closure",
  landscape_groups = 28L, landscape_chunks = nrow(ledger),
  landscape_pairs = sum(as.integer(ledger$pair_count)),
  measured_resource_chunks = nrow(measured_ledger),
  conservative_upper_bound_chunks = nrow(adopted_ledger),
  order_164_upper_bounds_not_measurements = TRUE,
  landscape_recomputation_records = 0L, retries = progress$retries,
  full_landscape_closure_bound = TRUE,
  canonical_resource_interpretation_published = TRUE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L,
  next_gate = "separate_label_closed_comparison_prefreeze",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(audit_chain, file.path(output_dir, "mv08zk-audit-chain.csv"))
atomic_csv(group_resource, file.path(output_dir, "mv08zk-group-resource.csv"))
atomic_csv(resource, file.path(output_dir, "mv08zk-resource-interpretation.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zk-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zk-decision.csv"))
atomic_text(c(
  "# MV8-ZK landscape recovery-provenance closure", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " checks pass; recovery provenance is bound to MV8-ZG's 152,744/152,744-pair closure."),
  "",
  paste0(
    "MV8-ZG remains the immutable independent scientific closure. MV8-ZK is its ",
    "mandatory resource-provenance companion: 627 chunk receipts contain measured ",
    "parent telemetry, while order 164 uses the prospectively frozen 3,600-second ",
    "and 4-GiB caps only as conservative upper bounds, not measurements."
  ),
  "",
  paste0(
    "All aggregate resource caps pass under the conservative interpretation. ",
    "No landscape was recomputed, no retry occurred, no private identity is ",
    "published, and comparisons plus every biological downstream stage remain closed."
  )
), file.path(output_dir, "MV08ZK_LANDSCAPE_RECOVERY_PROVENANCE_CLOSURE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zk-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zk-artifact-manifest.csv"))
cat("MV8-ZK recovery provenance passed ", sum(validation$passed), "/",
    nrow(validation), "; landscapes=152744/152744; retries=0\n", sep = "")
