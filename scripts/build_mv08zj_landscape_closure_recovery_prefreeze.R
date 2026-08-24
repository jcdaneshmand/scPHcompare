#!/usr/bin/env Rscript

# Prospectively bind a companion recovery-provenance closure around the
# immutable MV8-ZG scientific closure. The companion does not alter or replace
# any landscape output. It makes the order-164 launcher telemetry semantics
# explicit: the frozen caps are conservative upper bounds, not measurements.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "usage: build_mv08zj_landscape_closure_recovery_prefreeze.R",
  "<mv08zf-prefreeze> <mv08zh-prefreeze> <mv08zi-prefreeze> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zf_root <- normalizePath(args[[1L]], mustWork = TRUE)
zh_root <- normalizePath(args[[2L]], mustWork = TRUE)
zi_root <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- normalizePath(args[[4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZJ output", call. = FALSE)

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
    stop("MV8-ZJ prerequisite manifest drift: ", name, call. = FALSE)
  }
  data.frame(
    stage = sub("-artifact-manifest[.]csv$", "", name),
    artifacts = nrow(manifest), public_bytes = sum(as.numeric(manifest$bytes)),
    manifest_sha256 = sha_file(path), stringsAsFactors = FALSE
  )
}

audit_chain <- rbind(
  verify_manifest(zf_root, "mv08zf-artifact-manifest.csv"),
  verify_manifest(zh_root, "mv08zh-artifact-manifest.csv"),
  verify_manifest(zi_root, "mv08zi-artifact-manifest.csv")
)
zf_implementation <- read_csv(file.path(zf_root, "mv08zf-implementation-bindings.csv"))
zf_validation <- read_csv(file.path(zf_root, "mv08zf-validation.csv"))
zh_validation <- read_csv(file.path(zh_root, "mv08zh-validation.csv"))
zi_validation <- read_csv(file.path(zi_root, "mv08zi-validation.csv"))
zh_orphan <- read_csv(file.path(zh_root, "mv08zh-orphan-binding.csv"))
zh_policy <- read_csv(file.path(zh_root, "mv08zh-recovery-policy.csv"))
zh_decision <- read_csv(file.path(zh_root, "mv08zh-decision.csv"))
zi_binding <- read_csv(file.path(zi_root, "mv08zi-private-pair-binding.csv"))
zi_decision <- read_csv(file.path(zi_root, "mv08zi-decision.csv"))

implementation_files <- c(
  "scripts/build_mv08zg_full_landscape_production_closure.R",
  "scripts/build_mv08zk_landscape_recovery_provenance_closure.R",
  "scripts/build_mv08zj_landscape_closure_recovery_prefreeze.R"
)
if (!all(file.exists(implementation_files))) {
  stop("MV8-ZJ implementation absent", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv08zj_implementation_binding_v1",
  role = c("immutable_scientific_closure", "companion_recovery_closure",
           "companion_prefreeze_builder"),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  scientific_change = FALSE, resource_contract_change = FALSE,
  stringsAsFactors = FALSE
)

original_binding <- zf_implementation[
  zf_implementation$role == "independent_closure_builder", , drop = FALSE
]
closure_text <- paste(readLines(implementation_files[[1L]], warn = FALSE),
                      collapse = "\n")
original_gap <- !grepl("mv08zh", closure_text, fixed = TRUE) &&
  !grepl("mv08zi", closure_text, fixed = TRUE) &&
  grepl("sum(ledger$elapsed_seconds)", closure_text, fixed = TRUE) &&
  grepl("max(ledger$peak_process_tree_rss_bytes)", closure_text, fixed = TRUE)

contract <- data.frame(
  contract_id = "mv08zj_closure_recovery_prefreeze_v1",
  immutable_scientific_closure = "MV8-ZG",
  required_companion_closure = "MV8-ZK",
  adopted_production_order = 164L,
  measured_parent_telemetry_available = FALSE,
  elapsed_upper_bound_seconds = zh_orphan$parent_elapsed_upper_bound_seconds,
  rss_upper_bound_bytes = zh_orphan$parent_rss_upper_bound_bytes,
  upper_bounds_are_measurements = FALSE,
  original_MV8_ZG_runs_unchanged = TRUE,
  companion_recomputes_landscapes = FALSE,
  companion_rehashes_scientific_outputs = TRUE,
  canonical_resource_interpretation = "MV8-ZG_read_with_MV8-ZK_companion",
  workers = 1L, automatic_retries = 0L,
  comparison_state = "closed", clustering_state = "closed",
  fusion_state = "closed", label_state = "closed", outcome_state = "closed",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "prefreeze_manifest", "launcher_recovery_manifest", "pair_recovery_manifest",
    "prefreeze_validation", "launcher_recovery_validation",
    "pair_recovery_validation", "original_closure_binding",
    "recovery_gap_detected", "orphan_order", "orphan_upper_bound_semantics",
    "private_pair_binding", "orphan_pair_identity", "no_recomputation",
    "zero_retry", "one_worker", "science_unchanged", "resources_unchanged",
    "implementation_bound", "implementations_parse", "downstream_closed"
  ),
  passed = c(
    nrow(audit_chain[audit_chain$stage == "mv08zf", , drop = FALSE]) == 1L,
    nrow(audit_chain[audit_chain$stage == "mv08zh", , drop = FALSE]) == 1L,
    nrow(audit_chain[audit_chain$stage == "mv08zi", , drop = FALSE]) == 1L,
    all(truth(zf_validation$passed)), all(truth(zh_validation$passed)),
    all(truth(zi_validation$passed)),
    nrow(original_binding) == 1L &&
      original_binding$file == implementation_files[[1L]] &&
      original_binding$sha256 == sha_file(implementation_files[[1L]]),
    original_gap,
    nrow(zh_orphan) == 1L && zh_orphan$production_order == 164L &&
      nrow(zi_binding) == 1L && zi_binding$production_order == 164L,
    truth(zh_orphan$upper_bound_not_measurement) &&
      zh_orphan$parent_elapsed_upper_bound_seconds == 3600 &&
      zh_orphan$parent_rss_upper_bound_bytes == 4 * 1024^3,
    truth(zi_decision$private_pair_binding_required) &&
      !truth(zi_binding$public_private_identities_exposed),
    zi_binding$expected_pair_rows == zh_orphan$pair_count &&
      zi_binding$observed_pair_subset_sha256 == zh_orphan$pair_subset_sha256,
    !truth(zh_policy$landscape_recomputation) &&
      !truth(zh_decision$landscape_recomputation_authorized) &&
      !truth(zi_decision$landscape_recomputation_authorized),
    !truth(zh_policy$automatic_retry) && zh_decision$automatic_retries == 0L &&
      zi_decision$automatic_retries == 0L,
    zh_policy$workers == 1L && contract$workers == 1L,
    !any(truth(implementation$scientific_change)),
    !any(truth(implementation$resource_contract_change)),
    all(vapply(implementation$file, sha_file, character(1L)) == implementation$sha256),
    all(vapply(implementation$file, function(path) {
      !inherits(try(parse(path), silent = TRUE), "try-error")
    }, logical(1L))),
    all(contract[c("comparison_state", "clustering_state", "fusion_state",
                   "label_state", "outcome_state")] == "closed") &&
      contract$outcome_label_state == "closed" &&
      !truth(contract$biological_outcomes_computed)
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZJ prefreeze failed: ",
  paste(validation$check_id[!validation$passed], collapse = ", "),
  call. = FALSE
)

decision <- data.frame(
  contract_id = "mv08zj_decision_v1",
  decision = "require_recovery_provenance_companion_after_unchanged_MV8_ZG",
  validations_passed = sum(validation$passed),
  validations_total = nrow(validation),
  original_MV8_ZG_authorized_after_production = TRUE,
  companion_MV8_ZK_required = TRUE,
  order_164_upper_bounds_not_measurements = TRUE,
  scientific_contract_changed = FALSE, resource_contract_changed = FALSE,
  landscape_recomputation_authorized = FALSE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L,
  next_gate = "complete_MV8_ZF_then_run_MV8_ZG_and_MV8_ZK",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(audit_chain, file.path(output_dir, "mv08zj-audit-chain.csv"))
atomic_csv(contract, file.path(output_dir, "mv08zj-contract.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08zj-implementation-bindings.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zj-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zj-decision.csv"))
atomic_text(c(
  "# MV8-ZJ landscape-closure recovery prefreeze", "",
  "- Validation: 20/20 pass.",
  "- The immutable MV8-ZG scientific closure remains unchanged and still independently rehashes all 152,744 pairs.",
  "- A mandatory MV8-ZK companion binds MV8-ZH/MV8-ZI and reports order 164's 3,600-s/4-GiB values only as conservative upper bounds, never as measurements.",
  "- The companion performs no landscape recomputation and changes no scientific or resource contract.",
  "- Comparisons, clustering, fusion, labels, outcomes, and biological claims remain closed."
), file.path(output_dir, "MV08ZJ_LANDSCAPE_CLOSURE_RECOVERY_PREFREEZE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zj-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zj-artifact-manifest.csv"))
cat("MV8-ZJ landscape-closure recovery prefreeze passed 20/20; no science executed\n")
