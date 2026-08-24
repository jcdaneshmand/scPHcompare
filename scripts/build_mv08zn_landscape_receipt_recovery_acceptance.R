#!/usr/bin/env Rscript

# Accept the completed MV8-ZL receipt recovery after its final exact-double
# equality guard rejected a harmless CSV round-trip difference. This builder is
# read-only with respect to production and authorizes no landscape computation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv08zn_landscape_receipt_recovery_acceptance.R",
  "<mv08zl-root> <mv08zf-root> <production-private> <production-public>",
  "<output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zl_root <- normalizePath(args[[1L]], mustWork = TRUE)
zf_root <- normalizePath(args[[2L]], mustWork = TRUE)
private_root <- normalizePath(args[[3L]], mustWork = TRUE)
public_root <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZN output", call. = FALSE)

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
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(files, sha_file, character(1L)) == manifest$sha256)) {
    stop("MV8-ZN prerequisite manifest drift: ", name, call. = FALSE)
  }
  data.frame(stage = sub("-artifact-manifest[.]csv$", "", name),
             artifacts = nrow(manifest),
             manifest_sha256 = sha_file(file.path(root, name)),
             stringsAsFactors = FALSE)
}

audit_chain <- rbind(
  verify_manifest(zl_root, "mv08zl-artifact-manifest.csv"),
  verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")
)
zl_snapshot <- read_csv(file.path(zl_root, "mv08zl-stopped-snapshot.csv"))
zl_binding <- read_csv(file.path(zl_root, "mv08zl-order-280-binding.csv"))
zl_policy <- read_csv(file.path(zl_root, "mv08zl-recovery-policy.csv"))
zl_implementation <- read_csv(file.path(zl_root, "mv08zl-implementation-bindings.csv"))
zl_decision <- read_csv(file.path(zl_root, "mv08zl-decision.csv"))
zf_contract <- read_csv(file.path(zf_root, "mv08zf-contract.csv"))
ledger <- read_csv(file.path(public_root, "mv08zf-resource-ledger.csv"))
completed <- read_csv(file.path(public_root, "mv08zf-chunk-completions.csv"))
progress <- read_csv(file.path(public_root, "mv08zf-progress.csv"))
order_row <- ledger[ledger$production_order == 280L, , drop = FALSE]
completion_row <- completed[completed$production_order == 280L, , drop = FALSE]
preserved_prefix <- file.path(
  private_root, "recovery", "mv08zl", "mv08zf-chunk-completions-prefix-279.csv"
)
private_chunk <- file.path(private_root, "production", "group_10", "chunk_001")
private_paths <- c(
  distance = file.path(private_chunk, "distances.csv"),
  status = file.path(private_chunk, "status.csv"),
  stdout = file.path(private_root, "logs", "chunk_0280.stdout"),
  stderr = file.path(private_root, "logs", "chunk_0280.stderr")
)
partials <- c(
  list.files(public_root, pattern = "[.]partial$", full.names = TRUE),
  list.files(private_root, pattern = "(__partial__|[.]partial$)", recursive = TRUE,
             full.names = TRUE, all.files = TRUE)
)
runner_lines <- suppressWarnings(system2(
  "pgrep", c("-f", "[r]un_mv08zf_full_landscape_production[.]R"),
  stdout = TRUE, stderr = FALSE
))
executor_binding <- zl_implementation[
  zl_implementation$role == "receipt_recovery_executor", , drop = FALSE
]
aggregate_ledger <- sum(as.numeric(ledger$elapsed_seconds))
aggregate_progress <- as.numeric(progress$aggregate_child_seconds)
aggregate_difference <- aggregate_progress - aggregate_ledger
tolerance <- 1e-9

stop_record <- data.frame(
  contract_id = "mv08zn_validation_stop_v1",
  executor_exit_status = 1L,
  stop_stage = "post_publication_final_validation",
  failure_class = "csv_roundtrip_exact_double_equality_only",
  completion_promotion_finished = TRUE,
  progress_refresh_finished = TRUE,
  aggregate_progress_seconds = aggregate_progress,
  aggregate_ledger_seconds = aggregate_ledger,
  absolute_difference_seconds = abs(aggregate_difference),
  acceptance_tolerance_seconds = tolerance,
  tolerance_passed = abs(aggregate_difference) <= tolerance,
  scientific_failure = FALSE, resource_failure = FALSE,
  publication_failure = FALSE, stringsAsFactors = FALSE
)

downstream_columns <- c(
  "comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
  "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs"
)
validation <- data.frame(
  check_id = c(
    "audit_chain", "MV8_ZL_authorization", "MV8_ZL_executor_binding",
    "ledger_280", "completion_280", "progress_280", "strict_order",
    "ledger_completion_parity", "order_280_binding", "order_280_private_hashes",
    "preserved_prefix", "zero_partials", "runner_absent",
    "exact_equality_reproduced", "tolerance_passed", "difference_bounded",
    "measured_resource_preserved", "child_caps", "aggregate_cap",
    "private_cap", "one_worker", "zero_retry", "zero_recomputation",
    "no_ledger_rewrite", "WSL_only_monitoring", "downstream_firewall",
    "labels_outcomes_closed", "resume_exact_281"
  ),
  passed = c(
    nrow(audit_chain) == 2L,
    truth(zl_decision$receipt_promotion_authorized) &&
      zl_decision$adopted_production_order == 280L,
    nrow(executor_binding) == 1L && file.exists(executor_binding$file) &&
      sha_file(executor_binding$file) == executor_binding$sha256,
    nrow(ledger) == 280L && all(ledger$disposition == "completed") &&
      all(ledger$exit_status == 0L),
    nrow(completed) == 280L,
    nrow(progress) == 1L && progress$state == "running" &&
      progress$completed_chunks == 280L && progress$completed_pairs == 68884L,
    identical(as.integer(ledger$production_order), 1:280) &&
      identical(as.integer(completed$production_order), 1:280),
    all(ledger$distances_sha256 == completed$distances_sha256) &&
      all(ledger$status_sha256 == completed$status_sha256),
    nrow(order_row) == 1L && nrow(completion_row) == 1L &&
      order_row$distances_sha256 == zl_binding$distances_sha256 &&
      completion_row$status_sha256 == zl_binding$status_sha256,
    all(file.exists(private_paths)) &&
      sha_file(private_paths[["distance"]]) == zl_binding$distances_sha256 &&
      sha_file(private_paths[["status"]]) == zl_binding$status_sha256 &&
      sha_file(private_paths[["stdout"]]) == zl_binding$stdout_sha256 &&
      sha_file(private_paths[["stderr"]]) == zl_binding$stderr_sha256,
    file.exists(preserved_prefix) &&
      sha_file(preserved_prefix) == zl_snapshot$completion_sha256 &&
      nrow(read_csv(preserved_prefix)) == 279L,
    length(partials) == 0L, length(runner_lines) == 0L,
    aggregate_progress != aggregate_ledger,
    truth(stop_record$tolerance_passed),
    abs(aggregate_difference) <= tolerance && abs(aggregate_difference) < 1e-9,
    truth(zl_binding$telemetry_is_measured) &&
      order_row$elapsed_seconds == zl_binding$elapsed_seconds &&
      order_row$peak_process_tree_rss_bytes == zl_binding$peak_process_tree_rss_bytes,
    all(ledger$elapsed_seconds <= ledger$elapsed_cap_seconds) &&
      all(ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes),
    aggregate_ledger <= zf_contract$aggregate_elapsed_cap_seconds,
    progress$private_bytes <= zf_contract$private_storage_cap_bytes,
    all(ledger$workers == 1L) && all(completed$workers == 1L),
    all(ledger$retries == 0L) && all(completed$retries == 0L),
    !truth(zl_policy$landscape_recomputation), !truth(zl_policy$ledger_rewrite),
    truth(zl_policy$WSL_only_publication_and_monitoring),
    all(as.numeric(unlist(progress[downstream_columns], use.names = FALSE)) == 0),
    progress$outcome_label_state == "closed" &&
      !truth(progress$biological_outcomes_computed),
    zl_decision$resume_at_production_order == 281L
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZN acceptance failed: ",
  paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE
)

decision <- data.frame(
  contract_id = "mv08zn_decision_v1",
  decision = "accept_completed_order_280_receipt_recovery_and_resume_at_281",
  validations_passed = sum(validation$passed), validations_total = nrow(validation),
  receipt_recovery_complete = TRUE, executor_rerun_authorized = FALSE,
  resume_at_production_order = 281L,
  floating_tolerance_seconds = tolerance,
  landscape_recomputation_records = 0L, retry_records = 0L,
  scientific_contract_changed = FALSE, resource_contract_changed = FALSE,
  WSL_only_monitoring_required = TRUE, companion_MV8_ZM_required = TRUE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(audit_chain, file.path(output_dir, "mv08zn-audit-chain.csv"))
atomic_csv(stop_record, file.path(output_dir, "mv08zn-validation-stop.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zn-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zn-decision.csv"))
atomic_text(c(
  "# MV8-ZN landscape receipt-recovery acceptance", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " checks pass; the completed order-280 receipt recovery is accepted."),
  "",
  paste0(
    "The committed MV8-ZL executor completed both authorized publications, then ",
    "exited 1 only because CSV round-trip parsing changed the aggregate double by ",
    "", format(abs(aggregate_difference), scientific = TRUE), " seconds and its ",
    "final guard required exact equality. The frozen 1e-9-second acceptance ",
    "tolerance passes."
  ),
  "",
  "Ledger, completion, and progress now close a strict 280-row prefix; the original 279-row completion is privately preserved; zero partials, retries, recomputations, or downstream jobs exist. Do not rerun the executor. Resume unchanged MV8-ZF only at order 281 using WSL-only monitoring."
), file.path(output_dir, "MV08ZN_LANDSCAPE_RECEIPT_RECOVERY_ACCEPTANCE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zn-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zn-artifact-manifest.csv"))
cat("MV8-ZN accepted completed receipt recovery 28/28; resume_at=281\n")
