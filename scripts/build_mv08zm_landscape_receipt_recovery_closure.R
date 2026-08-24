#!/usr/bin/env Rscript

# Bind the MV8-ZL order-280 receipt-publication recovery to completed MV8-ZG/ZK
# closure. Outputs are aggregate-only and no landscape is recomputed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: build_mv08zm_landscape_receipt_recovery_closure.R",
  "<mv08zg-root> <mv08zk-root> <mv08zl-root> <mv08zf-root>",
  "<production-public> <production-private> <output-dir>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zg_root <- normalizePath(args[[1L]], mustWork = TRUE)
zk_root <- normalizePath(args[[2L]], mustWork = TRUE)
zl_root <- normalizePath(args[[3L]], mustWork = TRUE)
zf_root <- normalizePath(args[[4L]], mustWork = TRUE)
public_root <- normalizePath(args[[5L]], mustWork = TRUE)
private_root <- normalizePath(args[[6L]], mustWork = TRUE)
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-ZM output", call. = FALSE)

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
    stop("MV8-ZM prerequisite manifest drift: ", name, call. = FALSE)
  }
  data.frame(stage = sub("-artifact-manifest[.]csv$", "", name),
             artifacts = nrow(manifest), manifest_sha256 = sha_file(path),
             stringsAsFactors = FALSE)
}

audit_chain <- rbind(
  verify_manifest(zg_root, "mv08zg-artifact-manifest.csv"),
  verify_manifest(zk_root, "mv08zk-artifact-manifest.csv"),
  verify_manifest(zl_root, "mv08zl-artifact-manifest.csv"),
  verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")
)
zg_validation <- read_csv(file.path(zg_root, "mv08zg-validation.csv"))
zg_decision <- read_csv(file.path(zg_root, "mv08zg-decision.csv"))
zk_validation <- read_csv(file.path(zk_root, "mv08zk-validation.csv"))
zk_decision <- read_csv(file.path(zk_root, "mv08zk-decision.csv"))
zl_validation <- read_csv(file.path(zl_root, "mv08zl-validation.csv"))
zl_snapshot <- read_csv(file.path(zl_root, "mv08zl-stopped-snapshot.csv"))
zl_binding <- read_csv(file.path(zl_root, "mv08zl-order-280-binding.csv"))
zl_decision <- read_csv(file.path(zl_root, "mv08zl-decision.csv"))
ledger <- read_csv(file.path(public_root, "mv08zf-resource-ledger.csv"))
completed <- read_csv(file.path(public_root, "mv08zf-chunk-completions.csv"))
progress <- read_csv(file.path(public_root, "mv08zf-progress.csv"))
preserved_prefix <- file.path(
  private_root, "recovery", "mv08zl", "mv08zf-chunk-completions-prefix-279.csv"
)
order_row <- ledger[ledger$production_order == 280L, , drop = FALSE]
completion_row <- completed[completed$production_order == 280L, , drop = FALSE]
partials <- c(
  list.files(public_root, pattern = "[.]partial$", full.names = TRUE),
  list.files(private_root, pattern = "(__partial__|[.]partial$)", recursive = TRUE,
             full.names = TRUE, all.files = TRUE)
)
downstream_columns <- c(
  "comparison_jobs", "clustering_jobs", "fusion_jobs", "label_jobs",
  "outcome_jobs", "adoption_jobs", "manuscript_claim_jobs"
)

summary <- data.frame(
  contract_id = "mv08zm_receipt_recovery_summary_v1",
  recovery_order = 280L, preserved_prefix_rows = 279L,
  final_completion_rows = nrow(completed), final_ledger_rows = nrow(ledger),
  measured_elapsed_seconds = order_row$elapsed_seconds,
  measured_peak_process_tree_rss_bytes = order_row$peak_process_tree_rss_bytes,
  receipt_reconstructed = FALSE, ledger_rewritten = FALSE,
  landscape_recomputations = 0L, retries = order_row$retries,
  WSL_only_publication = TRUE, stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c(
    "audit_chain", "MV8_ZG_validation", "MV8_ZG_full_pair_closure",
    "MV8_ZK_validation", "MV8_ZK_recovery_semantics", "MV8_ZL_validation",
    "MV8_ZL_decision", "terminal_progress", "full_cardinality", "strict_order",
    "order_280_unique", "order_280_measured", "order_280_distance_binding",
    "order_280_status_binding", "preserved_prefix", "zero_partials",
    "one_worker", "zero_retry", "zero_recomputation", "ledger_not_rewritten",
    "execution_head", "downstream_firewall", "labels_outcomes_closed",
    "aggregate_only_schema", "next_gate_separate"
  ),
  passed = c(
    nrow(audit_chain) == 4L, all(truth(zg_validation$passed)),
    zg_decision$landscape_pairs == 152744L && truth(zg_decision$landscape_production_closed),
    all(truth(zk_validation$passed)), truth(zk_decision$full_landscape_closure_bound) &&
      truth(zk_decision$order_164_upper_bounds_not_measurements),
    all(truth(zl_validation$passed)), truth(zl_decision$receipt_promotion_authorized) &&
      zl_decision$adopted_production_order == 280L,
    progress$state == "landscape_production_complete_closure_pending" &&
      progress$completed_chunks == 628L && progress$completed_pairs == 152744L,
    nrow(ledger) == 628L && nrow(completed) == 628L,
    identical(as.integer(ledger$production_order), 1:628) &&
      identical(as.integer(completed$production_order), 1:628),
    nrow(order_row) == 1L && nrow(completion_row) == 1L,
    truth(zl_binding$telemetry_is_measured) &&
      order_row$elapsed_seconds == zl_binding$elapsed_seconds &&
      order_row$peak_process_tree_rss_bytes == zl_binding$peak_process_tree_rss_bytes,
    order_row$distances_sha256 == zl_binding$distances_sha256 &&
      completion_row$distances_sha256 == zl_binding$distances_sha256,
    order_row$status_sha256 == zl_binding$status_sha256 &&
      completion_row$status_sha256 == zl_binding$status_sha256,
    file.exists(preserved_prefix) &&
      sha_file(preserved_prefix) == zl_snapshot$completion_sha256 &&
      nrow(read_csv(preserved_prefix)) == 279L,
    length(partials) == 0L, all(ledger$workers == 1L) && all(completed$workers == 1L),
    all(ledger$retries == 0L) && all(completed$retries == 0L),
    summary$landscape_recomputations == 0L, !truth(summary$ledger_rewritten),
    length(unique(ledger$execution_head)) == 1L &&
      unique(ledger$execution_head) == progress$execution_head,
    all(as.numeric(unlist(progress[downstream_columns], use.names = FALSE)) == 0),
    progress$outcome_label_state == "closed" && !truth(progress$biological_outcomes_computed),
    !any(c("unit_id", "sample_id", "accession", "donor_id", "source_path",
           "private_path", "private_file") %in% c(names(audit_chain), names(summary))),
    zl_decision$companion_MV8_ZM_required &&
      zk_decision$next_gate == "separate_label_closed_comparison_prefreeze"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop(
  "MV8-ZM recovery closure failed: ",
  paste(validation$check_id[!validation$passed], collapse = ", "), call. = FALSE
)
decision <- data.frame(
  contract_id = "mv08zm_decision_v1",
  decision = "bind_order_280_receipt_recovery_to_full_landscape_closure",
  validations_passed = sum(validation$passed), validations_total = nrow(validation),
  landscape_pairs = 152744L, recovered_receipts = 1L,
  landscape_recomputation_records = 0L, retry_records = 0L,
  full_landscape_closure_bound = TRUE,
  comparison_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L,
  next_gate = "separate_label_closed_comparison_prefreeze",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE)
atomic_csv(audit_chain, file.path(output_dir, "mv08zm-audit-chain.csv"))
atomic_csv(summary, file.path(output_dir, "mv08zm-recovery-summary.csv"))
atomic_csv(validation, file.path(output_dir, "mv08zm-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08zm-decision.csv"))
atomic_text(c(
  "# MV8-ZM landscape receipt-recovery closure", "",
  paste0("**Result:** ", sum(validation$passed), "/", nrow(validation),
         " checks pass; order-280 publication recovery is bound to MV8-ZG/ZK full closure."),
  "",
  "Order 280 retained its original measured resource row and exact private scientific hashes. Its already-written completion receipt was promoted without reconstruction, landscape recomputation, retry, or ledger rewrite, and the prior 279-row prefix remains privately preserved.",
  "",
  "The final 628-chunk/152,744-pair corpus remains label-closed. Comparisons and every biological downstream stage require a separate prospective gate."
), file.path(output_dir, "MV08ZM_LANDSCAPE_RECEIPT_RECOVERY_CLOSURE.md"))
artifacts <- sort(setdiff(basename(list.files(output_dir, full.names = TRUE)),
                          "mv08zm-artifact-manifest.csv"))
manifest <- data.frame(
  artifact = artifacts,
  bytes = as.numeric(file.info(file.path(output_dir, artifacts))$size),
  sha256 = vapply(file.path(output_dir, artifacts), sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08zm-artifact-manifest.csv"))
cat("MV8-ZM receipt-recovery provenance passed ", sum(validation$passed), "/",
    nrow(validation), "; landscapes=152744/152744; retries=0\n", sep = "")
