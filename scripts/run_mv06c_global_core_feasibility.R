#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: run_mv06c_global_core_feasibility.R",
    "<repo-root> <d0-ledger.csv> <d0-cache-dir> <output-dir>"
  ), call. = FALSE)
}

repo <- normalizePath(args[[1L]], mustWork = TRUE)
ledger_path <- normalizePath(args[[2L]], mustWork = TRUE)
cache_dir <- normalizePath(args[[3L]], mustWork = TRUE)
output_dir <- args[[4L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
devtools::load_all(repo, quiet = TRUE)

file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
write_stable <- function(value, name) {
  path <- file.path(output_dir, name)
  utils::write.csv(value, path, row.names = FALSE, na = "", quote = TRUE)
  path
}
rss_now <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}

expected_ledger_sha <-
  "73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308"
if (!identical(file_sha(ledger_path), expected_ledger_sha)) {
  stop("MV6-C accepted D0 ledger identity mismatch.", call. = FALSE)
}
ledger <- utils::read.csv(
  ledger_path, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "sample_id", "seed", "normalization_cache_key", "payload_contract_id",
  "private_cache_file", "private_cache_size_bytes", "private_cache_sha256",
  "outcome_label_state", "biological_outcomes_computed"
)
if (!all(required %in% names(ledger)) || nrow(ledger) != 450L ||
    length(unique(ledger$sample_id)) != 90L ||
    !identical(sort(unique(as.integer(ledger$seed))), 20260805:20260809) ||
    any(table(ledger$seed) != 90L) ||
    anyDuplicated(paste(ledger$sample_id, ledger$seed, sep = "\r")) ||
    any(ledger$payload_contract_id != "mv05d0_sct_data_matrix_v1") ||
    any(ledger$outcome_label_state != "closed") ||
    any(as.logical(ledger$biological_outcomes_computed))) {
  stop("MV6-C D0 ledger violates its frozen axis or label boundary.",
       call. = FALSE)
}
ledger <- ledger[order(ledger$sample_id, ledger$seed, method = "radix"),,
                 drop = FALSE]

started <- proc.time()[["elapsed"]]
peak_rss <- rss_now()
common <- NULL
feature_union <- character()
verified_paths <- character(nrow(ledger))
for (index in seq_len(nrow(ledger))) {
  row <- ledger[index, , drop = FALSE]
  path <- file.path(cache_dir, row$private_cache_file)
  if (!file.exists(path) || !identical(file_sha(path), row$private_cache_sha256)) {
    stop("MV6-C D0 cache is missing or stale.", call. = FALSE)
  }
  record <- readRDS(path)
  mv05d0_validate_normalization_cache_record_v2(record)
  value <- mv05d0_sct_matrix_from_cache_v1(record)
  if (!identical(record$cache_key, row$normalization_cache_key) ||
      ncol(value) != 384L) {
    stop("MV6-C D0 cache identity or cell axis is stale.", call. = FALSE)
  }
  features <- rownames(value)
  common <- if (is.null(common)) features else intersect(common, features)
  feature_union <- union(feature_union, features)
  verified_paths[[index]] <- path
  peak_rss <- max(peak_rss, rss_now(), na.rm = TRUE)
  rm(record, value)
}
common <- sort(common, method = "radix")
feature_union <- sort(feature_union, method = "radix")
if (!length(common)) stop("MV6-C found no common feature axis.", call. = FALSE)

variances <- matrix(
  NA_real_, nrow = length(common), ncol = nrow(ledger),
  dimnames = list(common, paste(ledger$sample_id, ledger$seed, sep = "__"))
)
for (index in seq_len(nrow(ledger))) {
  path <- verified_paths[[index]]
  record <- readRDS(path)
  value <- mv05d0_sct_matrix_from_cache_v1(record)
  if (!all(common %in% rownames(value))) {
    stop("MV6-C common feature axis changed between passes.", call. = FALSE)
  }
  variances[, index] <- .mv03_row_variance(
    value[common, , drop = FALSE]
  )
  if (!identical(file_sha(path), ledger$private_cache_sha256[[index]])) {
    stop("MV6-C source changed between verification passes.", call. = FALSE)
  }
  peak_rss <- max(peak_rss, rss_now(), na.rm = TRUE)
  rm(record, value)
}

result <- mv06c_build_global_core_panel_v1(
  common, variances, as.integer(ledger$seed), panel_size = 500L
)
elapsed <- proc.time()[["elapsed"]] - started
source_bytes <- sum(as.numeric(ledger$private_cache_size_bytes))
resource <- data.frame(
  contract_id = "mv06c_inventory_resource_v1",
  elapsed_seconds = elapsed,
  peak_process_rss_bytes = peak_rss,
  unique_source_bytes = source_bytes,
  source_bytes_hashed = 2 * source_bytes,
  source_bytes_deserialized = 2 * source_bytes,
  output_bytes = NA_real_,
  elapsed_cap_seconds = 1800,
  rss_cap_bytes = 8 * 1024^3,
  elapsed_cap_pass = elapsed <= 1800,
  rss_cap_pass = is.finite(peak_rss) && peak_rss <= 8 * 1024^3,
  stringsAsFactors = FALSE
)
source_summary <- data.frame(
  contract_id = "mv06c_source_summary_v1",
  ledger_sha256 = expected_ledger_sha,
  source_set_sha256 = digest::digest(
    ledger[c("sample_id", "seed", "private_cache_sha256")],
    algo = "sha256", serialize = TRUE
  ),
  records_verified = nrow(ledger), biological_samples = 90L, seeds = 5L,
  cells_per_record = 384L,
  union_features = length(feature_union),
  common_features = length(common),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
panel <- result$panel[c(
  "panel_order", "feature_id", "gene", "category",
  "median_variance_rank", "minimum_variance", "finite_cache_count",
  "positive_cache_count", "selected_all_cache_core"
)]
panel$contract_id <- "mv06c_global_core_panel_v1"
panel$panel_sha256 <- result$panel_sha256
panel <- panel[c("contract_id", "panel_sha256", names(panel)[
  !names(panel) %in% c("contract_id", "panel_sha256")
])]
workload <- mv06c_future_workload_v1()
decision <- data.frame(
  contract_id = "mv06c_global_core_decision_v1",
  decision = result$decision,
  panel_sha256 = result$panel_sha256,
  selected_panel_size = nrow(panel),
  eligible_unique_canonical_genes =
    result$eligibility$eligible_unique_canonical_genes,
  eligibility_margin = result$eligibility$eligibility_margin,
  all_selected_present_finite_nonconstant =
    nrow(panel) == 500L && all(panel$finite_cache_count == 450L) &&
      all(panel$positive_cache_count == 450L),
  resource_caps_pass = resource$elapsed_cap_pass && resource$rss_cap_pass,
  cell_source_jobs_executed = 0L, gene_source_jobs_executed = 0L,
  pca_jobs_executed = 0L, ph_jobs_executed = 0L,
  landscape_jobs_executed = 0L, fusion_jobs_executed = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  next_allowed_scope = if (result$decision ==
    "go_bounded_matched_sct_profile") {
    "separately_prefrozen_bounded_fold_seed_profile"
  } else {
    "project_owner_contract_revision"
  },
  stringsAsFactors = FALSE
)
if (!decision$all_selected_present_finite_nonconstant ||
    !decision$resource_caps_pass) {
  decision$decision <- "stop_global_core_validation_or_resource_failure"
  decision$next_allowed_scope <- "project_owner_contract_revision"
}
contract <- data.frame(
  contract_id = "mv06c_global_core_matched_sct_feasibility_v1",
  prefreeze_commit = "b52c857",
  ledger_sha256 = expected_ledger_sha,
  variance_floor = .Machine$double.eps,
  requested_panel_size = 500L,
  technical_feature_rule = "technical_feature_rules_v1",
  harmonization_scope =
    "all_existing_samples_and_seeds_label_closed_transductive_technical",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

paths <- c(
  write_stable(contract, "mv06c-contract.csv"),
  write_stable(source_summary, "mv06c-source-summary.csv"),
  write_stable(result$eligibility, "mv06c-eligibility-summary.csv"),
  write_stable(panel, "mv06c-panel.csv"),
  write_stable(result$seed_stability, "mv06c-seed-stability.csv"),
  write_stable(workload, "mv06c-workload.csv"),
  write_stable(decision, "mv06c-decision.csv")
)
manifest <- data.frame(
  contract_id = "mv06c_artifact_manifest_v1",
  file = basename(paths), bytes = unname(file.info(paths)$size),
  sha256 = vapply(paths, file_sha, character(1L)),
  stringsAsFactors = FALSE
)
manifest_path <- write_stable(manifest, "mv06c-artifact-manifest.csv")
resource$output_bytes <- sum(file.info(c(paths, manifest_path))$size)
write_stable(resource, "mv06c-resource.csv")
message("MV6-C global-core feasibility complete: ", decision$decision)
