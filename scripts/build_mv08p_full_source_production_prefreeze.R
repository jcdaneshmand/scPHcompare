#!/usr/bin/env Rscript

# Freeze the remaining all-QC Pearson-residual source fits after MV8-O's
# bounded sentinel. This script only writes queue/decision metadata; it never
# opens a count matrix or computes PH, landscapes, clustering, or outcomes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: build_mv08p_full_source_production_prefreeze.R <output-dir>", call. = FALSE)
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-P output", call. = FALSE)
for (package in c("digest")) if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required", call. = FALSE)
dir.create(output_dir, recursive = TRUE)
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial"); utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(x, path) {
  partial <- paste0(path, ".partial"); writeLines(x, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}

source_queue <- utils::read.csv("docs/audits/mv08n-pearson-residual-migration-prefreeze-v1/mv08n-residual-source-queue.csv",
  check.names = FALSE, stringsAsFactors = FALSE)
sentinel <- utils::read.csv("docs/audits/mv08o-residual-source-sentinel-closure-v1/mv08o-source-sentinel-resource.csv",
  check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(source_queue) != 132L || sum(source_queue$authorization_state == "source_view_sentinel_authorized") != 3L ||
    nrow(sentinel) != 5L || !all(sentinel$disposition == "completed") ||
    !all(sentinel$elapsed_seconds <= sentinel$elapsed_cap_seconds & sentinel$peak_rss_bytes <= sentinel$rss_cap_bytes) ||
    !all(sentinel$all_geometry_valid)) stop("MV8-N/MV8-O prerequisite drift", call. = FALSE)
sentinel_sources <- source_queue[source_queue$authorization_state == "source_view_sentinel_authorized", , drop = FALSE]

remaining <- source_queue[source_queue$authorization_state == "closed_pending_sentinel", , drop = FALSE]
remaining <- remaining[order(remaining$fit_cells, remaining$unit_id), , drop = FALSE]
if (nrow(remaining) != 129L || any(remaining$fit_cells > max(sentinel_sources$fit_cells))) stop("MV8-P remaining-source range drift", call. = FALSE)
remaining$job_order <- seq_len(nrow(remaining))
remaining$authorization_state <- "authorized_after_mv08p_commit"
remaining$execution_policy <- "serial_one_worker_ascending_fit_cells_no_retry"
remaining$source_cache_state <- "private_cache_only"
remaining$ph_authorized <- FALSE
remaining$landscapes_authorized <- FALSE
remaining$comparisons_authorized <- FALSE
remaining$clustering_authorized <- FALSE
remaining$fusion_authorized <- FALSE
remaining$labels_authorized <- FALSE
remaining$outcomes_authorized <- FALSE

max_peak <- max(sentinel$peak_rss_bytes)
cap <- unique(sentinel$rss_cap_bytes)
contract <- data.frame(
  contract_id = "mv08p_full_source_production_prefreeze_v1",
  all_source_fits = nrow(source_queue), completed_sentinel_primary_fits = 3L,
  remaining_source_fits = nrow(remaining), remaining_internal_fits = sum(remaining$dataset_scope == "internal124"),
  remaining_external_fits = sum(remaining$dataset_scope == "external8"),
  maximum_sentinel_fit_cells = max(sentinel_sources$fit_cells), maximum_remaining_fit_cells = max(remaining$fit_cells),
  workers = 1L, retries = 0L, elapsed_cap_seconds = unique(sentinel$elapsed_cap_seconds),
  rss_cap_bytes = cap, observed_max_peak_rss_bytes = max_peak, observed_memory_margin_bytes = cap - max_peak,
  source_execution_state = "authorized_after_commit", topology_execution_state = "closed",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c("mv08n_source_universe", "mv08o_sentinel_complete", "sentinel_resource_caps", "remaining_fit_count",
    "maximum_remaining_bounded_by_sentinel", "serial_no_retry_policy", "topology_firewall", "outcome_firewall"),
  passed = c(nrow(source_queue) == 132L, nrow(sentinel) == 5L && all(sentinel$disposition == "completed"),
    all(sentinel$elapsed_seconds <= sentinel$elapsed_cap_seconds & sentinel$peak_rss_bytes <= sentinel$rss_cap_bytes),
    nrow(remaining) == 129L, max(remaining$fit_cells) <= max(sentinel_sources$fit_cells),
    all(remaining$workers == 1L & remaining$retries == 0L),
    !any(remaining$ph_authorized | remaining$landscapes_authorized | remaining$comparisons_authorized |
      remaining$clustering_authorized | remaining$fusion_authorized),
    all(remaining$outcome_label_state == "closed") && !any(remaining$biological_outcomes_computed)),
  evidence = c("132-source MV8-N universe", "three primary sentinel fits plus two deterministic repeats",
    "each sentinel attempt remained within 1,800 seconds and 12 GiB", "122 internal plus seven external remaining fits",
    "largest remaining source is no larger than the completed maximum sentinel", "one worker and zero automatic retries",
    "source construction only; PH through fusion remain closed", "labels and biological outcomes remain closed"), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-P prefreeze validation failed", call. = FALSE)
atomic_csv(contract, file.path(output_dir, "mv08p-contract.csv"))
atomic_csv(remaining, file.path(output_dir, "mv08p-remaining-source-queue.csv"))
atomic_csv(validation, file.path(output_dir, "mv08p-validation.csv"))
report <- c(
  "# MV8-P full Pearson-residual source-production prefreeze", "",
  "## Decision", "",
  "After this prefreeze is committed, exactly 129 remaining source fits may run serially in ascending frozen-QC cell count. The three completed MV8-O primary source fits are not rerun. Each new fit uses the MV8-N all-QC standard-SCTransform/Pearson-residual contract and the unchanged selected-384 axes.", "",
  "## Resource policy", "",
  "The policy remains one worker, zero automatic retries, 1,800 seconds, and 12 GiB per source process. The maximum completed sentinel used 12.708 GB, leaving 176.7 MB (about 1.37%) below the cap; the largest remaining source has 9,071 rather than 11,475 QC cells. Runs remain serial and stop safely on a resource/precondition failure.", "",
  "## Scope firewall", "",
  "This authorizes only source fitting and private cache construction. Persistence diagrams, landscapes, pairwise comparisons, clustering, fusion, labels, biological outcomes, manuscript claims, and default adoption remain closed.", "",
  "## Next gate", "",
  "A source-production closure must rehash every private cache and publish aggregate resource/geometry evidence before any PH or landscape execution is considered."
)
atomic_text(report, file.path(output_dir, "MV08P_FULL_SOURCE_PRODUCTION_PREFREEZE_2026-08-21.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08p-artifact-manifest.csv"]
manifest <- data.frame(artifact = basename(artifacts), bytes = file.info(artifacts)$size,
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE)
atomic_csv(manifest, file.path(output_dir, "mv08p-artifact-manifest.csv"))
cat("MV8-P prefreeze checks=", sum(validation$passed), "/", nrow(validation),
  "; 129 serial source fits authorized after commit; topology remains closed\n", sep = "")
