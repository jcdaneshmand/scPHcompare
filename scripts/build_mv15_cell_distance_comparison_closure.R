#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: build_mv15_cell_distance_comparison_closure.R <prefreeze>",
  "<mv14-private-root> <mv08zu-private-root> <private-output>",
  "<public-output> <gnu-time-stderr> <audit-output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:3], normalizePath, character(1L), mustWork = TRUE)
private <- normalizePath(args[[4L]], mustWork = TRUE)
public <- normalizePath(args[[5L]], mustWork = TRUE)
time_evidence <- normalizePath(args[[6L]], mustWork = TRUE)
output <- args[[7L]]
if (dir.exists(output)) stop("MV15 closure output exists", call. = FALSE)
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv15_cell_distance_comparison.R")
.mv08z_verify_manifest(prefreeze, "mv15-artifact-manifest.csv")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
contract <- readc(file.path(prefreeze, "mv15-contract.csv"))
bindings <- readc(file.path(prefreeze, "mv15-stack-bindings.csv"))
queue <- readc(file.path(prefreeze, "mv15-comparison-queue.csv"))
terminal <- readc(file.path(public, "mv15-terminal-receipt.csv"))
public_global <- readc(file.path(public, "mv15-global-summary.csv"))
public_neighbor <- readc(file.path(public, "mv15-neighbor-summary.csv"))
ledger <- readc(file.path(public, "mv15-resource-ledger.csv"))
time_lines <- readLines(time_evidence, warn = FALSE)
rss_line <- grep("Maximum resident set size [(]kbytes[)]:", time_lines,
                 value = TRUE)
exit_line <- grep("Exit status:", time_lines, value = TRUE)
peak_rss_bytes <- if (length(rss_line) == 1L) {
  as.numeric(sub(".*: *", "", rss_line)) * 1024
} else NA_real_
gnu_time_exit_status <- if (length(exit_line) == 1L) {
  as.integer(sub(".*: *", "", exit_line))
} else NA_integer_
if (nrow(contract) != 1L || nrow(bindings) != 28L || nrow(queue) != 36L ||
    nrow(terminal) != 1L || terminal$completion_state != "complete" ||
    nrow(public_global) != 36L || nrow(public_neighbor) != 42L ||
    nrow(ledger) != 36L) stop("MV15 closure cardinality drift")
loaded <- lapply(seq_len(nrow(bindings)), function(i) {
  result <- mv15_read_bound_stack_v1(
    bindings[i, , drop = FALSE], roots[[1L]], roots[[2L]]
  )
  if (result$payload_set_sha256 != bindings$payload_set_sha256[[i]] ||
      result$pair_axis_sha256 != bindings$pair_axis_sha256[[i]]) {
    stop("MV15 closure source rehash drift at stack ", i)
  }
  result
})
recompute <- vector("list", 36L)
rehash <- vector("list", 36L)
numeric_delta <- function(saved, fresh) {
  shared <- intersect(names(saved)[vapply(saved, is.numeric, logical(1L))],
                      names(fresh)[vapply(fresh, is.numeric, logical(1L))])
  max(abs(as.numeric(saved[1L, shared]) - as.numeric(fresh[1L, shared])))
}
for (i in 1:36) {
  row <- queue[i, , drop = FALSE]
  left <- loaded[[as.integer(row$left_stack_order)]]
  right <- loaded[[as.integer(row$right_stack_order)]]
  fresh <- mv15_compare_distance_pairs_v1(
    left$pairs, right$pairs, row$comparison_id,
    as.integer(strsplit(row$neighbor_k, ";", fixed = TRUE)[[1L]])
  )
  job <- file.path(private, "jobs", sprintf("job_%02d", i))
  paths <- file.path(job, c("summary.csv", "neighbor-summary.csv",
                            "neighbor.csv", "pair-axis.csv", "status.csv"))
  if (!all(file.exists(paths))) stop("MV15 closure private job missing")
  saved <- readc(paths[[1L]])
  saved_neighbor_summary <- readc(paths[[2L]])
  saved_neighbor <- readc(paths[[3L]])
  saved_axis <- readc(paths[[4L]])
  status <- readc(paths[[5L]])
  global_delta <- numeric_delta(saved, fresh$summary)
  neighbor_summary_delta <- max(abs(
    saved_neighbor_summary$mean_neighbor_jaccard -
      fresh$neighbor_summary$mean_neighbor_jaccard,
    saved_neighbor_summary$median_neighbor_jaccard -
      fresh$neighbor_summary$median_neighbor_jaccard,
    saved_neighbor_summary$p10_neighbor_jaccard -
      fresh$neighbor_summary$p10_neighbor_jaccard
  ))
  neighbor_delta <- max(abs(saved_neighbor$neighbor_jaccard -
                              fresh$neighbor$neighbor_jaccard))
  axis_identical <- identical(saved_axis$first_unit_id,
                              fresh$pair_axis$first_unit_id) &&
    identical(saved_axis$second_unit_id, fresh$pair_axis$second_unit_id)
  public_row <- public_global[public_global$comparison_order == i,
                              , drop = FALSE]
  public_neighbor_rows <- public_neighbor[
    public_neighbor$comparison_order == i, , drop = FALSE
  ]
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      global_delta > 1e-12 || neighbor_summary_delta > 1e-12 ||
      neighbor_delta > 1e-12 || !axis_identical || nrow(public_row) != 1L ||
      nrow(public_neighbor_rows) != nrow(saved_neighbor_summary) ||
      numeric_delta(public_row, saved) > 1e-14 ||
      max(abs(public_neighbor_rows$mean_neighbor_jaccard -
                saved_neighbor_summary$mean_neighbor_jaccard)) > 1e-14) {
    stop("MV15 independent reconstruction drift at job ", i)
  }
  rehash[[i]] <- data.frame(
    contract_id = "mv15_private_rehash_v1", comparison_order = i,
    summary_sha256 = sha(paths[[1L]]),
    neighbor_summary_sha256 = sha(paths[[2L]]),
    neighbor_sha256 = sha(paths[[3L]]),
    pair_axis_sha256 = sha(paths[[4L]]), status_sha256 = sha(paths[[5L]]),
    independently_rehashed = TRUE, stringsAsFactors = FALSE
  )
  recompute[[i]] <- data.frame(
    contract_id = "mv15_recomputation_v1", comparison_order = i,
    contrast_family = row$contrast_family, dataset_scope = row$dataset_scope,
    homology_dimension = row$homology_dimension,
    maximum_global_difference = global_delta,
    maximum_neighbor_summary_difference = neighbor_summary_delta,
    maximum_unit_neighbor_difference = neighbor_delta,
    pair_axis_identical = axis_identical,
    independently_recomputed = TRUE, stringsAsFactors = FALSE
  )
}
rehash <- do.call(rbind, rehash)
recompute <- do.call(rbind, recompute)
resource <- data.frame(
  contract_id = "mv15_resource_closure_v1", comparisons = 36L,
  aggregate_job_seconds = terminal$aggregate_job_seconds,
  total_elapsed_seconds = terminal$total_elapsed_seconds,
  elapsed_cap_seconds = contract$elapsed_cap_seconds,
  peak_rss_bytes = peak_rss_bytes, rss_cap_bytes = contract$rss_cap_bytes,
  GNU_time_exit_status = gnu_time_exit_status,
  private_bytes = terminal$private_bytes,
  private_storage_cap_bytes = contract$private_storage_cap_bytes,
  public_bytes = sum(file.info(list.files(public, full.names = TRUE))$size),
  public_storage_cap_bytes = contract$public_storage_cap_bytes,
  workers = terminal$workers, retries = terminal$retries,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv15_closure_validation_v1",
  check_id = c(
    "prefreeze_manifest", "terminal_complete", "stack_cardinality",
    "comparison_cardinality", "global_cardinality", "neighbor_cardinality",
    "source_rehash", "private_rehash", "fresh_recomputation",
    "pair_axis_identity", "numeric_tolerance", "cell_seed_complete",
    "cell_panel_complete", "cell_gene_complete", "H0_H1_separate",
    "external_k2_k3", "internal_k10", "resource_caps",
    "one_worker_zero_retry", "public_aggregate_only",
    "label_outcome_firewall", "downstream_firewall"
  ),
  passed = c(
    TRUE, terminal$completion_state == "complete", nrow(bindings) == 28L,
    nrow(recompute) == 36L, nrow(public_global) == 36L,
    nrow(public_neighbor) == 42L,
    all(vapply(seq_along(loaded), function(i)
      loaded[[i]]$payload_set_sha256 == bindings$payload_set_sha256[[i]],
      logical(1L))),
    all(rehash$independently_rehashed), all(recompute$independently_recomputed),
    all(recompute$pair_axis_identical),
    max(recompute$maximum_global_difference,
        recompute$maximum_neighbor_summary_difference,
        recompute$maximum_unit_neighbor_difference) <= 1e-12,
    sum(recompute$contrast_family == "cell_seed_stability") == 20L,
    sum(recompute$contrast_family == "cell_panel_sensitivity") == 2L,
    sum(recompute$contrast_family == "cell_gene_view_agreement") == 14L,
    sum(recompute$homology_dimension == "H0") == 18L &&
      sum(recompute$homology_dimension == "H1") == 18L,
    all(sort(unique(public_neighbor$k[
      public_neighbor$dataset_scope == "external8"])) == c(2L, 3L)),
    all(unique(public_neighbor$k[
      public_neighbor$dataset_scope == "internal124"]) == 10L),
    resource$total_elapsed_seconds <= resource$elapsed_cap_seconds &&
      resource$peak_rss_bytes <= resource$rss_cap_bytes &&
      resource$GNU_time_exit_status == 0L &&
      resource$private_bytes <= resource$private_storage_cap_bytes &&
      resource$public_bytes <= resource$public_storage_cap_bytes,
    resource$workers == 1L && resource$retries == 0L,
    !any(grepl("unit_id|pair_key", names(public_global))) &&
      !any(grepl("unit_id|pair_key", names(public_neighbor))),
    !terminal$labels_used && !terminal$outcomes_used,
    terminal$clustering_jobs == 0L && terminal$fusion_jobs == 0L &&
      terminal$inference_jobs == 0L && terminal$biological_claim_jobs == 0L &&
      terminal$manuscript_claim_jobs == 0L
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV15 independent closure failed")
decision <- data.frame(
  contract_id = "mv15_closure_decision_v1",
  distance_comparisons_independently_closed = TRUE,
  values_not_interpreted_by_closure = TRUE,
  descriptive_synthesis_prefreeze_eligible_next = TRUE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  biological_claims_authorized = FALSE, manuscript_claims_authorized = FALSE,
  next_action = "prospectively_freeze_threshold_free_descriptive_synthesis",
  stringsAsFactors = FALSE
)
atomic(rehash, file.path(output, "mv15-private-rehash.csv"))
atomic(recompute, file.path(output, "mv15-recomputation.csv"))
atomic(resource, file.path(output, "mv15-resource-closure.csv"))
atomic(validation, file.path(output, "mv15-validation.csv"))
atomic(decision, file.path(output, "mv15-decision.csv"))
writeLines(c(
  "# MV15 all-QC cell-distance comparison closure", "",
  "All 36 prospectively frozen comparisons were independently reloaded,",
  "rehashed, and recomputed. H0/H1 remain separate, unit identities remain",
  "private, and no result was interpreted by this closure.", "",
  "A separate threshold-free descriptive synthesis may now be prefrozen.",
  "Clustering, fusion, labels, outcomes, biology, and claims remain closed."
), file.path(output, "MV15_CELL_DISTANCE_COMPARISON_CLOSURE_2026-08-26.md"))
files <- sort(setdiff(list.files(output), "mv15-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv15_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv15-artifact-manifest.csv"))
message("Closed MV15 comparisons; checks=", nrow(validation))
