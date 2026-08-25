#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: build_mv08zz_comparison_closure.R <prefreeze> <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <private-output>",
  "<public-output> <audit-output>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
roots <- vapply(args[2:4], normalizePath, character(1L), mustWork = TRUE)
private_root <- normalizePath(args[[5L]], mustWork = TRUE)
public_root <- normalizePath(args[[6L]], mustWork = TRUE)
output <- args[[7L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV8-ZZ closure")
}
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
.mv08z_verify_manifest(prefreeze, "mv08zy-artifact-manifest.csv")
readc <- .mv08z_read_csv
sha <- .mv08z_sha256_file
atomic <- .mv08z_atomic_csv
contract <- readc(file.path(prefreeze, "mv08zy-contract.csv"))
queue <- readc(file.path(prefreeze, "mv08zy-comparison-queue.csv"))
catalog <- readc(file.path(prefreeze, "mv08zy-stack-bindings.csv"))
terminal <- readc(file.path(public_root, "mv08zy-terminal-receipt.csv"))
ledger <- readc(file.path(public_root, "mv08zy-resource-ledger.csv"))
completed <- readc(file.path(public_root, "mv08zy-completions.csv"))
public_summary <- readc(file.path(public_root, "mv08zy-comparison-summary.csv"))
if (nrow(contract) != 1L || nrow(queue) != 40L || nrow(catalog) != 38L ||
    nrow(terminal) != 1L || terminal$completion_state != "complete" ||
    nrow(ledger) != 40L || nrow(completed) != 40L ||
    nrow(public_summary) != 40L) stop("MV8-ZZ cardinality drift")

rehash <- vector("list", 40L)
recompute <- vector("list", 40L)
for (i in 1:40) {
  row <- queue[i, , drop = FALSE]
  left_binding <- catalog[catalog$catalog_order == row$left_catalog_order,
                          , drop = FALSE]
  right_binding <- catalog[catalog$catalog_order == row$right_catalog_order,
                           , drop = FALSE]
  left <- mv08zy_read_distance_stack_v1(left_binding, roots[[1L]], roots[[2L]],
                                        roots[[3L]])
  right <- mv08zy_read_distance_stack_v1(right_binding, roots[[1L]], roots[[2L]],
                                         roots[[3L]])
  root <- file.path(private_root, "jobs", sprintf("job_%02d", i))
  paths <- file.path(root, c("summary.csv", "neighbor.csv", "pair-axis.csv",
                             "status.csv"))
  if (!all(file.exists(paths))) stop("MV8-ZZ private job missing")
  saved <- readc(paths[[1L]])
  saved_neighbor <- readc(paths[[2L]])
  saved_axis <- readc(paths[[3L]])
  status <- readc(paths[[4L]])
  fresh <- mv08zy_compare_distance_pairs_v1(left$pairs, right$pairs,
                                             row$comparison_id)
  numeric_names <- names(fresh$summary)[vapply(fresh$summary, is.numeric,
                                               logical(1L))]
  shared_numeric <- intersect(numeric_names, names(saved))
  numeric_delta <- max(abs(as.numeric(saved[1, shared_numeric]) -
                           as.numeric(fresh$summary[1, shared_numeric])))
  neighbor_delta <- max(abs(saved_neighbor$neighbor_jaccard -
                            fresh$neighbor$neighbor_jaccard))
  valid <- nrow(status) == 1L && status$completion_state == "complete" &&
    status$execution_head == contract$execution_head &&
    left$payload_set_sha256 == left_binding$payload_set_sha256 &&
    right$payload_set_sha256 == right_binding$payload_set_sha256 &&
    left$pair_axis_sha256 == row$pair_axis_sha256 &&
    right$pair_axis_sha256 == row$pair_axis_sha256 &&
    identical(saved_axis, fresh$pair_axis) &&
    identical(saved_neighbor[c("comparison_id", "unit_id", "k")],
              fresh$neighbor[c("comparison_id", "unit_id", "k")]) &&
    numeric_delta <= 1e-12 && neighbor_delta <= 1e-12 &&
    sha(paths[[1L]]) == completed$summary_sha256[[i]] &&
    sha(paths[[2L]]) == completed$neighbor_sha256[[i]] &&
    sha(paths[[3L]]) == completed$pair_axis_payload_sha256[[i]] &&
    sha(paths[[4L]]) == completed$status_sha256[[i]]
  if (!valid) stop("MV8-ZZ independent reconstruction drift at ", i)
  rehash[[i]] <- data.frame(
    contract_id = "mv08zz_private_rehash_v1", comparison_order = i,
    summary_bytes = as.numeric(file.info(paths[[1L]])$size),
    summary_sha256 = sha(paths[[1L]]),
    neighbor_bytes = as.numeric(file.info(paths[[2L]])$size),
    neighbor_sha256 = sha(paths[[2L]]),
    pair_axis_bytes = as.numeric(file.info(paths[[3L]])$size),
    pair_axis_sha256 = sha(paths[[3L]]), status_sha256 = sha(paths[[4L]]),
    independently_rehashed = TRUE, stringsAsFactors = FALSE
  )
  recompute[[i]] <- data.frame(
    contract_id = "mv08zz_recomputation_v1", comparison_order = i,
    dataset_scope = row$dataset_scope, contrast_id = row$contrast_id,
    seed = row$seed, homology_dimension = row$homology_dimension,
    units = fresh$summary$units, unordered_pairs = fresh$summary$unordered_pairs,
    maximum_summary_numeric_difference = numeric_delta,
    maximum_neighbor_difference = neighbor_delta,
    pair_axis_sha256 = row$pair_axis_sha256,
    independently_recomputed = TRUE, stringsAsFactors = FALSE
  )
}
rehash <- do.call(rbind, rehash)
recompute <- do.call(rbind, recompute)
resource <- data.frame(
  contract_id = "mv08zz_resource_summary_v1", jobs = 40L,
  pair_alignments = sum(completed$unordered_pairs),
  aggregate_child_seconds = sum(ledger$elapsed_seconds),
  aggregate_elapsed_cap_seconds = contract$aggregate_elapsed_cap_seconds,
  peak_process_tree_rss_bytes = max(ledger$peak_process_tree_rss_bytes),
  child_rss_cap_bytes = contract$child_rss_cap_bytes,
  private_bytes = terminal$private_bytes,
  private_storage_cap_bytes = contract$private_storage_cap_bytes,
  elapsed_cap_passed = sum(ledger$elapsed_seconds) <=
    contract$aggregate_elapsed_cap_seconds,
  rss_cap_passed = all(ledger$peak_process_tree_rss_bytes <=
                         contract$child_rss_cap_bytes),
  storage_cap_passed = terminal$private_bytes <=
    contract$private_storage_cap_bytes, workers = 1L, retries = 0L,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv08zz_validation_v1",
  check_id = c("prefreeze_manifest", "terminal_complete", "job_cardinality",
               "source_rehash", "pair_axis_identity", "fresh_recomputation",
               "H0_H1_separate", "public_aggregate_only", "resource_caps",
               "one_worker_zero_retry", "label_outcome_firewall",
               "clustering_fusion_firewall"),
  passed = c(TRUE, terminal$completion_state == "complete",
             nrow(recompute) == 40L, all(rehash$independently_rehashed),
             all(nzchar(recompute$pair_axis_sha256)),
             all(recompute$independently_recomputed),
             sum(queue$homology_dimension == "H0") == 20L &&
               sum(queue$homology_dimension == "H1") == 20L,
             !any(grepl("unit_id|pair_key", names(public_summary))),
             all(resource[c("elapsed_cap_passed", "rss_cap_passed",
                            "storage_cap_passed")]),
             resource$workers == 1L && resource$retries == 0L,
             terminal$outcome_label_state == "closed" &&
               !terminal$biological_outcomes_computed,
             terminal$clustering_jobs == 0L && terminal$fusion_jobs == 0L),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-ZZ validation failed")
decision <- data.frame(
  contract_id = "mv08zz_decision_v1",
  decision = "close_label_closed_distance_comparisons",
  comparisons = 40L, interpretation =
    "descriptive_sensitivity_evidence_no_equivalence_or_biological_claim",
  next_stage = "separate_label_closed_robustness_synthesis_prefreeze",
  clustering_state = "closed", fusion_state = "closed",
  labels_outcomes_state = "closed", manuscript_claims_state = "closed",
  stringsAsFactors = FALSE
)
atomic(recompute, file.path(output, "mv08zz-recomputation-summary.csv"))
atomic(rehash, file.path(output, "mv08zz-private-rehash.csv"))
atomic(resource, file.path(output, "mv08zz-resource-summary.csv"))
atomic(validation, file.path(output, "mv08zz-validation.csv"))
atomic(decision, file.path(output, "mv08zz-decision.csv"))
writeLines(c(
  "# MV8-ZZ label-closed distance-comparison closure",
  "", "All 40 gene-topology comparisons were independently reloaded,",
  "rehashed, and recomputed. H0/H1 remain separate. Results are descriptive;",
  "labels, outcomes, clustering, fusion, adoption, and claims remain closed."
), file.path(output, "MV08ZZ_DISTANCE_COMPARISON_CLOSURE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv08zz-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv08zz_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv08zz-artifact-manifest.csv"))
message("Closed MV8-ZZ comparisons; checks=", nrow(validation))
