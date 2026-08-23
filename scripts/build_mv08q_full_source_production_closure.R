#!/usr/bin/env Rscript

# Independently rehash and close MV8-P's 129 private source caches. The three
# MV8-O primary source fits remain bound by their prior public closure. This
# builder publishes no matrices, barcodes, private paths, PH, or outcomes.

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(4L, 7L))) stop(paste(
  "usage: build_mv08q_full_source_production_closure.R <mv08p-audit-dir>",
  "<original-private-dir> <original-public-dir>",
  "[<mv08pr-prefreeze-dir> <recovery-private-dir> <recovery-public-dir>]",
  "<closure-output-dir>"), call. = FALSE)
has_recovery <- length(args) == 7L
audit_dir <- normalizePath(args[[1L]], mustWork = TRUE)
private_dir <- normalizePath(args[[2L]], mustWork = TRUE)
run_dir <- normalizePath(args[[3L]], mustWork = TRUE)
if (has_recovery) {
  recovery_audit_dir <- normalizePath(args[[4L]], mustWork = TRUE)
  recovery_private_dir <- normalizePath(args[[5L]], mustWork = TRUE)
  recovery_run_dir <- normalizePath(args[[6L]], mustWork = TRUE)
}
output_dir <- normalizePath(args[[if (has_recovery) 7L else 4L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-Q closure output", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required", call. = FALSE)
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

queue <- utils::read.csv(file.path(audit_dir, "mv08p-remaining-source-queue.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
original_resource <- utils::read.csv(file.path(run_dir, "mv08p-source-production-resource.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
original_progress <- utils::read.csv(file.path(run_dir, "mv08p-source-production-progress.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
sentinel <- utils::read.csv(
  "docs/audits/mv08o-residual-source-sentinel-closure-v1/mv08o-source-sentinel-resource.csv",
  check.names = FALSE, stringsAsFactors = FALSE)
private_dirs <- rep(private_dir, 129L)
recovery_summary <- NULL
if (has_recovery) {
  recovery_contract <- utils::read.csv(file.path(recovery_audit_dir, "mv08pr-contract.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  recovery_evidence <- utils::read.csv(file.path(recovery_audit_dir, "mv08pr-original-evidence.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  recovery_resource <- utils::read.csv(file.path(recovery_run_dir, "mv08pr-source-production-resource.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  recovery_progress <- utils::read.csv(file.path(recovery_run_dir, "mv08pr-source-production-progress.csv"),
    check.names = FALSE, stringsAsFactors = FALSE)
  stopped_paths <- c(file.path(run_dir, "mv08p-source-production-resource.csv"),
    file.path(run_dir, "mv08p-source-production-progress.csv"),
    file.path(private_dir, "logs", paste0(queue$unit_id[[124L]], "__primary-stderr.txt")),
    file.path(private_dir, "logs", paste0(queue$unit_id[[124L]], "__primary-stdout.txt")))
  stopped_evidence_valid <- nrow(recovery_evidence) == 4L && all(file.exists(stopped_paths)) &&
    identical(as.numeric(file.info(stopped_paths)$size), as.numeric(recovery_evidence$bytes)) &&
    identical(tolower(unname(vapply(stopped_paths, sha_file, character(1L)))),
      tolower(recovery_evidence$sha256))
  if (nrow(recovery_contract) != 1L || !stopped_evidence_valid ||
      nrow(original_resource) != 124L ||
      !identical(as.integer(original_resource$job_order), seq_len(124L)) ||
      !all(original_resource$disposition[1:123] == "completed") ||
      original_resource$disposition[[124L]] != "child_failed" ||
      nrow(original_progress) != 1L || original_progress$state != "stopped" ||
      original_progress$completed_jobs != 123L ||
      nrow(recovery_resource) != 6L ||
      !identical(as.integer(recovery_resource$job_order), 124:129) ||
      !all(recovery_resource$disposition == "completed") ||
      !identical(as.integer(recovery_resource$attempt_number), c(2L, rep(1L, 5L))) ||
      !all(recovery_resource$future_globals_max_size_bytes == 2 * 1024^3) ||
      nrow(recovery_progress) != 1L ||
      recovery_progress$state != "source_production_complete_closure_pending" ||
      recovery_progress$completed_jobs != 6L || recovery_progress$overall_completed_jobs != 129L) {
    stop("MV8-Q recovery prerequisite drift", call. = FALSE)
  }
  if (!all(names(original_resource) %in% names(recovery_resource))) {
    stop("MV8-Q recovery resource schema drift", call. = FALSE)
  }
  resource <- rbind(original_resource[1:123, , drop = FALSE],
    recovery_resource[, names(original_resource), drop = FALSE])
  resource$execution_source <- c(rep("mv08p_original_v1", 123L), rep("mv08pr_overlay_v1", 6L))
  resource$attempt_number <- c(rep(1L, 123L), as.integer(recovery_resource$attempt_number))
  private_dirs[124:129] <- recovery_private_dir
  progress_complete <- TRUE
  recovery_summary <- data.frame(
    contract_id = "mv08q_recovery_summary_v1", original_completed_jobs = 123L,
    original_failed_job_order = 124L, recovery_completed_jobs = 6L,
    explicit_retry_jobs = 1L, future_globals_max_size_bytes = 2 * 1024^3,
    original_resource_sha256 = sha_file(stopped_paths[[1L]]),
    original_progress_sha256 = sha_file(stopped_paths[[2L]]),
    original_failed_stderr_sha256 = sha_file(stopped_paths[[3L]]),
    original_failed_stdout_sha256 = sha_file(stopped_paths[[4L]]),
    original_evidence_preserved = TRUE, topology_execution_state = "closed",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
} else {
  resource <- original_resource
  progress_complete <- nrow(original_progress) == 1L &&
    original_progress$state == "source_production_complete_closure_pending" &&
    original_progress$completed_jobs == 129L
  resource$execution_source <- "mv08p_original_v1"
  resource$attempt_number <- 1L
}
if (nrow(queue) != 129L || nrow(resource) != 129L ||
    !identical(as.integer(queue$job_order), as.integer(resource$job_order)) ||
    !identical(as.character(queue$unit_id), as.character(resource$unit_id)) ||
    !progress_complete || nrow(sentinel) != 5L) {
  stop("MV8-Q prerequisite ledger drift", call. = FALSE)
}

audits <- vector("list", nrow(queue)); rehash <- vector("list", nrow(queue))
for (i in seq_len(nrow(queue))) {
  cache_path <- file.path(private_dirs[[i]], "cache", queue$output_file[[i]])
  audit_path <- file.path(private_dirs[[i]], "worker-audit", paste0(queue$unit_id[[i]], "__primary.csv"))
  if (!all(file.exists(c(cache_path, audit_path)))) stop("MV8-Q private artifact absent at order ", i, call. = FALSE)
  cache_sha <- sha_file(cache_path); audit_sha <- sha_file(audit_path)
  if (!identical(tolower(cache_sha), tolower(resource$cache_sha256[[i]])) ||
      !identical(tolower(audit_sha), tolower(resource$worker_audit_sha256[[i]]))) {
    stop("MV8-Q private artifact hash drift at order ", i, call. = FALSE)
  }
  audit <- utils::read.csv(audit_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(audit) != resource$worker_rows[[i]]) stop("MV8-Q audit row drift at order ", i, call. = FALSE)
  audits[[i]] <- audit
  rehash[[i]] <- data.frame(
    contract_id = "mv08q_source_cache_rehash_v1", job_order = i,
    unit_id = queue$unit_id[[i]], dataset_scope = queue$dataset_scope[[i]],
    fit_cells = queue$fit_cells[[i]], cache_bytes = file.info(cache_path)$size,
    cache_sha256 = cache_sha, worker_audit_sha256 = audit_sha,
    private_cache_rehashed = TRUE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
}
views <- do.call(rbind, audits)
cache_rehash <- do.call(rbind, rehash)
required <- views$dataset_scope == "internal124" | views$panel_id == "common475" |
  (views$panel_id == "exact500" &
    views$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
external_diagnostic <- views$dataset_scope == "external8" & views$panel_id == "exact500" &
  views$representation_id == "sct_data_all_qc_fit_selected384"
external_selection <- unique(views[views$dataset_scope == "external8",
  c("unit_id", "seed", "selected_cells", "selected_cell_sha256")])
external_selection <- external_selection[order(external_selection$unit_id), , drop = FALSE]
external_selection$contract_id <- "mv08q_external_selected_axis_v1"
external_selection$selection_state <- "frozen_by_mv08p_source_preflight"
external_selection$outcome_label_state <- "closed"
external_selection$biological_outcomes_computed <- FALSE

validation <- data.frame(
  check_id = c("remaining_queue_complete", "all_132_sources_covered", "resource_caps",
    "stderr_policy", "private_caches_rehashed", "worker_row_contract", "required_geometry",
    "external_exact500_data_diagnostic", "external_selection_frozen", "source_identity",
    "topology_firewall", "outcome_firewall"),
  passed = c(
    nrow(resource) == 129L && all(resource$disposition == "completed"),
    nrow(resource) + length(unique(sentinel$unit_id[sentinel$run_role == "primary"])) == 132L,
    all(resource$elapsed_seconds <= resource$elapsed_cap_seconds & resource$peak_rss_bytes <= resource$rss_cap_bytes),
    all(resource$stderr_class %in% c("empty", "known_glmGamPoi_native_fallback")),
    nrow(cache_rehash) == 129L && all(cache_rehash$private_cache_rehashed),
    all(resource$worker_rows == ifelse(resource$dataset_scope == "internal124", 10L, 4L)),
    all(views$values_finite[required] & views$zero_variance_gene_count[required] == 0L &
      views$correlation_chord_valid[required]),
    sum(external_diagnostic) == 7L,
    nrow(external_selection) == 7L && all(external_selection$selected_cells == 384L) &&
      all(grepl("^[0-9a-f]{64}$", external_selection$selected_cell_sha256)),
    identical(as.character(resource$unit_id), as.character(queue$unit_id)),
    !any(resource$persistence_computed | resource$landscapes_computed |
      resource$clustering_computed | resource$fusion_computed),
    all(resource$outcome_label_state == "closed") && !any(resource$biological_outcomes_computed)),
  evidence = c("all 129 MV8-P jobs completed", "129 new plus three MV8-O primary source fits",
    "every job within 1,800 seconds and 12 GiB", "only empty or documented glmGamPoi-native fallback stderr",
    "all private cache/audit hashes independently recomputed", "122 internal ten-row plus seven external four-row audits",
    "all PH-eligible correlation-chord geometries valid", "seven SCT-data exact500 views remain diagnostic-only",
    "seven external 384-cell axes frozen by deterministic preflight", "resource order matches the frozen queue",
    "no PH, landscapes, clustering, or fusion", "labels and biological outcomes remained closed"),
  stringsAsFactors = FALSE)
if (has_recovery) {
  validation <- rbind(validation, data.frame(
    check_id = c("stopped_run_evidence_preserved", "bounded_recovery_overlay",
      "single_explicit_retry"),
    passed = c(recovery_summary$original_evidence_preserved,
      nrow(recovery_resource) == 6L && all(recovery_resource$disposition == "completed") &&
        all(recovery_resource$future_globals_max_size_bytes == 2 * 1024^3),
      identical(as.integer(recovery_resource$attempt_number), c(2L, rep(1L, 5L)))),
    evidence = c("stopped v1 ledger, progress, and failed child logs match the recovery prefreeze",
      "jobs 124-129 completed in the fresh overlay with bounded 2-GiB future exports",
      "only job 124 received one explicit second attempt; jobs 125-129 were first attempts"),
    stringsAsFactors = FALSE))
}
if (!all(validation$passed)) stop("MV8-Q source-production closure validation failed", call. = FALSE)

atomic_csv(resource, file.path(output_dir, "mv08q-source-production-resource.csv"))
atomic_csv(cache_rehash, file.path(output_dir, "mv08q-source-cache-rehash.csv"))
atomic_csv(views, file.path(output_dir, "mv08q-source-view-summary.csv"))
atomic_csv(external_selection, file.path(output_dir, "mv08q-external-selected-axis.csv"))
atomic_csv(validation, file.path(output_dir, "mv08q-validation.csv"))
if (has_recovery) {
  atomic_csv(recovery_summary, file.path(output_dir, "mv08q-recovery-summary.csv"))
}
report <- c(
  "# MV8-Q full Pearson-residual source-production closure", "",
  "## Result", "",
  if (has_recovery)
    "The preserved MV8-P run contributed 123 accepted fits and the bounded MV8-PR overlay contributed jobs 124-129. Together with the three MV8-O primary source fits, the full 132-source representation layer is covered. Every private cache and worker audit was independently rehashed."
  else
    "All 129 MV8-P source fits completed within the frozen serial resource policy. Together with the three MV8-O primary source fits, the full 132-source representation layer is covered. Every private cache and worker audit was independently rehashed.", "",
  "## Scientific boundary", "",
  "This closes representation construction only. Required internal exact500 and external common475/exact500-residual geometries pass. External SCT-data exact500 remains diagnostic-only. No persistence diagrams, landscapes, comparisons, clustering, fusion, labels, outcomes, claims, or default adoption were produced.", "",
  "## Next gate", "",
  "A separate topology-production prefreeze must bind the PH queue, backend fallback, exact all-active-level landscape contract, resource staging, and comparison firewall before any diagram is computed."
)
atomic_text(report, file.path(output_dir, "MV08Q_FULL_SOURCE_PRODUCTION_CLOSURE_2026-08-23.md"))
artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08q-artifact-manifest.csv"]
manifest <- data.frame(artifact = basename(artifacts), bytes = file.info(artifacts)$size,
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE)
atomic_csv(manifest, file.path(output_dir, "mv08q-artifact-manifest.csv"))
cat("MV8-Q closure checks=", sum(validation$passed), "/", nrow(validation),
  "; 132/132 sources covered; topology remains closed\n", sep = "")
