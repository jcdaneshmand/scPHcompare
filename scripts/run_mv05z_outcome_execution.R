#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-Z outcome execution.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste("usage: run_mv05z_outcome_execution.R EXTERNAL_METADATA_PATH",
             "PRIVATE_ROOT AUDIT_DIR PREDICTION_RANKING_GZ"), call. = FALSE)
}
metadata_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
private_root <- args[[2L]]; audit_dir <- args[[3L]]
ranking_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
unit_root <- file.path(private_root, "units")
dir.create(unit_root, recursive = TRUE, showWarnings = FALSE)

readc <- function(path, ...) read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE, ...)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_rds <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(temporary)) stop("Stale MV5-Z outcome temporary artifact.", call. = FALSE)
  saveRDS(value, temporary, version = 3, compress = "xz")
  if (!file.rename(temporary, path)) {
    unlink(temporary); stop("MV5-Z outcome atomic rename failed.", call. = FALSE)
  }
}
atomic_csv <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  write_provenance_csv(value, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary); stop("MV5-Z outcome atomic CSV rename failed.", call. = FALSE)
  }
}

expected_head <- "c16f2b2"
head <- trimws(system2("git", c("rev-parse", "--short=7", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-Z outcomes require committed prediction lock ", expected_head,
       "; observed ", head, ".", call. = FALSE)
}

y_freeze_path <- "docs/audits/mv05y-source-freeze-2026-08-11.csv"
y_queue_path <- "docs/audits/mv05y-evaluation-queue-2026-08-11.csv"
prediction_lock_path <- "docs/audits/mv05z-prediction-lock-2026-08-11.csv"
prediction_completion_path <- "docs/audits/mv05z-ranking-unit-completion-2026-08-11.csv"
prediction_validation_path <-
  "docs/audits/mv05z-prediction-independent-validation-2026-08-11.csv"
baseline_sct_path <- "docs/audits/mv05e-query-endpoints-2026-08-08.csv"
baseline_integrated_path <- "docs/audits/mv05k-query-endpoints-2026-08-10.csv"
required_paths <- c(y_freeze_path, y_queue_path, prediction_lock_path,
                    prediction_completion_path, prediction_validation_path,
                    baseline_sct_path, baseline_integrated_path, ranking_path)
if (any(!file.exists(required_paths))) stop("MV5-Z outcome source missing.", call. = FALSE)

y_freeze <- readc(y_freeze_path); queue <- readc(y_queue_path)
prediction_lock <- readc(prediction_lock_path)
prediction_completion <- readc(prediction_completion_path)
prediction_validation <- readc(prediction_validation_path)
if (nrow(y_freeze) != 178L || nrow(queue) != 150L ||
    nrow(prediction_lock) != 1L || nrow(prediction_completion) != 150L ||
    !all(prediction_validation$passed) ||
    !prediction_lock$label_open_authorized ||
    prediction_lock$labels_opened || prediction_lock$outcomes_computed ||
    sha(ranking_path) != prediction_lock$ranking_sha256 ||
    any(prediction_completion$labels_opened) ||
    any(prediction_completion$outcomes_computed)) {
  stop("MV5-Z durable prediction lock failed before label access.", call. = FALSE)
}

# Revalidate every MV5-Y source before parsing any label or accepted endpoint.
source_checks <- lapply(seq_len(nrow(y_freeze)), function(index) {
  row <- y_freeze[index, , drop = FALSE]
  locator <- row$artifact_locator[[1L]]
  path <- if (row$source_id[[1L]] == "external_label_source") metadata_path
  else if (startsWith(locator, "ignored_tmp:")) {
    base <- sub("^ignored_tmp:", "", locator)
    if (grepl("^pc20_private_group_manifest_", row$source_id[[1L]]))
      file.path(base, "artifact_manifest.csv") else base
  } else locator
  actual <- sha(path)
  data.frame(
    contract_id = "mv05z_outcome_source_validation_v1",
    source_id = row$source_id, artifact_locator = locator,
    expected_sha256 = row$sha256, actual_sha256 = actual,
    matched = identical(actual, row$sha256[[1L]]),
    validated_before_label_access = TRUE,
    stringsAsFactors = FALSE)
})
source_checks <- do.call(rbind, source_checks)
if (!all(source_checks$matched)) stop("MV5-Z source drift before labels.", call. = FALSE)

expected_sct_sha <- y_freeze$sha256[y_freeze$source_id == "e_query_endpoints"]
expected_integrated_sha <- y_freeze$sha256[y_freeze$source_id == "k_query_endpoints"]
if (length(expected_sct_sha) != 1L || length(expected_integrated_sha) != 1L ||
    sha(baseline_sct_path) != expected_sct_sha ||
    sha(baseline_integrated_path) != expected_integrated_sha) {
  stop("MV5-Z accepted baseline endpoint hash drifted.", call. = FALSE)
}

implementation_paths <- c(
  engine = "R/mv05z_outcome_execution.R",
  mv05y_engine = "R/mv05y_robustness_outcome_prefreeze.R",
  specification = "docs/specifications/MV05Z_PC20_RETRIEVAL_ROBUSTNESS_EXECUTION_SPECIFICATION_V1.md",
  outcome_runner = "scripts/run_mv05z_outcome_execution.R",
  prediction_lock = prediction_lock_path,
  prediction_validation = prediction_validation_path,
  estimand_registry = "docs/audits/mv05y-estimand-registry-2026-08-11.csv")
implementation <- data.frame(
  contract_id = "mv05z_outcome_implementation_freeze_v1",
  source_id = names(implementation_paths),
  artifact_locator = unname(implementation_paths),
  sha256 = vapply(implementation_paths, sha, character(1L)),
  bytes = as.numeric(file.info(implementation_paths)$size),
  accepted_head = expected_head, stringsAsFactors = FALSE)
implementation$implementation_freeze_sha256 <- .mv05z_digest(
  paste(implementation$artifact_locator, implementation$sha256, sep = "\r"))
implementation_id <- unique(implementation$implementation_freeze_sha256)

# Only now, after the durable lock and 178 fresh source checks, parse sample,
# study, and tissue. Approach and every other label column are skipped.
header <- readc(metadata_path, nrows = 0L)
required_labels <- c("orig.ident", "SRA", "Tissue.x")
if (!all(required_labels %in% names(header))) stop("MV5-Z label schema drifted.", call. = FALSE)
classes <- rep("NULL", ncol(header)); names(classes) <- names(header)
classes[required_labels] <- "character"
raw <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE,
                colClasses = classes)
labels <- data.frame(
  sample_id = trimws(raw$orig.ident), study = trimws(raw$SRA),
  tissue = tolower(trimws(raw$Tissue.x)), stringsAsFactors = FALSE)
labels <- labels[labels$tissue %in% .mv05e_eligible_tissues, , drop = FALSE]
if (nrow(labels) != 90L || length(unique(labels$study)) != 15L ||
    anyNA(labels) || anyDuplicated(labels$sample_id) ||
    any(tapply(labels$tissue, labels$study, function(x) length(unique(x))) != 1L)) {
  stop("MV5-Z exact tissue-label join design drifted.", call. = FALSE)
}
label_provenance <- data.frame(
  contract_id = "mv05z_label_source_provenance_v1",
  source_file = basename(metadata_path), source_sha256 = sha(metadata_path),
  source_rows = 124L, candidate_samples = 90L, candidate_studies = 15L,
  tissue_classes = 5L, parsed_label_columns = "Tissue.x_only",
  approach_or_other_label_columns_parsed = FALSE,
  label_open_boundary = "after_committed_mv05z_prediction_lock_c16f2b2",
  labels_used_for_fit_or_ranking = FALSE, stringsAsFactors = FALSE)

rankings <- readc(ranking_path)
if (nrow(rankings) != 282800L || any(rankings$labels_opened) ||
    any(rankings$outcomes_computed) || any(!rankings$prediction_locked) ||
    !setequal(unique(c(rankings$query_sample_id, rankings$training_sample_id)),
              labels$sample_id)) {
  stop("MV5-Z ranking content drifted after label open.", call. = FALSE)
}

process <- ps::ps_handle(); results <- vector("list", nrow(queue))
completion <- vector("list", nrow(queue))
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  stem <- sub("^mv05y_eval_v1:", "", unit$evaluation_unit_id[[1L]])
  artifact <- file.path(unit_root, paste0(stem, ".rds"))
  status_path <- file.path(unit_root, paste0(stem, ".status.rds"))
  if (xor(file.exists(artifact), file.exists(status_path))) {
    stop("MV5-Z partial outcome artifact/status pair: ", stem, call. = FALSE)
  }
  reused <- FALSE
  if (file.exists(artifact)) {
    status <- readRDS(status_path); outcome <- readRDS(artifact)
    if (!identical(status$evaluation_unit_id, unit$evaluation_unit_id[[1L]]) ||
        !identical(status$prediction_lock_sha256, sha(prediction_lock_path)) ||
        !identical(status$implementation_freeze_sha256, implementation_id) ||
        !identical(status$artifact_sha256, sha(artifact)) ||
        !identical(status$state, "completed")) {
      stop("MV5-Z immutable outcome resume failed: ", stem, call. = FALSE)
    }
    reused <- TRUE
  } else {
    started <- proc.time()[["elapsed"]]
    group_rankings <- rankings[
      rankings$robustness_group_id == unit$robustness_group_id, , drop = FALSE]
    outcome <- mv05z_evaluate_group_v1(group_rankings, labels)
    elapsed <- proc.time()[["elapsed"]] - started
    if (elapsed > 300) stop("MV5-Z outcome group elapsed cap exceeded.", call. = FALSE)
    atomic_rds(outcome, artifact)
    status <- list(
      contract_id = "mv05z_outcome_unit_status_v1",
      evaluation_unit_id = unit$evaluation_unit_id[[1L]], state = "completed",
      artifact_sha256 = sha(artifact),
      artifact_bytes = unname(file.info(artifact)$size),
      elapsed_seconds = elapsed,
      process_rss_bytes = unname(ps::ps_memory_info(process)[["rss"]]),
      prediction_lock_sha256 = sha(prediction_lock_path),
      implementation_freeze_sha256 = implementation_id,
      completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
    atomic_rds(status, status_path)
  }
  results[[index]] <- outcome
  completion[[index]] <- data.frame(
    contract_id = "mv05z_outcome_unit_completion_v1",
    evaluation_unit_id = unit$evaluation_unit_id,
    evaluation_order = unit$evaluation_order,
    robustness_group_id = unit$robustness_group_id,
    fold_id = unit$fold_id, seed = unit$seed,
    representation = unit$representation,
    state = status$state, query_method_rows = nrow(outcome),
    artifact_sha256 = status$artifact_sha256,
    artifact_bytes = status$artifact_bytes,
    labels_opened_after_prediction_lock = TRUE,
    outcomes_computed = TRUE, evaluation_executed = TRUE,
    method_selection_executed = FALSE, stringsAsFactors = FALSE)
}
completion <- do.call(rbind, completion); pc20 <- do.call(rbind, results)
pc20 <- pc20[order(pc20$representation, pc20$family_id, pc20$fold_id,
                    pc20$seed, pc20$query_sample_id, method = "radix"), ]
rownames(pc20) <- NULL
if (nrow(completion) != 150L || nrow(pc20) != 3600L ||
    !all(completion$state == "completed") ||
    any(pc20$endpoint_status != "estimable")) {
  stop("MV5-Z outcome group completion failed.", call. = FALSE)
}

# Accepted outcome files are opened only after every PC20 outcome unit is
# immutable. They are paired, not refitted or reranked.
baseline_sct <- readc(baseline_sct_path)
baseline_integrated <- readc(baseline_integrated_path)
paired <- mv05z_pair_baseline_v1(pc20, baseline_sct, baseline_integrated)
query_estimands <- mv05z_query_estimands_v1(paired)
summaries <- mv05z_summarize_estimands_v1(query_estimands)
inference <- mv05z_block_inference_v1(summaries$sample)

long_endpoints <- rbind(
  data.frame(
    contract_id = "mv05z_pc20_query_endpoint_v1",
    pc20[c("fold_id", "held_out_study", "seed", "representation", "family_id",
           "method_id", "query_sample_id", "query_tissue", "training_samples",
           "training_studies", "endpoint_status")],
    endpoint_id = .mv05z_endpoints[[1L]], value = pc20$reciprocal_rank,
    stringsAsFactors = FALSE),
  data.frame(
    contract_id = "mv05z_pc20_query_endpoint_v1",
    pc20[c("fold_id", "held_out_study", "seed", "representation", "family_id",
           "method_id", "query_sample_id", "query_tissue", "training_samples",
           "training_studies", "endpoint_status")],
    endpoint_id = .mv05z_endpoints[[2L]], value = as.numeric(pc20$one_nn_correct),
    stringsAsFactors = FALSE))
long_endpoints <- long_endpoints[order(
  match(long_endpoints$endpoint_id, .mv05z_endpoints),
  long_endpoints$representation, long_endpoints$family_id,
  long_endpoints$fold_id, long_endpoints$seed,
  long_endpoints$query_sample_id, method = "radix"), ]
rownames(long_endpoints) <- NULL
if (nrow(long_endpoints) != 7200L) stop("MV5-Z long endpoint completion failed.", call. = FALSE)

production <- data.frame(
  contract_id = "mv05z_production_summary_v1",
  prediction_groups = 150L, ranking_rows = 282800L,
  outcome_groups = nrow(completion), query_method_rows = nrow(pc20),
  long_query_endpoint_rows = nrow(long_endpoints),
  paired_query_method_rows = nrow(paired),
  query_estimand_rows = nrow(query_estimands),
  sample_estimand_rows = nrow(summaries$sample),
  tissue_estimand_rows = nrow(summaries$tissue),
  macro_estimands = nrow(summaries$macro), intervals = nrow(inference$intervals),
  primary_tests = nrow(inference$primary),
  bootstrap_replicates = 2000L, randomization_replicates = 9999L,
  representations = 2L, methods = 4L, endpoints = 2L, samples = 90L,
  studies = 15L, tissues = 5L, seeds = 5L,
  clustering_executed = FALSE, other_configurations_executed = FALSE,
  refit_or_rerank_executed = FALSE, method_selection_executed = FALSE,
  equivalence_claim_authorized = FALSE,
  outcomes_computed = TRUE, evaluation_executed = TRUE,
  stringsAsFactors = FALSE)
resource <- data.frame(
  contract_id = "mv05z_outcome_resource_summary_v1",
  outcome_groups = 150L,
  aggregate_group_seconds = sum(vapply(seq_len(nrow(completion)), function(i) {
    stem <- sub("^mv05y_eval_v1:", "", completion$evaluation_unit_id[[i]])
    readRDS(file.path(unit_root, paste0(stem, ".status.rds")))$elapsed_seconds
  }, numeric(1L))),
  maximum_group_seconds = max(vapply(seq_len(nrow(completion)), function(i) {
    stem <- sub("^mv05y_eval_v1:", "", completion$evaluation_unit_id[[i]])
    readRDS(file.path(unit_root, paste0(stem, ".status.rds")))$elapsed_seconds
  }, numeric(1L))),
  maximum_process_rss_bytes = max(vapply(seq_len(nrow(completion)), function(i) {
    stem <- sub("^mv05y_eval_v1:", "", completion$evaluation_unit_id[[i]])
    readRDS(file.path(unit_root, paste0(stem, ".status.rds")))$process_rss_bytes
  }, numeric(1L))),
  private_artifact_bytes = sum(completion$artifact_bytes),
  worker_limit = 1L, aggregate_worker_hour_cap = 2,
  process_tree_rss_cap_bytes = 4294967296,
  stringsAsFactors = FALSE)
if (resource$aggregate_group_seconds > 7200 ||
    resource$maximum_process_rss_bytes > 4294967296) {
  stop("MV5-Z outcome resource cap failed.", call. = FALSE)
}

public <- list(
  "mv05z-outcome-source-validation-2026-08-11.csv" = source_checks,
  "mv05z-outcome-implementation-freeze-2026-08-11.csv" = implementation,
  "mv05z-label-source-provenance-2026-08-11.csv" = label_provenance,
  "mv05z-outcome-unit-completion-2026-08-11.csv" = completion,
  "mv05z-pc20-query-method-outcomes-2026-08-11.csv" = pc20,
  "mv05z-pc20-query-endpoints-2026-08-11.csv" = long_endpoints,
  "mv05z-paired-query-method-endpoints-2026-08-11.csv" = paired,
  "mv05z-query-estimands-2026-08-11.csv" = query_estimands,
  "mv05z-sample-estimands-2026-08-11.csv" = summaries$sample,
  "mv05z-tissue-estimands-2026-08-11.csv" = summaries$tissue,
  "mv05z-macro-estimands-2026-08-11.csv" = summaries$macro,
  "mv05z-estimand-intervals-2026-08-11.csv" = inference$intervals,
  "mv05z-primary-contrasts-2026-08-11.csv" = inference$primary,
  "mv05z-bootstrap-audit-2026-08-11.csv" = inference$bootstrap_audit,
  "mv05z-randomization-audit-2026-08-11.csv" = inference$randomization_audit,
  "mv05z-production-summary-2026-08-11.csv" = production,
  "mv05z-outcome-resource-summary-2026-08-11.csv" = resource)
for (name in names(public)) atomic_csv(public[[name]], file.path(audit_dir, name))
atomic_rds(list(
  bootstrap_counts = inference$bootstrap_counts,
  bootstrap_estimates = inference$bootstrap_estimates,
  sign_matrix = inference$sign_matrix),
  file.path(private_root, "inference-matrices.rds"))
message("MV5-Z outcomes passed: groups=150 query_methods=3600 endpoints=7200 estimands=24 primary_tests=4")
