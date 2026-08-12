#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-AL prediction locking.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste("usage: build_mv05al_prediction_lock.R MV05AJ_RESULT_ROOT",
             "EXTERNAL_METADATA_PATH PRIVATE_ROOT AUDIT_DIR EXPECTED_HEAD"),
       call. = FALSE)
}
x_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
metadata_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
private_root <- args[[3L]]; audit_dir <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
unit_root <- file.path(private_root, "units")
dir.create(unit_root, recursive = TRUE, showWarnings = FALSE)

file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
atomic_rds <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(temporary)) stop("Stale MV5-AL temporary artifact.", call. = FALSE)
  saveRDS(value, temporary, version = 3, compress = "xz")
  if (!file.rename(temporary, path)) {
    unlink(temporary); stop("MV5-AL atomic rename failed.", call. = FALSE)
  }
}
atomic_csv <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  write_provenance_csv(value, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary); stop("MV5-AL atomic CSV rename failed.", call. = FALSE)
  }
}
safe <- function(value) gsub("[^A-Za-z0-9._-]", "_", value)

if (!grepl("^[0-9a-f]{40}$", expected_head)) {
  stop("EXPECTED_HEAD must be a full 40-character Git commit identity.",
       call. = FALSE)
}
head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-AL prediction lock requires engine HEAD ", expected_head,
       "; observed ", head, ".", call. = FALSE)
}

y_freeze_path <- "docs/audits/mv05ak-source-freeze-2026-08-11.csv"
y_queue_path <- "docs/audits/mv05ak-evaluation-queue-2026-08-11.csv"
y_axis_path <- "docs/audits/mv05ak-prediction-axis-compatibility-2026-08-11.csv"
for (path in c(y_freeze_path, y_queue_path, y_axis_path)) {
  if (!file.exists(path)) stop("Missing accepted MV5-AK source: ", path, call. = FALSE)
}
y_freeze <- read_csv(y_freeze_path)
y_queue <- read_csv(y_queue_path)
axis <- read_csv(y_axis_path)
if (nrow(y_freeze) != 196L || nrow(y_queue) != 150L || nrow(axis) != 8L ||
    !all(axis$fold_seed_query_reference_axis_exact) ||
    any(y_queue$configuration_id != "nested_cells_256_pc30_euclidean_v1") ||
    any(!y_queue$coordinate_source_identity_exact) ||
    any(!y_queue$nested192_prefix_identity_exact) ||
    any(!y_queue$nested256_subset_384_identity_exact) ||
    any(y_queue$execution_authorized) || any(y_queue$ranking_executed) ||
    any(y_queue$outcomes_computed) || any(y_queue$evaluation_executed)) {
  stop("Accepted MV5-AK prefreeze axes or unopened state drifted.", call. = FALSE)
}

source_checks <- lapply(seq_len(nrow(y_freeze)), function(index) {
  row <- y_freeze[index, , drop = FALSE]
  locator <- row$artifact_locator[[1L]]
  path <- if (row$source_id[[1L]] == "external_label_source") {
    metadata_path
  } else if (startsWith(locator, "ignored_tmp:")) {
    base <- sub("^ignored_tmp:", "", locator)
    if (grepl("^nested256_private_group_manifest_", row$source_id[[1L]]))
      file.path(base, "artifact_manifest.csv") else base
  } else locator
  if (!file.exists(path)) stop("MV5-AL bound source missing: ", locator, call. = FALSE)
  actual <- file_sha(path)
  data.frame(
    contract_id = "mv05al_source_validation_v1",
    source_id = row$source_id, artifact_locator = locator,
    expected_sha256 = row$sha256, actual_sha256 = actual,
    matched = identical(actual, row$sha256[[1L]]),
    opened_as_outcome_source = FALSE,
    labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
})
source_checks <- do.call(rbind, source_checks)
if (!all(source_checks$matched)) stop("MV5-AL frozen source hash drift.", call. = FALSE)

implementation_paths <- c(
  engine = "R/mv05al_outcome_execution.R",
  mv05ak_engine = "R/mv05ak_nested256_outcome_prefreeze.R",
  specification = "docs/specifications/MV05AL_NESTED256_RETRIEVAL_ROBUSTNESS_EXECUTION_SPECIFICATION_V1.md",
  ranking_runner = "scripts/build_mv05al_prediction_lock.R",
  ranking_validator = "scripts/validate_mv05al_prediction_lock.R",
  tests = "tests/testthat/test-mv05al-outcome-execution.R")
implementation <- data.frame(
  contract_id = "mv05al_prediction_implementation_freeze_v1",
  source_id = names(implementation_paths),
  artifact_locator = unname(implementation_paths),
  sha256 = vapply(implementation_paths, file_sha, character(1L)),
  bytes = as.numeric(file.info(implementation_paths)$size),
  accepted_head = expected_head, labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE)
implementation$implementation_freeze_sha256 <- .mv05al_digest(
  paste(implementation$artifact_locator, implementation$sha256, sep = "\r"))
implementation_id <- unique(implementation$implementation_freeze_sha256)
y_source_id <- unique(y_freeze$source_freeze_sha256)

process <- ps::ps_handle()
results <- vector("list", nrow(y_queue)); completion <- vector("list", nrow(y_queue))
for (index in seq_len(nrow(y_queue))) {
  unit <- y_queue[index, , drop = FALSE]
  stem <- sub("^mv05ak_eval_v1:", "", unit$evaluation_unit_id[[1L]])
  artifact_path <- file.path(unit_root, paste0(stem, ".rds"))
  status_path <- file.path(unit_root, paste0(stem, ".status.rds"))
  both <- file.exists(artifact_path) && file.exists(status_path)
  if (xor(file.exists(artifact_path), file.exists(status_path))) {
    stop("MV5-AL partial ranking artifact/status pair: ", stem, call. = FALSE)
  }
  reused <- FALSE
  if (both) {
    status <- readRDS(status_path); ranking <- readRDS(artifact_path)
    if (!identical(status$evaluation_unit_id, unit$evaluation_unit_id[[1L]]) ||
        !identical(status$mv05ak_source_freeze_sha256, y_source_id) ||
        !identical(status$implementation_freeze_sha256, implementation_id) ||
        !identical(status$artifact_sha256, file_sha(artifact_path)) ||
        !identical(status$state, "completed")) {
      stop("MV5-AL immutable ranking resume failed: ", stem, call. = FALSE)
    }
    reused <- TRUE
  } else {
    started <- proc.time()[["elapsed"]]
    group_dir <- file.path(x_root, safe(unit$robustness_group_id[[1L]]))
    methods <- read_csv(file.path(group_dir, "method_rows.csv"))
    ranking <- mv05al_rank_group_v1(methods, unit)
    elapsed <- proc.time()[["elapsed"]] - started
    if (elapsed > 300) stop("MV5-AL ranking group elapsed cap exceeded.", call. = FALSE)
    atomic_rds(ranking, artifact_path)
    status <- list(
      contract_id = "mv05al_ranking_unit_status_v1",
      evaluation_unit_id = unit$evaluation_unit_id[[1L]], state = "completed",
      artifact_sha256 = file_sha(artifact_path),
      artifact_bytes = unname(file.info(artifact_path)$size),
      elapsed_seconds = elapsed,
      process_rss_bytes = unname(ps::ps_memory_info(process)[["rss"]]),
      mv05ak_source_freeze_sha256 = y_source_id,
      implementation_freeze_sha256 = implementation_id,
      completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
    atomic_rds(status, status_path)
  }
  results[[index]] <- ranking
  completion[[index]] <- data.frame(
    contract_id = "mv05al_ranking_unit_completion_v1",
    evaluation_unit_id = unit$evaluation_unit_id,
    evaluation_order = unit$evaluation_order,
    robustness_group_id = unit$robustness_group_id,
    fold_id = unit$fold_id, seed = unit$seed,
    representation = unit$representation,
    state = status$state, ranking_rows = nrow(ranking),
    artifact_sha256 = status$artifact_sha256,
    artifact_bytes = status$artifact_bytes,
    elapsed_seconds = status$elapsed_seconds,
    process_rss_bytes = status$process_rss_bytes,
    reused = reused, labels_opened = FALSE, outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
completion <- do.call(rbind, completion)
rankings <- do.call(rbind, results)
rankings <- rankings[order(rankings$representation, rankings$fold_id,
                            rankings$seed, rankings$method_id,
                            rankings$query_sample_id, rankings$neighbor_rank,
                            method = "radix"), , drop = FALSE]
rownames(rankings) <- NULL
if (nrow(rankings) != 282800L || anyDuplicated(rankings$ranking_id) ||
    any(rankings$labels_opened) || any(rankings$outcomes_computed) ||
    any(!rankings$prediction_locked) || nrow(completion) != 150L ||
    !all(completion$state == "completed") ||
    max(completion$process_rss_bytes) > 4294967296 ||
    sum(completion$elapsed_seconds[!completion$reused]) > 7200) {
  stop("MV5-AL prediction ranking completion/resource contract failed.", call. = FALSE)
}

ranking_path <- file.path(audit_dir, "mv05al-nested256-rankings-2026-08-11.csv.gz")
temporary <- paste0(ranking_path, ".tmp-", Sys.getpid())
connection <- gzfile(temporary, open = "wt", compression = 9)
utils::write.csv(rankings, connection, row.names = FALSE, na = "")
close(connection)
if (!file.rename(temporary, ranking_path)) stop("Atomic ranking publication failed.", call. = FALSE)
ranking_sha <- file_sha(ranking_path)
prediction_lock <- data.frame(
  contract_id = "mv05al_prediction_lock_v1",
  accepted_mv05ak_commit = "9bcf650a1f4b5910403cc32ecad41c5bccc4a9d1",
  engine_commit = expected_head,
  mv05ak_source_freeze_sha256 = y_source_id,
  implementation_freeze_sha256 = implementation_id,
  ranking_sha256 = ranking_sha, ranking_rows = nrow(rankings),
  ranking_groups = nrow(completion), completed_groups = sum(completion$state == "completed"),
  failed_groups = sum(completion$state != "completed"), methods = 4L,
  representations = 2L, seeds = 5L, folds = 15L,
  ranking_tie_policy =
    "ascending_distance_then_canonical_training_sample_id_radix_v1",
  label_open_authorized = TRUE,
  prediction_lock_status = "passed_before_label_open",
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
resource <- data.frame(
  contract_id = "mv05al_prediction_resource_summary_v1",
  groups = nrow(completion), executed_groups = sum(!completion$reused),
  reused_groups = sum(completion$reused), ranking_rows = nrow(rankings),
  aggregate_group_seconds = sum(completion$elapsed_seconds[!completion$reused]),
  maximum_group_seconds = max(completion$elapsed_seconds),
  maximum_process_rss_bytes = max(completion$process_rss_bytes),
  private_artifact_bytes = sum(completion$artifact_bytes),
  public_ranking_bytes = as.numeric(file.info(ranking_path)$size),
  worker_limit = 1L, labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv05al_prediction_lock_decision_v1",
  all_196_sources_valid = all(source_checks$matched),
  all_150_groups_complete = all(completion$state == "completed"),
  ranking_rows = nrow(rankings), canonical_ties_valid = TRUE,
  label_open_authorized = TRUE,
  authorization_scope = "mv05al_nested256_tissue_retrieval_robustness_only",
  clustering_authorized = FALSE, other_configuration_authorized = FALSE,
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE)

atomic_csv(source_checks, file.path(audit_dir,
  "mv05al-prediction-source-validation-2026-08-11.csv"))
atomic_csv(implementation, file.path(audit_dir,
  "mv05al-prediction-implementation-freeze-2026-08-11.csv"))
atomic_csv(completion, file.path(audit_dir,
  "mv05al-ranking-unit-completion-2026-08-11.csv"))
atomic_csv(prediction_lock, file.path(audit_dir,
  "mv05al-prediction-lock-2026-08-11.csv"))
atomic_csv(resource, file.path(audit_dir,
  "mv05al-prediction-resource-summary-2026-08-11.csv"))
atomic_csv(decision, file.path(audit_dir,
  "mv05al-prediction-lock-decision-2026-08-11.csv"))
message("MV5-AL prediction lock passed: groups=150 rankings=282800 labels=0 outcomes=0")
