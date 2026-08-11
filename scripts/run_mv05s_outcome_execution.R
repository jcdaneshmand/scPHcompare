#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-S execution.", call. = FALSE)
  }
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "usage: run_mv05s_outcome_execution.R EXTERNAL_METADATA_PATH",
    "PRIVATE_ROOT AUDIT_DIR EXECUTION_SOURCE_FREEZE"), call. = FALSE)
}
metadata_path <- normalizePath(args[[1L]], mustWork = TRUE)
private_root <- args[[2L]]
audit_dir <- args[[3L]]
execution_freeze_path <- args[[4L]]
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
unit_dir <- file.path(private_root, "units")
dir.create(unit_dir, recursive = TRUE, showWarnings = FALSE)

file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
read_public <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
atomic_rds <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(temporary)) stop("MV5-S stale temporary RDS exists.",
                                   call. = FALSE)
  saveRDS(value, temporary, version = 3, compress = "xz")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("MV5-S atomic RDS rename failed.", call. = FALSE)
  }
  invisible(path)
}
atomic_csv <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  write_provenance_csv(value, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("MV5-S atomic CSV rename failed.", call. = FALSE)
  }
  invisible(path)
}

mv05r_freeze_path <- "docs/audits/mv05r-source-freeze-2026-08-10.csv"
queue_path <- "docs/audits/mv05r-evaluation-queue-2026-08-10.csv"
selected_path <- "docs/audits/mv05q-selected-training-partitions-2026-08-10.csv.gz"
heldout_path <- "docs/audits/mv05q-heldout-assignments-2026-08-10.csv.gz"
for (path in c(mv05r_freeze_path, queue_path, selected_path, heldout_path,
               execution_freeze_path)) {
  if (!file.exists(path)) stop("MV5-S required source is missing: ", path,
                               call. = FALSE)
}

mv05r_freeze <- read_public(mv05r_freeze_path)
if (nrow(mv05r_freeze) != 18L ||
    length(unique(mv05r_freeze$source_freeze_sha256)) != 1L) {
  stop("MV5-S MV5-R source freeze is malformed.", call. = FALSE)
}
source_checks <- lapply(seq_len(nrow(mv05r_freeze)), function(index) {
  row <- mv05r_freeze[index, , drop = FALSE]
  path <- if (row$external[[1L]]) metadata_path else row$artifact_locator[[1L]]
  actual <- file_sha(path)
  data.frame(
    source_id = row$source_id[[1L]], artifact_locator = row$artifact_locator[[1L]],
    expected_sha256 = row$sha256[[1L]], actual_sha256 = actual,
    matched = identical(actual, row$sha256[[1L]]),
    external = row$external[[1L]], stringsAsFactors = FALSE)
})
source_checks <- do.call(rbind, source_checks)
if (!all(source_checks$matched)) stop("MV5-S MV5-R source hash drift.",
                                      call. = FALSE)

execution_freeze <- read_public(execution_freeze_path)
if (!nrow(execution_freeze) || any(execution_freeze$external) ||
    length(unique(execution_freeze$execution_freeze_sha256)) != 1L) {
  stop("MV5-S execution source freeze is malformed.", call. = FALSE)
}
execution_checks <- lapply(seq_len(nrow(execution_freeze)), function(index) {
  row <- execution_freeze[index, , drop = FALSE]
  actual <- file_sha(row$artifact_locator[[1L]])
  data.frame(
    source_id = row$source_id[[1L]], artifact_locator = row$artifact_locator[[1L]],
    expected_sha256 = row$sha256[[1L]], actual_sha256 = actual,
    matched = identical(actual, row$sha256[[1L]]), external = FALSE,
    stringsAsFactors = FALSE)
})
execution_checks <- do.call(rbind, execution_checks)
if (!all(execution_checks$matched)) stop("MV5-S execution source hash drift.",
                                         call. = FALSE)
source_checks <- rbind(source_checks, execution_checks)

queue <- read_public(queue_path)
mv05r_validate_evaluation_queue_v1(queue)
if (nrow(queue) != 2400L || any(queue$outcomes_computed) ||
    any(queue$evaluation_executed) || any(queue$method_selection_executed) ||
    !identical(unique(queue$source_freeze_sha256),
               unique(mv05r_freeze$source_freeze_sha256))) {
  stop("MV5-S queue identity or unopened state drifted.", call. = FALSE)
}

selected <- read_public(selected_path)
heldout <- read_public(heldout_path)
if (nrow(selected) != 126000L || nrow(heldout) != 9000L ||
    any(selected$biological_outcomes_computed) ||
    any(heldout$biological_outcomes_computed) ||
    !setequal(queue$analysis_group_id, selected$analysis_group_id) ||
    !setequal(queue$analysis_group_id, heldout$analysis_group_id)) {
  stop("MV5-S MV5-Q partition/assignment axes drifted.", call. = FALSE)
}

# The label file is opened only after every frozen source and unopened queue
# identity above has passed.
raw_labels <- utils::read.csv(metadata_path, stringsAsFactors = FALSE,
                              check.names = TRUE)
required_labels <- c("orig.ident", "SRA", "Tissue.x", "Approach.x")
if (!all(required_labels %in% names(raw_labels))) {
  stop("MV5-S external label schema drifted.", call. = FALSE)
}
labels <- data.frame(
  sample_id = trimws(as.character(raw_labels$orig.ident)),
  study = trimws(as.character(raw_labels$SRA)),
  tissue = tolower(trimws(as.character(raw_labels$Tissue.x))),
  approach = trimws(as.character(raw_labels$Approach.x)),
  stringsAsFactors = FALSE)
labels <- labels[labels$tissue %in% c("bone marrow", "colon", "liver", "pbmc",
                                      "testis"), , drop = FALSE]
expected_samples <- sort(unique(c(selected$sample_id, heldout$query_sample_id)),
                         method = "radix")
if (nrow(labels) != 90L || length(unique(labels$study)) != 15L ||
    anyNA(labels) || anyDuplicated(labels$sample_id) ||
    !identical(sort(labels$sample_id, method = "radix"), expected_samples)) {
  stop("MV5-S exact 90-sample label join failed.", call. = FALSE)
}
study_tissue <- unique(labels[c("study", "tissue")])
approaches_per_study <- tapply(labels$approach, labels$study,
                               function(x) length(unique(x)))
if (anyDuplicated(study_tissue$study) ||
    sum(approaches_per_study > 1L) != 3L) {
  stop("MV5-S label block structure drifted.", call. = FALSE)
}

mv05r_source_id <- unique(mv05r_freeze$source_freeze_sha256)
execution_source_id <- unique(execution_freeze$execution_freeze_sha256)
process <- ps::ps_handle()
statuses <- vector("list", nrow(queue))
results <- vector("list", nrow(queue))

for (index in seq_len(nrow(queue))) {
  queue_row <- queue[index, , drop = FALSE]
  unit_stem <- sub("^mv05r_eval_v1:", "", queue_row$evaluation_unit_id[[1L]])
  artifact_path <- file.path(unit_dir, paste0(unit_stem, ".rds"))
  status_path <- file.path(unit_dir, paste0(unit_stem, ".status.rds"))
  artifact_exists <- file.exists(artifact_path)
  status_exists <- file.exists(status_path)
  reused <- FALSE
  if (xor(artifact_exists, status_exists)) {
    stop("MV5-S partial artifact/status pair detected for ", unit_stem,
         call. = FALSE)
  }
  if (artifact_exists && status_exists) {
    status <- readRDS(status_path)
    if (!identical(status$evaluation_unit_id,
                   queue_row$evaluation_unit_id[[1L]]) ||
        !identical(status$mv05r_source_freeze_sha256, mv05r_source_id) ||
        !identical(status$execution_freeze_sha256, execution_source_id) ||
        !identical(status$artifact_sha256, file_sha(artifact_path)) ||
        !identical(status$state, "completed")) {
      stop("MV5-S immutable resume validation failed for ", unit_stem,
           call. = FALSE)
    }
    result <- readRDS(artifact_path)
    reused <- TRUE
  } else {
    started <- proc.time()[["elapsed"]]
    result <- if (queue_row$evaluation_scope[[1L]] ==
                  "overlapping_training_partition_alignment") {
      mv05s_evaluate_training_unit_v1(queue_row, selected, labels)
    } else {
      mv05s_evaluate_heldout_unit_v1(queue_row, selected, heldout, labels)
    }
    elapsed <- proc.time()[["elapsed"]] - started
    if (elapsed > 300) stop("MV5-S per-unit elapsed cap exceeded.", call. = FALSE)
    atomic_rds(result, artifact_path)
    status <- list(
      contract_id = "mv05s_unit_status_v1",
      evaluation_unit_id = queue_row$evaluation_unit_id[[1L]],
      state = "completed", artifact_sha256 = file_sha(artifact_path),
      artifact_bytes = unname(file.info(artifact_path)$size),
      elapsed_seconds = elapsed,
      process_rss_bytes = unname(ps::ps_memory_info(process)[["rss"]]),
      mv05r_source_freeze_sha256 = mv05r_source_id,
      execution_freeze_sha256 = execution_source_id,
      completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
    atomic_rds(status, status_path)
  }
  statuses[[index]] <- data.frame(
    evaluation_unit_id = status$evaluation_unit_id,
    execution_order = queue_row$execution_order[[1L]],
    analysis_group_id = queue_row$analysis_group_id[[1L]],
    algorithm_id = queue_row$algorithm_id[[1L]],
    endpoint_id = queue_row$endpoint_id[[1L]],
    state = status$state, artifact_sha256 = status$artifact_sha256,
    artifact_bytes = status$artifact_bytes,
    elapsed_seconds = status$elapsed_seconds,
    process_rss_bytes = status$process_rss_bytes,
    reused = reused, outcomes_computed = TRUE, evaluation_executed = TRUE,
    method_selection_executed = FALSE, stringsAsFactors = FALSE)
  results[[index]] <- result
}

statuses <- do.call(rbind, statuses)
if (nrow(statuses) != 2400L || anyDuplicated(statuses$evaluation_unit_id) ||
    !all(statuses$state == "completed") ||
    sum(statuses$elapsed_seconds[!statuses$reused]) > 7200 ||
    max(statuses$process_rss_bytes) > 4294967296) {
  stop("MV5-S completion or resource envelope failed.", call. = FALSE)
}

is_training <- queue$evaluation_scope == "overlapping_training_partition_alignment"
training_seed <- do.call(rbind, lapply(results[is_training], `[[`, "seed"))
training_summary <- do.call(rbind, lapply(results[is_training], `[[`, "summary"))
heldout_seed <- do.call(rbind, lapply(results[!is_training], `[[`, "seed"))
heldout_summary <- do.call(rbind, lapply(results[!is_training], `[[`, "summary"))
heldout_private <- do.call(rbind, lapply(results[!is_training], `[[`, "private"))

context_columns <- c("representation", "distance_id", "algorithm_id",
                     "algorithm_role", "endpoint_id", "label_axis")
training_context <- do.call(rbind, lapply(
  split(training_summary, interaction(training_summary[context_columns],
                                      drop = TRUE, lex.order = TRUE)),
  function(part) {
    values <- part$seed_mean[part$status == "completed"]
    data.frame(
      representation = part$representation[[1L]],
      distance_id = part$distance_id[[1L]],
      algorithm_id = part$algorithm_id[[1L]],
      algorithm_role = part$algorithm_role[[1L]],
      endpoint_id = part$endpoint_id[[1L]],
      label_axis = part$label_axis[[1L]],
      fold_mean = mean(values), fold_median = stats::median(values),
      fold_minimum = min(values), fold_maximum = max(values),
      completed_folds = length(values), expected_folds = 15L,
      inference_scope = "descriptive_overlapping_training_folds",
      p_value_computed = FALSE, method_selection_executed = FALSE,
      stringsAsFactors = FALSE)
  }))

sample_context_columns <- c("representation", "distance_id", "algorithm_id",
                            "algorithm_role", "endpoint_id", "label_axis",
                            "query_sample_id", "heldout_study", "true_label")
heldout_sample <- do.call(rbind, lapply(
  split(heldout_private, interaction(heldout_private[sample_context_columns],
                                     drop = TRUE, lex.order = TRUE)),
  function(part) {
    if (nrow(part) != 5L || length(unique(part$seed)) != 5L) {
      stop("MV5-S held-out sample seed aggregation is incomplete.",
           call. = FALSE)
    }
    data.frame(
      representation = part$representation[[1L]],
      distance_id = part$distance_id[[1L]],
      algorithm_id = part$algorithm_id[[1L]],
      algorithm_role = part$algorithm_role[[1L]],
      endpoint_id = part$endpoint_id[[1L]],
      label_axis = part$label_axis[[1L]],
      sample_id = part$query_sample_id[[1L]], study = part$heldout_study[[1L]],
      true_label = part$true_label[[1L]], seed_mean_correct = mean(part$correct),
      completed_seeds = 5L, stringsAsFactors = FALSE)
  }))

tissue_counts <- mv05s_bootstrap_counts_v1(labels, "tissue", 2000L, 20260810L)
approach_counts <- mv05s_bootstrap_counts_v1(labels, "approach", 2000L, 20260810L)
heldout_macro <- list()
bootstrap_audit <- list()
heldout_contexts <- split(
  heldout_sample,
  interaction(heldout_sample[context_columns], drop = TRUE, lex.order = TRUE))
for (index in seq_along(heldout_contexts)) {
  part <- heldout_contexts[[index]]
  axis <- part$label_axis[[1L]]
  if (nrow(part) != 90L || anyDuplicated(part$sample_id)) {
    stop("MV5-S held-out 90-sample context is incomplete.", call. = FALSE)
  }
  names(part)[names(part) == "true_label"] <- axis
  score <- mv05s_balanced_accuracy_v1(part$seed_mean_correct, part[[axis]])
  counts <- if (axis == "tissue") tissue_counts else approach_counts
  bootstrap <- mv05s_apply_bootstrap_v1(part, axis, counts)
  heldout_macro[[index]] <- data.frame(
    representation = part$representation[[1L]],
    distance_id = part$distance_id[[1L]],
    algorithm_id = part$algorithm_id[[1L]],
    algorithm_role = part$algorithm_role[[1L]],
    endpoint_id = part$endpoint_id[[1L]], label_axis = axis,
    estimate = score$estimate, ci_lower = bootstrap$ci_lower,
    ci_upper = bootstrap$ci_upper,
    samples = nrow(part), studies = length(unique(part$study)),
    label_classes = length(score$class_values), completed_seeds = 5L,
    bootstrap_replicates = bootstrap$requested_replicates,
    estimable_bootstrap_replicates = bootstrap$estimable_replicates,
    bootstrap_seed = 20260810L, interval_type = "percentile_type7",
    status = bootstrap$status, p_value_computed = FALSE,
    method_selection_executed = FALSE, stringsAsFactors = FALSE)
  bootstrap_audit[[index]] <- data.frame(
    representation = part$representation[[1L]],
    distance_id = part$distance_id[[1L]],
    algorithm_id = part$algorithm_id[[1L]],
    endpoint_id = part$endpoint_id[[1L]], label_axis = axis,
    resampling_unit = if (axis == "tissue")
      "tissue_stratified_study_block" else "global_study_block_preserve_mixed",
    block_count_matrix_sha256 = .mv05s_digest(counts),
    replicate_estimates_sha256 = .mv05s_digest(bootstrap$estimates),
    bootstrap_replicates = bootstrap$requested_replicates,
    estimable_replicates = bootstrap$estimable_replicates,
    seeds_treated_as_independent = FALSE, stringsAsFactors = FALSE)
}
heldout_macro <- do.call(rbind, heldout_macro)
bootstrap_audit <- do.call(rbind, bootstrap_audit)

public_objects <- list(
  "mv05s-source-validation-2026-08-10.csv" = source_checks,
  "mv05s-unit-completion-2026-08-10.csv" = statuses,
  "mv05s-training-seed-outcomes-2026-08-10.csv" = training_seed,
  "mv05s-training-unit-summary-2026-08-10.csv" = training_summary,
  "mv05s-training-context-summary-2026-08-10.csv" = training_context,
  "mv05s-heldout-seed-outcomes-2026-08-10.csv" = heldout_seed,
  "mv05s-heldout-unit-summary-2026-08-10.csv" = heldout_summary,
  "mv05s-heldout-macro-outcomes-2026-08-10.csv" = heldout_macro,
  "mv05s-bootstrap-audit-2026-08-10.csv" = bootstrap_audit)

prohibited_columns <- c("true_label", "predicted_label")
label_values <- unique(c(labels$tissue, labels$approach))
for (name in names(public_objects)) {
  value <- public_objects[[name]]
  if (any(names(value) %in% prohibited_columns) ||
      any(vapply(value, function(column) {
        any(tolower(as.character(column)) %in% tolower(label_values))
      }, logical(1L)))) {
    stop("MV5-S public label safety failed for ", name, call. = FALSE)
  }
  atomic_csv(value, file.path(audit_dir, name))
}

resource <- data.frame(
  contract_id = "mv05s_resource_summary_v1",
  evaluation_units = nrow(statuses), reused_units = sum(statuses$reused),
  executed_units = sum(!statuses$reused),
  aggregate_unit_seconds = sum(statuses$elapsed_seconds[!statuses$reused]),
  maximum_unit_seconds = max(statuses$elapsed_seconds),
  maximum_process_rss_bytes = max(statuses$process_rss_bytes),
  private_unit_artifact_bytes = sum(statuses$artifact_bytes),
  public_output_bytes = sum(file.info(file.path(audit_dir,
    names(public_objects)))$size), worker_limit = 1L,
  outcomes_computed = TRUE, evaluation_executed = TRUE,
  method_selection_executed = FALSE, stringsAsFactors = FALSE)
if (resource$public_output_bytes > 1073741824) {
  stop("MV5-S public output storage cap exceeded.", call. = FALSE)
}
atomic_csv(resource, file.path(audit_dir,
                               "mv05s-resource-summary-2026-08-10.csv"))
atomic_rds(heldout_sample, file.path(private_root,
                                     "heldout-sample-summary-private.rds"))

message("MV5-S execution passed: units=", nrow(statuses),
        " executed=", sum(!statuses$reused), " reused=", sum(statuses$reused),
        " training_seed_rows=", nrow(training_seed),
        " heldout_seed_rows=", nrow(heldout_seed),
        " heldout_macro_rows=", nrow(heldout_macro),
        " p_values=0 method_selection=0")
