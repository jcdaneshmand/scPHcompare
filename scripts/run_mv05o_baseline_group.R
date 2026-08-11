#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 15L) {
  stop(paste(
    "usage: run_mv05o_baseline_group.R GROUP_REQUESTS GROUP_QUEUE BASELINE_QUEUE",
    "D1_RESOURCE D1_CACHE_ROOT MEAN_AUDIT MEAN_ROOT G_MANIFEST G_RESOURCE",
    "G_RESULT_ROOT OUTPUT_DIR STATUS_DIR MV05N_SOURCE IMPLEMENTATION_SHA",
    "SOURCE_FREEZE_SHA"
  ), call. = FALSE)
}
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05d5_retrieval_inputs.R")
source("R/mv05f_integration_gate.R")
source("R/mv05h_integrated_ph_production.R")
source(args[[13L]])
read_public <- function(path) utils::read.csv(normalizePath(path, mustWork = TRUE),
                                               stringsAsFactors = FALSE,
                                               check.names = FALSE)
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
write_atomic <- function(value, path) {
  temporary <- paste0(path, ".partial.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write_provenance_csv(value, temporary)
  if (file.exists(path)) stop("Refusing to overwrite ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Atomic rename failed.", call. = FALSE)
}
implementation_sha <- args[[14L]]
source_freeze_sha <- args[[15L]]
if (file_sha("scripts/run_mv05o_baseline_group.R") != implementation_sha ||
    !grepl("^[0-9a-f]{64}$", source_freeze_sha)) {
  stop("MV5-O baseline implementation/source hash is stale.", call. = FALSE)
}
requests <- read_public(args[[1L]])
groups <- read_public(args[[2L]])
baseline_queue <- read_public(args[[3L]])
group_ids <- unique(requests$production_group_id)
if (length(group_ids) != 1L || any(requests$outcome_label_state != "closed") ||
    any(as.logical(requests$biological_outcomes_computed)) ||
    any(c("tissue", "approach", "class", "label", "outcome") %in%
          tolower(names(requests)))) {
  stop("MV5-O baseline group request boundary is invalid.", call. = FALSE)
}
group <- groups[groups$production_group_id == group_ids, , drop = FALSE]
units <- baseline_queue[baseline_queue$production_group_id == group_ids, , drop = FALSE]
if (nrow(group) != 1L || nrow(units) != if (group$representation == "sct_whole") 2L else 1L) {
  stop("MV5-O baseline queue does not match its production group.", call. = FALSE)
}
pairs <- unique(requests[requests$homology_dimension == "H0",
                         c("pair_request_id", "pair_ordinal", "first_sample_id",
                           "second_sample_id"), drop = FALSE])
if (nrow(pairs) != group$unordered_training_pairs ||
    anyDuplicated(pairs[c("first_sample_id", "second_sample_id")])) {
  stop("MV5-O baseline pair set is incomplete.", call. = FALSE)
}
d1 <- read_public(args[[4L]])
d1_root <- normalizePath(args[[5L]], mustWork = TRUE)
mean_audit <- read_public(args[[6L]])
mean_root <- normalizePath(args[[7L]], mustWork = TRUE)
g_manifest <- read_public(args[[8L]])
g_resource <- read_public(args[[9L]])
g_root <- normalizePath(args[[10L]], mustWork = TRUE)
seed <- as.integer(group$seed)
held_out_study <- sub("^large_loso_v1:", "", group$fold_id)
d1_row <- d1[d1$held_out_study == held_out_study & d1$seed == seed, ]
if (nrow(d1_row) != 1L || d1_row$exit_status != 0L) {
  stop("MV5-O D1 fold source is not accepted.", call. = FALSE)
}
fold_path <- file.path(d1_root, d1_row$private_cache_file)
if (!file.exists(fold_path) || file_sha(fold_path) != d1_row$private_cache_sha256) {
  stop("MV5-O D1 fold source is missing or stale.", call. = FALSE)
}
fold <- readRDS(fold_path)
mv05d1_validate_cell_fold_record_v1(fold)
mean_row <- mean_audit[mean_audit$seed == seed, ]
if (nrow(mean_row) != 1L) stop("MV5-O mean profile is unavailable.", call. = FALSE)
mean_path <- file.path(mean_root, mean_row$private_file)
if (!file.exists(mean_path) || file_sha(mean_path) != mean_row$private_file_sha256) {
  stop("MV5-O mean profile is missing or stale.", call. = FALSE)
}
mean_bundle <- readRDS(mean_path)
mv05d5_validate_mean_profile_bundle_v1(mean_bundle)
vectors <- mv05d5_pseudobulk_vectors_v1(mean_bundle$payload$profiles, fold)
if (group$representation == "sct_whole") {
  views <- fold$payload$cell_views
  energy_source <- fold$cache_key
} else {
  g_row <- g_manifest[g_manifest$held_out_study == held_out_study &
                        g_manifest$seed == seed, ]
  g_metric <- g_resource[g_resource$held_out_study == held_out_study &
                           g_resource$seed == seed, ]
  safe_group <- safe_name(g_row$group_id)
  g_path <- file.path(g_root, safe_group, paste0(safe_group, ".rds"))
  if (nrow(g_row) != 1L || nrow(g_metric) != 1L || !file.exists(g_path) ||
      file_sha(g_path) != g_metric$private_result_sha256) {
    stop("MV5-O integrated coordinates are missing or stale.", call. = FALSE)
  }
  record <- readRDS(g_path)
  mv05f_validate_group_record_v1(record)
  ids <- sort(names(record$payload$coordinates), method = "radix")
  views <- lapply(ids, function(sample_id) mv05h_new_integrated_cell_view_v1(
    record$payload$coordinates[[sample_id]], sample_id, fold$identity$fold_id,
    fold$identity$fit_scope_id, seed, record$cache_key,
    record$payload$coordinate_set_sha256
  ))
  names(views) <- ids
  energy_source <- record$cache_key
}
dir.create(args[[11L]], recursive = TRUE, showWarnings = FALSE)
dir.create(args[[12L]], recursive = TRUE, showWarnings = FALSE)
for (index in seq_len(nrow(units))) {
  unit <- units[index, , drop = FALSE]
  stem <- safe_name(unit$baseline_unit_id)
  output_path <- file.path(args[[11L]], paste0(stem, ".csv"))
  status_path <- file.path(args[[12L]], paste0(stem, "__status.csv"))
  if (file.exists(output_path) || file.exists(status_path)) {
    if (!file.exists(output_path) || !file.exists(status_path)) {
      stop("MV5-O baseline has a partial output/status pair.", call. = FALSE)
    }
    status <- read_public(status_path)
    if (nrow(status) != 1L || status$status != "completed" ||
        status$output_sha256 != file_sha(output_path) ||
        status$implementation_sha256 != implementation_sha ||
        status$source_freeze_sha256 != source_freeze_sha) {
      stop("MV5-O baseline resume validation failed.", call. = FALSE)
    }
    next
  }
  started <- proc.time()[["elapsed"]]
  if (unit$baseline_method == "cell_distribution_energy_v1") {
    values <- mv05n_training_energy_pairs_v1(views, pairs)
    source_identity <- energy_source
  } else {
    values <- mv05n_training_pseudobulk_pairs_v1(vectors, pairs)
    source_identity <- mean_bundle$cache_key
  }
  elapsed <- proc.time()[["elapsed"]] - started
  if (elapsed > 900 || any(!is.finite(values$distance)) ||
      any(values$distance < 0)) {
    stop("MV5-O baseline numerical/resource guard failed.", call. = FALSE)
  }
  output <- data.frame(
    contract_id = "mv05o_baseline_group_output_v1",
    baseline_unit_id = unit$baseline_unit_id,
    production_group_id = group$production_group_id,
    fold_id = group$fold_id, seed = seed,
    representation = group$representation,
    method_id = unit$baseline_method,
    pair_request_id = values$pair_request_id,
    first_sample_id = values$first_sample_id,
    second_sample_id = values$second_sample_id,
    distance = values$distance, source_identity = source_identity,
    implementation_sha256 = implementation_sha,
    source_freeze_sha256 = source_freeze_sha,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, stringsAsFactors = FALSE
  )
  write_atomic(output, output_path)
  status <- data.frame(
    contract_id = "mv05o_baseline_group_status_v1",
    baseline_unit_id = unit$baseline_unit_id, status = "completed",
    pair_rows = nrow(output), elapsed_seconds = elapsed,
    output_sha256 = file_sha(output_path),
    output_size_bytes = unname(file.info(output_path)$size),
    implementation_sha256 = implementation_sha,
    source_freeze_sha256 = source_freeze_sha,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, stringsAsFactors = FALSE
  )
  write_atomic(status, status_path)
}
message("Completed or reused MV5-O baselines for ", group$production_group_id, ".")
