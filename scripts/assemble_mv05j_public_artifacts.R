#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: assemble_mv05j_public_artifacts.R BUNDLE_DIR AUDIT_DIR ",
    "METRICS_CSV OUTPUT_DIR BUILD_DIR", call. = FALSE
  )
}
bundle_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
audit_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
metrics_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[4L]]
build_dir <- args[[5L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/mv05d5_retrieval_inputs.R")
source("R/mv05j_integrated_retrieval_inputs.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
metrics <- utils::read.csv(
  metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(metrics) != 75L || any(metrics$disposition != "completed")) {
  stop("MV5-J public assembly requires 75 accepted production groups.",
       call. = FALSE)
}

pair_rows <- list()
completion_rows <- list()
index_rows <- list()
scale_rows <- list()
for (index in seq_len(nrow(metrics))) {
  metric <- metrics[index, , drop = FALSE]
  stem <- paste0(gsub("[^A-Za-z0-9_.-]", "_", metric$held_out_study),
                 "__", metric$seed)
  bundle_path <- file.path(bundle_dir, paste0(stem, "__retrieval.rds"))
  audit_path <- file.path(audit_dir, paste0(stem, "__audit.csv"))
  bundle <- readRDS(bundle_path)
  mv05j_validate_group_bundle_v1(bundle)
  audit <- utils::read.csv(
    audit_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  pair_rows[[index]] <- bundle$payload$pairs
  completion_rows[[index]] <- bundle$payload$completion
  index_rows[[index]] <- data.frame(
    contract_id = "mv05j_public_group_index_v1",
    group_order = metric$group_order,
    group_id = bundle$identity$group_id, fold_id = bundle$identity$fold_id,
    fit_scope_id = bundle$identity$fit_scope_id,
    seed = bundle$identity$seed,
    query_samples = audit$query_samples,
    training_samples = audit$training_samples,
    biological_pairs = audit$biological_pairs,
    retrieval_rows = audit$retrieval_rows,
    completed_methods = audit$completed_methods,
    failed_methods = audit$failed_methods,
    bundle_cache_key = bundle$cache_key,
    payload_sha256 = bundle$payload_sha256,
    private_bundle_sha256 = file_sha(bundle_path),
    fold_cache_key = bundle$identity$fold_cache_key,
    mean_profile_cache_key = bundle$identity$mean_profile_cache_key,
    mv05g_group_cache_key = bundle$identity$mv05g_group_cache_key,
    mv05g_payload_sha256 = bundle$identity$mv05g_payload_sha256,
    mv05g_coordinate_set_sha256 =
      bundle$identity$mv05g_coordinate_set_sha256,
    mv05g_private_file_sha256 =
      bundle$identity$mv05g_private_file_sha256,
    mv05i_group_sha256 = bundle$identity$mv05i_group_sha256,
    implementation_sha256 = bundle$identity$implementation_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  scale_rows[[index]] <- data.frame(
    contract_id = "mv05j_component_scale_disposition_v1",
    group_id = bundle$identity$group_id, fold_id = bundle$identity$fold_id,
    fit_scope_id = bundle$identity$fit_scope_id,
    seed = bundle$identity$seed,
    method_id = c(
      "integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1",
      "integrated_cell_landscape_h0_h1_raw_euclidean_v1"
    ),
    primary_component = c(TRUE, TRUE, FALSE),
    within_training_topology_pairs_available = 0L,
    held_out_query_pairs_used_for_scale = 0L,
    scale_value = NA_real_,
    scale_disposition = c(
      "not_required_for_within_method_rank",
      "not_required_for_within_method_rank",
      "not_fitted_training_pair_scope_absent"
    ),
    resulting_distance = c(
      "raw_h0", "raw_h1", "raw_euclidean_secondary_only"
    ),
    additional_topology_jobs_executed = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
pairs <- do.call(rbind, pair_rows)
completion <- do.call(rbind, completion_rows)
group_index <- do.call(rbind, index_rows)
scale_disposition <- do.call(rbind, scale_rows)
registry <- mv05j_method_registry_v1()
registry$contract_id <- "mv05j_method_registry_v1"
registry$outcome_label_state <- "closed"
registry$biological_outcomes_computed <- FALSE
registry <- registry[, c(
  "contract_id", setdiff(names(registry), "contract_id")
), drop = FALSE]

if (nrow(pairs) != 176750L || anyDuplicated(pairs$pair_id) ||
    nrow(completion) != 375L || any(completion$status != "completed") ||
    nrow(group_index) != 75L || nrow(scale_disposition) != 225L ||
    any(c("tissue", "approach") %in% names(pairs)) ||
    any(pairs$outcome_label_state != "closed") ||
    any(as.logical(pairs$biological_outcomes_computed))) {
  stop("MV5-J public assembly violates count or label invariants.",
       call. = FALSE)
}
pairs <- pairs[order(
  pairs$fold_id, pairs$seed, pairs$method_id, pairs$query_sample_id,
  pairs$neighbor_rank, method = "radix"
), , drop = FALSE]
completion <- completion[order(
  completion$fold_id, completion$seed, completion$method_id,
  method = "radix"
), , drop = FALSE]
group_index <- group_index[order(group_index$group_order), , drop = FALSE]
scale_disposition <- scale_disposition[order(
  scale_disposition$fold_id, scale_disposition$seed,
  scale_disposition$method_id, method = "radix"
), , drop = FALSE]
rownames(pairs) <- rownames(completion) <- rownames(group_index) <-
  rownames(scale_disposition) <- NULL

write_plain <- function(data, name) {
  path <- file.path(output_dir, name)
  if (file.exists(path)) {
    stop("Refusing to overwrite existing public artifact: ", name,
         call. = FALSE)
  }
  write_provenance_csv(data, path)
  path
}
write_gzip <- function(data, name) {
  final <- file.path(output_dir, name)
  if (file.exists(final)) {
    stop("Refusing to overwrite existing public artifact: ", name,
         call. = FALSE)
  }
  csv <- file.path(build_dir, sub("[.]gz$", "", name))
  if (file.exists(csv) || file.exists(paste0(csv, ".gz"))) {
    stop("MV5-J build directory contains a stale compressed artifact.",
         call. = FALSE)
  }
  write_provenance_csv(data, csv)
  status <- system2("gzip", c("-n", "-f", csv))
  gz <- paste0(csv, ".gz")
  if (!identical(status, 0L) || !file.exists(gz) || !file.rename(gz, final)) {
    stop("Failed to deterministically compress a public MV5-J artifact.",
         call. = FALSE)
  }
  final
}

ranking_path <- write_gzip(
  pairs, "mv05j-integrated-cell-retrieval-rankings-2026-08-09.csv.gz"
)
completion_path <- write_plain(
  completion, "mv05j-method-completion-2026-08-09.csv"
)
index_path <- write_plain(
  group_index, "mv05j-group-bundle-index-2026-08-09.csv"
)
registry_path <- write_plain(
  registry, "mv05j-method-registry-2026-08-09.csv"
)
scale_path <- write_plain(
  scale_disposition,
  "mv05j-component-scale-disposition-2026-08-09.csv"
)
summary <- data.frame(
  contract_id = "mv05j_public_assembly_summary_v1",
  groups = nrow(group_index), seeds = length(unique(group_index$seed)),
  biological_pairs = nrow(pairs) / nrow(registry), methods = nrow(registry),
  retrieval_rows = nrow(pairs), method_completion_rows = nrow(completion),
  failed_method_groups = sum(completion$status != "completed"),
  component_scale_rows = nrow(scale_disposition),
  ranking_file_sha256 = file_sha(ranking_path),
  ranking_file_size_bytes = unname(file.info(ranking_path)$size),
  completion_file_sha256 = file_sha(completion_path),
  group_index_file_sha256 = file_sha(index_path),
  method_registry_file_sha256 = file_sha(registry_path),
  scale_disposition_file_sha256 = file_sha(scale_path),
  retrieval_evaluation_jobs_executed = 0L,
  clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
  gene_topology_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, held_out_scale_fit_jobs_executed = 0L,
  biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
write_plain(summary, "mv05j-public-assembly-summary-2026-08-09.csv")
message("Assembled deterministic public MV5-J retrieval artifacts.")
