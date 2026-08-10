#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop(
    "usage: validate_mv05j_production.R BUNDLE_DIR AUDIT_DIR ",
    "D1_RESOURCE_CSV FOLD_CACHE_DIR MEAN_AUDIT_CSV MEAN_DIR G_RESOURCE_CSV ",
    "G_RESULT_ROOT I_COMPONENT_GZ METRICS_CSV GROUP_VALIDATION_CSV ",
    "SUMMARY_CSV", call. = FALSE
  )
}
bundle_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
audit_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
d1_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
fold_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
mean_audit_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
mean_dir <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
g_resource_path <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
g_result_root <- normalizePath(args[[8L]], winslash = "/", mustWork = TRUE)
i_path <- normalizePath(args[[9L]], winslash = "/", mustWork = TRUE)
metrics_path <- normalizePath(args[[10L]], winslash = "/", mustWork = TRUE)
group_output <- args[[11L]]
summary_output <- args[[12L]]

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
source("R/mv05j_integrated_retrieval_inputs.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
key <- function(data) paste(
  data$query_sample_id, data$training_sample_id, sep = "\r"
)

metrics <- utils::read.csv(
  metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
d1 <- utils::read.csv(d1_path, stringsAsFactors = FALSE, check.names = FALSE)
mean_audit <- utils::read.csv(
  mean_audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
g_resource <- utils::read.csv(
  g_resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
i_components <- utils::read.csv(
  i_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(metrics) != 75L || any(metrics$disposition != "completed") ||
    any(metrics$failed_methods != 0L) || nrow(d1) != 75L ||
    nrow(mean_audit) != 5L || nrow(g_resource) != 75L ||
    nrow(i_components) != 35350L ||
    any(c("tissue", "approach") %in% names(i_components)) ||
    any(i_components$outcome_label_state != "closed") ||
    any(as.logical(i_components$biological_outcomes_computed))) {
  stop("MV5-J production inputs fail independent preflight.",
       call. = FALSE)
}

mean_bundles <- list()
group_rows <- list()
all_pair_ids <- character()
all_ranking_ids <- character()
all_bundle_keys <- character()
all_method_counts <- integer(nrow(mv05j_method_registry_v1()))
names(all_method_counts) <- mv05j_method_registry_v1()$method_id
for (index in seq_len(nrow(metrics))) {
  metric <- metrics[index, , drop = FALSE]
  stem <- paste0(safe_name(metric$held_out_study), "__", metric$seed)
  bundle_path <- file.path(bundle_dir, paste0(stem, "__retrieval.rds"))
  audit_path <- file.path(audit_dir, paste0(stem, "__audit.csv"))
  if (!file.exists(bundle_path) || !file.exists(audit_path)) {
    stop("A declared MV5-J group artifact is absent.", call. = FALSE)
  }
  bundle <- readRDS(bundle_path)
  mv05j_validate_group_bundle_v1(bundle)
  audit <- utils::read.csv(
    audit_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  fold_row <- d1[d1$fold_id == metric$fold_id & d1$seed == metric$seed,
                 , drop = FALSE]
  if (nrow(audit) != 1L || nrow(fold_row) != 1L ||
      !identical(file_sha(bundle_path), audit$output_file_sha256) ||
      !identical(bundle$cache_key, audit$bundle_cache_key) ||
      !identical(bundle$payload_sha256, audit$payload_sha256) ||
      !identical(bundle$identity$fold_cache_key, fold_row$fold_cache_key)) {
    stop("A group bundle, audit, or fold identity is inconsistent.",
         call. = FALSE)
  }
  fold_path <- file.path(fold_dir, fold_row$private_cache_file)
  if (!identical(file_sha(fold_path), fold_row$private_cache_sha256)) {
    stop("A fold cache hash differs during MV5-J validation.",
         call. = FALSE)
  }
  fold <- readRDS(fold_path)
  mv05d1_validate_cell_fold_record_v1(fold)
  seed_key <- as.character(metric$seed)
  if (is.null(mean_bundles[[seed_key]])) {
    mean_row <- mean_audit[mean_audit$seed == metric$seed, , drop = FALSE]
    mean_path <- file.path(mean_dir, mean_row$private_file)
    if (nrow(mean_row) != 1L ||
        !identical(file_sha(mean_path), mean_row$private_file_sha256)) {
      stop("A mean-profile bundle hash differs during validation.",
           call. = FALSE)
    }
    mean_bundles[[seed_key]] <- readRDS(mean_path)
    mv05d5_validate_mean_profile_bundle_v1(mean_bundles[[seed_key]])
  }
  mean_bundle <- mean_bundles[[seed_key]]
  pairs <- bundle$payload$pairs
  methods <- split(pairs, pairs$method_id)
  registry <- mv05j_method_registry_v1()
  if (!setequal(names(methods), registry$method_id) ||
      any(vapply(methods, nrow, integer(1L)) != audit$biological_pairs) ||
      any(vapply(methods, function(rows) {
        !setequal(key(rows), key(methods[[1L]]))
      }, logical(1L))) ||
      !identical(sort(unique(pairs$query_sample_id), method = "radix"),
                 fold$identity$query_ids) ||
      !identical(sort(unique(pairs$training_sample_id), method = "radix"),
                 fold$identity$training_ids) ||
      any(pairs$query_sample_id %in% fold$identity$training_ids) ||
      any(pairs$training_sample_id %in% fold$identity$query_ids) ||
      any(pairs$outcome_label_state != "closed") ||
      any(as.logical(pairs$biological_outcomes_computed)) ||
      any(c("tissue", "approach") %in% names(pairs))) {
    stop("A group violates method-axis, partition, or label closure.",
         call. = FALSE)
  }

  topology <- i_components[
    i_components$fold_id == metric$fold_id & i_components$seed == metric$seed,
                 , drop = FALSE]
  topology <- topology[order(
    topology$query_sample_id, topology$training_sample_id, method = "radix"
  ), , drop = FALSE]
  topology_difference <- numeric()
  for (method_id in c("integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1",
                      "integrated_cell_landscape_h0_h1_raw_euclidean_v1")) {
    observed <- methods[[method_id]]
    observed <- observed[match(key(topology), key(observed)), , drop = FALSE]
    expected <- switch(
      method_id, integrated_cell_landscape_h0_v1 = topology$h0_distance,
      integrated_cell_landscape_h1_v1 = topology$h1_distance,
      integrated_cell_landscape_h0_h1_raw_euclidean_v1 =
        sqrt(topology$h0_distance^2 + topology$h1_distance^2)
    )
    topology_difference[[method_id]] <- max(abs(observed$distance - expected))
  }

  oracle_indices <- unique(round(seq(
    1, audit$biological_pairs, length.out = 3L
  )))
  energy <- methods[["integrated_cell_distribution_energy_v1"]]
  energy <- energy[order(energy$query_sample_id, energy$training_sample_id,
                         method = "radix"), , drop = FALSE]
  g_metric <- g_resource[
    g_resource$held_out_study == metric$held_out_study &
      g_resource$seed == metric$seed, , drop = FALSE
  ]
  safe_group <- safe_name(g_metric$group_id)
  g_path <- file.path(g_result_root, safe_group, paste0(safe_group, ".rds"))
  if (nrow(g_metric) != 1L || !file.exists(g_path) ||
      !identical(file_sha(g_path), g_metric$private_result_sha256)) {
    stop("An MV5-G record is absent or hash-mismatched during validation.",
         call. = FALSE)
  }
  g_record <- readRDS(g_path)
  mv05f_validate_group_record_v1(g_record)
  if (!identical(g_record$cache_key, bundle$identity$mv05g_group_cache_key) ||
      !identical(g_record$payload$coordinate_set_sha256,
                 bundle$identity$mv05g_coordinate_set_sha256)) {
    stop("An MV5-G record differs from its MV5-J identity.", call. = FALSE)
  }
  energy_expected <- vapply(oracle_indices, function(pair_index) {
    .mv05_empirical_energy_distance(
      g_record$payload$coordinates[[energy$query_sample_id[[pair_index]]]],
      g_record$payload$coordinates[[energy$training_sample_id[[pair_index]]]]
    )
  }, numeric(1L))
  energy_difference <- max(abs(
    energy$distance[oracle_indices] - energy_expected
  ))

  profile_vector <- function(sample_id) {
    profile <- mean_bundle$payload$profiles[[sample_id]]
    panel <- fold$payload$panel
    present <- panel$feature_id %in% names(profile)
    result <- numeric(nrow(panel))
    result[present] <- (profile[panel$feature_id[present]] -
                          fold$payload$center[present]) /
      fold$payload$scale[present]
    result
  }
  pseudobulk <- methods[["pseudobulk_training_standardized_panel_v1"]]
  pseudobulk <- pseudobulk[order(
    pseudobulk$query_sample_id, pseudobulk$training_sample_id,
    method = "radix"
  ), , drop = FALSE]
  pseudobulk_expected <- vapply(oracle_indices, function(pair_index) {
    sqrt(sum((profile_vector(pseudobulk$query_sample_id[[pair_index]]) -
                profile_vector(pseudobulk$training_sample_id[[pair_index]]))^2))
  }, numeric(1L))
  pseudobulk_difference <- max(abs(
    pseudobulk$distance[oracle_indices] - pseudobulk_expected
  ))

  all_pair_ids <- c(all_pair_ids, pairs$pair_id)
  all_ranking_ids <- c(all_ranking_ids, pairs$ranking_id)
  all_bundle_keys <- c(all_bundle_keys, bundle$cache_key)
  counts <- table(factor(pairs$method_id, levels = names(all_method_counts)))
  all_method_counts <- all_method_counts + as.integer(counts)
  group_rows[[length(group_rows) + 1L]] <- data.frame(
    contract_id = "mv05j_independent_group_validation_v1",
    group_order = metric$group_order, group_id = bundle$identity$group_id,
    fold_id = metric$fold_id, seed = metric$seed,
    query_samples = length(fold$identity$query_ids),
    training_samples = length(fold$identity$training_ids),
    biological_pairs = audit$biological_pairs,
    retrieval_rows = nrow(pairs), methods = length(methods),
    h0_max_absolute_difference =
      topology_difference[["integrated_cell_landscape_h0_v1"]],
    h1_max_absolute_difference =
      topology_difference[["integrated_cell_landscape_h1_v1"]],
    combined_max_absolute_difference =
      topology_difference[["integrated_cell_landscape_h0_h1_raw_euclidean_v1"]],
    energy_oracle_pairs = length(oracle_indices),
    energy_max_absolute_difference = energy_difference,
    pseudobulk_oracle_pairs = length(oracle_indices),
    pseudobulk_max_absolute_difference = pseudobulk_difference,
    exact_distance_tie_rows = sum(pairs$distance_tied),
    bundle_cache_key = bundle$cache_key,
    bundle_file_sha256 = file_sha(bundle_path),
    partition_closed = TRUE, outcomes_absent = TRUE,
    retrieval_evaluation_jobs_executed = 0L,
    clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
    gene_topology_jobs_executed = 0L, fusion_jobs_executed = 0L,
    new_data_jobs_executed = 0L, held_out_scale_fit_jobs_executed = 0L,
    biological_outcomes_computed = FALSE,
    outcome_label_state = "closed", stringsAsFactors = FALSE
  )
}
groups <- do.call(rbind, group_rows)
numeric_checks <- c(
  groups$h0_max_absolute_difference, groups$h1_max_absolute_difference,
  groups$combined_max_absolute_difference,
  groups$energy_max_absolute_difference,
  groups$pseudobulk_max_absolute_difference
)
if (nrow(groups) != 75L || anyDuplicated(all_pair_ids) ||
    anyDuplicated(all_ranking_ids) ||
    anyDuplicated(all_bundle_keys) || length(all_pair_ids) != 176750L ||
    !all(all_method_counts == 35350L) || any(numeric_checks > 1e-12) ||
    sum(groups$biological_pairs) != 35350L ||
    sum(groups$retrieval_rows) != 176750L) {
  stop("MV5-J independent production validation failed.", call. = FALSE)
}
summary <- data.frame(
  contract_id = "mv05j_independent_validation_summary_v1",
  groups = nrow(groups), seeds = length(unique(groups$seed)),
  biological_pairs = sum(groups$biological_pairs),
  methods = length(all_method_counts), retrieval_rows = sum(groups$retrieval_rows),
  completed_method_groups = 75L * length(all_method_counts),
  failed_method_groups = 0L,
  maximum_numeric_difference = max(numeric_checks),
  exact_distance_tie_rows = sum(groups$exact_distance_tie_rows),
  unique_pair_ids = length(unique(all_pair_ids)),
  unique_ranking_ids = length(unique(all_ranking_ids)),
  unique_bundle_keys = length(unique(all_bundle_keys)),
  partition_closed_groups = sum(groups$partition_closed),
  outcomes_absent_groups = sum(groups$outcomes_absent),
  retrieval_evaluation_jobs_executed = 0L,
  clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
  gene_topology_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, held_out_scale_fit_jobs_executed = 0L,
  biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
write_provenance_csv(groups, group_output)
write_provenance_csv(summary, summary_output)
message("Independently validated all 75 MV5-J production groups.")
