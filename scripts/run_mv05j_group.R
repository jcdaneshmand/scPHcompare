#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop(
    "usage: run_mv05j_group.R D1_RESOURCE_CSV FOLD_CACHE_DIR ",
    "MEAN_AUDIT_CSV MEAN_DIR G_MANIFEST_CSV G_RESOURCE_CSV ",
    "G_RESULT_ROOT I_COMPONENT_GZ OUTPUT_RDS AUDIT_CSV HELD_OUT_STUDY SEED",
    call. = FALSE
  )
}
d1_resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mean_audit_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
mean_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
g_manifest_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
g_resource_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
g_result_root <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
i_path <- normalizePath(args[[8L]], winslash = "/", mustWork = TRUE)
output_path <- args[[9L]]
audit_path <- args[[10L]]
held_out_study <- args[[11L]]
seed <- as.integer(args[[12L]])
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)

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
implementation_files <- c(
  "R/mv05_resource_safe_execution.R", "R/mv05_benchmark_execution.R",
  "R/mv05_inductive_mapping.R",
  "R/mv05d5_retrieval_inputs.R",
  "R/mv05f_integration_gate.R", "R/mv05h_integrated_ph_production.R",
  "R/mv05j_integrated_retrieval_inputs.R",
  "scripts/run_mv05j_group.R"
)
implementation_sha <- .mv05j_digest(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)), implementation_files
))

d1 <- utils::read.csv(
  d1_resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
row <- d1[d1$held_out_study == held_out_study & d1$seed == seed,
          , drop = FALSE]
if (nrow(row) != 1L || row$disposition != "built_atomic" ||
    row$exit_status != 0L || row$outcome_label_state != "closed" ||
    as.logical(row$biological_outcomes_computed) ||
    any(c("tissue", "approach") %in% names(d1))) {
  stop("Requested MV5-D1 fold is not an accepted label-closed input.",
       call. = FALSE)
}
fold_path <- file.path(fold_cache_dir, row$private_cache_file)
if (!file.exists(fold_path) ||
    !identical(file_sha(fold_path), row$private_cache_sha256)) {
  stop("MV5-D1 fold cache is missing or hash-mismatched.", call. = FALSE)
}
fold <- readRDS(fold_path)
mv05d1_validate_cell_fold_record_v1(fold)
if (!identical(fold$cache_key, row$fold_cache_key) ||
    !identical(fold$payload_sha256, row$payload_sha256)) {
  stop("MV5-D1 fold payload differs from its public manifest.", call. = FALSE)
}

mean_audit <- utils::read.csv(
  mean_audit_path, stringsAsFactors = FALSE, check.names = FALSE
)
mean_row <- mean_audit[mean_audit$seed == seed, , drop = FALSE]
if (nrow(mean_row) != 1L ||
    !mean_row$disposition %in% c("built_atomic", "reused_validated") ||
    mean_row$outcome_label_state != "closed" ||
    as.logical(mean_row$biological_outcomes_computed)) {
  stop("MV5-J mean-profile audit does not admit this seed.", call. = FALSE)
}
mean_path <- file.path(mean_dir, mean_row$private_file)
if (!file.exists(mean_path) ||
    !identical(file_sha(mean_path), mean_row$private_file_sha256)) {
  stop("MV5-J mean-profile bundle is missing or hash-mismatched.",
       call. = FALSE)
}
mean_bundle <- readRDS(mean_path)
mv05d5_validate_mean_profile_bundle_v1(mean_bundle)
if (!identical(mean_bundle$cache_key, mean_row$mean_profile_cache_key) ||
    !identical(mean_bundle$identity$normalization_cache_keys,
               fold$identity$normalization_cache_keys)) {
  stop("MV5-D5 mean profiles do not match the fold normalization inputs.",
       call. = FALSE)
}

g_manifest <- utils::read.csv(
  g_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
g_resource <- utils::read.csv(
  g_resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
g_row <- g_manifest[
  g_manifest$held_out_study == held_out_study & g_manifest$seed == seed,
  , drop = FALSE
]
g_metric <- g_resource[
  g_resource$held_out_study == held_out_study & g_resource$seed == seed,
  , drop = FALSE
]
if (nrow(g_row) != 1L || nrow(g_metric) != 1L ||
    g_row$outcome_label_state != "closed" ||
    as.logical(g_row$biological_outcomes_computed) ||
    g_metric$disposition != "completed" || g_metric$exit_status != 0L ||
    g_metric$outcome_label_state != "closed" ||
    as.logical(g_metric$biological_outcomes_computed) ||
    !identical(g_row$group_id, g_metric$group_id) ||
    any(c("tissue", "approach") %in% names(g_manifest)) ||
    any(c("tissue", "approach") %in% names(g_resource))) {
  stop("MV5-G public evidence does not admit this integrated group.",
       call. = FALSE)
}
safe_group <- gsub("[^A-Za-z0-9_.-]", "_", g_row$group_id)
g_path <- file.path(g_result_root, safe_group, paste0(safe_group, ".rds"))
if (!file.exists(g_path) ||
    !identical(file_sha(g_path), g_metric$private_result_sha256)) {
  stop("MV5-G private record is missing or hash-mismatched.", call. = FALSE)
}
g_record <- readRDS(g_path)
mv05f_validate_group_record_v1(g_record)
if (!identical(g_record$identity$fold_id, fold$identity$fold_id) ||
    !identical(g_record$identity$fit_scope_id, fold$identity$fit_scope_id) ||
    !identical(g_record$identity$seed, seed) ||
    !identical(g_record$identity$d1_fold_cache_key, fold$cache_key) ||
    !identical(g_record$payload$coordinate_set_sha256,
               g_metric$coordinate_set_sha256)) {
  stop("MV5-G record does not match the accepted fold or public resource.",
       call. = FALSE)
}

components <- utils::read.csv(
  i_path, stringsAsFactors = FALSE, check.names = FALSE
)
topology <- components[
  components$fold_id == fold$identity$fold_id & components$seed == seed,
  , drop = FALSE
]
query_ids <- fold$identity$query_ids
training_ids <- fold$identity$training_ids
expected_pairs <- length(query_ids) * length(training_ids)
if (nrow(topology) != expected_pairs ||
    !identical(sort(unique(topology$query_sample_id), method = "radix"),
               query_ids) ||
    !identical(sort(unique(topology$training_sample_id), method = "radix"),
               training_ids) ||
    any(!is.finite(topology$h0_distance)) ||
    any(!is.finite(topology$h1_distance)) ||
    any(!is.finite(topology$combined_distance)) ||
    any(topology$outcome_label_state != "closed") ||
    any(as.logical(topology$biological_outcomes_computed)) ||
    any(c("tissue", "approach") %in% names(topology))) {
  stop("MV5-I component rows do not match the accepted fold scope.",
       call. = FALSE)
}
topology <- topology[order(
  topology$query_sample_id, topology$training_sample_id, method = "radix"
), , drop = FALSE]
i_group_sha <- .mv05j_digest(topology)
identity <- list(
  contract_id = "mv05j_group_identity_v1",
  group_id = paste0("mv05j_group__", held_out_study, "__", seed),
  mv05i_source_group_id = unique(topology$group_id),
  fold_id = fold$identity$fold_id,
  fit_scope_id = fold$identity$fit_scope_id, seed = seed,
  fold_cache_key = fold$cache_key,
  fold_payload_sha256 = fold$payload_sha256,
  mean_profile_cache_key = mean_bundle$cache_key,
  mean_profile_file_sha256 = mean_row$private_file_sha256,
  mv05g_group_cache_key = g_record$cache_key,
  mv05g_payload_sha256 = g_record$payload_sha256,
  mv05g_coordinate_set_sha256 = g_record$payload$coordinate_set_sha256,
  mv05g_private_file_sha256 = g_metric$private_result_sha256,
  mv05i_group_sha256 = i_group_sha,
  implementation_sha256 = implementation_sha
)
if (length(identity$mv05i_source_group_id) != 1L) {
  stop("MV5-I rows contain multiple group identities.", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
disposition <- "built_atomic"
if (file.exists(output_path)) {
  bundle <- readRDS(output_path)
  mv05j_validate_group_bundle_v1(bundle)
  observed <- bundle$identity
  observed$cache_key <- NULL
  if (!identical(observed, identity)) {
    stop("Existing MV5-J group is stale; refusing overwrite.",
         call. = FALSE)
  }
  disposition <- "reused_validated"
} else {
  registry <- mv05j_method_registry_v1()
  add_method <- function(method_id, distances, source_pair_keys) {
    method <- registry[registry$method_id == method_id, , drop = FALSE]
    if (nrow(method) != 1L || nrow(distances) != expected_pairs ||
        length(source_pair_keys) != expected_pairs) {
      stop("MV5-J method assembly received an invalid pair table.",
           call. = FALSE)
    }
    data.frame(
      contract_id = "mv05j_retrieval_pair_v1",
      group_id = identity$group_id, fold_id = identity$fold_id,
      fit_scope_id = identity$fit_scope_id, seed = seed,
      representation = method$coordinate_representation,
      view_id = if (method_id ==
        "pseudobulk_training_standardized_panel_v1"
      ) "pseudobulk_context_v1" else "cell_topology_v1",
      method_id = method_id, method_role = method$role,
      distance_policy = method$distance_policy,
      scale_policy = method$scale_policy,
      primary_retrieval = method$primary_retrieval,
      query_sample_id = distances$query_sample_id,
      training_sample_id = distances$training_sample_id,
      distance = distances$distance,
      source_pair_key = source_pair_keys,
      source_identity = source_pair_keys,
      fold_cache_key = fold$cache_key,
      mean_profile_cache_key = if (
        method_id == "pseudobulk_training_standardized_panel_v1"
      ) mean_bundle$cache_key else "not_applicable",
      mv05i_group_sha256 = if (grepl("^integrated_cell_landscape", method_id)) {
        i_group_sha
      } else "not_applicable",
      mv05g_group_cache_key = if (grepl("^integrated_cell", method_id)) {
        g_record$cache_key
      } else "not_applicable",
      prediction_payload = "ranked_training_sample_ids_only_no_labels",
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
  base <- topology[, c("query_sample_id", "training_sample_id")]
  h0 <- transform(base, distance = topology$h0_distance)
  h1 <- transform(base, distance = topology$h1_distance)
  combined <- transform(base, distance = topology$combined_distance)
  expected_combined <- sqrt(topology$h0_distance^2 + topology$h1_distance^2)
  if (max(abs(expected_combined - topology$combined_distance)) > 1e-12) {
    stop("MV5-I raw composite is inconsistent with its H0/H1 components.",
         call. = FALSE)
  }
  coordinates <- g_record$payload$coordinates
  ids <- sort(c(query_ids, training_ids), method = "radix")
  if (!identical(sort(names(coordinates), method = "radix"), ids)) {
    stop("MV5-G coordinates do not match the frozen sample partition.",
         call. = FALSE)
  }
  integrated_views <- lapply(ids, function(sample_id) {
    mv05h_new_integrated_cell_view_v1(
      coordinates[[sample_id]], sample_id, fold$identity$fold_id,
      fold$identity$fit_scope_id, seed, g_record$cache_key,
      g_record$payload$coordinate_set_sha256
    )
  })
  names(integrated_views) <- ids
  energy <- mv05j_energy_pairs_v1(
    integrated_views, query_ids, training_ids
  )
  vectors <- mv05j_pseudobulk_vectors_v1(
    mean_bundle$payload$profiles, fold
  )
  pseudobulk <- mv05j_pseudobulk_pairs_v1(
    vectors, query_ids, training_ids
  )
  pair_key <- paste(base$query_sample_id, base$training_sample_id, sep = "\r")
  order_like_topology <- function(values) {
    key <- paste(values$query_sample_id, values$training_sample_id, sep = "\r")
    values[match(pair_key, key), , drop = FALSE]
  }
  energy <- order_like_topology(energy)
  pseudobulk <- order_like_topology(pseudobulk)
  view_keys <- vapply(integrated_views, `[[`, character(1L),
                      "cache_key")
  energy_keys <- vapply(seq_len(expected_pairs), function(index) paste0(
    "mv05j_energy_pair_v1:", .mv05j_digest(list(
      query = view_keys[[base$query_sample_id[[index]]]],
      training = view_keys[[base$training_sample_id[[index]]]]
    ))
  ), character(1L))
  normalization_keys <- fold$identity$normalization_cache_keys
  pseudobulk_keys <- vapply(seq_len(expected_pairs), function(index) paste0(
    "mv05j_pseudobulk_pair_v1:", .mv05j_digest(list(
      query = normalization_keys[[base$query_sample_id[[index]]]],
      training = normalization_keys[[base$training_sample_id[[index]]]],
      standardization_id = fold$payload$standardization_id,
      panel = fold$payload$panel
    ))
  ), character(1L))
  pairs <- do.call(rbind, list(
    add_method("integrated_cell_landscape_h0_v1", h0, topology$h0_pair_request_id),
    add_method("integrated_cell_landscape_h1_v1", h1, topology$h1_pair_request_id),
    add_method(
      "integrated_cell_landscape_h0_h1_raw_euclidean_v1", combined,
      paste0("mv05j_raw_combined_pair_v1:", vapply(
        seq_len(expected_pairs), function(index) .mv05j_digest(list(
          h0 = topology$h0_pair_request_id[[index]],
          h1 = topology$h1_pair_request_id[[index]]
        )), character(1L)
      ))
    ),
    add_method(
      "integrated_cell_distribution_energy_v1", energy, energy_keys
    ),
    add_method(
      "pseudobulk_training_standardized_panel_v1", pseudobulk,
      pseudobulk_keys
    )
  ))
  pairs$pair_id <- paste0("mv05j_retrieval_pair_v1:", vapply(
    seq_len(nrow(pairs)), function(index) .mv05j_digest(list(
      group_id = pairs$group_id[[index]],
      method_id = pairs$method_id[[index]],
      query_sample_id = pairs$query_sample_id[[index]],
      training_sample_id = pairs$training_sample_id[[index]],
      source_pair_key = pairs$source_pair_key[[index]]
    )), character(1L)
  ))
  completion <- data.frame(
    contract_id = "mv05j_method_completion_v1",
    group_id = identity$group_id, fold_id = identity$fold_id, seed = seed,
    method_id = registry$method_id, status = "completed",
    expected_pair_rows = expected_pairs, completed_pair_rows = expected_pairs,
    failure_code = "none", failure_detail = "",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  bundle <- mv05j_new_group_bundle_v1(identity, pairs, completion)
  mv05j_validate_group_bundle_v1(bundle)
  partial <- tempfile(pattern = basename(output_path),
                      tmpdir = dirname(output_path))
  saveRDS(bundle, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output_path)) {
    unlink(partial)
    stop("Failed to atomically publish the MV5-J group bundle.",
         call. = FALSE)
  }
}
elapsed <- proc.time()[["elapsed"]] - started
audit <- data.frame(
  contract_id = "mv05j_group_audit_v1",
  group_id = identity$group_id, fold_id = identity$fold_id,
  fit_scope_id = identity$fit_scope_id, held_out_study = held_out_study,
  seed = seed, query_samples = length(query_ids),
  training_samples = length(training_ids),
  biological_pairs = expected_pairs, method_count = 5L,
  retrieval_rows = nrow(bundle$payload$pairs),
  completed_methods = sum(bundle$payload$completion$status == "completed"),
  failed_methods = sum(bundle$payload$completion$status != "completed"),
  exact_distance_ties = sum(bundle$payload$pairs$distance_tied),
  bundle_cache_key = bundle$cache_key,
  payload_sha256 = bundle$payload_sha256,
  fold_cache_key = fold$cache_key,
  mean_profile_cache_key = mean_bundle$cache_key,
  mv05g_group_cache_key = g_record$cache_key,
  mv05g_coordinate_set_sha256 = g_record$payload$coordinate_set_sha256,
  mv05g_private_file_sha256 = g_metric$private_result_sha256,
  mv05i_group_sha256 = i_group_sha,
  implementation_sha256 = implementation_sha,
  output_file = basename(output_path),
  output_file_sha256 = file_sha(output_path),
  output_size_bytes = unname(file.info(output_path)$size),
  disposition = disposition, operation_seconds = elapsed,
  retrieval_input_methods_executed = 5L,
  retrieval_evaluation_jobs_executed = 0L,
  clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
  gene_topology_jobs_executed = 0L, fusion_jobs_executed = 0L,
  new_data_jobs_executed = 0L, held_out_scale_fit_jobs_executed = 0L,
  biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
write_provenance_csv(audit, audit_path)
message("Completed MV5-J group: ", held_out_study, " / ", seed)
