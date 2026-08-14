#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(
    "usage: run_mv05d5_group.R D1_RESOURCE_CSV FOLD_CACHE_DIR ",
    "MEAN_AUDIT_CSV MEAN_DIR D4_COMPONENT_GZ OUTPUT_RDS AUDIT_CSV ",
    "HELD_OUT_STUDY SEED", call. = FALSE
  )
}
d1_resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mean_audit_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
mean_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
d4_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output_path <- args[[6L]]
audit_path <- args[[7L]]
held_out_study <- args[[8L]]
seed <- as.integer(args[[9L]])
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d5_retrieval_inputs.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
implementation_files <- c(
  "R/mv05_resource_safe_execution.R", "R/mv05d5_retrieval_inputs.R",
  "scripts/run_mv05d5_group.R"
)
implementation_sha <- .mv05d5_digest(stats::setNames(
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
  stop("MV5-D5 mean-profile audit does not admit this seed.", call. = FALSE)
}
mean_path <- file.path(mean_dir, mean_row$private_file)
if (!file.exists(mean_path) ||
    !identical(file_sha(mean_path), mean_row$private_file_sha256)) {
  stop("MV5-D5 mean-profile bundle is missing or hash-mismatched.",
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

d4 <- utils::read.csv(
  d4_path, stringsAsFactors = FALSE, check.names = FALSE
)
topology <- d4[d4$fold_id == fold$identity$fold_id & d4$seed == seed,
               , drop = FALSE]
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
  stop("MV5-D4 component rows do not match the accepted fold scope.",
       call. = FALSE)
}
topology <- topology[order(
  topology$query_sample_id, topology$training_sample_id, method = "radix"
), , drop = FALSE]
d4_group_sha <- .mv05d5_digest(topology)
identity <- list(
  contract_id = "mv05d5_group_identity_v1",
  group_id = unique(topology$group_id), fold_id = fold$identity$fold_id,
  fit_scope_id = fold$identity$fit_scope_id, seed = seed,
  fold_cache_key = fold$cache_key,
  fold_payload_sha256 = fold$payload_sha256,
  mean_profile_cache_key = mean_bundle$cache_key,
  mv05d4_group_sha256 = d4_group_sha,
  implementation_sha256 = implementation_sha
)
if (length(identity$group_id) != 1L) {
  stop("MV5-D4 rows contain multiple group identities.", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
disposition <- "built_atomic"
if (file.exists(output_path)) {
  bundle <- readRDS(output_path)
  mv05d5_validate_group_bundle_v1(bundle)
  observed <- bundle$identity
  observed$cache_key <- NULL
  if (!identical(observed, identity)) {
    stop("Existing MV5-D5 group is stale; refusing overwrite.",
         call. = FALSE)
  }
  disposition <- "reused_validated"
} else {
  registry <- mv05d5_method_registry_v1()
  add_method <- function(method_id, distances, source_pair_keys) {
    method <- registry[registry$method_id == method_id, , drop = FALSE]
    if (nrow(method) != 1L || nrow(distances) != expected_pairs ||
        length(source_pair_keys) != expected_pairs) {
      stop("MV5-D5 method assembly received an invalid pair table.",
           call. = FALSE)
    }
    data.frame(
      contract_id = "mv05d5_retrieval_pair_v1",
      group_id = identity$group_id, fold_id = identity$fold_id,
      fit_scope_id = identity$fit_scope_id, seed = seed,
      representation = "sct_fold", view_id = "cell_topology_v1",
      method_id = method_id, method_role = method$role,
      distance_policy = method$distance_policy,
      scale_policy = method$scale_policy,
      primary_retrieval = method$primary_retrieval,
      query_sample_id = distances$query_sample_id,
      training_sample_id = distances$training_sample_id,
      distance = distances$distance,
      source_pair_key = source_pair_keys,
      fold_cache_key = fold$cache_key,
      mean_profile_cache_key = if (
        method_id == "pseudobulk_shared_panel_euclidean_v1"
      ) mean_bundle$cache_key else "not_applicable",
      mv05d4_group_sha256 = if (grepl("^cell_landscape", method_id)) {
        d4_group_sha
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
  energy <- mv05d5_energy_pairs_v1(
    fold$payload$cell_views, query_ids, training_ids
  )
  vectors <- mv05d5_pseudobulk_vectors_v1(
    mean_bundle$payload$profiles, fold
  )
  pseudobulk <- mv05d5_pseudobulk_pairs_v1(
    vectors, query_ids, training_ids
  )
  pair_key <- paste(base$query_sample_id, base$training_sample_id, sep = "\r")
  order_like_topology <- function(values) {
    key <- paste(values$query_sample_id, values$training_sample_id, sep = "\r")
    values[match(pair_key, key), , drop = FALSE]
  }
  energy <- order_like_topology(energy)
  pseudobulk <- order_like_topology(pseudobulk)
  view_keys <- vapply(fold$payload$cell_views, `[[`, character(1L),
                      "cache_key")
  energy_keys <- vapply(seq_len(expected_pairs), function(index) paste0(
    "mv05d5_energy_pair_v1:", .mv05d5_digest(list(
      query = view_keys[[base$query_sample_id[[index]]]],
      training = view_keys[[base$training_sample_id[[index]]]]
    ))
  ), character(1L))
  normalization_keys <- fold$identity$normalization_cache_keys
  pseudobulk_keys <- vapply(seq_len(expected_pairs), function(index) paste0(
    "mv05d5_pseudobulk_pair_v1:", .mv05d5_digest(list(
      query = normalization_keys[[base$query_sample_id[[index]]]],
      training = normalization_keys[[base$training_sample_id[[index]]]],
      standardization_id = fold$payload$standardization_id,
      panel = fold$payload$panel
    ))
  ), character(1L))
  pairs <- do.call(rbind, list(
    add_method("cell_landscape_h0_v1", h0, topology$h0_pair_request_id),
    add_method("cell_landscape_h1_v1", h1, topology$h1_pair_request_id),
    add_method(
      "cell_landscape_h0_h1_raw_euclidean_v1", combined,
      paste0("mv05d5_raw_combined_pair_v1:", vapply(
        seq_len(expected_pairs), function(index) .mv05d5_digest(list(
          h0 = topology$h0_pair_request_id[[index]],
          h1 = topology$h1_pair_request_id[[index]]
        )), character(1L)
      ))
    ),
    add_method(
      "cell_distribution_energy_shared_pca_v1", energy, energy_keys
    ),
    add_method(
      "pseudobulk_shared_panel_euclidean_v1", pseudobulk,
      pseudobulk_keys
    )
  ))
  pairs$pair_id <- paste0("mv05d5_retrieval_pair_v1:", vapply(
    seq_len(nrow(pairs)), function(index) .mv05d5_digest(list(
      group_id = pairs$group_id[[index]],
      method_id = pairs$method_id[[index]],
      query_sample_id = pairs$query_sample_id[[index]],
      training_sample_id = pairs$training_sample_id[[index]],
      source_pair_key = pairs$source_pair_key[[index]]
    )), character(1L)
  ))
  failures <- data.frame(
    contract_id = "mv05d5_method_completion_v1",
    group_id = identity$group_id, fold_id = identity$fold_id, seed = seed,
    method_id = registry$method_id, status = "completed",
    expected_pair_rows = expected_pairs, completed_pair_rows = expected_pairs,
    failure_code = "none", failure_detail = "",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  bundle <- mv05d5_new_group_bundle_v1(identity, pairs, failures)
  mv05d5_validate_group_bundle_v1(bundle)
  partial <- tempfile(pattern = basename(output_path),
                      tmpdir = dirname(output_path))
  saveRDS(bundle, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, output_path)) {
    unlink(partial)
    stop("Failed to atomically publish the MV5-D5 group bundle.",
         call. = FALSE)
  }
}
elapsed <- proc.time()[["elapsed"]] - started
audit <- data.frame(
  contract_id = "mv05d5_group_audit_v1",
  group_id = identity$group_id, fold_id = identity$fold_id,
  fit_scope_id = identity$fit_scope_id, held_out_study = held_out_study,
  seed = seed, query_samples = length(query_ids),
  training_samples = length(training_ids),
  biological_pairs = expected_pairs, method_count = 5L,
  retrieval_rows = nrow(bundle$payload$pairs),
  completed_methods = sum(bundle$payload$failures$status == "completed"),
  failed_methods = sum(bundle$payload$failures$status != "completed"),
  exact_distance_ties = sum(bundle$payload$pairs$distance_tied),
  bundle_cache_key = bundle$cache_key,
  payload_sha256 = bundle$payload_sha256,
  fold_cache_key = fold$cache_key,
  mean_profile_cache_key = mean_bundle$cache_key,
  mv05d4_group_sha256 = d4_group_sha,
  implementation_sha256 = implementation_sha,
  output_file = basename(output_path),
  output_file_sha256 = file_sha(output_path),
  output_size_bytes = unname(file.info(output_path)$size),
  disposition = disposition, operation_seconds = elapsed,
  clustering_jobs_executed = 0L, integration_jobs_executed = 0L,
  gene_view_jobs_executed = 0L, biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
write_provenance_csv(audit, audit_path)
message("Completed MV5-D5 group: ", held_out_study, " / ", seed)
