#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(paste(
    "usage: run_mv05n_baseline_admission.R ADMISSION_REQUESTS D1_RESOURCE",
    "D1_CACHE_ROOT MEAN_AUDIT MEAN_ROOT G_MANIFEST G_RESOURCE G_RESULT_ROOT",
    "PAIR_OUTPUT RESOURCE_OUTPUT IMPLEMENTATION_SOURCE"
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
source(args[[11L]])

read_public <- function(path) {
  utils::read.csv(normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
                  check.names = FALSE)
}
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)

requests <- read_public(args[[1L]])
d1 <- read_public(args[[2L]])
d1_root <- normalizePath(args[[3L]], mustWork = TRUE)
mean_audit <- read_public(args[[4L]])
mean_root <- normalizePath(args[[5L]], mustWork = TRUE)
g_manifest <- read_public(args[[6L]])
g_resource <- read_public(args[[7L]])
g_root <- normalizePath(args[[8L]], mustWork = TRUE)

if (nrow(requests) != 384L ||
    any(c("tissue", "approach", "class", "label", "outcome") %in%
          tolower(names(requests))) ||
    any(requests$outcome_label_state != "closed") ||
    any(as.logical(requests$biological_outcomes_computed))) {
  stop("MV5-N baseline admission request boundary is invalid.", call. = FALSE)
}
pair_requests <- unique(requests[c(
  "profile", "fold_id", "seed", "representation", "first_sample_id",
  "second_sample_id", "pair_ordinal"
)])
pair_requests <- pair_requests[order(
  pair_requests$profile, pair_requests$representation,
  pair_requests$pair_ordinal, method = "radix"
), ]
if (nrow(pair_requests) != 192L ||
    any(table(pair_requests$profile, pair_requests$representation) != 32L)) {
  stop("MV5-N baseline admission does not resolve to 192 unique pairs.",
       call. = FALSE)
}

seed <- 20260805L
mean_row <- mean_audit[mean_audit$seed == seed, , drop = FALSE]
if (nrow(mean_row) != 1L ||
    !mean_row$disposition %in% c("built_atomic", "reused_validated")) {
  stop("MV5-N mean-profile source is not accepted.", call. = FALSE)
}
mean_path <- file.path(mean_root, mean_row$private_file)
if (!file.exists(mean_path) || file_sha(mean_path) != mean_row$private_file_sha256) {
  stop("MV5-N mean-profile source is missing or stale.", call. = FALSE)
}
mean_bundle <- readRDS(mean_path)
mv05d5_validate_mean_profile_bundle_v1(mean_bundle)

pair_rows <- list()
resource_rows <- list()
pair_cursor <- resource_cursor <- 0L
for (profile in c("minimum", "representative", "maximum")) {
  profile_pairs <- pair_requests[pair_requests$profile == profile, ]
  fold_id <- unique(profile_pairs$fold_id)
  held_out_study <- sub("^large_loso_v1:", "", fold_id)
  d1_row <- d1[d1$held_out_study == held_out_study & d1$seed == seed, ]
  if (nrow(d1_row) != 1L || d1_row$exit_status != 0L ||
      d1_row$outcome_label_state != "closed" ||
      as.logical(d1_row$biological_outcomes_computed)) {
    stop("MV5-N fold source is not accepted.", call. = FALSE)
  }
  fold_path <- file.path(d1_root, d1_row$private_cache_file)
  if (!file.exists(fold_path) || file_sha(fold_path) != d1_row$private_cache_sha256) {
    stop("MV5-N fold source is missing or stale.", call. = FALSE)
  }
  fold <- readRDS(fold_path)
  mv05d1_validate_cell_fold_record_v1(fold)
  vectors <- mv05d5_pseudobulk_vectors_v1(mean_bundle$payload$profiles, fold)

  for (representation in c("sct_whole", "inductive_integrated")) {
    selected_pairs <- profile_pairs[
      profile_pairs$representation == representation,
      c("first_sample_id", "second_sample_id", "pair_ordinal"), drop = FALSE
    ]
    if (representation == "sct_whole") {
      views <- fold$payload$cell_views
      source_identity <- fold$cache_key
    } else {
      g_row <- g_manifest[g_manifest$held_out_study == held_out_study &
                            g_manifest$seed == seed, , drop = FALSE]
      g_metric <- g_resource[g_resource$held_out_study == held_out_study &
                               g_resource$seed == seed, , drop = FALSE]
      safe_group <- gsub("[^A-Za-z0-9_.-]", "_", g_row$group_id)
      g_path <- file.path(g_root, safe_group, paste0(safe_group, ".rds"))
      if (nrow(g_row) != 1L || nrow(g_metric) != 1L || !file.exists(g_path) ||
          file_sha(g_path) != g_metric$private_result_sha256) {
        stop("MV5-N integrated coordinate source is missing or stale.",
             call. = FALSE)
      }
      g_record <- readRDS(g_path)
      mv05f_validate_group_record_v1(g_record)
      ids <- sort(names(g_record$payload$coordinates), method = "radix")
      views <- lapply(ids, function(sample_id) {
        mv05h_new_integrated_cell_view_v1(
          g_record$payload$coordinates[[sample_id]], sample_id,
          fold$identity$fold_id, fold$identity$fit_scope_id, seed,
          g_record$cache_key, g_record$payload$coordinate_set_sha256
        )
      })
      names(views) <- ids
      source_identity <- g_record$cache_key
    }

    energy_started <- proc.time()[["elapsed"]]
    energy <- mv05n_training_energy_pairs_v1(views, selected_pairs)
    energy_seconds <- proc.time()[["elapsed"]] - energy_started
    pseudobulk_started <- proc.time()[["elapsed"]]
    pseudobulk <- mv05n_training_pseudobulk_pairs_v1(vectors, selected_pairs)
    pseudobulk_seconds <- proc.time()[["elapsed"]] - pseudobulk_started

    add_rows <- function(values, method_id, source_id) {
      data.frame(
        contract_id = "mv05n_baseline_admission_pair_v1", profile = profile,
        fold_id = fold_id, seed = seed, representation = representation,
        method_id = method_id, pair_ordinal = values$pair_ordinal,
        first_sample_id = values$first_sample_id,
        second_sample_id = values$second_sample_id, distance = values$distance,
        source_identity = source_id, outcome_label_state = "closed",
        biological_outcomes_computed = FALSE, clustering_jobs_executed = 0L,
        stringsAsFactors = FALSE
      )
    }
    pair_cursor <- pair_cursor + 1L
    pair_rows[[pair_cursor]] <- add_rows(
      energy, "cell_distribution_energy_v1", source_identity
    )
    pair_cursor <- pair_cursor + 1L
    pair_rows[[pair_cursor]] <- add_rows(
      pseudobulk, "pseudobulk_training_standardized_panel_v1",
      mean_bundle$cache_key
    )
    for (method in c("energy", "pseudobulk")) {
      resource_cursor <- resource_cursor + 1L
      seconds <- if (method == "energy") energy_seconds else pseudobulk_seconds
      resource_rows[[resource_cursor]] <- data.frame(
        contract_id = "mv05n_baseline_admission_resource_v1",
        profile = profile, fold_id = fold_id, seed = seed,
        representation = representation,
        method_id = if (method == "energy") "cell_distribution_energy_v1" else
          "pseudobulk_training_standardized_panel_v1",
        pair_rows = nrow(selected_pairs), operation_seconds = seconds,
        operation_seconds_per_pair = seconds / nrow(selected_pairs),
        elapsed_cap_seconds = 900, elapsed_cap_passed = seconds <= 900,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        clustering_jobs_executed = 0L, stringsAsFactors = FALSE
      )
    }
  }
}

pairs <- do.call(rbind, pair_rows)
resources <- do.call(rbind, resource_rows)
pairs <- pairs[order(pairs$profile, pairs$representation, pairs$method_id,
                     pairs$pair_ordinal, method = "radix"), ]
resources <- resources[order(resources$profile, resources$representation,
                             resources$method_id, method = "radix"), ]
if (nrow(pairs) != 384L || nrow(resources) != 12L ||
    any(!is.finite(pairs$distance)) || any(pairs$distance < 0) ||
    any(!resources$elapsed_cap_passed)) {
  stop("MV5-N baseline admission failed numerical or resource checks.",
       call. = FALSE)
}
if (file.exists(args[[9L]]) || file.exists(args[[10L]])) {
  stop("Refusing to overwrite MV5-N baseline admission artifacts.",
       call. = FALSE)
}
write_provenance_csv(pairs, args[[9L]])
write_provenance_csv(resources, args[[10L]])
message("Completed 384 MV5-N baseline admission rows.")
