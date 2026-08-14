#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: stage_mv05d5_mean_profiles.R SCT_RESOURCE_CSV SCT_CACHE_DIR ",
    "OUTPUT_DIR AUDIT_CSV", call. = FALSE
  )
}
resource_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cache_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3L]]
audit_path <- args[[4L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
  "scripts/stage_mv05d5_mean_profiles.R"
)
implementation_sha <- .mv05d5_digest(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)), implementation_files
))
manifest_sha <- file_sha(resource_path)
resources <- utils::read.csv(
  resource_path, stringsAsFactors = FALSE, check.names = FALSE
)
required <- c(
  "sample_id", "seed", "normalization_cache_key", "private_cache_file",
  "private_cache_sha256", "disposition", "exit_status",
  "outcome_label_state", "biological_outcomes_computed"
)
if (nrow(resources) != 450L || !all(required %in% names(resources)) ||
    anyDuplicated(paste(resources$sample_id, resources$seed, sep = "\r")) ||
    !identical(sort(unique(resources$seed)), 20260805:20260809) ||
    any(resources$disposition != "built_atomic") ||
    any(resources$exit_status != 0L) ||
    any(resources$outcome_label_state != "closed") ||
    any(as.logical(resources$biological_outcomes_computed)) ||
    any(c("tissue", "approach") %in% names(resources))) {
  stop("MV5-D0 resource manifest violates the MV5-D5 staging contract.",
       call. = FALSE)
}

rows <- list()
for (seed in sort(unique(resources$seed))) {
  group <- resources[resources$seed == seed, , drop = FALSE]
  group <- group[order(group$sample_id, method = "radix"), , drop = FALSE]
  path <- file.path(output_dir, paste0(seed, "__mean_profiles.rds"))
  started <- proc.time()[["elapsed"]]
  disposition <- "built_atomic"
  if (file.exists(path)) {
    bundle <- readRDS(path)
    mv05d5_validate_mean_profile_bundle_v1(bundle)
    if (!identical(bundle$identity$normalization_cache_keys,
                   stats::setNames(group$normalization_cache_key,
                                   group$sample_id)) ||
        !identical(bundle$identity$source_manifest_sha256, manifest_sha) ||
        !identical(bundle$identity$implementation_sha256,
                   implementation_sha)) {
      stop("Existing MV5-D5 mean-profile bundle is stale; refusing overwrite.",
           call. = FALSE)
    }
    disposition <- "reused_validated"
  } else {
    profiles <- vector("list", nrow(group))
    names(profiles) <- group$sample_id
    for (index in seq_len(nrow(group))) {
      cache_path <- file.path(cache_dir, group$private_cache_file[[index]])
      if (!file.exists(cache_path) ||
          !identical(file_sha(cache_path),
                     group$private_cache_sha256[[index]])) {
        stop("An MV5-D0 cache file is absent or hash-mismatched.",
             call. = FALSE)
      }
      record <- readRDS(cache_path)
      mv05d0_validate_normalization_cache_record_v2(record)
      if (!identical(record$cache_key,
                     group$normalization_cache_key[[index]]) ||
          !identical(record$identity$sample_id, group$sample_id[[index]]) ||
          record$identity$seed != seed) {
        stop("An MV5-D0 cache identity differs from its public manifest.",
             call. = FALSE)
      }
      matrix <- mv05d0_sct_matrix_from_cache_v1(record)
      profile <- Matrix::rowMeans(matrix)
      names(profile) <- rownames(matrix)
      profiles[[index]] <- profile
    }
    bundle <- mv05d5_new_mean_profile_bundle_v1(
      seed, profiles,
      stats::setNames(group$normalization_cache_key, group$sample_id),
      manifest_sha, implementation_sha
    )
    mv05d5_validate_mean_profile_bundle_v1(bundle)
    partial <- tempfile(pattern = basename(path), tmpdir = output_dir)
    saveRDS(bundle, partial, compress = FALSE, version = 3)
    if (!file.rename(partial, path)) {
      unlink(partial)
      stop("Failed to atomically publish an MV5-D5 mean-profile bundle.",
           call. = FALSE)
    }
  }
  elapsed <- proc.time()[["elapsed"]] - started
  rows[[length(rows) + 1L]] <- data.frame(
    contract_id = "mv05d5_mean_profile_staging_audit_v1",
    seed = seed, sample_count = nrow(group),
    normalization_cache_count = length(
      bundle$identity$normalization_cache_keys
    ),
    mean_profile_cache_key = bundle$cache_key,
    payload_sha256 = bundle$payload_sha256,
    source_manifest_sha256 = manifest_sha,
    implementation_sha256 = implementation_sha,
    private_file = basename(path), private_file_sha256 = file_sha(path),
    private_file_size_bytes = unname(file.info(path)$size),
    disposition = disposition, operation_seconds = elapsed,
    ph_jobs_executed = 0L, landscape_jobs_executed = 0L,
    distance_jobs_executed = 0L, clustering_jobs_executed = 0L,
    integration_jobs_executed = 0L, gene_view_jobs_executed = 0L,
    biological_outcomes_computed = FALSE, outcome_label_state = "closed",
    stringsAsFactors = FALSE
  )
  write_provenance_csv(do.call(rbind, rows), audit_path)
}
message("Staged five immutable MV5-D5 mean-profile bundles.")
