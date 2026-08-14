#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    "usage: run_mv06d_source_entry.R SENTINEL_CSV PANEL_CSV CANDIDATE_CSV ",
    "RESOURCE_CSV CACHE_DIR FOLD_ID SEED OUTPUT_RDS", call. = FALSE
  )
}
sentinel_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
panel_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
resource_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
cache_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
fold_id <- args[[6L]]
seed <- as.integer(args[[7L]])
output_path <- args[[8L]]

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv06d_matched_profile.R")

file_sha <- mv06d_file_sha256_v1
implementation_files <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv03_pilot.R",
  "R/mv05_resource_safe_execution.R", "R/mv06d_matched_profile.R",
  "scripts/run_mv06d_source_entry.R"
)
implementation_sha <- digest::digest(stats::setNames(vapply(
  implementation_files, file_sha, character(1L)), implementation_files),
  algo = "sha256", serialize = TRUE)

sentinels <- utils::read.csv(sentinel_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
panel_public <- utils::read.csv(panel_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
candidate <- utils::read.csv(candidate_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
resources <- utils::read.csv(resource_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
selected <- sentinels[sentinels$fold_id == fold_id & sentinels$seed == seed,
                      , drop = FALSE]
if (nrow(selected) != 2L || !setequal(selected$role,
                                     c("held_out", "training"))) {
  stop("Requested MV6-D fold sentinel is not unique.", call. = FALSE)
}
panel <- panel_public[order(panel_public$panel_order),
                      c("feature_id", "gene"), drop = FALSE]
panel_sha <- unique(panel_public$panel_sha256)
if (nrow(panel) != 500L || length(panel_sha) != 1L ||
    !identical(panel_public$panel_order, seq_len(500L)) ||
    any(!as.logical(panel_public$selected_all_cache_core))) {
  stop("MV6-D panel differs from the accepted global core.", call. = FALSE)
}
candidate <- candidate[order(candidate$sample_id, method = "radix"),
                       , drop = FALSE]
study <- selected$held_out_study[[1L]]
training_ids <- candidate$sample_id[candidate$study != study]
query_ids <- candidate$sample_id[candidate$study == study]
seed_resources <- resources[resources$seed == seed, , drop = FALSE]
seed_resources <- seed_resources[order(seed_resources$sample_id,
                                       method = "radix"), , drop = FALSE]
if (nrow(seed_resources) != 90L ||
    !identical(seed_resources$sample_id, candidate$sample_id)) {
  stop("MV6-D seed cache axis differs from the candidate manifest.",
       call. = FALSE)
}
cache_paths <- file.path(cache_dir, seed_resources$private_cache_file)
if (!all(file.exists(cache_paths)) ||
    !identical(unname(vapply(cache_paths, file_sha, character(1L))),
               tolower(unname(seed_resources$private_cache_sha256)))) {
  stop("MV6-D source caches are absent or differ from accepted hashes.",
       call. = FALSE)
}
records <- lapply(cache_paths, readRDS)
names(records) <- seed_resources$sample_id
invisible(lapply(records, mv05d0_validate_normalization_cache_record_v2))
if (!identical(unname(vapply(records, `[[`, character(1L), "cache_key")),
               unname(seed_resources$normalization_cache_key)) ||
    !identical(unname(vapply(records, function(record) {
      record$identity$selected_cell_sha256
    }, character(1L))), unname(seed_resources$selected_cell_sha256))) {
  stop("MV6-D cache keys or selected-cell identities are stale.",
       call. = FALSE)
}
matrices <- lapply(records, mv05d0_sct_matrix_from_cache_v1)
normalization_keys <- stats::setNames(seed_resources$normalization_cache_key,
                                      seed_resources$sample_id)
prepared <- mv06d_prepare_matched_sources_v1(
  matrices, panel, training_ids, fold_id, selected$fit_scope_id[[1L]], seed,
  normalization_keys
)
pca_model <- fit_cell_topology_pca(prepared$sources[training_ids],
                                   n_components = 30L, pca_seed = seed)
views <- lapply(c("held_out", "training"), function(role) {
  sample_id <- selected$sample_id[selected$role == role]
  source_object <- prepared$sources[[sample_id]]
  coordinates <- t(source_object$matrix) %*% pca_model$rotation
  list(
    cell_topology_v1 = construct_frozen_cell_topology_view(
      source_object, coordinates, "mv06d_training_fitted_pca_v1",
      pca_model$cache_key
    ),
    gene_topology_v1 = construct_gene_topology_view(source_object)
  )
})
names(views) <- c("held_out", "training")
identity <- mv06d_source_identity_v1(
  selected, training_ids, query_ids, normalization_keys,
  panel_sha, file_sha(panel_path), file_sha(resource_path),
  file_sha(candidate_path),
  "50379f98cd4927c5c8cb19dbd9ca8ecc7b7b3a9af2e04eb9c8358ecb0b722c6d",
  implementation_sha
)
record <- mv06d_new_source_record_v1(identity, panel, prepared, pca_model,
                                     views)
stem <- basename(output_path)
path <- output_path
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(path)) {
  existing <- readRDS(path)
  mv06d_validate_source_record_v1(existing)
  if (!identical(existing$cache_key, record$cache_key)) {
    stop("Existing MV6-D source result has a stale identity; refusing overwrite.",
         call. = FALSE)
  }
  message("Reused validated MV6-D source bundle: ", stem)
} else {
  partial <- tempfile(pattern = stem, tmpdir = dirname(path))
  saveRDS(record, partial, compress = FALSE, version = 3)
  if (!file.rename(partial, path)) {
    unlink(partial)
    stop("Failed to atomically publish MV6-D source bundle.", call. = FALSE)
  }
  message("Built MV6-D source bundle: ", stem)
}
