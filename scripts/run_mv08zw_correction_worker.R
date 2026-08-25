#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: run_mv08zw_correction_worker.R <prefreeze> <ph-private-root>",
  "<rust-library> <H0|H1> <private-output-root> <execution-head>"
), call. = FALSE)
for (package in c("digest")) if (!requireNamespace(package, quietly = TRUE)) {
  stop(package, " required", call. = FALSE)
}

prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
ph_root <- normalizePath(args[[2L]], mustWork = TRUE)
rust_library <- normalizePath(args[[3L]], mustWork = TRUE)
dimension <- args[[4L]]
output_root <- normalizePath(args[[5L]], mustWork = FALSE)
execution_head <- tolower(trimws(args[[6L]]))
if (!dimension %in% c("H0", "H1") ||
    !grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV8-ZW worker arguments are invalid", call. = FALSE)
}

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv08z_landscape_production.R")

.mv08z_verify_manifest(prefreeze, "mv08zw-artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(prefreeze, "mv08zw-contract.csv"))
queue <- .mv08z_read_csv(file.path(prefreeze, "mv08zw-correction-queue.csv"))
implementation <- .mv08z_read_csv(file.path(
  prefreeze, "mv08zw-implementation-bindings.csv"))
source_freeze <- .mv08z_read_csv(file.path(prefreeze, "mv08zw-source-freeze.csv"))
row <- queue[queue$homology_dimension == dimension, , drop = FALSE]
if (nrow(contract) != 1L || nrow(row) != 1L ||
    execution_head != contract$execution_head ||
    .mv08z_sha256_file(rust_library) != contract$rust_library_sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, .mv08z_sha256_file, character(1L)) ==
           implementation$sha256) ||
    !all(file.exists(source_freeze$artifact)) ||
    !all(vapply(source_freeze$artifact, .mv08z_sha256_file, character(1L)) ==
           source_freeze$sha256)) {
  stop("MV8-ZW worker binding drift", call. = FALSE)
}

ph_queue <- .mv08z_read_csv(
  "docs/audits/mv08r-full-topology-production-prefreeze-v1/mv08r-ph-queue.csv")
ph_rehash <- .mv08z_read_csv(
  "docs/audits/mv08t-ph-sentinel-closure-v1/mv08t-private-artifact-rehash.csv")
wanted <- ph_queue[
  ph_queue$dataset_scope == "external8" &
    ph_queue$representation_id == "sct_data_selected384_fit_same_axis" &
    ph_queue$panel_id == "common475" &
    ph_queue$seed == 20260805L &
    ph_queue$view_kind == "gene_topology_v1", , drop = FALSE
]
rehash <- ph_rehash[ph_rehash$artifact_type == "ph" &
                      ph_rehash$artifact_id %in% wanted$job_id, , drop = FALSE]
files <- list.files(file.path(ph_root, "ph"), pattern = "[.]rds$",
                    full.names = TRUE)
records <- lapply(files, readRDS)
job_ids <- vapply(records, function(value) value$identity$job_id, character(1L))
selected <- match(wanted$job_id, job_ids)
if (nrow(wanted) != 8L || nrow(rehash) != 8L || anyNA(selected) ||
    anyDuplicated(job_ids) || !all(rehash$independently_validated)) {
  stop("MV8-ZW closed PH selection drift", call. = FALSE)
}
records <- records[selected]
files <- files[selected]
for (index in seq_along(records)) {
  record <- records[[index]]
  mv08s_validate_ph_record_v1(record)
  receipt <- rehash[rehash$artifact_id == record$identity$job_id, , drop = FALSE]
  if (nrow(receipt) != 1L ||
      .mv08z_sha256_file(files[[index]]) != receipt$primary_sha256 ||
      as.numeric(file.info(files[[index]])$size) != receipt$primary_bytes ||
      record$identity$dataset_scope != "external8" ||
      record$identity$representation_id !=
        "sct_data_selected384_fit_same_axis" ||
      record$identity$panel_id != "common475" ||
      record$identity$view_id != "gene_topology_v1") {
    stop("MV8-ZW PH record identity drift", call. = FALSE)
  }
}
ordering <- order(vapply(records, function(value) value$identity$unit_id,
                         character(1L)), method = "radix")
records <- records[ordering]
bindings <- data.frame(
  axis_order = seq_along(records),
  job_id = vapply(records, function(value) value$identity$job_id, character(1L)),
  unit_id = vapply(records, function(value) value$identity$unit_id, character(1L)),
  diagram_sha256 = vapply(records, function(value)
    value$topology_result$provenance$diagram_sha256, character(1L)),
  ph_cache_key = vapply(records, function(value) value$cache_key, character(1L)),
  stringsAsFactors = FALSE
)
group_id <- row$source_group_id
pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(bindings), group_id)
if (nrow(pairs) != 28L) stop("MV8-ZW pair-axis drift", call. = FALSE)
intervals <- lapply(records, .mv08z_finite_intervals,
                    homology_dimension = dimension)
dimension_number <- as.integer(sub("H", "", dimension, fixed = TRUE))
started <- proc.time()[["elapsed"]]
results <- lapply(seq_len(nrow(pairs)), function(index) {
  pair <- pairs[index, , drop = FALSE]
  first <- as.integer(pair$first_axis_order)
  second <- as.integer(pair$second_axis_order)
  value <- landscape_rust_prototype_dimension(
    intervals[[first]], intervals[[second]], dimension_number,
    library = rust_library
  )
  expected_depth <- max(.mv08z_active_depth(intervals[[first]]),
                        .mv08z_active_depth(intervals[[second]]))
  if (!isTRUE(value$rust_used) || value$status != 0L ||
      value$engine_version != 2L || !is.finite(value$squared_distance) ||
      value$squared_distance < 0 || value$active_levels != expected_depth) {
    stop("MV8-ZW engine-v2 calculation failed at pair ", index, call. = FALSE)
  }
  data.frame(
    contract_id = "mv08zw_landscape_distance_v1",
    execution_head = execution_head, group_id = group_id,
    homology_dimension = dimension, pair_ordinal = pair$pair_ordinal,
    pair_identity_sha256 = pair$pair_identity_sha256,
    first_unit_id = pair$first_unit_id, second_unit_id = pair$second_unit_id,
    first_diagram_sha256 = pair$first_diagram_sha256,
    second_diagram_sha256 = pair$second_diagram_sha256,
    first_ph_cache_key = records[[first]]$cache_key,
    second_ph_cache_key = records[[second]]$cache_key,
    squared_distance = value$squared_distance,
    distance = sqrt(value$squared_distance),
    active_levels = as.integer(value$active_levels),
    event_segments = as.integer(value$event_segments),
    first_finite_intervals = nrow(intervals[[first]]),
    second_finite_intervals = nrow(intervals[[second]]),
    exact = TRUE, all_active_levels = TRUE, essential_h0_excluded = TRUE,
    grid_points = 0L, level_cap_applied = FALSE,
    engine_id = "rust_scph_landscape_kernel_v2",
    rust_library_sha256 = contract$rust_library_sha256,
    workers = 1L, retries = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, comparison_jobs = 0L,
    clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
    outcome_jobs = 0L, stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, results)
elapsed <- proc.time()[["elapsed"]] - started
final_dir <- file.path(output_root, tolower(dimension))
if (dir.exists(final_dir)) stop("MV8-ZW output already exists", call. = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = paste0(tolower(dimension), "__partial__"),
                    tmpdir = output_root)
dir.create(partial)
.mv08z_atomic_csv(results, file.path(partial, "distances.csv"))
status <- data.frame(
  contract_id = "mv08zw_group_status_v1", execution_head = execution_head,
  group_id = group_id, homology_dimension = dimension,
  completion_state = "complete", units = 8L, pair_count = 28L,
  pair_axis_sha256 = .mv08z_sha256_text(results$pair_identity_sha256),
  elapsed_seconds = elapsed,
  distances_bytes = as.numeric(file.info(file.path(partial, "distances.csv"))$size),
  distances_sha256 = .mv08z_sha256_file(file.path(partial, "distances.csv")),
  rust_library_sha256 = contract$rust_library_sha256,
  scientific_engine_version = 2L, exact = TRUE,
  all_active_levels = TRUE, essential_h0_excluded = TRUE,
  grid_points = 0L, level_cap_applied = FALSE, workers = 1L, retries = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
.mv08z_atomic_csv(status, file.path(partial, "status.csv"))
if (!file.rename(partial, final_dir)) {
  stop("MV8-ZW atomic group promotion failed", call. = FALSE)
}
message("Completed MV8-ZW ", dimension, "; pairs=28")
