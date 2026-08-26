#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "usage: run_mv14_cell_landscape_chunk.R <prefreeze> <private-bindings>",
  "<mv13-private-groups> <rust-library> <group-order> <chunk-order>",
  "<private-production-root> <execution-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
binding_path <- normalizePath(args[[2L]], mustWork = TRUE)
group_root <- normalizePath(args[[3L]], mustWork = TRUE)
rust_library <- normalizePath(args[[4L]], mustWork = TRUE)
group_order <- as.integer(args[[5L]])
chunk_order <- as.integer(args[[6L]])
production_root <- normalizePath(args[[7L]], mustWork = FALSE)
execution_head <- tolower(trimws(args[[8L]]))
if (is.na(group_order) || is.na(chunk_order) ||
    !grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV14 chunk arguments are invalid.", call. = FALSE)
}

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv13_allqc_cell_topology.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv14_cell_landscape.R")

.mv14_verify_manifest(prefreeze, "mv14-artifact-manifest.csv")
contract <- .mv14_read_csv(file.path(prefreeze, "mv14-contract.csv"))
groups <- .mv14_read_csv(file.path(prefreeze, "mv14-group-queue.csv"))
queue <- .mv14_read_csv(file.path(prefreeze, "mv14-production-queue.csv"))
inputs <- .mv14_read_csv(file.path(prefreeze, "mv14-input-bindings.csv"))
implementation <- .mv14_read_csv(file.path(
  prefreeze, "mv14-implementation-bindings.csv"
))
decision <- .mv14_read_csv(file.path(prefreeze, "mv14-decision.csv"))
if (nrow(contract) != 1L || contract$execution_head != execution_head ||
    nrow(decision) != 1L ||
    !.mv14_truth(decision$production_authorized_after_commit) ||
    .mv14_sha256_file(binding_path) !=
      inputs$sha256[inputs$role == "private_axis_bindings"] ||
    .mv14_sha256_file(rust_library) != contract$rust_library_sha256 ||
    !all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, .mv14_sha256_file, character(1L)) ==
           implementation$sha256)) {
  stop("MV14 chunk implementation/input binding drift.", call. = FALSE)
}
group <- groups[groups$group_order == group_order, , drop = FALSE]
chunk <- queue[queue$group_order == group_order &
                 queue$chunk_order == chunk_order, , drop = FALSE]
bindings <- .mv14_read_csv(binding_path)
bindings <- bindings[bindings$group_order == group_order, , drop = FALSE]
if (nrow(group) != 1L || nrow(chunk) != 1L ||
    nrow(bindings) != group$units ||
    chunk$authorization_state != "authorized_after_mv14_prefreeze_commit") {
  stop("MV14 chunk is not prospectively authorized.", call. = FALSE)
}
pairs <- .mv14_group_pairs(bindings, group$group_id)
pairs <- pairs[pairs$pair_ordinal >= chunk$pair_start &
                 pairs$pair_ordinal <= chunk$pair_end, , drop = FALSE]
if (nrow(pairs) != chunk$pair_count ||
    .mv14_sha256_text(pairs$pair_identity_sha256) != chunk$pair_subset_sha256) {
  stop("MV14 deterministic pair subset drift.", call. = FALSE)
}

artifact_path <- file.path(group_root, group$artifact_file)
if (!file.exists(artifact_path) ||
    as.numeric(file.info(artifact_path)$size) != group$artifact_bytes ||
    .mv14_sha256_file(artifact_path) != group$artifact_sha256) {
  stop("MV14 source group artifact drift.", call. = FALSE)
}
artifact <- readRDS(artifact_path)
mv13_validate_cell_group_v1(artifact)
needed <- sort(unique(c(pairs$first_axis_order, pairs$second_axis_order)))
intervals <- vector("list", nrow(bindings))
cache_keys <- character(nrow(bindings))
for (axis in needed) {
  binding <- bindings[bindings$axis_order == axis, , drop = FALSE]
  record <- artifact$records[[binding$record_index]]
  if (nrow(binding) != 1L || record$unit_id != binding$unit_id ||
      record$result$provenance$diagram_sha256 != binding$diagram_sha256) {
    stop("MV14 record identity drift at axis ", axis, call. = FALSE)
  }
  value <- .mv14_intervals(record, group$homology_dimension)
  if (nrow(value) != binding$finite_intervals ||
      .mv14_active_depth(value) != binding$active_depth) {
    stop("MV14 interval inventory drift at axis ", axis, call. = FALSE)
  }
  intervals[[axis]] <- value
  cache_keys[[axis]] <- record$result$cache_key
}

result <- vector("list", nrow(pairs))
dimension_number <- match(group$homology_dimension, c("H0", "H1")) - 1L
started <- proc.time()[["elapsed"]]
for (index in seq_len(nrow(pairs))) {
  pair <- pairs[index, , drop = FALSE]
  first <- as.integer(pair$first_axis_order)
  second <- as.integer(pair$second_axis_order)
  value <- landscape_rust_prototype_dimension(
    intervals[[first]], intervals[[second]], dimension_number,
    library = rust_library
  )
  expected_depth <- max(.mv14_active_depth(intervals[[first]]),
                        .mv14_active_depth(intervals[[second]]))
  if (!isTRUE(value$rust_used) || value$status != 0L ||
      value$engine_version != 2L || !is.finite(value$squared_distance) ||
      value$squared_distance < 0 || value$active_levels != expected_depth) {
    stop("MV14 Rust calculation failed at pair ", pair$pair_ordinal,
         call. = FALSE)
  }
  result[[index]] <- data.frame(
    contract_id = "mv14_cell_landscape_distance_v1",
    execution_head = execution_head, scientific_engine_version = 2L,
    group_order = group_order, group_id = group$group_id,
    chunk_order = chunk_order, pair_ordinal = pair$pair_ordinal,
    pair_identity_sha256 = pair$pair_identity_sha256,
    first_axis_order = first, second_axis_order = second,
    first_unit_id = pair$first_unit_id, second_unit_id = pair$second_unit_id,
    first_result_cache_key = cache_keys[[first]],
    second_result_cache_key = cache_keys[[second]],
    first_diagram_sha256 = pair$first_diagram_sha256,
    second_diagram_sha256 = pair$second_diagram_sha256,
    homology_dimension = group$homology_dimension,
    squared_distance = value$squared_distance,
    distance = sqrt(value$squared_distance),
    active_levels = value$active_levels, event_segments = value$event_segments,
    first_finite_intervals = value$first_finite_intervals,
    second_finite_intervals = value$second_finite_intervals,
    exact = TRUE, all_active_levels = TRUE, grid_points = 0L,
    level_cap_applied = FALSE, rust_status = value$status,
    rust_engine_version = value$engine_version,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
    label_jobs = 0L, outcome_jobs = 0L, manuscript_claim_jobs = 0L,
    stringsAsFactors = FALSE
  )
}
result <- do.call(rbind, result)
elapsed <- proc.time()[["elapsed"]] - started
group_dir <- file.path(production_root, .mv14_safe_group(group_order))
final_dir <- file.path(group_dir, .mv14_safe_chunk(chunk_order))
if (dir.exists(final_dir)) stop("MV14 chunk output already exists.", call. = FALSE)
dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = "mv14_partial_", tmpdir = group_dir)
dir.create(partial)
.mv14_atomic_csv(result, file.path(partial, "distances.csv"))
status <- data.frame(
  contract_id = "mv14_cell_landscape_chunk_status_v1",
  execution_head = execution_head, scientific_engine_version = 2L,
  group_order = group_order, group_id = group$group_id,
  chunk_order = chunk_order, completion_state = "complete",
  pair_start = chunk$pair_start, pair_end = chunk$pair_end,
  pair_count = nrow(result), pair_subset_sha256 = chunk$pair_subset_sha256,
  elapsed_seconds = elapsed,
  rust_library_sha256 = contract$rust_library_sha256,
  distances_bytes = as.numeric(file.info(file.path(partial, "distances.csv"))$size),
  distances_sha256 = .mv14_sha256_file(file.path(partial, "distances.csv")),
  workers = 1L, retries = 0L, fallback_used = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, manuscript_claim_jobs = 0L,
  stringsAsFactors = FALSE
)
.mv14_atomic_csv(status, file.path(partial, "status.csv"))
if (!file.rename(partial, final_dir)) {
  stop("MV14 failed to atomically publish a completed chunk.", call. = FALSE)
}
message("Completed MV14 ", .mv14_safe_group(group_order), "/",
        .mv14_safe_chunk(chunk_order), "; pairs=", nrow(result))
