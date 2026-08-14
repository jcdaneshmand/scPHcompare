#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 11L) {
  stop("usage: run_mv06f_group.R QUEUE CONTRACT SOURCES CANDIDATE FOLDS ",
       "RESOURCES PANEL CACHE_DIR RUST_LIBRARY GROUP_ID OUTPUT_ROOT",
       call. = FALSE)
}
input_paths <- vapply(args[1:9], normalizePath, character(1L), winslash = "/",
                      mustWork = TRUE)
group_id <- args[[10L]]
output_root <- args[[11L]]

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv06d_matched_profile.R")
source("R/mv06f_production.R")
source("R/landscape_rust_prototype.R")

queue <- utils::read.csv(input_paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
contract <- utils::read.csv(input_paths[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
sources <- utils::read.csv(input_paths[[3L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
candidate <- utils::read.csv(input_paths[[4L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
folds <- utils::read.csv(input_paths[[5L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
resources <- utils::read.csv(input_paths[[6L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
panel_public <- utils::read.csv(input_paths[[7L]], stringsAsFactors = FALSE,
                                check.names = FALSE)
cache_dir <- input_paths[[8L]]
rust_library <- input_paths[[9L]]

mv06f_validate_queue_v1(queue)
mv06f_validate_prefreeze_inputs_v1(candidate, folds, resources, panel_public)
row <- queue[queue$group_id == group_id, , drop = FALSE]
if (nrow(row) != 1L || nrow(contract) != 1L ||
    contract$contract_id != "mv06f_production_prefreeze_v1" ||
    contract$queue_root_sha256 != mv06f_queue_root_v1(queue) ||
    !all(file.exists(sources$path)) ||
    !identical(tolower(unname(vapply(
      sources$path, .mv06f_sha256, character(1L)
    ))), tolower(unname(sources$sha256))) ||
    .mv06f_sha256(rust_library) != contract$rust_library_sha256) {
  stop("MV6-F group preflight identity validation failed.", call. = FALSE)
}
fold <- folds[folds$fold_id == row$fold_id & folds$seed == row$seed,
              , drop = FALSE]
if (nrow(fold) != 1L || fold$fit_scope_id != row$fit_scope_id ||
    fold$held_out_study != row$held_out_study) {
  stop("MV6-F group differs from the fold plan.", call. = FALSE)
}

safe_id <- gsub("[^A-Za-z0-9_.-]", "_", group_id)
final_dir <- file.path(output_root, safe_id)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
partial_pattern <- paste0("^", safe_id, "\\.partial\\.")
partials <- list.files(output_root, pattern = partial_pattern, full.names = TRUE)
if (length(partials)) {
  stop("MV6-F found a prior partial group; quarantine is required.",
       call. = FALSE)
}
if (dir.exists(final_dir)) {
  mv06f_validate_group_directory_v1(
    final_dir, row, contract$queue_root_sha256,
    contract$implementation_root_sha256, contract$rust_library_sha256
  )
  message("Reused validated MV6-F group: ", group_id)
  quit(status = 0L, save = "no")
}
partial_dir <- paste0(final_dir, ".partial.", Sys.getpid())
if (!dir.create(partial_dir, recursive = FALSE, showWarnings = FALSE)) {
  stop("MV6-F could not create its private partial group directory.",
       call. = FALSE)
}
started <- proc.time()[["elapsed"]]

candidate <- candidate[order(candidate$sample_id, method = "radix"),
                       , drop = FALSE]
training_ids <- candidate$sample_id[candidate$study != row$held_out_study]
query_ids <- candidate$sample_id[candidate$study == row$held_out_study]
seed_resources <- resources[resources$seed == row$seed, , drop = FALSE]
seed_resources <- seed_resources[order(seed_resources$sample_id,
                                       method = "radix"), , drop = FALSE]
if (!identical(seed_resources$sample_id, candidate$sample_id) ||
    length(training_ids) != row$training_samples ||
    length(query_ids) != row$held_out_samples) {
  stop("MV6-F group sample roles differ from the queue.", call. = FALSE)
}
cache_paths <- file.path(cache_dir, seed_resources$private_cache_file)
if (!all(file.exists(cache_paths)) ||
    !identical(tolower(unname(vapply(
      cache_paths, .mv06f_sha256, character(1L)
    ))), tolower(unname(seed_resources$private_cache_sha256)))) {
  stop("MV6-F source caches are absent or differ from accepted hashes.",
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
  stop("MV6-F cache keys or selected-cell identities are stale.",
       call. = FALSE)
}
matrices <- lapply(records, mv05d0_sct_matrix_from_cache_v1)
normalization_keys <- stats::setNames(seed_resources$normalization_cache_key,
                                      seed_resources$sample_id)
panel <- panel_public[order(panel_public$panel_order),
                      c("feature_id", "gene"), drop = FALSE]
prepared <- mv06d_prepare_matched_sources_v1(
  matrices, panel, training_ids, row$fold_id, row$fit_scope_id, row$seed,
  normalization_keys
)
rm(records, matrices)
invisible(gc())
pca_model <- fit_cell_topology_pca(prepared$sources[training_ids],
                                   n_components = 30L, pca_seed = row$seed)
group_identity <- list(
  contract_id = "mv06f_source_group_identity_v1", group_id = group_id,
  queue_root_sha256 = contract$queue_root_sha256,
  implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  panel_scientific_sha256 =
    "7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8",
  fold_id = row$fold_id, fit_scope_id = row$fit_scope_id,
  held_out_study = row$held_out_study, seed = as.integer(row$seed),
  training_ids = training_ids, query_ids = query_ids,
  normalization_keys = normalization_keys,
  pca_cache_key = pca_model$cache_key, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE
)
group_key <- paste0("mv06f_source_group_v1:", digest::digest(
  group_identity, algo = "sha256", serialize = TRUE
))

ph_records <- list()
manifest_rows <- list()
intervals <- list()
for (sample_id in candidate$sample_id) {
  role <- if (sample_id %in% query_ids) "held_out" else "training"
  source_object <- prepared$sources[[sample_id]]
  coordinates <- t(source_object$matrix) %*% pca_model$rotation
  views <- list(
    cell_topology_v1 = construct_frozen_cell_topology_view(
      source_object, coordinates, "mv06f_training_fitted_pca_v1",
      pca_model$cache_key
    ),
    gene_topology_v1 = construct_gene_topology_view(source_object)
  )
  for (view_id in names(views)) {
    result <- run_topology_view_ph(views[[view_id]], max_dim = 1L,
                                   threshold = -1, field = 2L)
    record <- mv06f_new_ph_record_v1(
      group_key, sample_id, role, views[[view_id]], result
    )
    record_name <- paste(sample_id, view_id, sep = "\r")
    ph_records[[record_name]] <- record
    h0 <- mv06f_finite_intervals_v1(record, "H0")
    h1 <- mv06f_finite_intervals_v1(record, "H1")
    intervals[[paste(sample_id, view_id, "H0", sep = "\r")]] <- h0
    intervals[[paste(sample_id, view_id, "H1", sep = "\r")]] <- h1
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      contract_id = "mv06f_diagram_manifest_v1", group_id = group_id,
      sample_id = sample_id, role = role, view_id = view_id,
      ph_cache_key = record$cache_key,
      diagram_sha256 = result$provenance$diagram_sha256,
      finite_h0_intervals = nrow(h0), finite_h1_intervals = nrow(h1),
      essential_h0_intervals = 1L, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }
  rm(views, coordinates, source_object)
  invisible(gc(FALSE))
}
manifest <- do.call(rbind, manifest_rows)
manifest <- manifest[order(manifest$sample_id, manifest$view_id,
                           method = "radix"), , drop = FALSE]
names(ph_records) <- names(ph_records)
ph_records <- ph_records[order(names(ph_records), method = "radix")]

distance_rows <- vector("list", row$landscape_component_rows)
index <- 0L
for (query_id in sort(query_ids, method = "radix")) {
  for (training_id in sort(training_ids, method = "radix")) {
    for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
      for (dimension in c("H0", "H1")) {
        first <- intervals[[paste(query_id, view_id, dimension, sep = "\r")]]
        second <- intervals[[paste(training_id, view_id, dimension, sep = "\r")]]
        result <- landscape_rust_prototype_dimension(
          first, second, as.integer(sub("H", "", dimension, fixed = TRUE)),
          library = rust_library
        )
        if (!isTRUE(result$rust_used) || result$status != 0L ||
            !is.finite(result$squared_distance) || result$squared_distance < 0) {
          stop("MV6-F Rust landscape calculation failed closed.",
               call. = FALSE)
        }
        index <- index + 1L
        distance_rows[[index]] <- data.frame(
          contract_id = "mv06f_exact_landscape_distance_v1",
          engine_id = "rust_scph_landscape_kernel_v1", group_id = group_id,
          pair_id = mv06f_pair_id_v1(
            group_id, query_id, training_id, view_id, dimension
          ), fold_id = row$fold_id, seed = as.integer(row$seed),
          query_sample_id = query_id, training_sample_id = training_id,
          view_id = view_id, homology_dimension = dimension,
          first_finite_intervals = result$first_finite_intervals,
          second_finite_intervals = result$second_finite_intervals,
          squared_distance = result$squared_distance,
          distance = sqrt(result$squared_distance),
          active_levels = result$active_levels,
          event_segments = result$event_segments,
          exact = TRUE, all_active_levels = TRUE, level_cap_applied = FALSE,
          rust_status = result$status, rust_engine_version = result$engine_version,
          outcome_label_state = "closed", biological_outcomes_computed = FALSE,
          fusion_jobs = 0L, clustering_jobs = 0L, outcome_jobs = 0L,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
distances <- do.call(rbind, distance_rows)
if (index != row$landscape_component_rows ||
    anyDuplicated(distances$pair_id)) {
  stop("MV6-F generated an incomplete or duplicate group pair axis.",
       call. = FALSE)
}

elapsed <- proc.time()[["elapsed"]] - started
status_lines <- if (file.exists("/proc/self/status")) {
  readLines("/proc/self/status", warn = FALSE)
} else character()
peak_line <- grep("^VmHWM:", status_lines, value = TRUE)
peak_rss <- if (length(peak_line) == 1L) {
  as.numeric(gsub("[^0-9]", "", peak_line)) * 1024
} else NA_real_
metrics <- data.frame(
  contract_id = "mv06f_group_metrics_v1", group_id = group_id,
  elapsed_seconds = elapsed, peak_self_rss_bytes = peak_rss,
  ph_jobs = length(ph_records), diagram_dimension_records = 2L * length(ph_records),
  biological_pairs = row$biological_pairs,
  landscape_component_rows = nrow(distances), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
artifact_paths <- file.path(partial_dir, c(
  "diagrams.rds", "diagram-manifest.csv", "distances.csv", "metrics.csv"
))
saveRDS(ph_records, artifact_paths[[1L]], compress = FALSE, version = 3)
utils::write.csv(manifest, artifact_paths[[2L]], row.names = FALSE, na = "")
utils::write.csv(distances, artifact_paths[[3L]], row.names = FALSE, na = "")
utils::write.csv(metrics, artifact_paths[[4L]], row.names = FALSE, na = "")
artifact_hashes <- vapply(artifact_paths, .mv06f_sha256, character(1L))
artifact_bytes <- unname(file.info(artifact_paths)$size)
status <- data.frame(
  contract_id = "mv06f_group_status_v1", group_id = group_id,
  queue_root_sha256 = contract$queue_root_sha256,
  implementation_root_sha256 = contract$implementation_root_sha256,
  rust_library_sha256 = contract$rust_library_sha256,
  completion_state = "complete",
  diagrams_sha256 = artifact_hashes[[1L]], diagrams_bytes = artifact_bytes[[1L]],
  diagram_manifest_sha256 = artifact_hashes[[2L]],
  diagram_manifest_bytes = artifact_bytes[[2L]],
  distances_sha256 = artifact_hashes[[3L]], distances_bytes = artifact_bytes[[3L]],
  metrics_sha256 = artifact_hashes[[4L]], metrics_bytes = artifact_bytes[[4L]],
  ph_jobs = length(ph_records), diagram_dimension_records = 2L * length(ph_records),
  biological_pairs = row$biological_pairs,
  landscape_component_rows = nrow(distances), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_jobs = 0L,
  clustering_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
utils::write.csv(status, file.path(partial_dir, "status.csv"), row.names = FALSE,
                 na = "")
mv06f_validate_group_directory_v1(
  partial_dir, row, contract$queue_root_sha256,
  contract$implementation_root_sha256, contract$rust_library_sha256
)
if (!file.rename(partial_dir, final_dir)) {
  stop("MV6-F failed to atomically publish the completed group directory.",
       call. = FALSE)
}
message("Completed MV6-F group: ", group_id, "; ", nrow(distances),
        " component rows.")
