#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: run_mv06g_stage1_group.R QUEUE PARENT_CONTRACT SOURCE_GROUPS ",
       "LAUNCH SOURCES GROUP_ROOT RUST_LIBRARY OUTPUT_ROOT", call. = FALSE)
}
source("R/landscape_rust_prototype.R")
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
source("R/mv06g_stage1.R")
paths <- vapply(args[1:7], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
queue <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
source_groups <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
launch <- utils::read.csv(paths[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
sources <- utils::read.csv(paths[[5L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
source_group <- source_groups[source_groups$group_id == stage$group_id,
                              , drop = FALSE]
root <- normalizePath(paths[[6L]], winslash = "/", mustWork = TRUE)
rust <- paths[[7L]]
if (nrow(stage) != 1L || nrow(parent) != 1L || nrow(launch) != 1L ||
    nrow(source_group) != 1L || parent$queue_root_sha256 !=
      "f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d" ||
    launch$group_id != stage$group_id ||
    launch$parent_contract_sha256 != .mv06f_sha256(paths[[2L]]) ||
    parent$rust_library_sha256 != .mv06f_sha256(rust) ||
    !identical(sources$path, mv06g_stage1_source_paths_v1()) ||
    !all(file.exists(sources$path)) ||
    !identical(unname(vapply(sources$path, .mv06f_sha256, character(1L))),
               unname(sources$sha256)) ||
    launch$stage1_implementation_root_sha256 != digest::digest(
      stats::setNames(sources$sha256, sources$path), algo = "sha256",
      serialize = TRUE
    )) {
  stop("MV6-G stage-one execution identities are stale.", call. = FALSE)
}
safe <- gsub("[^A-Za-z0-9_.-]", "_", stage$group_id)
source_dir <- file.path(root, safe)
if (.mv06f_sha256(file.path(source_dir, "diagrams.rds")) !=
      source_group$diagrams_sha256 ||
    .mv06f_sha256(file.path(source_dir, "distances.csv")) !=
      source_group$distances_sha256) {
  stop("MV6-G stage-one accepted source group is stale.", call. = FALSE)
}
output_root <- args[[8L]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
final_dir <- file.path(output_root, safe)
partials <- list.files(output_root, pattern = paste0("^", safe, "\\.partial\\."),
                       full.names = TRUE)
if (length(partials)) {
  stop("MV6-G stage one found partial state; quarantine is required.",
       call. = FALSE)
}
if (dir.exists(final_dir)) {
  mv06g_validate_stage1_group_v1(final_dir, stage, parent, launch,
                                 source_group)
  message("Reused validated MV6-G stage-one group.")
  quit(status = 0L, save = "no")
}
partial_dir <- paste0(final_dir, ".partial.", Sys.getpid())
if (!dir.create(partial_dir, recursive = FALSE)) {
  stop("MV6-G stage one could not create an atomic partial directory.",
       call. = FALSE)
}
started <- proc.time()[["elapsed"]]
records <- readRDS(file.path(source_dir, "diagrams.rds"))
invisible(lapply(records, mv06f_validate_ph_record_v1))
training_ids <- sort(unique(vapply(records, function(record) {
  if (record$identity$role == "training") record$identity$sample_id else NA_character_
}, character(1L))), method = "radix", na.last = NA)
held_out_ids <- sort(unique(vapply(records, function(record) {
  if (record$identity$role == "held_out") record$identity$sample_id else NA_character_
}, character(1L))), method = "radix", na.last = NA)
if (length(training_ids) != 65L || length(held_out_ids) != 25L) {
  stop("MV6-G stage-one record roles differ from the frozen maximum group.",
       call. = FALSE)
}
record_map <- stats::setNames(records, names(records))
intervals <- list()
for (sample_id in training_ids) {
  for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
    record <- record_map[[paste(sample_id, view_id, sep = "\r")]]
    for (dimension in c("H0", "H1")) {
      intervals[[paste(sample_id, view_id, dimension, sep = "\r")]] <-
        mv06f_finite_intervals_v1(record, dimension)
    }
  }
}
pairs <- utils::combn(training_ids, 2L)
rows <- vector("list", ncol(pairs) * 4L)
index <- 0L
for (pair_index in seq_len(ncol(pairs))) {
  first_id <- pairs[1L, pair_index]
  second_id <- pairs[2L, pair_index]
  for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
    for (dimension in c("H0", "H1")) {
      first <- intervals[[paste(first_id, view_id, dimension, sep = "\r")]]
      second <- intervals[[paste(second_id, view_id, dimension, sep = "\r")]]
      result <- landscape_rust_prototype_dimension(
        first, second, as.integer(sub("H", "", dimension, fixed = TRUE)), rust
      )
      if (!isTRUE(result$rust_used) || result$status != 0L ||
          !is.finite(result$squared_distance) || result$squared_distance < 0) {
        stop("MV6-G stage-one Rust calculation failed closed.", call. = FALSE)
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        contract_id = "mv06g_training_landscape_distance_v1",
        group_id = stage$group_id, fold_id = stage$fold_id,
        seed = as.integer(stage$seed),
        pair_id = mv06g_training_pair_id_v1(
          stage$group_id, first_id, second_id, view_id, dimension
        ), first_sample_id = first_id, second_sample_id = second_id,
        view_id = view_id, homology_dimension = dimension,
        component_id = mv06g_component_id_v1(view_id, dimension),
        first_finite_intervals = result$first_finite_intervals,
        second_finite_intervals = result$second_finite_intervals,
        squared_distance = result$squared_distance,
        distance = sqrt(result$squared_distance),
        active_levels = result$active_levels,
        event_segments = result$event_segments, exact = TRUE,
        all_active_levels = TRUE, level_cap_applied = FALSE,
        rust_status = result$status,
        rust_engine_version = result$engine_version,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        fusion_evaluations = 0L, outcome_jobs = 0L,
        stringsAsFactors = FALSE
      )
    }
  }
}
training <- do.call(rbind, rows)
training <- training[order(training$component_id, training$first_sample_id,
                           training$second_sample_id, method = "radix"),
                     , drop = FALSE]
rownames(training) <- NULL
scales <- mv06g_fit_component_scales_v1(training)
query <- utils::read.csv(file.path(source_dir, "distances.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
rankings <- mv06g_build_rankings_v1(query, scales)
elapsed <- proc.time()[["elapsed"]] - started
status_lines <- if (file.exists("/proc/self/status")) {
  readLines("/proc/self/status", warn = FALSE)
} else character()
peak_line <- grep("^VmHWM:", status_lines, value = TRUE)
peak_self <- if (length(peak_line) == 1L) {
  as.numeric(gsub("[^0-9]", "", peak_line)) * 1024
} else NA_real_
metrics <- data.frame(
  contract_id = "mv06g_stage1_group_metrics_v1", group_id = stage$group_id,
  training_biological_pairs = ncol(pairs),
  training_component_rows = nrow(training), component_scales = nrow(scales),
  query_biological_pairs = nrow(query) / 4L,
  query_ranking_rows = nrow(rankings), elapsed_seconds = elapsed,
  peak_self_rss_bytes = peak_self, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_evaluations = 0L,
  outcome_jobs = 0L, stringsAsFactors = FALSE
)
artifacts <- file.path(partial_dir, c("training-distances.csv", "scales.csv",
                                     "rankings.csv", "metrics.csv"))
utils::write.csv(training, artifacts[[1L]], row.names = FALSE, na = "")
utils::write.csv(scales, artifacts[[2L]], row.names = FALSE, na = "")
utils::write.csv(rankings, artifacts[[3L]], row.names = FALSE, na = "")
utils::write.csv(metrics, artifacts[[4L]], row.names = FALSE, na = "")
hashes <- vapply(artifacts, .mv06f_sha256, character(1L))
bytes <- unname(file.info(artifacts)$size)
status <- data.frame(
  contract_id = "mv06g_stage1_group_status_v1", group_id = stage$group_id,
  completion_state = "complete",
  parent_contract_sha256 = launch$parent_contract_sha256,
  stage1_implementation_root_sha256 = launch$stage1_implementation_root_sha256,
  rust_library_sha256 = parent$rust_library_sha256,
  source_diagrams_sha256 = source_group$diagrams_sha256,
  source_distances_sha256 = source_group$distances_sha256,
  training_distances_sha256 = hashes[[1L]],
  training_distances_bytes = bytes[[1L]], scales_sha256 = hashes[[2L]],
  scales_bytes = bytes[[2L]], rankings_sha256 = hashes[[3L]],
  rankings_bytes = bytes[[3L]], metrics_sha256 = hashes[[4L]],
  metrics_bytes = bytes[[4L]], training_biological_pairs = ncol(pairs),
  training_component_rows = nrow(training), component_scales = nrow(scales),
  query_biological_pairs = nrow(query) / 4L,
  query_ranking_rows = nrow(rankings), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_evaluations = 0L,
  outcome_jobs = 0L, stringsAsFactors = FALSE
)
utils::write.csv(status, file.path(partial_dir, "status.csv"), row.names = FALSE,
                 na = "")
mv06g_validate_stage1_group_v1(partial_dir, stage, parent, launch,
                               source_group)
if (!file.rename(partial_dir, final_dir)) {
  stop("MV6-G stage one failed atomic publication.", call. = FALSE)
}
message("Completed MV6-G stage-one group: ", stage$group_id, "; ",
        nrow(training), " training rows and ", nrow(rankings),
        " ranking rows.")
