#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop("usage: run_mv06g_group.R QUEUE PARENT SOURCE_GROUPS POLICY SOURCES ",
       "GROUP_ROOT RUST_LIBRARY GROUP_ID OUTPUT_ROOT", call. = FALSE)
}
source("R/landscape_rust_prototype.R")
source("R/mv06f_production.R")
source("R/mv06g_fusion_prefreeze.R")
source("R/mv06g_stage1.R")
source("R/mv06g_production.R")
paths <- vapply(args[1:7], normalizePath, character(1L), winslash = "/",
                mustWork = TRUE)
queue <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
parent <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
source_groups <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                                 check.names = FALSE)
policy <- utils::read.csv(paths[[4L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
sources <- utils::read.csv(paths[[5L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
row <- queue[queue$group_id == args[[8L]], , drop = FALSE]
source_group <- source_groups[source_groups$group_id == args[[8L]],
                              , drop = FALSE]
if (nrow(row) != 1L || nrow(source_group) != 1L || nrow(parent) != 1L ||
    nrow(policy) != 1L || !identical(sources$path,
                                      mv06g_production_source_paths_v1()) ||
    !identical(unname(vapply(sources$path, .mv06f_sha256, character(1L))),
               unname(sources$sha256)) ||
    policy$production_implementation_root_sha256 != digest::digest(
      stats::setNames(sources$sha256, sources$path), algo = "sha256",
      serialize = TRUE
    ) || policy$parent_contract_sha256 != .mv06f_sha256(paths[[2L]]) ||
    policy$rust_library_sha256 != .mv06f_sha256(paths[[7L]])) {
  stop("MV6-G production execution identities are stale.", call. = FALSE)
}
safe <- gsub("[^A-Za-z0-9_.-]", "_", row$group_id)
source_dir <- file.path(paths[[6L]], safe)
if (.mv06f_sha256(file.path(source_dir, "diagrams.rds")) !=
      source_group$diagrams_sha256 ||
    .mv06f_sha256(file.path(source_dir, "distances.csv")) !=
      source_group$distances_sha256) {
  stop("MV6-G production source group is stale.", call. = FALSE)
}
output_root <- args[[9L]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
final_dir <- file.path(output_root, safe)
partials <- list.files(output_root, pattern = paste0("^", safe, "\\.partial\\."),
                       full.names = TRUE)
if (length(partials)) stop("MV6-G production partial state requires quarantine.",
                           call. = FALSE)
if (dir.exists(final_dir)) {
  mv06g_validate_production_group_v1(final_dir, row, parent, policy,
                                     source_group)
  message("Reused validated MV6-G production group.")
  quit(status = 0L, save = "no")
}
partial <- paste0(final_dir, ".partial.", Sys.getpid())
if (!dir.create(partial)) stop("MV6-G production partial directory failed.",
                               call. = FALSE)
started <- proc.time()[["elapsed"]]
records <- readRDS(file.path(source_dir, "diagrams.rds"))
record_map <- stats::setNames(records, names(records))
invisible(lapply(records, mv06f_validate_ph_record_v1))
training_ids <- sort(unique(vapply(records, function(record) {
  if (record$identity$role == "training") record$identity$sample_id else NA_character_
}, character(1L))), method = "radix", na.last = NA)
if (length(training_ids) != row$training_samples) {
  stop("MV6-G production training role count drifted.", call. = FALSE)
}
intervals <- list()
for (sample_id in training_ids) for (view_id in c("cell_topology_v1",
                                                   "gene_topology_v1")) {
  record <- record_map[[paste(sample_id, view_id, sep = "\r")]]
  for (dimension in c("H0", "H1")) intervals[[paste(
    sample_id, view_id, dimension, sep = "\r"
  )]] <- mv06f_finite_intervals_v1(record, dimension)
}
pairs <- utils::combn(training_ids, 2L)
rows <- vector("list", 4L * ncol(pairs)); index <- 0L
for (pair_index in seq_len(ncol(pairs))) {
  first_id <- pairs[1L, pair_index]; second_id <- pairs[2L, pair_index]
  for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
    for (dimension in c("H0", "H1")) {
      result <- landscape_rust_prototype_dimension(
        intervals[[paste(first_id, view_id, dimension, sep = "\r")]],
        intervals[[paste(second_id, view_id, dimension, sep = "\r")]],
        as.integer(sub("H", "", dimension, fixed = TRUE)), paths[[7L]]
      )
      if (!isTRUE(result$rust_used) || result$status != 0L) {
        stop("MV6-G production Rust calculation failed closed.", call. = FALSE)
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        contract_id = "mv06g_training_landscape_distance_v1",
        group_id = row$group_id, fold_id = row$fold_id,
        seed = as.integer(row$seed),
        pair_id = mv06g_training_pair_id_v1(row$group_id, first_id, second_id,
                                            view_id, dimension),
        first_sample_id = first_id, second_sample_id = second_id,
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
        fusion_evaluations = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
      )
    }
  }
}
training <- do.call(rbind, rows)
training <- training[order(training$component_id, training$first_sample_id,
                           training$second_sample_id, method = "radix"),
                     , drop = FALSE]
rownames(training) <- NULL
scales <- mv06g_fit_group_scales_v1(training, ncol(pairs))
query <- utils::read.csv(file.path(source_dir, "distances.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
rankings <- mv06g_build_group_rankings_v1(
  query, scales, row$training_samples, row$biological_pairs
)
elapsed <- proc.time()[["elapsed"]] - started
lines <- if (file.exists("/proc/self/status")) readLines(
  "/proc/self/status", warn = FALSE
) else character()
peak_line <- grep("^VmHWM:", lines, value = TRUE)
peak <- if (length(peak_line) == 1L) as.numeric(gsub(
  "[^0-9]", "", peak_line
)) * 1024 else NA_real_
metrics <- data.frame(
  contract_id = "mv06g_production_group_metrics_v1", group_id = row$group_id,
  elapsed_seconds = elapsed, peak_self_rss_bytes = peak,
  training_biological_pairs = ncol(pairs),
  training_component_rows = nrow(training), component_scales = 4L,
  query_biological_pairs = row$biological_pairs,
  query_ranking_rows = nrow(rankings), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_evaluations = 0L,
  outcome_jobs = 0L, stringsAsFactors = FALSE
)
artifacts <- file.path(partial, c("training-distances.csv", "scales.csv",
                                 "rankings.csv", "metrics.csv"))
utils::write.csv(training, artifacts[[1L]], row.names = FALSE, na = "")
utils::write.csv(scales, artifacts[[2L]], row.names = FALSE, na = "")
utils::write.csv(rankings, artifacts[[3L]], row.names = FALSE, na = "")
utils::write.csv(metrics, artifacts[[4L]], row.names = FALSE, na = "")
hashes <- vapply(artifacts, .mv06f_sha256, character(1L))
bytes <- unname(file.info(artifacts)$size)
status <- data.frame(
  contract_id = "mv06g_production_group_status_v1", group_id = row$group_id,
  completion_state = "complete",
  parent_contract_sha256 = policy$parent_contract_sha256,
  production_implementation_root_sha256 =
    policy$production_implementation_root_sha256,
  rust_library_sha256 = parent$rust_library_sha256,
  source_diagrams_sha256 = source_group$diagrams_sha256,
  source_distances_sha256 = source_group$distances_sha256,
  training_distances_sha256 = hashes[[1L]], training_distances_bytes = bytes[[1L]],
  scales_sha256 = hashes[[2L]], scales_bytes = bytes[[2L]],
  rankings_sha256 = hashes[[3L]], rankings_bytes = bytes[[3L]],
  metrics_sha256 = hashes[[4L]], metrics_bytes = bytes[[4L]],
  training_biological_pairs = ncol(pairs),
  training_component_rows = nrow(training), component_scales = 4L,
  query_biological_pairs = row$biological_pairs,
  query_ranking_rows = nrow(rankings), outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, fusion_evaluations = 0L,
  outcome_jobs = 0L, stringsAsFactors = FALSE
)
utils::write.csv(status, file.path(partial, "status.csv"), row.names = FALSE,
                 na = "")
mv06g_validate_production_group_v1(partial, row, parent, policy, source_group)
if (!file.rename(partial, final_dir)) stop("MV6-G production atomic publish failed.",
                                          call. = FALSE)
message("Completed MV6-G production group: ", row$group_id, ".")
