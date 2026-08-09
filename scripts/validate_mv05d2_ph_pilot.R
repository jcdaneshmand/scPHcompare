#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)
options(warn = 2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(
    "usage: validate_mv05d2_ph_pilot.R PILOT_CSV PRIMARY_METRICS_CSV ",
    "REPEAT_METRICS_CSV PRIMARY_RESULT_DIR REPEAT_RESULT_DIR ",
    "FOLD_CACHE_DIR DETERMINISM_CSV CROSS_ENGINE_CSV PROJECTION_CSV ",
    "PROFILE_CSV COMPLETION_CSV",
    call. = FALSE
  )
}
manifest_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
primary_metrics_path <- normalizePath(
  args[[2L]], winslash = "/", mustWork = TRUE
)
repeat_metrics_path <- normalizePath(
  args[[3L]], winslash = "/", mustWork = TRUE
)
primary_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
repeat_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
fold_cache_dir <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
determinism_path <- args[[7L]]
cross_engine_path <- args[[8L]]
projection_path <- args[[9L]]
profile_path <- args[[10L]]
completion_path <- args[[11L]]

for (package in c("digest", "ripserr", "TDA")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-D2 validation.", call. = FALSE)
  }
}
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05d2_ph_profiling.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
primary <- utils::read.csv(
  primary_metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
repeated_metrics <- utils::read.csv(
  repeat_metrics_path, stringsAsFactors = FALSE, check.names = FALSE
)
mv05d2_validate_resource_metrics_v1(primary, expected_jobs = 30L)
mv05d2_validate_resource_metrics_v1(repeated_metrics, expected_jobs = 5L)
expected_repeat <- manifest$job_id[as.logical(manifest$repeat_required)]
if (!setequal(repeated_metrics$job_id, expected_repeat) ||
    any(c("tissue", "approach") %in% names(manifest)) ||
    any(manifest$outcome_label_state != "closed") ||
    any(as.logical(manifest$biological_outcomes_computed))) {
  stop("MV5-D2 manifest or repeat axis violates the frozen contract.",
       call. = FALSE)
}

determinism_rows <- lapply(sort(expected_repeat, method = "radix"), function(id) {
  first_metric <- primary[primary$job_id == id, , drop = FALSE]
  second_metric <- repeated_metrics[
    repeated_metrics$job_id == id, , drop = FALSE
  ]
  first_path <- file.path(primary_dir, first_metric$result_file)
  second_path <- file.path(repeat_dir, second_metric$result_file)
  first <- readRDS(first_path)
  second <- readRDS(second_path)
  mv05d2_validate_ph_record_v1(first)
  mv05d2_validate_ph_record_v1(second)
  data.frame(
    contract_id = "mv05d2_exact_repeat_validation_v1",
    job_id = id,
    fold_id = first_metric$fold_id,
    seed = first_metric$seed,
    execution_role = first_metric$execution_role,
    mapping_stratum = first_metric$mapping_stratum,
    implementation_sha256 = first$identity$implementation_sha256,
    primary_diagram_sha256 = first_metric$diagram_sha256,
    repeat_diagram_sha256 = second_metric$diagram_sha256,
    diagram_sha256_identical = identical(
      first_metric$diagram_sha256, second_metric$diagram_sha256
    ),
    record_cache_key_identical = identical(first$cache_key, second$cache_key),
    record_object_identical = identical(first, second),
    primary_file_sha256 = file_sha(first_path),
    repeat_file_sha256 = file_sha(second_path),
    file_bytes_identical = identical(file_sha(first_path), file_sha(second_path)),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
determinism <- do.call(rbind, determinism_rows)
if (any(!as.matrix(determinism[c(
  "diagram_sha256_identical", "record_cache_key_identical",
  "record_object_identical", "file_bytes_identical"
)]))) {
  stop("MV5-D2 exact repeat determinism failed.", call. = FALSE)
}
if (length(unique(determinism$implementation_sha256)) != 1L) {
  stop("MV5-D2 repeated records have inconsistent implementation hashes.",
       call. = FALSE)
}

normalize_diagram <- function(diagram) {
  diagram <- as.matrix(diagram)
  storage.mode(diagram) <- "double"
  if (ncol(diagram) != 3L) {
    stop("Cross-engine diagram has an unexpected schema.", call. = FALSE)
  }
  colnames(diagram) <- c("dimension", "birth", "death")
  diagram[
    diagram[, "dimension"] %in% 0:1 & is.finite(diagram[, "death"]) &
      diagram[, "death"] > diagram[, "birth"], , drop = FALSE
  ]
}
remove_capped_essential_h0 <- function(diagram, point_count) {
  h0 <- which(diagram[, "dimension"] == 0)
  removed <- FALSE
  if (length(h0) == point_count) {
    capped <- h0[[which.max(diagram[h0, "death"])]]
    diagram <- diagram[-capped, , drop = FALSE]
    removed <- TRUE
  }
  list(diagram = diagram, removed = removed)
}
compare_dimension <- function(first, second, dimension, tolerance = 1e-6) {
  first <- first[first[, "dimension"] == dimension, c("birth", "death"),
                 drop = FALSE]
  second <- second[second[, "dimension"] == dimension, c("birth", "death"),
                   drop = FALSE]
  order_rows <- function(value) {
    value[order(value[, "birth"], value[, "death"], method = "radix"),
          , drop = FALSE]
  }
  first <- order_rows(first)
  second <- order_rows(second)
  maximum_error <- if (nrow(first) == nrow(second) && nrow(first) > 0L) {
    max(abs(first - second))
  } else if (nrow(first) == 0L && nrow(second) == 0L) {
    0
  } else {
    Inf
  }
  list(
    first_count = nrow(first), second_count = nrow(second),
    maximum_absolute_error = maximum_error,
    passed = is.finite(maximum_error) && maximum_error <= tolerance
  )
}

cross_rows <- list()
cross_index <- 0L
repeat_jobs <- manifest[manifest$job_id %in% expected_repeat, , drop = FALSE]
repeat_jobs <- repeat_jobs[order(repeat_jobs$job_id, method = "radix"),
                           , drop = FALSE]
for (index in seq_len(nrow(repeat_jobs))) {
  job <- repeat_jobs[index, , drop = FALSE]
  cache_path <- file.path(fold_cache_dir, job$fold_cache_file)
  if (!file.exists(cache_path) || file_sha(cache_path) != job$fold_cache_sha256) {
    stop("Cross-engine source cache differs from the pilot manifest.",
         call. = FALSE)
  }
  fold_record <- readRDS(cache_path)
  mv05d1_validate_cell_fold_record_v1(fold_record)
  view <- fold_record$payload$cell_views[[job$sample_id]]
  validate_topology_view(view)
  coordinates <- view$payload[seq_len(32L), , drop = FALSE]
  maximum_scale <- max(stats::dist(coordinates))
  ripser_normalized <- remove_capped_essential_h0(normalize_diagram(
    ripserr::vietoris_rips(
    dataset = coordinates, max_dim = 1L, threshold = -1, p = 2L,
    return_format = "mat"
  )), nrow(coordinates))
  gudhi_normalized <- remove_capped_essential_h0(normalize_diagram(TDA::ripsDiag(
    X = coordinates, maxdimension = 1L, maxscale = maximum_scale,
    dist = "euclidean", library = "GUDHI", location = FALSE,
    printProgress = FALSE
  )$diagram), nrow(coordinates))
  ripser_diagram <- ripser_normalized$diagram
  gudhi_diagram <- gudhi_normalized$diagram
  for (dimension in 0:1) {
    comparison <- compare_dimension(
      ripser_diagram, gudhi_diagram, dimension = dimension
    )
    cross_index <- cross_index + 1L
    cross_rows[[cross_index]] <- data.frame(
      contract_id = "mv05d2_reduced_cross_engine_validation_v1",
      job_id = job$job_id,
      fold_id = job$fold_id,
      seed = job$seed,
      execution_role = job$execution_role,
      mapping_stratum = job$mapping_stratum,
      deterministic_subset = "first_32_ordered_cells",
      subset_cells = 32L,
      coordinate_count = ncol(coordinates),
      homology_dimension = paste0("H", dimension),
      ripserr_intervals = comparison$first_count,
      gudhi_intervals = comparison$second_count,
      ripserr_capped_essential_h0_removed = ripser_normalized$removed,
      gudhi_capped_essential_h0_removed = gudhi_normalized$removed,
      maximum_absolute_error = comparison$maximum_absolute_error,
      tolerance = 1e-6,
      passed = comparison$passed,
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
cross_engine <- do.call(rbind, cross_rows)
if (any(!cross_engine$passed)) {
  failed <- cross_engine[!cross_engine$passed, , drop = FALSE]
  stop(
    "MV5-D2 reduced cross-engine validation failed: ",
    paste(
      paste0(
        failed$homology_dimension, " ", failed$ripserr_intervals, "/",
        failed$gudhi_intervals, " max_error=",
        failed$maximum_absolute_error
      ),
      collapse = "; "
    ),
    call. = FALSE
  )
}

projection <- mv05d2_project_full_ph_v1(primary, full_jobs = 6750L)
profile_groups <- split(
  primary,
  interaction(primary$execution_role, primary$mapping_stratum,
              drop = TRUE, lex.order = TRUE)
)
profile <- do.call(rbind, lapply(profile_groups, function(group) {
  data.frame(
    contract_id = "mv05d2_cell_ph_profile_summary_v1",
    execution_role = group$execution_role[[1L]],
    mapping_stratum = group$mapping_stratum[[1L]],
    jobs = nrow(group),
    h0_intervals_minimum = min(group$h0_intervals),
    h0_intervals_maximum = max(group$h0_intervals),
    h1_intervals_minimum = min(group$h1_intervals),
    h1_intervals_median = stats::median(group$h1_intervals),
    h1_intervals_maximum = max(group$h1_intervals),
    elapsed_seconds_median = stats::median(group$elapsed_seconds),
    elapsed_seconds_maximum = max(group$elapsed_seconds),
    peak_process_tree_rss_bytes = max(group$peak_process_tree_rss_bytes),
    result_bytes_sum = sum(group$result_size_bytes),
    h0_mst_maximum_absolute_error =
      max(group$h0_mst_maximum_absolute_error),
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))
rownames(profile) <- NULL

completion <- data.frame(
  contract_id = "mv05d2_cell_ph_pilot_completion_v1",
  implementation_sha256 = unique(determinism$implementation_sha256),
  manifest_sha256 = file_sha(manifest_path),
  primary_metrics_sha256 = file_sha(primary_metrics_path),
  repeat_metrics_sha256 = file_sha(repeat_metrics_path),
  pilot_jobs = nrow(primary),
  completed_jobs = sum(primary$disposition == "completed"),
  failed_jobs = sum(primary$disposition != "completed"),
  folds = length(unique(primary$fold_id)),
  seeds = length(unique(primary$seed)),
  held_out_jobs = sum(primary$execution_role == "held_out"),
  training_jobs = sum(primary$execution_role == "training"),
  mapped_held_out_jobs = sum(
    primary$execution_role == "held_out" &
      primary$mapping_stratum == "training_schema_mapped"
  ),
  unmapped_held_out_jobs = sum(
    primary$execution_role == "held_out" &
      primary$mapping_stratum == "no_missing_training_features"
  ),
  h0_intervals_total = sum(primary$h0_intervals),
  h1_intervals_total = sum(primary$h1_intervals),
  worker_seconds_sum = sum(primary$elapsed_seconds),
  elapsed_seconds_median = stats::median(primary$elapsed_seconds),
  elapsed_seconds_maximum = max(primary$elapsed_seconds),
  peak_process_tree_rss_bytes = max(primary$peak_process_tree_rss_bytes),
  private_result_bytes = sum(primary$result_size_bytes),
  h0_mst_checks_passed = sum(primary$h0_mst_oracle_passed),
  exact_repeat_jobs = nrow(determinism),
  exact_repeat_jobs_passed = sum(determinism$file_bytes_identical),
  cross_engine_checks = nrow(cross_engine),
  cross_engine_checks_passed = sum(cross_engine$passed),
  full_production_jobs_launched = 0L,
  landscape_jobs_executed = 0L,
  distance_jobs_executed = 0L,
  clustering_jobs_executed = 0L,
  integration_jobs_executed = 0L,
  gene_view_jobs_executed = 0L,
  outcome_label_state = "closed",
  biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

write_provenance_csv(determinism, determinism_path)
write_provenance_csv(cross_engine, cross_engine_path)
write_provenance_csv(projection, projection_path)
write_provenance_csv(profile, profile_path)
write_provenance_csv(completion, completion_path)
message("Validated bounded MV5-D2 PH pilot and wrote public-safe audits.")
