#!/usr/bin/env Rscript

args <- getOption("mv06e.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 10L) {
  stop("usage: validate_mv06e_acceleration.R REFERENCES PAIRS ",
       "PERSIM_A PERSIM_A_METRICS PERSIM_B PERSIM_B_METRICS ",
       "RUST_A RUST_A_METRICS RUST_B RUST_B_METRICS", call. = FALSE)
}
paths <- lapply(args, normalizePath, winslash = "/", mustWork = TRUE)
read_one <- function(index) utils::read.csv(
  paths[[index]], stringsAsFactors = FALSE, check.names = FALSE
)
references <- read_one(1L); pairs <- read_one(2L)
persim_a <- read_one(3L); persim_a_metrics <- read_one(4L)
persim_b <- read_one(5L); persim_b_metrics <- read_one(6L)
rust_a <- read_one(7L); rust_a_metrics <- read_one(8L)
rust_b <- read_one(9L); rust_b_metrics <- read_one(10L)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
checks <- list()
add <- function(category, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv06e_acceleration_validation_v1", category = category,
    passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE
  )
}

scientific_schema <- c(
  "row_type", "pair_id", "view_id", "homology_dimension",
  "first_diagram_id", "second_diagram_id", "first_finite_intervals",
  "second_finite_intervals", "squared_distance", "distance",
  "active_levels", "event_segments", "exact", "all_active_levels",
  "level_cap_applied", "outcome_label_state", "biological_outcomes_computed"
)
valid_frame <- function(value, engine) {
  is.data.frame(value) && nrow(value) == 240L &&
    all(scientific_schema %in% names(value)) &&
    identical(sort(as.integer(table(value$row_type))),
              sort(c(180L, 20L, 40L))) &&
    all(value$engine_id == engine) &&
    all(is.finite(value$squared_distance) & value$squared_distance >= 0) &&
    all(is.finite(value$distance) & value$distance >= 0) &&
    all(abs(value$distance^2 - value$squared_distance) <=
          1e-12 + 1e-12 * abs(value$squared_distance)) &&
    all(as.logical(value$exact)) && all(as.logical(value$all_active_levels)) &&
    !any(as.logical(value$level_cap_applied)) &&
    all(value$outcome_label_state == "closed") &&
    !any(as.logical(value$biological_outcomes_computed))
}
add("engine_output_contracts",
    valid_frame(persim_a, "persim_0.3.8_mv05d4_exact_critical_pairs_v1") &&
      valid_frame(rust_a, "rust_scph_landscape_kernel_v1"),
    "both candidates contain 180 primary, 20 reverse, and 40 self exact rows")

reference_check <- function(candidate) {
  primary <- candidate[candidate$row_type == "primary", ]
  matched <- match(references$pair_id, primary$pair_id)
  if (anyNA(matched)) return(rep(FALSE, nrow(references)))
  observed <- primary$squared_distance[matched]
  tolerance <- ifelse(
    as.logical(references$reference_exact),
    1e-10 + 1e-10 * abs(references$reference_squared_distance),
    references$achieved_absolute_error +
      100 * .Machine$double.eps * pmax(
        1, abs(references$reference_squared_distance)
      )
  )
  is.finite(observed) &
    abs(observed - references$reference_squared_distance) <= tolerance &
    primary$first_finite_intervals[matched] ==
      references$first_finite_intervals &
    primary$second_finite_intervals[matched] ==
      references$second_finite_intervals
}
persim_reference <- reference_check(persim_a)
rust_reference <- reference_check(rust_a)
add("canonical_r_references",
    length(persim_reference) == 20L && all(persim_reference) &&
      length(rust_reference) == 20L && all(rust_reference),
    "grouped Persim and Rust each pass all 20 exact/adaptive-certified R rows")

pa <- persim_a[persim_a$row_type == "primary", ]
ra <- rust_a[rust_a$row_type == "primary", ]
matched <- match(pa$pair_id, ra$pair_id)
cross_error <- abs(pa$squared_distance - ra$squared_distance[matched])
cross_tolerance <- 1e-10 + 1e-10 * abs(pa$squared_distance)
add("candidate_cross_equivalence",
    nrow(pa) == 180L && !anyNA(matched) && all(cross_error <= cross_tolerance) &&
      identical(pa$first_finite_intervals,
                ra$first_finite_intervals[matched]) &&
      identical(pa$second_finite_intervals,
                ra$second_finite_intervals[matched]),
    "all 180 throughput rows agree across exact engines within frozen tolerance")

reverse_check <- function(candidate) {
  primary <- candidate[candidate$row_type == "primary", ]
  reverse <- candidate[candidate$row_type == "reverse", ]
  matched <- match(reverse$pair_id, primary$pair_id)
  !anyNA(matched) && all(abs(reverse$squared_distance -
    primary$squared_distance[matched]) <= 1e-12 +
      1e-12 * abs(primary$squared_distance[matched])) &&
    identical(reverse$first_finite_intervals,
              primary$second_finite_intervals[matched]) &&
    identical(reverse$second_finite_intervals,
              primary$first_finite_intervals[matched])
}
add("reverse_symmetry", reverse_check(persim_a) && reverse_check(rust_a),
    "20/20 reversed reference rows pass for both engines")

self_check <- function(candidate) {
  self <- candidate[candidate$row_type == "self", ]
  nrow(self) == 40L && all(self$squared_distance == 0) &&
    all(self$distance == 0) && identical(self$first_finite_intervals,
                                         self$second_finite_intervals)
}
add("self_zero", self_check(persim_a) && self_check(rust_a),
    "40/40 diagram/dimension self rows are exactly zero for both engines")

deterministic_columns <- intersect(names(persim_a), names(persim_b))
persim_repeat <- identical(persim_a[deterministic_columns],
                           persim_b[deterministic_columns]) &&
  identical(file.info(paths[[3L]])$size, file.info(paths[[5L]])$size) &&
  identical(file_sha(paths[[3L]]), file_sha(paths[[5L]]))
deterministic_columns <- intersect(names(rust_a), names(rust_b))
rust_repeat <- identical(rust_a[deterministic_columns],
                         rust_b[deterministic_columns]) &&
  identical(file.info(paths[[7L]])$size, file.info(paths[[9L]])$size) &&
  identical(file_sha(paths[[7L]]), file_sha(paths[[9L]]))
add("byte_repeat", persim_repeat && rust_repeat,
    "both 240-row normalized scientific outputs repeat byte-identically")

metrics_valid <- function(value, engine) {
  nrow(value) == 1L && value$engine_id == engine &&
    value$diagram_dimensions_built %in% c(0L, 40L) &&
    value$primary_pair_rows == 180L && value$reverse_rows == 20L &&
    value$self_rows == 40L && value$total_seconds > 0 &&
    value$total_seconds <= 3600 && value$peak_process_rss_bytes <= 8 * 1024^3 &&
    value$outcome_label_state == "closed" &&
    !as.logical(value$biological_outcomes_computed)
}
add("resource_guards",
    metrics_valid(persim_a_metrics,
                  "persim_0.3.8_mv05d4_exact_critical_pairs_v1") &&
      metrics_valid(persim_b_metrics,
                  "persim_0.3.8_mv05d4_exact_critical_pairs_v1") &&
      metrics_valid(rust_a_metrics, "rust_scph_landscape_kernel_v1") &&
      metrics_valid(rust_b_metrics, "rust_scph_landscape_kernel_v1"),
    "both clean runs remain below the 3,600-second and 8-GiB guards")

project_one <- function(metric, engine) {
  if (grepl("persim", engine)) {
    build_unit <- metric$landscape_build_seconds / 40
    pair_unit <- metric$primary_pair_seconds / 180
    overhead <- max(0, metric$total_seconds - metric$landscape_build_seconds -
                      metric$primary_pair_seconds - metric$reverse_self_seconds)
    worker_seconds <- 27000 * build_unit + 141400 * pair_unit + 75 * overhead
  } else {
    build_unit <- 0
    pair_unit <- metric$primary_pair_seconds / 180
    overhead <- max(0, metric$total_seconds - metric$primary_pair_seconds -
                      metric$reverse_self_seconds)
    worker_seconds <- 141400 * pair_unit + 75 * overhead
  }
  data.frame(
    contract_id = "mv06e_landscape_production_projection_v1",
    engine_id = engine, run_id = NA_character_,
    diagram_dimension_build_seconds = build_unit,
    component_pair_seconds = pair_unit, group_overhead_seconds = overhead,
    projected_diagram_dimensions = 27000L,
    projected_component_pairs = 141400L,
    projected_groups = 75L, projected_worker_hours = worker_seconds / 3600,
    peak_process_rss_bytes = metric$peak_process_rss_bytes,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
projection <- rbind(
  transform(project_one(persim_a_metrics, persim_a_metrics$engine_id),
            run_id = "A"),
  transform(project_one(persim_b_metrics, persim_b_metrics$engine_id),
            run_id = "B"),
  transform(project_one(rust_a_metrics, rust_a_metrics$engine_id), run_id = "A"),
  transform(project_one(rust_b_metrics, rust_b_metrics$engine_id), run_id = "B")
)
persim_max <- max(projection$projected_worker_hours[
  grepl("persim", projection$engine_id)])
rust_max <- max(projection$projected_worker_hours[
  grepl("rust", projection$engine_id)])
speedup <- persim_max / rust_max
add("production_projection",
    nrow(projection) == 4L &&
      all(is.finite(projection$projected_worker_hours) &
            projection$projected_worker_hours > 0) &&
      all(projection$peak_process_rss_bytes <= 8 * 1024^3),
    "four finite cap-evaluable production projections were reconstructed")

if (!all(vapply(checks, function(value) value$passed[[1L]], logical(1L)))) {
  decision_id <- "stop_acceleration_candidate"
} else if (persim_max <= 60 && rust_max <= 60 && speedup >= 3) {
  decision_id <- "admit_both_with_rust_preferred_private"
} else if (persim_max <= 60 && rust_max <= 60) {
  decision_id <- "admit_grouped_persim"
} else if (rust_max <= 60) {
  decision_id <- "admit_explicit_rust_wsl"
} else if (persim_max <= 60) {
  decision_id <- "admit_grouped_persim"
} else {
  decision_id <- "revise_batch_design"
}

decision <- data.frame(
  contract_id = "mv06e_acceleration_decision_v1", decision = decision_id,
  persim_maximum_projected_worker_hours = persim_max,
  rust_maximum_projected_worker_hours = rust_max,
  rust_projection_speedup = speedup,
  persim_peak_rss_bytes = max(persim_a_metrics$peak_process_rss_bytes,
                              persim_b_metrics$peak_process_rss_bytes),
  rust_peak_rss_bytes = max(rust_a_metrics$peak_process_rss_bytes,
                            rust_b_metrics$peak_process_rss_bytes),
  rust_library_sha256 = rust_a_metrics$rust_library_sha256,
  full_production_authorized = FALSE, fusion_jobs_executed = 0L,
  clustering_jobs_executed = 0L, outcome_jobs_executed = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

validation <- do.call(rbind, checks)
output_dir <- dirname(paths[[3L]])
utils::write.csv(validation, file.path(output_dir,
  "mv06e-validation.csv"), row.names = FALSE, na = "")
utils::write.csv(projection, file.path(output_dir,
  "mv06e-projection.csv"), row.names = FALSE, na = "")
utils::write.csv(decision, file.path(output_dir,
  "mv06e-decision.csv"), row.names = FALSE, na = "")
if (!all(validation$passed)) {
  stop("MV6-E validation failed: ",
       paste(validation$category[!validation$passed], collapse = ", "),
       call. = FALSE)
}
message("Validated MV6-E: ", decision_id, "; Rust projection speedup ",
        format(speedup, digits = 6), "x.")
