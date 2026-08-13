args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: run_mv05bd_full_rust_equivalence.R <intervals.csv> ",
    "<references.csv> <rust-library> <output-dir>", call. = FALSE
  )
}

interval_path <- args[[1L]]
reference_path <- args[[2L]]
rust_library <- args[[3L]]
output_dir <- args[[4L]]

source("R/landscape_rust_prototype.R")

atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  utils::write.csv(value, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not atomically publish ", path, ".", call. = FALSE)
  }
}

intervals <- utils::read.csv(interval_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
references <- utils::read.csv(reference_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
references <- references[order(
  references$stratum_id, references$pair_order, references$dimension,
  references$pair_cache_key
), , drop = FALSE]
if (nrow(references) != 408L ||
    sum(as.logical(references$reference_exact)) != 318L ||
    sum(!as.logical(references$reference_exact)) != 90L) {
  stop("MV5-BD requires the complete frozen 318/90 reference corpus.",
       call. = FALSE)
}

intervals <- intervals[order(intervals$diagram_id, intervals$dimension,
                             intervals$interval_order), , drop = FALSE]
interval_keys <- paste(intervals$diagram_id, intervals$dimension, sep = "|")
interval_map <- split(intervals[, c("birth", "death"), drop = FALSE],
                      interval_keys, drop = TRUE)

interval_for <- function(diagram_id, dimension) {
  key <- paste(diagram_id, as.integer(sub("^H", "", dimension)), sep = "|")
  value <- interval_map[[key]]
  if (is.null(value)) {
    return(matrix(numeric(), nrow = 0L, ncol = 2L,
                  dimnames = list(NULL, c("birth", "death"))))
  }
  as.matrix(value)
}

call_kernel <- function(first_id, second_id, dimension) {
  landscape_rust_prototype_dimension(
    interval_for(first_id, dimension), interval_for(second_id, dimension),
    dimension = as.integer(sub("^H", "", dimension)), library = rust_library
  )
}

threshold_for <- function(row) {
  reference <- as.double(row$reference_squared_distance)
  if (isTRUE(as.logical(row$reference_exact))) {
    1e-10 + 1e-10 * abs(reference)
  } else {
    as.double(row$achieved_absolute_error_estimate) +
      100 * .Machine$double.eps * max(1, abs(reference))
  }
}

equivalence <- lapply(seq_len(nrow(references)), function(index) {
  row <- references[index, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  forward <- call_kernel(
    row$first_diagram_id, row$second_diagram_id, row$dimension
  )
  elapsed <- proc.time()[["elapsed"]] - started
  reverse <- call_kernel(
    row$second_diagram_id, row$first_diagram_id, row$dimension
  )
  reference <- as.double(row$reference_squared_distance)
  threshold <- threshold_for(row)
  error <- abs(forward$squared_distance - reference)
  tier <- if (isTRUE(as.logical(row$reference_exact))) "D" else "E"
  identity_payload <- paste(
    "mv05bd_result_v1", tier, row$stratum_id, row$pair_order,
    row$pair_cache_key, row$dimension,
    sprintf("%.17g", forward$squared_distance), forward$active_levels,
    forward$event_segments, sep = "|"
  )
  data.frame(
    contract_id = "mv05bd_private_equivalence_v1", tier = tier,
    corpus_order = index, stratum_id = row$stratum_id,
    pair_order = row$pair_order, pair_cache_key = row$pair_cache_key,
    dimension = row$dimension,
    reference_squared_distance = reference,
    candidate_squared_distance = forward$squared_distance,
    absolute_error = error, acceptance_threshold = threshold,
    equivalent = isTRUE(forward$rust_used) && error <= threshold,
    status = forward$status, engine_version = forward$engine_version,
    active_levels = forward$active_levels,
    event_segments = forward$event_segments,
    first_finite_intervals = forward$first_finite_intervals,
    expected_first_finite_intervals = row$first_finite_intervals,
    second_finite_intervals = forward$second_finite_intervals,
    expected_second_finite_intervals = row$second_finite_intervals,
    finite_counts_match =
      forward$first_finite_intervals == row$first_finite_intervals &&
      forward$second_finite_intervals == row$second_finite_intervals,
    reverse_bit_identical = identical(
      forward$squared_distance, reverse$squared_distance
    ),
    reverse_counts_swap =
      forward$first_finite_intervals == reverse$second_finite_intervals &&
      forward$second_finite_intervals == reverse$first_finite_intervals,
    reverse_diagnostics_match =
      forward$active_levels == reverse$active_levels &&
      forward$event_segments == reverse$event_segments,
    elapsed_seconds = elapsed,
    result_identity = paste0(
      "mv05bd_rust_result_v1:",
      digest::digest(identity_payload, algo = "sha256", serialize = FALSE)
    ),
    stringsAsFactors = FALSE
  )
})
equivalence <- do.call(rbind, equivalence)

diagram_dimensions <- unique(intervals[, c("diagram_id", "dimension"),
                                       drop = FALSE])
diagram_dimensions <- diagram_dimensions[order(
  diagram_dimensions$diagram_id, diagram_dimensions$dimension
), , drop = FALSE]
self_rows <- lapply(seq_len(nrow(diagram_dimensions)), function(index) {
  item <- diagram_dimensions[index, , drop = FALSE]
  dimension <- paste0("H", item$dimension)
  result <- call_kernel(item$diagram_id, item$diagram_id, dimension)
  data.frame(
    contract_id = "mv05bd_private_self_v1", self_order = index,
    diagram_id = item$diagram_id, dimension = dimension,
    squared_distance = result$squared_distance,
    exactly_zero = identical(result$squared_distance, 0),
    status = result$status, engine_version = result$engine_version,
    finite_intervals = result$first_finite_intervals,
    counts_match = result$first_finite_intervals ==
      result$second_finite_intervals,
    stringsAsFactors = FALSE
  )
})
self_rows <- do.call(rbind, self_rows)

environment <- data.frame(
  contract_id = "mv05bd_environment_v1",
  engine_id = "rust_exact_critical_pairs_segment_tree_v1",
  rust_toolchain = "1.97.1-x86_64-unknown-linux-gnu",
  engine_version = 1L, corpus_results = nrow(equivalence),
  exact_results = sum(equivalence$tier == "D"),
  adaptive_certified_results = sum(equivalence$tier == "E"),
  reverse_results = nrow(equivalence), self_results = nrow(self_rows),
  labels_opened = FALSE, outcomes_computed = FALSE,
  defaults_changed = FALSE, production_adoption_authorized = FALSE,
  stringsAsFactors = FALSE
)

atomic_csv(equivalence, file.path(output_dir, "private-equivalence.csv"))
atomic_csv(self_rows, file.path(output_dir, "private-self.csv"))
atomic_csv(environment, file.path(output_dir, "environment.csv"))

passed <- all(equivalence$equivalent) &&
  all(equivalence$finite_counts_match) &&
  all(equivalence$reverse_bit_identical) &&
  all(equivalence$reverse_counts_swap) &&
  all(equivalence$reverse_diagnostics_match) &&
  all(equivalence$status == 0L) &&
  all(equivalence$engine_version == 1L) &&
  nrow(self_rows) == 112L && all(self_rows$exactly_zero) &&
  all(self_rows$counts_match) && all(self_rows$status == 0L)
if (!passed) stop("MV5-BD full-equivalence gate failed.", call. = FALSE)

cat(
  "MV5-BD full equivalence: D ",
  sum(equivalence$equivalent[equivalence$tier == "D"]), "/318; E ",
  sum(equivalence$equivalent[equivalence$tier == "E"]),
  "/90; reverse ", sum(equivalence$reverse_bit_identical), "/408; self ",
  sum(self_rows$exactly_zero), "/112\n", sep = ""
)
