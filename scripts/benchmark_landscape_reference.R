#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste("Usage: benchmark_landscape_reference.R",
          "<historical-pd-rds> <output-dir>"),
    call. = FALSE
  )
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
historical_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_dir, "R/landscape_contract.R"), local = .GlobalEnv)
source(file.path(project_dir, "R/landscape_reference.R"), local = .GlobalEnv)

empty_diagram <- matrix(numeric(), nrow = 0L, ncol = 3L)
diagram <- function(...) matrix(c(...), ncol = 3L, byrow = TRUE)

timed_reference <- function(case, first, second, method, repetitions = 2L,
                            exact_max_intervals = 200L) {
  invisible(landscape_reference_distance(
    first, second, method = method,
    exact_max_intervals = exact_max_intervals,
    abs_tol = 1e-8, rel_tol = 1e-8
  ))
  runs <- vector("list", repetitions)
  elapsed <- numeric(repetitions)
  for (index in seq_len(repetitions)) {
    timing <- system.time({
      runs[[index]] <- landscape_reference_distance(
        first, second, method = method,
        exact_max_intervals = exact_max_intervals,
        abs_tol = 1e-8, rel_tol = 1e-8
      )
    })
    elapsed[[index]] <- unname(timing[["elapsed"]])
  }
  distances <- do.call(rbind, lapply(runs, `[[`, "distances"))
  data.frame(
    case = case,
    method = method,
    repetitions = repetitions,
    elapsed_median_seconds = stats::median(elapsed),
    elapsed_max_seconds = max(elapsed),
    result_bytes = as.numeric(utils::object.size(runs[[1]])),
    H0_distance = distances[1L, "H0"],
    H1_distance = distances[1L, "H1"],
    combined_distance = distances[1L, "combined"],
    maximum_repetition_delta = max(abs(
      sweep(distances, 2L, distances[1L, ], "-")
    )),
    H0_error_estimate = runs[[1]]$dimensions$H0$achieved_absolute_error_estimate,
    H1_error_estimate = runs[[1]]$dimensions$H1$achieved_absolute_error_estimate,
    stringsAsFactors = FALSE
  )
}

analytical_cases <- list(
  single_tent = list(diagram(0, 0, 2), empty_diagram),
  translated_tents = list(diagram(0, 0, 2), diagram(0, 3, 5)),
  overlapping_dimensions = list(
    diagram(0, 0, 4, 0, 1, 3, 1, 0, 2),
    diagram(0, 0.5, 3.5, 0, 1.25, 2.75, 1, 0.25, 1.75)
  ),
  narrow_feature = list(diagram(0, 0.499, 0.501), empty_diagram),
  deep_landscape = list(
    do.call(rbind, lapply(seq_len(25L), function(index) {
      c(0, index / 100, 3 - index / 100)
    })),
    empty_diagram
  )
)

rows <- list()
for (case in names(analytical_cases)) {
  inputs <- analytical_cases[[case]]
  rows[[paste(case, "exact")]] <- timed_reference(
    paste0("analytical_", case), inputs[[1]], inputs[[2]], "exact"
  )
  rows[[paste(case, "adaptive")]] <- timed_reference(
    paste0("analytical_", case), inputs[[1]], inputs[[2]], "adaptive"
  )
}

diagrams <- readRDS(historical_file)
if (!is.list(diagrams) || is.null(names(diagrams))) {
  stop("Historical persistence-diagram artifact must be a named list.",
       call. = FALSE)
}
finite_counts <- vapply(diagrams, function(pd) {
  sum(vapply(0:1, function(dimension) {
    nrow(finite_landscape_diagram(pd, dimension))
  }, integer(1)))
}, integer(1))
eligible <- names(finite_counts)[finite_counts > 0L & finite_counts <= 200L]
eligible <- eligible[order(finite_counts[eligible], eligible)]
if (length(eligible) < 2L) {
  stop("Need two complete existing diagrams with at most 200 finite intervals.",
       call. = FALSE)
}
existing_ids <- eligible[1:2]
rows$existing_full_exact <- timed_reference(
  "existing_full_small_pair", diagrams[[existing_ids[[1]]]],
  diagrams[[existing_ids[[2]]]], "exact"
)
rows$existing_full_adaptive <- timed_reference(
  "existing_full_small_pair", diagrams[[existing_ids[[1]]]],
  diagrams[[existing_ids[[2]]]], "adaptive"
)

quantile_subset <- function(pd, per_dimension = 25L) {
  pieces <- lapply(0:1, function(dimension) {
    values <- finite_landscape_diagram(pd, dimension)
    if (nrow(values) <= per_dimension) return(values)
    persistence <- values[, 3] - values[, 2]
    ordered <- values[order(persistence, values[, 2], values[, 3]), , drop = FALSE]
    indices <- unique(as.integer(round(seq(
      1, nrow(ordered), length.out = per_dimension
    ))))
    ordered[indices, , drop = FALSE]
  })
  do.call(rbind, pieces)
}
large_ids <- names(sort(finite_counts, decreasing = TRUE))[1:2]
large_subsets <- lapply(diagrams[large_ids], quantile_subset)
rows$existing_quantile_subset_exact <- timed_reference(
  "existing_large_quantile_subset_25_per_dimension",
  large_subsets[[1]], large_subsets[[2]], "exact"
)
rows$existing_quantile_subset_adaptive <- timed_reference(
  "existing_large_quantile_subset_25_per_dimension",
  large_subsets[[1]], large_subsets[[2]], "adaptive"
)

benchmark <- do.call(rbind, rows)
exact_lookup <- benchmark[benchmark$method == "exact",
                          c("case", "combined_distance")]
names(exact_lookup)[[2]] <- "exact_combined_distance"
benchmark <- merge(benchmark, exact_lookup, by = "case", all.x = TRUE,
                   sort = FALSE)
benchmark$absolute_distance_error_vs_exact <- abs(
  benchmark$combined_distance - benchmark$exact_combined_distance
)

manifest <- data.frame(
  historical_file = basename(historical_file),
  historical_file_md5 = unname(tools::md5sum(historical_file)),
  existing_small_first = existing_ids[[1]],
  existing_small_second = existing_ids[[2]],
  existing_large_first = large_ids[[1]],
  existing_large_second = large_ids[[2]],
  large_subset_policy = "25 persistence-quantile intervals per dimension",
  warmup_repetitions = 1L,
  timed_repetitions = 2L,
  exact_max_intervals = 200L,
  abs_tol = 1e-8,
  rel_tol = 1e-8,
  scientific_interpretation = FALSE,
  stringsAsFactors = FALSE
)

utils::write.csv(
  benchmark,
  file.path(output_dir, "landscape-reference-benchmark.csv"),
  row.names = FALSE
)
utils::write.csv(
  manifest,
  file.path(output_dir, "landscape-reference-benchmark-manifest.csv"),
  row.names = FALSE
)
writeLines(
  sub("[[:space:]]+$", "", capture.output(sessionInfo())),
  file.path(output_dir, "landscape-reference-session-info.txt")
)
