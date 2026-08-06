#!/usr/bin/env Rscript

options(warn = 2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    paste("Usage: audit_landscape_oracle_and_diagrams.R",
          "<historical-dir> <python-corpus-csv>",
          "<candidate-benchmark-csv> <boundary-benchmark-csv> <output-dir>"),
    call. = FALSE
  )
}
historical_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
python_corpus_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
candidate_benchmark_file <- normalizePath(
  args[[3]], winslash = "/", mustWork = TRUE
)
boundary_benchmark_file <- normalizePath(
  args[[4]], winslash = "/", mustWork = TRUE
)
output_dir <- args[[5]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_dir, "R", "landscape_contract.R"), local = .GlobalEnv)
source(file.path(project_dir, "R", "landscape_reference.R"), local = .GlobalEnv)

diagram <- function(h0 = matrix(numeric(), 0L, 2L),
                    h1 = matrix(numeric(), 0L, 2L)) {
  h0 <- matrix(h0, ncol = 2L)
  h1 <- matrix(h1, ncol = 2L)
  bind_dimension <- function(values, dimension) {
    if (!nrow(values)) return(matrix(numeric(), nrow = 0L, ncol = 3L))
    cbind(dimension, values)
  }
  rbind(bind_dimension(h0, 0), bind_dimension(h1, 1))
}
deep <- cbind(seq_len(12L) / 100, 3 - seq_len(12L) / 100)
fixtures <- list(
  single_tent = list(diagram(c(0, 2)), diagram()),
  translated_tents = list(diagram(c(0, 2)), diagram(c(3, 5))),
  overlapping_dimensions = list(
    diagram(rbind(c(0, 4), c(1, 3)), c(0, 2)),
    diagram(rbind(c(0.5, 3.5), c(1.25, 2.75)), c(0.25, 1.75))
  ),
  narrow_feature = list(diagram(c(0.499, 0.501)), diagram()),
  deep_landscape = list(diagram(deep), diagram()),
  both_dimensions = list(
    diagram(c(0, 2), c(0, 1)),
    diagram(c(0.25, 2.25), c(0.2, 0.9))
  )
)

r_rows <- do.call(rbind, lapply(names(fixtures), function(case) {
  exact <- landscape_reference_distance(
    fixtures[[case]][[1]], fixtures[[case]][[2]], "exact"
  )
  adaptive <- landscape_reference_distance(
    fixtures[[case]][[1]], fixtures[[case]][[2]], "adaptive",
    abs_tol = 1e-10, rel_tol = 1e-10
  )
  do.call(rbind, lapply(c("H0", "H1"), function(dimension) {
    data.frame(
      case = case,
      dimension = dimension,
      r_exact_distance = exact$distances[[dimension]],
      r_adaptive_distance = adaptive$distances[[dimension]],
      r_adaptive_error_estimate = adaptive$dimensions[[dimension]]$
        achieved_absolute_error_estimate,
      stringsAsFactors = FALSE
    )
  }))
}))
python_rows <- utils::read.csv(
  python_corpus_file, stringsAsFactors = FALSE, check.names = FALSE
)
persim <- python_rows[python_rows$engine == "persim_exact_0.3.8", ]
names(persim)[names(persim) == "distance"] <- "persim_exact_distance"
scipy <- python_rows[python_rows$engine == "scipy_quad_1.15.3", ]
names(scipy)[names(scipy) == "distance"] <- "scipy_quad_distance"
persim_corrected <- python_rows[
  python_rows$engine == "persim_0.3.8_corrected_pnorm",
]
names(persim_corrected)[names(persim_corrected) == "distance"] <-
  "persim_corrected_distance"
crosscheck <- merge(
  r_rows,
  persim[, c("case", "dimension", "persim_exact_distance")],
  by = c("case", "dimension"), all = TRUE
)
crosscheck <- merge(
  crosscheck,
  persim_corrected[, c("case", "dimension", "persim_corrected_distance")],
  by = c("case", "dimension"), all = TRUE
)
crosscheck <- merge(
  crosscheck,
  scipy[, c("case", "dimension", "scipy_quad_distance")],
  by = c("case", "dimension"), all = TRUE
)
crosscheck$r_exact_absolute_error_vs_persim <- abs(
  crosscheck$r_exact_distance - crosscheck$persim_exact_distance
)
crosscheck$r_adaptive_absolute_error_vs_persim <- abs(
  crosscheck$r_adaptive_distance - crosscheck$persim_exact_distance
)
crosscheck$r_exact_absolute_error_vs_scipy <- abs(
  crosscheck$r_exact_distance - crosscheck$scipy_quad_distance
)
crosscheck$persim_absolute_error_vs_scipy <- abs(
  crosscheck$persim_exact_distance - crosscheck$scipy_quad_distance
)
crosscheck$persim_corrected_absolute_error_vs_scipy <- abs(
  crosscheck$persim_corrected_distance - crosscheck$scipy_quad_distance
)
gudhi <- python_rows[python_rows$engine == "gudhi_grid_3.12.0", ]
gudhi <- merge(
  gudhi,
  persim[, c("case", "dimension", "persim_exact_distance")],
  by = c("case", "dimension"), all.x = TRUE
)
gudhi$absolute_error_vs_persim <- abs(
  gudhi$distance - gudhi$persim_exact_distance
)
gudhi$relative_error_vs_persim <- ifelse(
  gudhi$persim_exact_distance > 0,
  gudhi$absolute_error_vs_persim / gudhi$persim_exact_distance,
  ifelse(gudhi$absolute_error_vs_persim == 0, 0, Inf)
)

metadata_files <- file.path(
  historical_dir,
  c("joined_metadata_cellcounts.csv", "joined_metadata_cellcounts_bonemarrow.csv")
)
metadata_rows <- lapply(metadata_files[file.exists(metadata_files)], function(path) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
  data.frame(
    sample_id = values[["orig.ident"]],
    filtered_cells = suppressWarnings(as.numeric(
      values[["Number_of_Cells_After_Filtering"]]
    )),
    stringsAsFactors = FALSE
  )
})
metadata <- do.call(rbind, metadata_rows)
cell_lookup <- setNames(metadata$filtered_cells, metadata$sample_id)
pd_files <- list.files(
  historical_dir, pattern = "^PD_list.*[.][Rr][Dd][Ss]$", full.names = TRUE
)
threshold_files <- list.files(
  historical_dir, pattern = "^thresholds.*[.][Rr][Dd][Ss]$", full.names = TRUE
)
threshold_objects <- lapply(threshold_files, readRDS)
names(threshold_objects) <- basename(threshold_files)

detail_rows <- list()
summary_rows <- list()
for (pd_file in pd_files) {
  diagrams <- readRDS(pd_file)
  if (!is.list(diagrams) || is.null(names(diagrams))) next
  artifact_sha256 <- digest::digest(file = pd_file, algo = "sha256")
  artifact_rows <- lapply(names(diagrams), function(sample_id) {
    pd <- as.matrix(diagrams[[sample_id]])
    h0 <- if (ncol(pd) >= 3L) pd[pd[, 1] == 0, , drop = FALSE] else pd[0, ]
    h1 <- if (ncol(pd) >= 3L) pd[pd[, 1] == 1, , drop = FALSE] else pd[0, ]
    cells <- if (sample_id %in% names(cell_lookup)) {
      unname(cell_lookup[[sample_id]])
    } else NA_real_
    threshold_matches <- vapply(names(threshold_objects), function(source) {
      values <- threshold_objects[[source]]
      if (is.null(names(values)) || !sample_id %in% names(values)) return(NA_real_)
      as.numeric(values[[sample_id]])
    }, numeric(1))
    matched <- threshold_matches[is.finite(threshold_matches)]
    lower_bound <- nrow(h0) + 1L
    data.frame(
      artifact = basename(pd_file),
      artifact_sha256 = artifact_sha256,
      sample_id = sample_id,
      pd_rows = nrow(pd),
      h0_intervals = nrow(h0),
      h0_infinite = if (nrow(h0)) sum(!is.finite(h0[, 3])) else 0L,
      h1_intervals = nrow(h1),
      h1_infinite = if (nrow(h1)) sum(!is.finite(h1[, 3])) else 0L,
      inferred_point_count_lower_bound = lower_bound,
      filtered_cells = cells,
      inferred_points_per_cell = if (is.finite(cells) && cells > 0) {
        lower_bound / cells
      } else NA_real_,
      lower_bound_exceeds_filtered_cells = is.finite(cells) && lower_bound > cells,
      threshold_candidate_count = length(matched),
      threshold_candidate_min = if (length(matched)) min(matched) else NA_real_,
      threshold_candidate_max = if (length(matched)) max(matched) else NA_real_,
      code_input_orientation = "feature_by_cell_passed_without_transpose",
      intended_orientation = "cell_points_by_expression_dimensions",
      scientific_reuse_eligibility = "ineligible_orientation_conflict",
      allowed_use = "performance_stress_only",
      stringsAsFactors = FALSE
    )
  })
  artifact <- do.call(rbind, artifact_rows)
  detail_rows[[basename(pd_file)]] <- artifact
  summary_rows[[basename(pd_file)]] <- data.frame(
    artifact = basename(pd_file),
    artifact_sha256 = artifact$artifact_sha256[[1]],
    samples = nrow(artifact),
    samples_with_cell_metadata = sum(is.finite(artifact$filtered_cells)),
    fraction_point_lower_bound_exceeds_cells = mean(
      artifact$lower_bound_exceeds_filtered_cells[is.finite(artifact$filtered_cells)]
    ),
    median_inferred_points_per_cell = stats::median(
      artifact$inferred_points_per_cell, na.rm = TRUE
    ),
    minimum_h0_intervals = min(artifact$h0_intervals),
    maximum_h0_intervals = max(artifact$h0_intervals),
    sample_identifier_provenance = if (all(grepl("^[0-9]+$", artifact$sample_id))) {
      "positional_only"
    } else "stable_sample_ids",
    threshold_provenance = "candidate_files_exist_but_not_bound_to_artifact",
    scientific_reuse_eligibility = "ineligible_orientation_conflict",
    allowed_use = "performance_stress_only",
    stringsAsFactors = FALSE
  )
}
detail <- do.call(rbind, detail_rows)
summary <- do.call(rbind, summary_rows)

summarize_benchmark <- function(path) {
  benchmark <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  benchmark_groups <- split(
    benchmark,
    interaction(benchmark$candidate, benchmark$intervals, drop = TRUE)
  )
  result <- do.call(rbind, lapply(benchmark_groups, function(rows) {
    completed <- rows[rows$status == "completed", , drop = FALSE]
    data.frame(
      candidate = rows$candidate[[1]],
      implementation = rows$implementation[[1]],
      intervals = rows$intervals[[1]],
      maximum_active_levels = rows$maximum_active_levels[[1]],
      repetitions = nrow(rows),
      completed_repetitions = nrow(completed),
      median_end_to_end_seconds = if (nrow(completed)) {
        stats::median(completed$elapsed_seconds)
      } else NA_real_,
      median_kernel_seconds = if (nrow(completed)) {
        stats::median(completed$kernel_elapsed_seconds)
      } else NA_real_,
      maximum_peak_process_tree_rss_bytes = if (nrow(completed)) {
        max(completed$peak_process_tree_rss_bytes)
      } else NA_real_,
      median_distance = if (nrow(completed)) {
        stats::median(completed$distance)
      } else NA_real_,
      maximum_repetition_distance_delta = if (nrow(completed)) {
        diff(range(completed$distance))
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  scipy_reference <- result[
    result$candidate == "scipy_quad", c("intervals", "median_distance")
  ]
  names(scipy_reference)[[2]] <- "scipy_quad_distance"
  result <- merge(result, scipy_reference, by = "intervals", all.x = TRUE)
  result$absolute_distance_error_vs_scipy <- abs(
    result$median_distance - result$scipy_quad_distance
  )
  result$relative_distance_error_vs_scipy <- ifelse(
    result$scipy_quad_distance > 0,
    result$absolute_distance_error_vs_scipy / result$scipy_quad_distance,
    0
  )
  result
}
benchmark_summary <- summarize_benchmark(candidate_benchmark_file)
boundary_benchmark_summary <- summarize_benchmark(boundary_benchmark_file)
cell_counts <- metadata$filtered_cells[is.finite(metadata$filtered_cells) &
  metadata$filtered_cells > 0]
cell_quantiles <- stats::quantile(
  cell_counts, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE
)
cell_scale_summary <- data.frame(
  samples = length(cell_counts),
  minimum_filtered_cells = cell_quantiles[[1]],
  q25_filtered_cells = cell_quantiles[[2]],
  median_filtered_cells = cell_quantiles[[3]],
  q75_filtered_cells = cell_quantiles[[4]],
  maximum_filtered_cells = cell_quantiles[[5]],
  stringsAsFactors = FALSE
)

utils::write.csv(r_rows, file.path(output_dir, "r-landscape-oracle-corpus.csv"),
                 row.names = FALSE)
utils::write.csv(crosscheck,
                 file.path(output_dir, "persim-r-exact-crosscheck.csv"),
                 row.names = FALSE)
utils::write.csv(gudhi,
                 file.path(output_dir, "gudhi-grid-crosscheck.csv"),
                 row.names = FALSE)
utils::write.csv(detail,
                 file.path(output_dir, "historical-diagram-eligibility-detail.csv"),
                 row.names = FALSE)
utils::write.csv(summary,
                 file.path(output_dir, "historical-diagram-eligibility-summary.csv"),
                 row.names = FALSE)
utils::write.csv(
  benchmark_summary,
  file.path(output_dir, "candidate-scaling-benchmark-summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  boundary_benchmark_summary,
  file.path(output_dir, "candidate-production-boundary-summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  cell_scale_summary,
  file.path(output_dir, "corrected-cell-scale-summary.csv"),
  row.names = FALSE
)
writeLines(
  sub("[[:space:]]+$", "", capture.output(sessionInfo())),
  file.path(output_dir, "r-session-info.txt")
)
