#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: audit_legacy_landscape_artifacts.R <legacy-dir> <output.csv>",
    call. = FALSE
  )
}
legacy_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_file <- args[[2]]

cases <- data.frame(
  analysis = c(
    "large_integrated", "large_raw", "large_sct_individual", "large_sct_whole",
    "bone_integrated", "bone_raw", "bone_sct_individual", "bone_sct_whole"
  ),
  landscape_file = c(
    "landscape_list_integrated.Rds",
    "landscape_list_Raw.Rds",
    "landscape_list_sctInd.Rds",
    "landscape_list_sctWhole.Rds",
    "landscape_list_integrated_bonemarrow.Rds",
    "landscape_list_Raw_bonemarrow.Rds",
    "landscape_list_sctInd_bonemarrow.Rds",
    "landscape_list_sctWhole_bonemarrow.Rds"
  ),
  distance_file = c(
    "landscape_l2_distance_matrix_integrated.Rds",
    "landscape_l2_distance_matrix_Raw.Rds",
    "landscape_l2_distance_matrix_sctInd.Rds",
    "landscape_l2_distance_matrix_sctWhole.Rds",
    "landscape_l2_distance_matrix_integrated_bonemarrow.Rds",
    "landscape_l2_distance_matrix_Raw_bonemarrow.Rds",
    "landscape_l2_distance_matrix_sctInd_bonemarrow.Rds",
    "landscape_l2_distance_matrix_sctWhole_bonemarrow.Rds"
  ),
  pd_file = c(
    "PD_list_dim1_th-1_integrated.Rds",
    "PD_list_after_retries_unintegrated_Raw.rds",
    "PD_list_after_retries_unintegrated_sctInd.rds",
    "PD_list_after_retries_unintegrated_sctWhole.rds",
    "PD_list_dim1_th-1_integrated_bonemarrow.Rds",
    "PD_list_dim1_th-1_unintegrated_Raw_bonemarrow.Rds",
    "PD_list_dim1_th-1_unintegrated_sctInd_bonemarrow.Rds",
    "PD_list_dim1_th-1_unintegrated_sctWhole_bonemarrow.Rds"
  ),
  stringsAsFactors = FALSE
)

profile_case <- function(case) {
  paths <- file.path(
    legacy_dir,
    c(case$landscape_file, case$distance_file, case$pd_file)
  )
  if (!all(file.exists(paths))) {
    stop("Missing legacy inputs for ", case$analysis, ": ",
         paste(paths[!file.exists(paths)], collapse = ", "), call. = FALSE)
  }
  landscapes <- readRDS(paths[[1]])
  saved_distance <- as.matrix(readRDS(paths[[2]]))
  diagrams <- readRDS(paths[[3]])
  if (length(landscapes) != length(diagrams) || length(landscapes) < 2L) {
    stop("Landscape/diagram sample counts do not reconcile for ", case$analysis,
         call. = FALSE)
  }

  first_shape <- dim(landscapes[[1]]$dim0)
  shapes <- unique(vapply(landscapes, function(x) {
    paste(dim(x$dim0), collapse = "x")
  }, character(1)))
  if (length(shapes) != 1L) {
    stop("Mixed landscape matrix shapes in ", case$analysis, call. = FALSE)
  }
  orientation <- if (first_shape[[1]] > first_shape[[2]]) {
    "time_by_level"
  } else {
    "level_by_time"
  }
  canonical_first <- if (orientation == "time_by_level") {
    landscapes[[1]]$dim0
  } else {
    t(landscapes[[1]]$dim0)
  }
  correct_l2_curve <- sqrt(rowSums(canonical_first ^ 2))
  legacy_curve <- rowMeans(landscapes[[1]]$dim0)
  legacy_interpolated <- stats::approx(
    seq(0, 1, length.out = length(legacy_curve)),
    legacy_curve,
    xout = seq(0, 1, length.out = length(correct_l2_curve)),
    rule = 2
  )$y

  pairs <- utils::combn(seq_along(landscapes), 2L)
  distances <- apply(pairs, 2L, function(indices) {
    first <- landscapes[[indices[[1]]]]
    second <- landscapes[[indices[[2]]]]
    d0 <- sqrt(sum((first$dim0 - second$dim0) ^ 2))
    d1 <- sqrt(sum((first$dim1 - second$dim1) ^ 2))
    total <- d0 ^ 2 + d1 ^ 2
    c(d0 = d0, d1 = d1, h1_fraction = if (total == 0) NA_real_ else d1 ^ 2 / total)
  })
  direct_pair12 <- sqrt(
    sum((landscapes[[1]]$dim0 - landscapes[[2]]$dim0) ^ 2) +
      sum((landscapes[[1]]$dim1 - landscapes[[2]]$dim1) ^ 2)
  )
  max_deaths <- vapply(diagrams, function(pd) {
    finite_deaths <- pd[, 3][is.finite(pd[, 3])]
    max(finite_deaths, na.rm = TRUE)
  }, numeric(1))
  h1_quantiles <- stats::quantile(
    distances["h1_fraction", ], c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE
  )

  data.frame(
    analysis = case$analysis,
    samples = length(landscapes),
    landscape_shape = shapes,
    inferred_orientation = orientation,
    finite_death_max_min = min(max_deaths),
    finite_death_max_median = stats::median(max_deaths),
    finite_death_max_max = max(max_deaths),
    endpoint_ratio = max(max_deaths) / min(max_deaths),
    saved_pair12 = saved_distance[1, 2],
    direct_frobenius_pair12 = direct_pair12,
    saved_direct_pair12_error = abs(saved_distance[1, 2] - direct_pair12),
    legacy_curve_vs_l2_correlation = stats::cor(
      legacy_interpolated, correct_l2_curve
    ),
    median_h0_distance = stats::median(distances["d0", ], na.rm = TRUE),
    median_h1_distance = stats::median(distances["d1", ], na.rm = TRUE),
    h1_energy_fraction_min = h1_quantiles[[1]],
    h1_energy_fraction_q25 = h1_quantiles[[2]],
    h1_energy_fraction_median = h1_quantiles[[3]],
    h1_energy_fraction_q75 = h1_quantiles[[4]],
    h1_energy_fraction_max = h1_quantiles[[5]],
    stringsAsFactors = FALSE
  )
}

rows <- lapply(seq_len(nrow(cases)), function(i) {
  profile_case(cases[i, , drop = FALSE])
})
profile <- do.call(rbind, rows)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(profile, output_file, row.names = FALSE, na = "")
print(profile)
cat("Wrote:", normalizePath(output_file, winslash = "/", mustWork = TRUE), "\n")
