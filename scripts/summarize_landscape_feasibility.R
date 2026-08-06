#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: summarize_landscape_feasibility.R <profile-dir>", call. = FALSE)
}
profile_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
exact <- utils::read.csv(
  file.path(profile_dir, "landscape-exact-full-convergence.csv"),
  stringsAsFactors = FALSE
)
exact_recommendations <- utils::read.csv(
  file.path(profile_dir, "landscape-exact-full-recommendations.csv"),
  stringsAsFactors = FALSE
)
grid <- utils::read.csv(
  file.path(profile_dir, "landscape-grid-energy-fidelity-summary.csv"),
  stringsAsFactors = FALSE
)

strata <- unique(exact$stratum)
rows <- lapply(strata, function(stratum) {
  exact_group <- exact[exact$stratum == stratum, , drop = FALSE]
  exact_reference <- exact_group[
    exact_group$h0_evaluated_levels == exact_group$h0_exact_full_levels &
      exact_group$h1_evaluated_levels == exact_group$h1_exact_full_levels,
    , drop = FALSE
  ]
  exact_reference <- exact_reference[nrow(exact_reference), , drop = FALSE]
  first_level <- exact_group[exact_group$requested_levels == 1L, , drop = FALSE]
  k25 <- exact_group[exact_group$requested_levels == 25L, , drop = FALSE]
  grid_1000 <- grid[
    grid$stratum == stratum & grid$grid_points == 1000L, , drop = FALSE
  ]
  level_recommendation <- exact_recommendations[
    exact_recommendations$stratum == stratum, , drop = FALSE
  ]
  tested_grid_passes <- nrow(grid_1000) == 2L &&
    all(grid_1000$maximum_absolute_relative_energy_error <= 0.01)
  data.frame(
    stratum = stratum,
    h0_exact_full_levels = exact_reference$h0_exact_full_levels,
    h1_exact_full_levels = exact_reference$h1_exact_full_levels,
    provisional_level_cap_on_250_grid =
      level_recommendation$provisional_levels_on_250_grid,
    k1_spearman_vs_250_grid_full = first_level$combined_spearman,
    k1_relative_error_vs_250_grid_full = first_level$combined_relative_error,
    k25_spearman_vs_250_grid_full = k25$combined_spearman,
    k25_relative_error_vs_250_grid_full = k25$combined_relative_error,
    exact_full_h1_energy_median_on_250_grid = exact_reference$h1_energy_median,
    h0_1000_grid_max_energy_error = grid_1000$maximum_absolute_relative_energy_error[
      grid_1000$dimension == "H0"
    ],
    h1_1000_grid_max_energy_error = grid_1000$maximum_absolute_relative_energy_error[
      grid_1000$dimension == "H1"
    ],
    tested_1000_grid_passes_one_percent_worst_case = tested_grid_passes,
    primary_level_policy = "all_levels_no_universal_cap",
    primary_grid_policy = "error_controlled_adaptive_or_exact_integration",
    compatibility_sensitivity = "k1_common_grid",
    activation_status = "not_ready_grid_method_required",
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, rows)
summary <- summary[order(summary$stratum), , drop = FALSE]
utils::write.csv(
  summary, file.path(profile_dir, "landscape-feasibility-decision-summary.csv"),
  row.names = FALSE
)
print(summary, row.names = FALSE)
