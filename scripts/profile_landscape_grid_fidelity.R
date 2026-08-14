#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: profile_landscape_grid_fidelity.R",
      "<legacy-dir> <first-pass-dir> <output-dir>"
    ),
    call. = FALSE
  )
}
legacy_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
first_pass_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_dir, "R/landscape_contract.R"), local = .GlobalEnv)
source(file.path(project_dir, "R/landscape_convergence.R"), local = .GlobalEnv)

cases <- data.frame(
  stratum = c(
    "large_integrated", "large_raw", "large_sct_individual",
    "large_sct_whole", "bone_integrated", "bone_raw",
    "bone_sct_individual", "bone_sct_whole"
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
selected <- utils::read.csv(
  file.path(first_pass_dir, "landscape-diagnostic-sample-manifest.csv"),
  stringsAsFactors = FALSE
)
sampled <- utils::read.csv(
  file.path(first_pass_dir, "landscape-energy-capture.csv"),
  stringsAsFactors = FALSE
)
sampled <- sampled[!duplicated(sampled[c(
  "stratum", "sample_id", "dimension", "grid_points"
)]), c("stratum", "sample_id", "dimension", "grid_points", "total_energy")]
names(sampled)[names(sampled) == "total_energy"] <- "sampled_total_energy"

rows <- list()
row_index <- 0L
for (case_index in seq_len(nrow(cases))) {
  case <- cases[case_index, , drop = FALSE]
  message("Grid-fidelity audit: ", case$stratum)
  diagrams <- readRDS(file.path(legacy_dir, case$pd_file))
  selected_ids <- selected$sample_id[selected$stratum == case$stratum]
  for (dimension in 0:1) {
    dimension_name <- paste0("H", dimension)
    for (grid_points in c(250L, 500L, 1000L)) {
      grid <- derive_common_landscape_grid(
        diagrams, grid_points = grid_points, dimension = dimension
      )
      grid_step <- stats::median(diff(grid))
      for (sample_id in selected_ids) {
        finite_pd <- finite_landscape_diagram(
          diagrams[[sample_id]], dimension
        )
        persistence <- finite_pd[, 3] - finite_pd[, 2]
        interval_energy <- persistence ^ 3 / 12
        exact_total <- landscape_exact_total_energy(
          diagrams[[sample_id]], dimension
        )
        sampled_match <- sampled[
          sampled$stratum == case$stratum &
            sampled$sample_id == sample_id &
            sampled$dimension == dimension_name &
            sampled$grid_points == grid_points,
          "sampled_total_energy"
        ]
        if (length(sampled_match) != 1L) {
          stop("Sampled energy did not reconcile for ", case$stratum, "/",
               sample_id, "/", dimension_name, "/", grid_points,
               call. = FALSE)
        }
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          stratum = case$stratum,
          sample_id = sample_id,
          dimension = dimension_name,
          grid_points = grid_points,
          grid_minimum = min(grid),
          grid_maximum = max(grid),
          grid_step = grid_step,
          finite_intervals = length(persistence),
          intervals_narrower_than_step = sum(persistence < grid_step),
          narrow_interval_fraction = if (length(persistence)) {
            mean(persistence < grid_step)
          } else 0,
          narrow_exact_energy_fraction = if (exact_total > 0) {
            sum(interval_energy[persistence < grid_step]) / exact_total
          } else 0,
          exact_total_energy = exact_total,
          sampled_total_energy = sampled_match,
          signed_relative_energy_error = if (exact_total > 0) {
            (sampled_match - exact_total) / exact_total
          } else 0,
          absolute_relative_energy_error = if (exact_total > 0) {
            abs(sampled_match - exact_total) / exact_total
          } else 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
detail <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(
  split(detail, list(detail$stratum, detail$dimension, detail$grid_points),
        drop = TRUE),
  function(group) {
    data.frame(
      stratum = group$stratum[[1]],
      dimension = group$dimension[[1]],
      grid_points = group$grid_points[[1]],
      samples = nrow(group),
      grid_step = group$grid_step[[1]],
      median_absolute_relative_energy_error = stats::median(
        group$absolute_relative_energy_error
      ),
      maximum_absolute_relative_energy_error = max(
        group$absolute_relative_energy_error
      ),
      median_narrow_interval_fraction = stats::median(
        group$narrow_interval_fraction
      ),
      maximum_narrow_interval_fraction = max(group$narrow_interval_fraction),
      median_narrow_exact_energy_fraction = stats::median(
        group$narrow_exact_energy_fraction
      ),
      maximum_narrow_exact_energy_fraction = max(
        group$narrow_exact_energy_fraction
      ),
      stringsAsFactors = FALSE
    )
  }
))
summary <- summary[order(summary$stratum, summary$dimension,
                         summary$grid_points), , drop = FALSE]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  detail, file.path(output_dir, "landscape-grid-energy-fidelity-detail.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary,
  file.path(output_dir, "landscape-grid-energy-fidelity-summary.csv"),
  row.names = FALSE
)
print(summary, row.names = FALSE)
