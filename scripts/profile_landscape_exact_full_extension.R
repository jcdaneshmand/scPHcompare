#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage: profile_landscape_exact_full_extension.R",
      "<legacy-dir> <first-pass-dir> <output-dir>"
    ),
    call. = FALSE
  )
}
legacy_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
first_pass_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("The locked R environment must contain Rcpp.", call. = FALSE)
}
source(file.path(project_dir, "R/landscape_contract.R"), local = .GlobalEnv)
source(file.path(project_dir, "R/landscape_convergence.R"), local = .GlobalEnv)
Rcpp::sourceCpp(file.path(project_dir, "scripts/landscape_topk_engine.cpp"))

validation_diagram <- rbind(
  c(0, 0, 4), c(0, 1, 3), c(0, 1.5, 2.5),
  c(1, 0, 2), c(1, 0.5, 1.5)
)
validation_grid <- seq(0, 4, length.out = 101L)
for (validation_dimension in 0:1) {
  accelerated <- top_k_landscape_profile_cpp(
    validation_diagram, validation_dimension, validation_grid, 3L
  )$values
  reference <- compute_landscape_values(
    validation_diagram, validation_dimension, validation_grid, 1:3
  )
  if (max(abs(accelerated - reference)) > 1e-10) {
    stop("The compiled landscape engine disagrees with TDA.", call. = FALSE)
  }
}

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
selected_manifest <- utils::read.csv(
  file.path(first_pass_dir, "landscape-diagnostic-sample-manifest.csv"),
  stringsAsFactors = FALSE
)

tail_tolerance <- 1e-4
rank_tolerance <- 0.999
relative_error_tolerance <- 0.01
clustering_tolerance <- 0.95
base_levels <- c(1L, 5L, 10L, 25L, 100L, 250L, 500L, 1000L,
                 2000L, 5000L, 10000L, 20000L)

current_rss <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  unname(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}

profile_case <- function(case) {
  message("Exact-full extension: ", case$stratum)
  diagrams <- readRDS(file.path(legacy_dir, case$pd_file))
  selected_ids <- selected_manifest$sample_id[
    selected_manifest$stratum == case$stratum
  ]
  if (!all(selected_ids %in% names(diagrams))) {
    stop("Selected diagram IDs do not reconcile for ", case$stratum,
         call. = FALSE)
  }
  selected <- diagrams[selected_ids]
  full_levels <- vapply(0:1, function(dimension) {
    max(vapply(selected, landscape_interval_depth, integer(1),
               dimension = dimension))
  }, integer(1))
  names(full_levels) <- c("H0", "H1")
  full_levels[full_levels < 1L] <- 1L
  maximum_full <- max(full_levels)
  requested_levels <- sort(unique(c(
    base_levels[base_levels < maximum_full], maximum_full
  )))

  energy_rows <- list()
  runtime_rows <- list()
  distance_store <- list()
  for (dimension in 0:1) {
    dimension_name <- paste0("H", dimension)
    grid <- derive_common_landscape_grid(
      diagrams, grid_points = 250L, dimension = dimension
    )
    evaluated_levels <- full_levels[[dimension_name]]
    rss_before <- current_rss()
    started <- proc.time()[["elapsed"]]
    dimension_profiles <- lapply(selected, function(pd) {
      top_k_landscape_profile_cpp(
        as.matrix(pd), dimension, grid, evaluated_levels
      )
    })
    names(dimension_profiles) <- selected_ids
    elapsed <- proc.time()[["elapsed"]] - started
    rss_after <- current_rss()
    values_list <- lapply(dimension_profiles, `[[`, "values")
    object_bytes <- sum(vapply(dimension_profiles, object.size, numeric(1)))
    dimension_requested <- sort(unique(pmin(requested_levels, evaluated_levels)))
    energy_rows[[dimension_name]] <- do.call(rbind, lapply(
      seq_along(dimension_profiles), function(i) {
        capture <- landscape_energy_capture(
          dimension_profiles[[i]]$values,
          dimension_profiles[[i]]$total_squared,
          grid, dimension_requested
        )
        capture$stratum <- case$stratum
        capture$sample_id <- selected_ids[[i]]
        capture$dimension <- dimension_name
        capture$grid_points <- 250L
        capture
      }
    ))
    distance_store[[dimension_name]] <- lapply(
      dimension_requested, function(levels) {
        landscape_distance_matrix_chunked(
          values_list, grid, levels, level_chunk = 250L
        )
      }
    )
    names(distance_store[[dimension_name]]) <- as.character(dimension_requested)
    runtime_rows[[dimension_name]] <- data.frame(
      stratum = case$stratum,
      dimension = dimension_name,
      grid_points = 250L,
      exact_full_levels = evaluated_levels,
      selected_samples = length(selected_ids),
      evaluation_elapsed_seconds = elapsed,
      profile_object_bytes = object_bytes,
      rss_before_bytes = rss_before,
      rss_after_evaluation_bytes = rss_after,
      stringsAsFactors = FALSE
    )
    rm(dimension_profiles, values_list)
    gc(verbose = FALSE)
  }

  energy <- do.call(rbind, energy_rows)
  reference_h0 <- distance_store$H0[[as.character(full_levels[["H0"]])]]
  reference_h1 <- distance_store$H1[[as.character(full_levels[["H1"]])]]
  reference <- combine_landscape_distance_matrices(reference_h0, reference_h1)
  rows <- list()
  clustering_rows <- list()
  for (i in seq_along(requested_levels)) {
    requested <- requested_levels[[i]]
    h0_levels <- min(requested, full_levels[["H0"]])
    h1_levels <- min(requested, full_levels[["H1"]])
    h0 <- distance_store$H0[[as.character(h0_levels)]]
    h1 <- distance_store$H1[[as.character(h1_levels)]]
    combined <- combine_landscape_distance_matrices(h0, h1)
    stability <- landscape_distance_stability(combined, reference)
    clustering <- landscape_clustering_stability(
      combined, reference, cluster_counts = c(2L, 5L, 8L),
      methods = c("average", "ward.D2")
    )
    h0_energy <- energy[
      energy$dimension == "H0" & energy$levels == h0_levels, , drop = FALSE
    ]
    h1_energy <- energy[
      energy$dimension == "H1" & energy$levels == h1_levels, , drop = FALSE
    ]
    contribution <- landscape_h1_energy_summary(h0, h1)
    rows[[i]] <- data.frame(
      stratum = case$stratum,
      grid_points = 250L,
      requested_levels = requested,
      h0_evaluated_levels = h0_levels,
      h1_evaluated_levels = h1_levels,
      h0_exact_full_levels = full_levels[["H0"]],
      h1_exact_full_levels = full_levels[["H1"]],
      combined_spearman = stability$spearman,
      combined_relative_error = stability$relative_frobenius_error,
      h0_max_tail_fraction = max(h0_energy$tail_energy_fraction),
      h1_max_tail_fraction = max(h1_energy$tail_energy_fraction),
      minimum_cluster_ari = min(clustering$adjusted_rand_index),
      median_cluster_ari = stats::median(clustering$adjusted_rand_index),
      minimum_cophenetic_correlation = min(
        clustering$cophenetic_correlation
      ),
      contribution,
      stringsAsFactors = FALSE
    )
    clustering$stratum <- case$stratum
    clustering$requested_levels <- requested
    clustering_rows[[i]] <- clustering
  }
  convergence <- do.call(rbind, rows)
  eligible <- convergence[
    convergence$combined_spearman >= rank_tolerance &
      convergence$combined_relative_error <= relative_error_tolerance &
      convergence$minimum_cluster_ari >= clustering_tolerance &
      convergence$h0_max_tail_fraction <= tail_tolerance &
      convergence$h1_max_tail_fraction <= tail_tolerance,
    , drop = FALSE
  ]
  recommended <- if (nrow(eligible)) min(eligible$requested_levels) else NA_integer_
  recommendation <- data.frame(
    stratum = case$stratum,
    samples_profiled = length(selected_ids),
    h0_exact_full_levels = full_levels[["H0"]],
    h1_exact_full_levels = full_levels[["H1"]],
    provisional_levels_on_250_grid = recommended,
    status = if (is.na(recommended)) {
      "no_level_candidate_met_250_grid_thresholds"
    } else {
      "level_convergence_only_grid_fidelity_required"
    },
    stringsAsFactors = FALSE
  )
  list(
    energy = energy,
    runtime = do.call(rbind, runtime_rows),
    convergence = convergence,
    clustering = do.call(rbind, clustering_rows),
    recommendation = recommendation
  )
}

results <- lapply(seq_len(nrow(cases)), function(i) {
  profile_case(cases[i, , drop = FALSE])
})
bind_result <- function(name) do.call(rbind, lapply(results, `[[`, name))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
outputs <- list(
  "landscape-exact-full-energy.csv" = bind_result("energy"),
  "landscape-exact-full-runtime-memory.csv" = bind_result("runtime"),
  "landscape-exact-full-convergence.csv" = bind_result("convergence"),
  "landscape-exact-full-clustering-detail.csv" = bind_result("clustering"),
  "landscape-exact-full-recommendations.csv" = bind_result("recommendation")
)
for (filename in names(outputs)) {
  utils::write.csv(
    outputs[[filename]], file.path(output_dir, filename), row.names = FALSE,
    na = ""
  )
}
print(outputs[["landscape-exact-full-recommendations.csv"]])
cat("Wrote exact-full extension to:",
    normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
