#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste(
      "Usage: profile_landscape_convergence.R",
      "<legacy-dir> <output-dir> [max-reference-k]"
    ),
    call. = FALSE
  )
}

legacy_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
max_reference_k <- if (length(args) >= 3L) as.integer(args[[3]]) else 1000L
if (is.na(max_reference_k) || max_reference_k < 25L) {
  stop("max-reference-k must be an integer of at least 25.", call. = FALSE)
}

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c(
  "R/landscape_contract.R", "R/landscape_convergence.R",
  "inst/extdata/inputs/metadata_MultiTissueAnalysis.csv"
)
if (!all(file.exists(file.path(project_dir, required_project_files)))) {
  stop("Run this script from the scPHcompare repository root.", call. = FALSE)
}
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("The locked R environment must contain Rcpp.", call. = FALSE)
}

source(file.path(project_dir, "R/landscape_contract.R"), local = .GlobalEnv)
source(file.path(project_dir, "R/landscape_convergence.R"), local = .GlobalEnv)

Rcpp::cppFunction(
  plugins = "cpp11",
  code = '
  Rcpp::List top_k_landscape_profile(
      Rcpp::NumericMatrix pd, int dimension, Rcpp::NumericVector grid,
      int max_k) {
    if (pd.ncol() < 3 || max_k < 1) {
      Rcpp::stop("pd must have three columns and max_k must be positive.");
    }
    std::vector<double> births;
    std::vector<double> deaths;
    births.reserve(pd.nrow());
    deaths.reserve(pd.nrow());
    for (int i = 0; i < pd.nrow(); ++i) {
      const double dim = pd(i, 0);
      const double birth = pd(i, 1);
      const double death = pd(i, 2);
      if (R_finite(dim) && static_cast<int>(dim) == dimension &&
          R_finite(birth) && R_finite(death) && birth < death) {
        births.push_back(birth);
        deaths.push_back(death);
      }
    }
    Rcpp::NumericMatrix values(grid.size(), max_k);
    Rcpp::NumericVector total_squared(grid.size());
    Rcpp::IntegerVector active_levels(grid.size());
    std::vector<double> tents;
    tents.reserve(births.size());
    for (int g = 0; g < grid.size(); ++g) {
      tents.clear();
      double total = 0.0;
      const double t = grid[g];
      for (std::size_t i = 0; i < births.size(); ++i) {
        const double height = std::min(t - births[i], deaths[i] - t);
        if (height > 0.0) {
          tents.push_back(height);
          total += height * height;
        }
      }
      active_levels[g] = static_cast<int>(tents.size());
      total_squared[g] = total;
      if (tents.size() > static_cast<std::size_t>(max_k)) {
        std::nth_element(
          tents.begin(), tents.begin() + max_k, tents.end(),
          std::greater<double>()
        );
        tents.resize(max_k);
      }
      std::sort(tents.begin(), tents.end(), std::greater<double>());
      for (std::size_t k = 0; k < tents.size(); ++k) values(g, k) = tents[k];
    }
    return Rcpp::List::create(
      Rcpp::Named("values") = values,
      Rcpp::Named("total_squared") = total_squared,
      Rcpp::Named("active_levels") = active_levels,
      Rcpp::Named("finite_intervals") = static_cast<int>(births.size())
    );
  }'
)

verify_engine <- function() {
  diagram <- rbind(
    c(0, 0, 4), c(0, 1, 3), c(0, 1.5, 2.5),
    c(1, 0, 2), c(1, 0.5, 1.5)
  )
  grid <- seq(0, 4, length.out = 101L)
  rows <- lapply(0:1, function(dimension) {
    evaluated <- top_k_landscape_profile(diagram, dimension, grid, 3L)
    expected <- compute_landscape_values(
      diagram, dimension, grid, levels = 1:3
    )
    data.frame(
      dimension = paste0("H", dimension),
      grid_points = length(grid),
      levels = 3L,
      maximum_absolute_error = max(abs(evaluated$values - expected)),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (any(result$maximum_absolute_error > 1e-10)) {
    stop("Accelerated evaluator disagrees with TDA analytical fixture.",
         call. = FALSE)
  }
  result
}

cases <- data.frame(
  stratum = c(
    "large_integrated", "large_raw", "large_sct_individual",
    "large_sct_whole", "bone_integrated", "bone_raw",
    "bone_sct_individual", "bone_sct_whole"
  ),
  dataset = c(rep("large", 4L), rep("bone", 4L)),
  representation = rep(
    c("integrated", "raw", "sct_individual", "sct_whole"), 2L
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

candidate_grid_points <- c(250L, 500L, 1000L)
candidate_levels <- c(1L, 5L, 10L, 25L, 100L, 250L, 500L, 1000L)
tail_tolerance <- 1e-4
rank_tolerance <- 0.999
relative_error_tolerance <- 0.01
clustering_tolerance <- 0.95

metadata <- utils::read.csv(
  file.path(project_dir, "inst/extdata/inputs/metadata_MultiTissueAnalysis.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
metadata$diagram_id <- paste(metadata$SRA, metadata[["Sample Name"]], sep = "_")

evenly_spaced <- function(values, count) {
  values <- sort(unique(values))
  if (length(values) <= count) return(values)
  values[unique(as.integer(round(seq(1, length(values), length.out = count))))]
}

select_diagnostic_samples <- function(diagrams, dataset) {
  ids <- names(diagrams)
  if (is.null(ids)) stop("Historical diagram lists must be named.", call. = FALSE)
  if (dataset == "bone") return(sort(ids))
  matched <- metadata[metadata$diagram_id %in% ids, , drop = FALSE]
  selected <- unlist(lapply(split(matched$diagram_id, matched$Tissue), function(x) {
    evenly_spaced(x, 2L)
  }), use.names = FALSE)
  selected <- sort(unique(selected))
  if (length(selected) != 16L) {
    stop("Expected two matched diagnostic samples from each of eight tissues.",
         call. = FALSE)
  }
  selected
}

current_rss <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  unname(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}

profile_dimension_grid <- function(diagrams, selected_ids, dimension, grid,
                                   evaluated_levels, stratum) {
  profiles <- vector("list", length(selected_ids))
  names(profiles) <- selected_ids
  sample_runtime <- numeric(length(selected_ids))
  rss_before <- current_rss()
  for (i in seq_along(selected_ids)) {
    started <- proc.time()[["elapsed"]]
    profiles[[i]] <- top_k_landscape_profile(
      as.matrix(diagrams[[selected_ids[[i]]]]), dimension, grid,
      evaluated_levels
    )
    sample_runtime[[i]] <- proc.time()[["elapsed"]] - started
  }
  rss_after <- current_rss()
  values_list <- lapply(profiles, `[[`, "values")
  profile_bytes <- sum(vapply(profiles, object.size, numeric(1)))
  available_levels <- candidate_levels[candidate_levels <= evaluated_levels]
  if (!evaluated_levels %in% available_levels) {
    available_levels <- sort(unique(c(available_levels, evaluated_levels)))
  }
  energy <- do.call(rbind, lapply(seq_along(profiles), function(i) {
    capture <- landscape_energy_capture(
      profiles[[i]]$values, profiles[[i]]$total_squared, grid,
      available_levels
    )
    capture$stratum <- stratum
    capture$sample_id <- selected_ids[[i]]
    capture$dimension <- paste0("H", dimension)
    capture$grid_points <- length(grid)
    capture
  }))
  distances <- lapply(available_levels, function(level) {
    landscape_distance_matrix_from_values(values_list, grid, level)
  })
  names(distances) <- as.character(available_levels)
  runtime <- data.frame(
    stratum = stratum,
    dimension = paste0("H", dimension),
    grid_points = length(grid),
    evaluated_levels = evaluated_levels,
    selected_samples = length(selected_ids),
    total_elapsed_seconds = sum(sample_runtime),
    median_sample_seconds = stats::median(sample_runtime),
    maximum_sample_seconds = max(sample_runtime),
    profile_object_bytes = profile_bytes,
    rss_before_bytes = rss_before,
    rss_after_evaluation_bytes = rss_after,
    stringsAsFactors = FALSE
  )
  list(distances = distances, energy = energy, runtime = runtime)
}

summarize_case <- function(case) {
  message("Profiling ", case$stratum, " ...")
  diagrams <- readRDS(file.path(legacy_dir, case$pd_file))
  if (!is.list(diagrams) || length(diagrams) < 2L) {
    stop("Invalid diagram list for ", case$stratum, call. = FALSE)
  }
  selected_ids <- select_diagnostic_samples(diagrams, case$dataset)
  inventory <- landscape_interval_inventory(diagrams, case$stratum)
  selected_manifest <- data.frame(
    stratum = case$stratum,
    dataset = case$dataset,
    representation = case$representation,
    sample_id = selected_ids,
    selection = if (case$dataset == "large") {
      "deterministic_two_per_tissue"
    } else {
      "all_samples"
    },
    tissue = if (case$dataset == "large") {
      metadata$Tissue[match(selected_ids, metadata$diagram_id)]
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )

  maximum_levels <- vapply(0:1, function(dimension) {
    selected_inventory <- inventory[
      inventory$sample_id %in% selected_ids &
        inventory$dimension == paste0("H", dimension), , drop = FALSE
    ]
    min(max_reference_k, max(selected_inventory$maximum_active_levels))
  }, integer(1))
  names(maximum_levels) <- c("H0", "H1")
  maximum_levels[maximum_levels < 1L] <- 1L

  profiles <- list()
  energy_rows <- list()
  runtime_rows <- list()
  for (grid_points in candidate_grid_points) {
    for (dimension in 0:1) {
      grid <- derive_common_landscape_grid(
        diagrams, grid_points = grid_points, dimension = dimension
      )
      key <- paste0("G", grid_points, "_H", dimension)
      profiles[[key]] <- profile_dimension_grid(
        diagrams, selected_ids, dimension, grid,
        maximum_levels[[paste0("H", dimension)]], case$stratum
      )
      energy_rows[[key]] <- profiles[[key]]$energy
      runtime_rows[[key]] <- profiles[[key]]$runtime
      gc(verbose = FALSE)
    }
  }
  energy <- do.call(rbind, energy_rows)
  runtime <- do.call(rbind, runtime_rows)

  reference_h0 <- profiles[["G1000_H0"]]$distances[[
    as.character(maximum_levels[["H0"]])
  ]]
  reference_h1 <- profiles[["G1000_H1"]]$distances[[
    as.character(maximum_levels[["H1"]])
  ]]
  reference_combined <- combine_landscape_distance_matrices(
    reference_h0, reference_h1
  )

  convergence_rows <- list()
  clustering_rows <- list()
  row_index <- 0L
  cluster_index <- 0L
  for (grid_points in candidate_grid_points) {
    for (requested_levels in candidate_levels) {
      h0_levels <- min(requested_levels, maximum_levels[["H0"]])
      h1_levels <- min(requested_levels, maximum_levels[["H1"]])
      h0_key <- paste0("G", grid_points, "_H0")
      h1_key <- paste0("G", grid_points, "_H1")
      if (!as.character(h0_levels) %in% names(profiles[[h0_key]]$distances) ||
          !as.character(h1_levels) %in% names(profiles[[h1_key]]$distances)) {
        next
      }
      h0 <- profiles[[h0_key]]$distances[[as.character(h0_levels)]]
      h1 <- profiles[[h1_key]]$distances[[as.character(h1_levels)]]
      combined <- combine_landscape_distance_matrices(h0, h1)
      stability_h0 <- landscape_distance_stability(h0, reference_h0)
      stability_h1 <- landscape_distance_stability(h1, reference_h1)
      stability_combined <- landscape_distance_stability(
        combined, reference_combined
      )
      clustering <- landscape_clustering_stability(
        combined, reference_combined,
        cluster_counts = c(2L, 5L, 8L),
        methods = c("average", "ward.D2")
      )
      contribution <- landscape_h1_energy_summary(h0, h1)
      h0_energy <- energy[
        energy$grid_points == grid_points & energy$dimension == "H0" &
          energy$levels == h0_levels, , drop = FALSE
      ]
      h1_energy <- energy[
        energy$grid_points == grid_points & energy$dimension == "H1" &
          energy$levels == h1_levels, , drop = FALSE
      ]
      row_index <- row_index + 1L
      convergence_rows[[row_index]] <- data.frame(
        stratum = case$stratum,
        grid_points = grid_points,
        requested_levels = requested_levels,
        h0_evaluated_levels = h0_levels,
        h1_evaluated_levels = h1_levels,
        h0_reference_levels = maximum_levels[["H0"]],
        h1_reference_levels = maximum_levels[["H1"]],
        h0_spearman = stability_h0$spearman,
        h1_spearman = stability_h1$spearman,
        combined_spearman = stability_combined$spearman,
        h0_relative_error = stability_h0$relative_frobenius_error,
        h1_relative_error = stability_h1$relative_frobenius_error,
        combined_relative_error = stability_combined$relative_frobenius_error,
        h0_max_tail_fraction = max(h0_energy$tail_energy_fraction),
        h0_median_tail_fraction = stats::median(h0_energy$tail_energy_fraction),
        h1_max_tail_fraction = max(h1_energy$tail_energy_fraction),
        h1_median_tail_fraction = stats::median(h1_energy$tail_energy_fraction),
        minimum_cluster_ari = min(clustering$adjusted_rand_index),
        median_cluster_ari = stats::median(clustering$adjusted_rand_index),
        minimum_cophenetic_correlation = min(
          clustering$cophenetic_correlation
        ),
        contribution,
        stringsAsFactors = FALSE
      )
      cluster_index <- cluster_index + 1L
      clustering$stratum <- case$stratum
      clustering$grid_points <- grid_points
      clustering$requested_levels <- requested_levels
      clustering_rows[[cluster_index]] <- clustering
    }
  }
  convergence <- do.call(rbind, convergence_rows)
  clustering <- do.call(rbind, clustering_rows)

  reference_row <- convergence[
    convergence$grid_points == 1000L &
      convergence$requested_levels == 1000L &
      convergence$h0_evaluated_levels == maximum_levels[["H0"]] &
      convergence$h1_evaluated_levels == maximum_levels[["H1"]],
    , drop = FALSE
  ]
  reference_certified <- nrow(reference_row) == 1L &&
    reference_row$h0_max_tail_fraction <= tail_tolerance &&
    reference_row$h1_max_tail_fraction <= tail_tolerance

  level_candidates <- convergence[
    convergence$grid_points == 1000L &
      convergence$combined_spearman >= rank_tolerance &
      convergence$combined_relative_error <= relative_error_tolerance &
      convergence$minimum_cluster_ari >= clustering_tolerance &
      convergence$h0_max_tail_fraction <= tail_tolerance &
      convergence$h1_max_tail_fraction <= tail_tolerance,
    , drop = FALSE
  ]
  chosen_levels <- if (reference_certified && nrow(level_candidates)) {
    min(level_candidates$requested_levels)
  } else {
    NA_integer_
  }
  grid_level <- if (is.na(chosen_levels)) max(candidate_levels) else chosen_levels
  grid_candidates <- convergence[
    convergence$requested_levels == grid_level &
      convergence$combined_spearman >= rank_tolerance &
      convergence$combined_relative_error <= relative_error_tolerance &
      convergence$minimum_cluster_ari >= clustering_tolerance,
    , drop = FALSE
  ]
  chosen_grid <- if (nrow(grid_candidates)) {
    min(grid_candidates$grid_points)
  } else {
    NA_integer_
  }
  recommendation <- data.frame(
    stratum = case$stratum,
    samples_total = length(diagrams),
    samples_profiled = length(selected_ids),
    h0_reference_levels = maximum_levels[["H0"]],
    h1_reference_levels = maximum_levels[["H1"]],
    h0_reference_max_tail_fraction = reference_row$h0_max_tail_fraction,
    h1_reference_max_tail_fraction = reference_row$h1_max_tail_fraction,
    reference_certified_effectively_full = reference_certified,
    recommended_levels = chosen_levels,
    recommended_grid_points = chosen_grid,
    status = if (reference_certified && !is.na(chosen_levels) &&
        !is.na(chosen_grid)) {
      "methodologically_converged_on_diagnostic_samples"
    } else if (!reference_certified) {
      "reference_cap_not_effectively_full"
    } else {
      "candidate_thresholds_not_met"
    },
    limitation = paste(
      "Diagnostic sensitivity only; diagram eligibility remains conditional",
      "on observation-unit and filtration audits."
    ),
    stringsAsFactors = FALSE
  )

  list(
    inventory = inventory,
    selected = selected_manifest,
    energy = energy,
    runtime = runtime,
    convergence = convergence,
    clustering = clustering,
    recommendation = recommendation
  )
}

engine_validation <- verify_engine()
results <- lapply(seq_len(nrow(cases)), function(i) {
  summarize_case(cases[i, , drop = FALSE])
})

bind_result <- function(name) {
  do.call(rbind, lapply(results, `[[`, name))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
outputs <- list(
  "landscape-engine-validation.csv" = engine_validation,
  "landscape-interval-level-inventory.csv" = bind_result("inventory"),
  "landscape-diagnostic-sample-manifest.csv" = bind_result("selected"),
  "landscape-energy-capture.csv" = bind_result("energy"),
  "landscape-runtime-memory-profile.csv" = bind_result("runtime"),
  "landscape-distance-clustering-convergence.csv" = bind_result("convergence"),
  "landscape-clustering-stability-detail.csv" = bind_result("clustering"),
  "landscape-stratum-recommendations.csv" = bind_result("recommendation")
)
for (filename in names(outputs)) {
  utils::write.csv(
    outputs[[filename]], file.path(output_dir, filename), row.names = FALSE,
    na = ""
  )
}

input_manifest <- data.frame(
  stratum = cases$stratum,
  pd_file = cases$pd_file,
  sha256 = vapply(file.path(legacy_dir, cases$pd_file), function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1)),
  max_reference_k = max_reference_k,
  candidate_grid_points = paste(candidate_grid_points, collapse = ";"),
  candidate_levels = paste(candidate_levels, collapse = ";"),
  tail_tolerance = tail_tolerance,
  rank_tolerance = rank_tolerance,
  relative_error_tolerance = relative_error_tolerance,
  clustering_tolerance = clustering_tolerance,
  stringsAsFactors = FALSE
)
utils::write.csv(
  input_manifest, file.path(output_dir, "landscape-profile-input-manifest.csv"),
  row.names = FALSE
)

print(outputs[["landscape-stratum-recommendations.csv"]])
cat("Wrote convergence profile to:",
    normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
