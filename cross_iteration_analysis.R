# cross_iteration_analysis.R
# Description: Functions for comparing results (Betti curves, cluster metrics) across different data processing iterations.

# ---------------------------
# Libraries 
# ---------------------------
# Include all libraries needed by the functions in this file
library(Seurat) # Needed for Seurat object access within helpers
library(dplyr)
library(ggplot2)
library(mclust) # adjustedRandIndex
library(aricode) # NMI
library(entropy) # variation_of_information
library(parallel) # mclapply
library(reshape2)
library(gridExtra) # grid.arrange
library(digest) # Hashing for caching
library(transport) # Wasserstein distance
library(TDA) # Betti curve computation might use this or TDAstats
library(TDAstats) # Betti curve computation might use this or TDA
library(cowplot) # plot_grid or similar might be used implicitly
library(svglite) # For saving SVG plots

# ---------------------------
# General Helper Functions (Consider moving to a separate utils.R file)
# ---------------------------

# Simple logging function.
log_message <- function(msg, task = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  prefix <- if (!is.null(task)) sprintf("[%s] Task: %s -", timestamp, task) else sprintf("[%s]", timestamp)
  cat(paste(prefix, msg, "\n"))
}

# --- Helper: Directory Creation ---
ensure_directory <- function(path) {
  if (length(path) == 0 || is.na(path) || path == "") {
    stop("Invalid directory path provided.")
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# --- Helper: Normalization with Division-by-Zero Check ---
normalize_curve <- function(curve, grid_points) {
  tot <- sum(curve, na.rm = TRUE) # Added na.rm=TRUE
  if (tot == 0 || is.na(tot)) rep(0, grid_points) else curve / tot
}

# --- Helper: Integrated Difference (Trapezoidal Rule) ---
integrated_diff <- function(curve_diff, tau_vals) {
  if (length(curve_diff) < 2 || length(tau_vals) < 2) return(0)
  # Ensure lengths match for diff()
  if (length(curve_diff) != length(tau_vals)) {
    log_message("Interpolating curve_diff to match tau_vals length in integrated_diff.")
    # Simple linear interpolation
    curve_diff <- approx(x = seq(0, 1, length.out = length(curve_diff)),
                          y = curve_diff,
                          xout = seq(0, 1, length.out = length(tau_vals)),
                          rule = 2)$y
  }
  delta_tau <- diff(tau_vals)
  # Ensure vectors used for trapezoidal rule have compatible lengths
  term1 <- curve_diff[-length(curve_diff)] ^ 2
  term2 <- curve_diff[-1] ^ 2
  if (length(delta_tau) != length(term1)) {
    log_message(paste("Length mismatch in integrated_diff: delta_tau=", length(delta_tau), "terms=", length(term1)))
    # Fallback or error handling needed - here, return NA
    return(NA)
  }
  sqrt(sum(delta_tau * (term1 + term2) / 2, na.rm = TRUE))
}


# ---------------------------
# Betti Curve Analysis Helpers
# ---------------------------

# --- Curve Caching ---
generate_cache_filename <- function(dataset_name, hash, results_folder) {
  # Define cache path relative to the provided results_folder
  cache_base <- file.path(results_folder, "plots", "betti_plots", "betti_cache")
  path <- file.path(cache_base, dataset_name, paste0("cache_", hash, ".rds"))
  ensure_directory(dirname(path))
  path
}

load_curve_cache <- function(dataset_name, hash, results_folder) {
  cache_file <- generate_cache_filename(dataset_name, hash, results_folder)
  if (file.exists(cache_file)) {
    log_message(paste("Loaded cached curve from:", cache_file))
    return(readRDS(cache_file))
  }
  NULL
}

save_curve_cache <- function(curve, dataset_name, hash, results_folder) {
  cache_file <- generate_cache_filename(dataset_name, hash, results_folder)
  ensure_directory(dirname(cache_file)) # Ensure directory exists before saving
  tryCatch({
    saveRDS(curve, cache_file)
    log_message(paste("Saved cached curve to:", cache_file))
  }, error = function(e) {
    log_message(paste("ERROR saving cache file:", cache_file, "-", e$message))
  })
}

generate_hash <- function(pd, dimension, dataset_name, base_sigma, grid_points, tau_max) {
  # Make sure pd is consistent for hashing (e.g., convert to data.frame)
  pd_hashable <- tryCatch(as.data.frame(pd), error = function(e) NULL)
  if (is.null(pd_hashable)) pd_hashable <- "invalid_pd_structure"

  digest::digest(list(pd_hashable, dimension, dataset_name, base_sigma, grid_points, tau_max), algo = "sha256")
}


# --- Compute Betti Curve (Gaussian Kernel) ---
compute_betti_curve <- function(pd, dimension, grid_points, tau_max, base_sigma) {
  pd <- as.data.frame(pd)
  # Ensure correct column names or indices if using indices
  # Assuming columns are: dimension, birth, death
  valid_pd <- pd[pd[, 1] == dimension & pd[, 2] < pd[, 3] & is.finite(pd[, 2]) & is.finite(pd[, 3]),]
  if (nrow(valid_pd) == 0) return(rep(0, grid_points))

  tau_grid <- seq(min(pd[, 2], na.rm = TRUE), tau_max, length.out = grid_points) # Adjust grid start?

  gaussian <- function(x, mu, sigma) {
    # Ensure sigma is positive and non-zero
    sigma <- max(sigma, 1e-9)
    exp(-((x - mu) ^ 2) / (2 * sigma ^ 2)) / (sqrt(2 * pi) * sigma)
  }

  curve_vals <- sapply(tau_grid, function(tau) {
    sum(sapply(1:nrow(valid_pd), function(i) {
      b <- valid_pd[i, 2]
      d <- valid_pd[i, 3]
      persistence <- d - b
      # Ensure persistence is positive
      if (persistence <= 0) return(0)

      interval_sigma <- max(base_sigma, persistence / 2)

      # Calculate Gaussian contribution, ensure it's finite
      val <- persistence * gaussian(tau, (b + d) / 2, interval_sigma)
      ifelse(is.finite(val), val, 0)
    }))
  })

  # Ensure the final curve is finite
  curve_vals[!is.finite(curve_vals)] <- 0
  return(curve_vals)
}

# --- Get or Compute Curve (Handles Caching) ---
get_or_compute_curve <- function(pd, dimension, dataset_name, base_sigma, grid_points, tau_max, results_folder) {
  # Generate hash based on relevant parameters
  hash <- generate_hash(pd, dimension, dataset_name, base_sigma, grid_points, tau_max)

  # Try loading from cache
  cached_curve <- load_curve_cache(dataset_name, hash, results_folder)
  if (!is.null(cached_curve) && length(cached_curve) == grid_points && all(is.finite(cached_curve))) {
    return(cached_curve) # Return cached & normalized curve
  }

  log_message("Cache miss or invalid cache. Computing curve.")
  # Compute the raw curve
  curve <- compute_betti_curve(pd, dimension, grid_points, tau_max, base_sigma)

  # Normalize the computed curve
  norm_curve <- normalize_curve(curve, grid_points)

  # Save the *normalized* curve to cache
  if (all(is.finite(norm_curve))) {
    save_curve_cache(norm_curve, dataset_name, hash, results_folder)
  } else {
    log_message("Warning: Computed curve contains non-finite values after normalization. Not caching.")
    norm_curve <- rep(0, grid_points) # Fallback to zero curve
  }

  return(norm_curve)
}

# --- Bootstrapping for PD Curves ---
bootstrap_curve <- function(pd_subset, dimension, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, seed = 42, results_folder, num_cores = 8) {
  # (Keep the definition from your script)
  log_message(paste("Starting bootstrapping for Betti curve in dimension", dimension, "with", n_bootstrap, "samples."))
  individual_pd_curves <- lapply(seq_along(pd_subset), function(pd_idx) {
    pd <- pd_subset[[pd_idx]]
    if (is.matrix(pd) && nrow(pd) > 0) {
      curve <- get_or_compute_curve(
         pd = pd, dimension = dimension, dataset_name = dataset_name, base_sigma = base_sigma,
         grid_points = grid_points, tau_max = tau_max, results_folder = results_folder
       )
      # Normalization happens within get_or_compute_curve now
      return(curve)
    } else {
      return(rep(0, grid_points))
    }
  })
  names(individual_pd_curves) <- paste0("pd_", seq_along(pd_subset))

  bootstrap_curves <- parallel::mclapply(1:n_bootstrap, function(rep) {
    # Use parallel::mclapply
    set.seed(seed + rep)
    sampled_pds <- pd_subset[sample(seq_along(pd_subset), replace = TRUE)]
    aggregated_curve <- numeric(grid_points)
    for (pd in sampled_pds) {
      if (is.matrix(pd) && nrow(pd) > 0) {
        # Fetch normalized curve
        curve <- get_or_compute_curve(
           pd = pd, dimension = dimension, dataset_name = dataset_name, base_sigma = base_sigma,
           grid_points = grid_points, tau_max = tau_max, results_folder = results_folder
         )
        aggregated_curve <- aggregated_curve + curve
      }
    }
    # Normalize the aggregated curve from the bootstrap sample
    normalize_curve(aggregated_curve, grid_points)
  }, mc.cores = num_cores) # Use parallelization

  valid_curves <- bootstrap_curves[!sapply(bootstrap_curves, is.null)]
  if (length(valid_curves) == 0) {
    return(list(
       mean = rep(0, grid_points),
       lower = rep(0, grid_points),
       upper = rep(0, grid_points),
       individual_pd_curves = individual_pd_curves
     ))
  }

  # Combine results into a matrix for faster calculation
  bootstrap_matrix <- do.call(cbind, valid_curves)

  aggregated_results <- list(
     mean = rowMeans(bootstrap_matrix, na.rm = TRUE),
     lower = apply(bootstrap_matrix, 1, quantile, probs = 0.025, na.rm = TRUE),
     upper = apply(bootstrap_matrix, 1, quantile, probs = 0.975, na.rm = TRUE),
     individual_pd_curves = individual_pd_curves # Store the *normalized* individual curves
   )

  log_message(paste("Completed bootstrapping for Betti curve in dimension", dimension))
  return(aggregated_results)
}

# --- Euler Curve Computation ---
compute_euler_curve <- function(pd_subset, dimensions, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, results_folder) {
  # (Keep the definition from your script)
  log_message("Starting total Euler curve computation with bootstrapping.")
  # Get individual normalized curves for each dimension first
  individual_betti_curves_by_dim <- lapply(dimensions, function(d) {
    lapply(seq_along(pd_subset), function(pd_idx) {
      pd <- pd_subset[[pd_idx]]
      if (is.matrix(pd) && nrow(pd) > 0) {
        get_or_compute_curve(pd, d, dataset_name, base_sigma, grid_points, tau_max, results_folder)
      } else {
        rep(0, grid_points)
      }
    })
  })

  # Compute individual Euler curves by combining Betti curves with signs
  individual_euler_curves <- lapply(seq_along(pd_subset), function(pd_idx) {
    euler_curve <- rep(0, grid_points)
    for (dim_idx in seq_along(dimensions)) {
      euler_curve <- euler_curve + ((-1) ^ dimensions[dim_idx]) * individual_betti_curves_by_dim[[dim_idx]][[pd_idx]]
    }
    # Normalize the final Euler curve for this PD
    normalize_curve(euler_curve, grid_points)
  })
  names(individual_euler_curves) <- paste0("pd_", seq_along(pd_subset))

  # Bootstrap the already computed *individual normalized Euler curves*
  bootstrap_euler_curves <- parallel::mclapply(1:n_bootstrap, function(rep) {
    set.seed(42 + rep) # Ensure reproducibility
    sampled_indices <- sample(seq_along(individual_euler_curves), replace = TRUE)
    sampled_curves <- individual_euler_curves[sampled_indices]

    # Check if sampled_curves is not empty and contains valid vectors
    valid_sampled_curves <- Filter(function(x) is.numeric(x) && length(x) == grid_points, sampled_curves)
    if (length(valid_sampled_curves) == 0) return(rep(0, grid_points))

    # Aggregate (sum) the sampled normalized Euler curves
    aggregated_curve <- rowSums(do.call(cbind, valid_sampled_curves), na.rm = TRUE)

    # Normalize the aggregated bootstrap replicate curve
    normalize_curve(aggregated_curve, grid_points)
  }, mc.cores = detectCores() - 1) # Use parallelization

  valid_bootstrap_curves <- bootstrap_euler_curves[!sapply(bootstrap_euler_curves, is.null)]
  if (length(valid_bootstrap_curves) == 0) {
    return(list(mean = rep(0, grid_points), lower = rep(0, grid_points), upper = rep(0, grid_points),
                   individual_euler_curves = individual_euler_curves))
  }

  bootstrap_matrix <- do.call(cbind, valid_bootstrap_curves)

  log_message("Successfully computed the total Euler curve.")
  return(list(
     mean = rowMeans(bootstrap_matrix, na.rm = TRUE),
     lower = apply(bootstrap_matrix, 1, quantile, probs = 0.025, na.rm = TRUE),
     upper = apply(bootstrap_matrix, 1, quantile, probs = 0.975, na.rm = TRUE),
     individual_euler_curves = individual_euler_curves # Store normalized individual Euler curves
   ))
}

# --- Check Betti-Euler Consistency ---
# (Keep the definition from your script)
# No changes needed here, it uses the mean curves.

# --- Null Stats Computation ---
# (Keep the definition from your script, ensuring it uses the updated helpers)
compute_null_stats <- function(pd_list, metadata_list, tau_vals, dimensions,
                               curve_type = "Euler", dataset_name, base_sigma,
                               grid_points, tau_max, results_folder, verbose = TRUE,
                               num_cores = 8, landscape_list = NULL) {
  # Added landscape_list

  log_message <- function(msg) {
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
  }

  log_message(paste("Computing Null Stats for", curve_type))

  integrated_diff <- function(curve_diff, tau_vals) {
    # Local helper
    delta_tau <- diff(tau_vals)
    integrated_sq <- sum(delta_tau * (curve_diff[-length(curve_diff)] ^ 2 + curve_diff[-1] ^ 2) / 2, na.rm = TRUE)
    sqrt(integrated_sq)
  }

  combined_meta <- do.call(rbind, metadata_list)
  # Use ignore.case=TRUE for robustness
  bs_cols <- grep("^random_group_bootstrap_", colnames(combined_meta), value = TRUE, ignore.case = TRUE)
  log_message(paste("Found", length(bs_cols), "random bootstrap columns."))
  if (length(bs_cols) == 0) {
    log_message("No random bootstrap columns found in combined metadata. Cannot compute null stats.")
    return(NULL)
  }

  # Parallel processing over bootstrap columns
  results <- parallel::mclapply(bs_cols, function(col) {
    group1_ids <- combined_meta$orig.ident.Iter[which(combined_meta[[col]] == 1)]
    group2_ids <- combined_meta$orig.ident.Iter[which(combined_meta[[col]] == 2)]

    if (length(group1_ids) == 0 || length(group2_ids) == 0) return(NULL)

    pd_group1 <- pd_list[names(pd_list) %in% group1_ids]
    pd_group2 <- pd_list[names(pd_list) %in% group2_ids]

    # Betti/Euler Calculation
    curve1 <- compute_euler_curve(pd_group1, dimensions, "group1", dataset_name, grid_points, tau_max, base_sigma, n_bootstrap = 50, results_folder = results_folder)
    curve2 <- compute_euler_curve(pd_group2, dimensions, "group2", dataset_name, grid_points, tau_max, base_sigma, n_bootstrap = 50, results_folder = results_folder)

    diff_val_euler <- integrated_diff(curve1$mean - curve2$mean, tau_vals)
    denom_euler <- sd(c(curve1$mean, curve2$mean), na.rm = TRUE)
    eff_size_euler <- if (denom_euler == 0) 0 else diff_val_euler / denom_euler

    # KS calculation using individual curves
    inds1_euler <- do.call(cbind, curve1$individual_euler_curves)
    inds2_euler <- do.call(cbind, curve2$individual_euler_curves)
    ks_stats_euler <- NA # Default if calculation fails
    if (is.matrix(inds1_euler) && is.matrix(inds2_euler) && ncol(inds1_euler) > 0 && ncol(inds2_euler)) {
      ks_stats_list <- lapply(1:nrow(inds1_euler), function(i) {
        res <- tryCatch(ks.test(inds1_euler[i,], inds2_euler[i,]), error = function(e) NULL)
        if (is.null(res)) NA else as.numeric(res$statistic)
      })
      ks_stats_euler <- max(unlist(ks_stats_list), na.rm = TRUE)
      if (!is.finite(ks_stats_euler)) ks_stats_euler <- NA
    }

    # Landscape Calculation (if landscape_list is provided)
    eff_size_landscape <- NA
    ks_stat_landscape <- NA
    if (!is.null(landscape_list)) {
      landscape_group1 <- landscape_list[names(landscape_list) %in% group1_ids]
      landscape_group2 <- landscape_list[names(landscape_list) %in% group2_ids]

      if (length(landscape_group1) > 0 && length(landscape_group2) > 0) {
        grid_land <- seq(0, 1, length.out = grid_points) # Grid for landscape
        curve1_land <- compute_aggregated_landscape_curve(landscape_group1, grid = grid_land)
        curve2_land <- compute_aggregated_landscape_curve(landscape_group2, grid = grid_land)

        diff_val_land <- integrated_diff(curve1_land - curve2_land, grid_land)
        denom_land <- sd(c(curve1_land, curve2_land), na.rm = TRUE)
        eff_size_landscape <- if (denom_land == 0) 0 else diff_val_land / denom_land

        # KS for landscape requires individual landscape curves (if available & computed)
        # This part assumes `bootstrap_aggregated_landscape_curve` was run and stored individual curves
        # Let's compute KS based on the aggregated curves for now as a proxy
        if (sum(curve1_land) > 0 && sum(curve2_land) > 0) {
          cdf1 <- cumsum(curve1_land / sum(curve1_land))
          cdf2 <- cumsum(curve2_land / sum(curve2_land))
          ks_stat_landscape <- max(abs(cdf1 - cdf2), na.rm = TRUE)
        }
      }
    }

    list(
           euler_eff_size = eff_size_euler,
           euler_ks_stat = ks_stats_euler,
           landscape_eff_size = eff_size_landscape,
           landscape_ks_stat = ks_stat_landscape
       )
  }, mc.cores = num_cores)

  results <- results[!sapply(results, is.null)] # Filter NULLs

  # Aggregate stats
  null_effects_euler <- sapply(results, `[[`, "euler_eff_size")
  null_ks_euler <- sapply(results, `[[`, "euler_ks_stat")
  null_effects_landscape <- sapply(results, `[[`, "landscape_eff_size")
  null_ks_landscape <- sapply(results, `[[`, "landscape_ks_stat")

  summarize_stats <- function(vals) {
    if (length(vals) > 0 && any(is.finite(vals))) {
      list(mean = mean(vals, na.rm = T), sd = sd(vals, na.rm = T), quantiles = quantile(vals, probs = c(0.025, 0.5, 0.975), na.rm = T))
    } else { NULL }
  }

  list(
       euler_effect_size = summarize_stats(null_effects_euler),
       euler_ks = summarize_stats(null_ks_euler),
       landscape_effect_size = summarize_stats(null_effects_landscape),
       landscape_ks = summarize_stats(null_ks_landscape)
   )
}


# --- Compute Null Stats for Individual Landscapes via Bootstrapping ---
# (Keep definition from your script)
# No changes needed.

# --- Compute Single Aggregated Landscape Curve ---
# (Keep definition from your script)
# No changes needed.

# --- Compute Aggregated Landscape Curve for a Group ---
# (Keep definition from your script)
# No changes needed.

# --- Compute Null Stats for Aggregated Landscape Curves ---
# (Keep definition from your script)
# No changes needed.

# --- Bootstrap Aggregated Landscape Curve ---
# (Keep definition from your script)
# No changes needed.

# --- Analyze and Plot Curves (Pairwise) ---
# (Keep definition from your script)
# Make sure it uses the passed null_effect_stats and null_ks_stats correctly
# Ensure plot saving uses ensure_directory()

# --- Perform Pairwise Comparisons (Main Logic) ---
# (Keep definition from your script)
# Make sure it calls the corrected analyze_and_plot_curves

# --- Plot Aggregated Curves ---
# (Keep definition from your script)
# Make sure it uses ensure_directory()

# --- Main Function for Cross-Iteration Betti Comparison ---
# (Keep the definition of compute_and_compare_betti_iterations)
# Ensure it calls the corrected helper functions defined above.
# Make sure it handles the landscape_list input and calls landscape comparisons.

# ---------------------------
# Cluster Metric Comparison Helpers & Main Function
# ---------------------------

# --- Metric Functions (Keep definitions) ---
# variation_of_information(), jaccard_index(), purity_single(), 
# safe_normalize(), safe_metric(), metric_ari(), metric_nmi(), etc.
# (Keep these definitions from your script)

# --- Helper: Global Linked Cluster Mapping ---
# (Keep definition from your script)

# --- Helper: Process Linked Clusters ---
# (Keep definition from your script)

# --- Main Function for Cross-Iteration Cluster Metric Comparison ---
# (Keep definition of compare_all_iterations_with_mapped_methods_random_baseline_including_purity_pvals)
# Ensure it calls the correct helper functions defined above.
# Ensure plot saving uses ensure_directory().

# ---------------------------
# Example Usage Section (Optional - keep for testing/documentation)
# ---------------------------

# # --- Example Setup ---
# # This would typically be done in a main calling script
# results_folder <- "results_march25_cross_iter_test" 
# config_paths <- setup_directories(results_folder) 
# plots_folder <- config_paths$plots_folder 
# plot_folders <- config_paths$plot_folders
# num_cores <- 4 # Example core count
# 
# # --- Load Final Seurat Objects (replace with actual paths) ---
# # These should be the outputs saved by 03_run_analysis.R
# data_iterations_final <- list(
#   list(
#     name = "Raw",
#     seurat_obj = readRDS(file.path(results_folder, "../seurat_objects/raw_seurat_object_final.rds")),
#     pd_list = file.path(results_folder, "../preprocessed/raw_pd_list_aligned.rds"), # Path to aligned list
#     landscape_list = file.path(results_folder, "../preprocessed/raw_landscape_list_aligned.rds") # Path to aligned list
#   ),
#   list(
#     name = "SCT_Individual",
#     seurat_obj = readRDS(file.path(results_folder, "../seurat_objects/sct_individual_seurat_object_final.rds")),
#     pd_list = file.path(results_folder, "../preprocessed/sctInd_pd_list_aligned.rds"),
#     landscape_list = file.path(results_folder, "../preprocessed/sctInd_landscape_list_aligned.rds")
#   ),
#   list(
#     name = "Integrated",
#     seurat_obj = readRDS(file.path(results_folder, "../seurat_objects/integrated_seurat_object_final.rds")),
#     pd_list = file.path(results_folder, "../preprocessed/integrated_pd_list_aligned.rds"),
#     landscape_list = file.path(results_folder, "../preprocessed/integrated_landscape_list_aligned.rds")
#   )
#   # ... add other iterations ...
# )
# 
# # --- Example Call for Betti Curve Comparison ---
# # Load PD lists and metadata lists required by the function
# all_pd_list <- list()
# all_metadata_list <- list()
# all_landscape_list <- list() # Initialize
# 
# for (iter in data_iterations_final) {
#     pd <- readRDS(iter$pd_list)
#     meta <- iter$seurat_obj@meta.data
#     meta$orig.ident.Iter <- rownames(meta) # Assuming rownames are unique identifiers like "Raw_Sample1"
#     all_pd_list <- c(all_pd_list, pd)
#     all_metadata_list <- c(all_metadata_list, list(meta))
#     
#     # Load landscape list if path exists
#     if (!is.null(iter$landscape_list) && file.exists(iter$landscape_list)) {
#         land <- readRDS(iter$landscape_list)
#         # Ensure names match pd_list for consistency if needed, might need adjustment
#         names(land) <- names(pd) 
#         all_landscape_list <- c(all_landscape_list, land)
#     } else {
#         log_message(paste("Landscape list not found or specified for", iter$name))
#         # Handle cases where landscape list might be missing for some iterations
#     }
# }
# # Ensure landscape list has same names as pd list if used together
# if(length(all_landscape_list) > 0 && length(all_landscape_list) == length(all_pd_list)) {
#    names(all_landscape_list) <- names(all_pd_list)
# } else if (length(all_landscape_list) > 0) {
#    log_message("Warning: Landscape list length mismatch with PD list. Landscape comparisons might be affected.")
#    # Set to NULL if inconsistent to avoid errors in downstream functions expecting alignment
#    all_landscape_list <- NULL 
# }
# 
# 
# cross_iteration_comparison_with_betti(
#   pd_list = all_pd_list,
#   metadata_list = all_metadata_list,
#   landscape_list = if(length(all_landscape_list) > 0) all_landscape_list else NULL, # Pass combined list
#   group_by_col = "Tissue", # Example grouping
#   dataset_name = "VDT_Cross_Iteration_Tissue", # Descriptive name
#   output_folder = file.path(results_folder, "cross_iteration_comparisons"),
#   num_cores = num_cores,
#   # Keep other parameters like dimensions, grid_points etc.
#   grid_points = 500, 
#   num_permutations = 1000, # Reduce for testing if needed
#   bootstrap_samples = 100, # Reduce for testing
#   base_sigma = 1,
#   verbose = TRUE
# )
# 
# # --- Example Setup & Call for Cluster Metric Comparison ---
# # Define method mappings based on columns present in the *final* Seurat objects
# method_mapping_list_example <- list(
#   "seurat_cluster" = setNames(paste0("seurat_cluster_", tolower(sapply(data_iterations_final, `[[`, "name"))), sapply(data_iterations_final, `[[`, "name")),
#   "kmeans_cluster_sra" = setNames(paste0("kmeans_cluster_", tolower(sapply(data_iterations_final, `[[`, "name")), "_sra"), sapply(data_iterations_final, `[[`, "name"))
#   # ... add other unified methods based on the actual columns present ...
# )
# 
# results_cluster_metrics <- compare_all_iterations_with_mapped_methods_random_baseline_including_purity_pvals (
#   data_iterations = data_iterations_final, # Use the list containing loaded Seurat objects
#   method_mapping_list = method_mapping_list_example,
#   results_folder    = file.path(results_folder, "cross_iteration_comparisons"),
#   dataset_name = "VDT_Cross_Iteration_Metrics",
#   n_cores = num_cores
# )

log_message("Cross-iteration analysis script loaded.")