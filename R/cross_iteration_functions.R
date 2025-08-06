# ---------------------------
# Cross-Iteration Betti Curve and Cluster Metric Comparisons
# ---------------------------
# Add this to the top of your script with other library calls

#' Log a timestamped message
#'
#' Prepends a timestamp to a message and prints it to the console. Useful for
#' lightweight progress reporting throughout long-running computations.
#'
#' @param msg Character string to log.
#' @return Invisibly returns `NULL` after printing the message.
#' @examples
#' log_message("Processing started")
#' @keywords internal
log_message <- function(msg) {
  message(sprintf("[%s] %s", Sys.time(), msg))
}

# --- Helper: Directory Creation ---
#' Ensure that a directory exists
#'
#' Creates the directory path recursively if it does not already exist.
#'
#' @param path Character path to check or create.
#' @return Invisibly returns `NULL` after ensuring the directory exists.
#' @examples
#' ensure_directory(tempdir())
#' @keywords internal
ensure_directory <- function(path) {
  if (length(path) == 0 || is.na(path) || path == "") {
    stop("Invalid directory path provided.")
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# --- Helper: Normalization with Division-by-Zero Check ---
#' Normalize a curve to unit sum
#'
#' Scales a numeric vector so that it sums to one. If the total is zero, a
#' vector of zeros is returned to avoid division by zero.
#'
#' @param curve Numeric vector representing the curve values.
#' @param grid_points Integer indicating the length of the curve.
#' @return Numeric vector of length `grid_points` that sums to one.
#' @examples
#' normalize_curve(c(1, 2, 3), grid_points = 3)
#' @keywords internal
normalize_curve <- function(curve, grid_points) {
  tot <- sum(curve)
  if (tot == 0) rep(0, grid_points) else curve / tot
}

# --- Curve Caching ---
#' Generate a cache file name for a Betti curve
#'
#' Constructs the file path used to store a cached Betti curve and ensures that
#' the directory exists.
#'
#' @param dataset_name Name of the dataset as a string.
#' @param hash Hash string representing the curve parameters.
#' @param results_folder Base directory where results are stored.
#' @return Character string containing the full cache file path.
#' @examples
#' generate_cache_filename("dataset", "abc123", tempdir())
#' @keywords internal
generate_cache_filename <- function(dataset_name, hash, results_folder) {
  path <- file.path(results_folder, "plots", "betti_plots", "betti_cache", dataset_name, paste0("cache_", hash, ".rds"))
  ensure_directory(dirname(path))
  path
}

# --- Compute Betti Curve ---
#' Compute a smoothed Betti curve
#'
#' Converts a persistence diagram into a smoothed Betti curve using Gaussian
#' kernels whose width depends on the persistence of each feature.
#'
#' @param pd Persistence diagram as a matrix or data frame with columns
#'   `dimension`, `birth`, and `death`.
#' @param dimension Integer dimension for which to compute the Betti curve.
#' @param grid_points Number of evaluation points along the filtration scale.
#' @param tau_max Maximum filtration value.
#' @param base_sigma Minimum standard deviation for the Gaussian kernel.
#'
#' @return Numeric vector of length `grid_points` representing the Betti curve.
#' @examples
#' pd <- matrix(c(0, 0.1, 0.4), ncol = 3)
#' compute_betti_curve(pd, dimension = 0, grid_points = 10, tau_max = 1, base_sigma = 0.05)
#' @keywords internal
compute_betti_curve <- function(pd, dimension, grid_points, tau_max, base_sigma) {
  pd <- as.data.frame(pd)
  valid_pd <- pd[pd[, 1] == dimension & pd[, 2] < pd[, 3],]
  if (nrow(valid_pd) == 0) return(rep(0, grid_points))
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  gaussian <- function(x, mu, sigma) {
    exp(-((x - mu) ^ 2) / (2 * sigma ^ 2)) / (sqrt(2 * pi) * sigma)
  }
  sapply(tau_grid, function(tau) {
    sum(sapply(1:nrow(valid_pd), function(i) {
      b <- valid_pd[i, 2]
      d <- valid_pd[i, 3]
      persistence <- d - b
      interval_sigma <- max(base_sigma, persistence / 2)
      persistence * gaussian(tau, (b + d) / 2, interval_sigma)
    }))
  })
}

# --- Check Betti-Euler Consistency ---
#' Verify consistency between Betti and Euler curves
#'
#' Compares an Euler curve with the alternating sum of Betti curves across
#' dimensions to check whether they satisfy the Euler characteristic identity.
#'
#' @param euler_curve List containing at least a `mean` component representing
#'   the Euler curve.
#' @param betti_curves List of lists where each element has a `mean` component
#'   for the corresponding Betti curve.
#' @param dimensions Integer vector of homology dimensions represented in the
#'   Betti curves.
#' @param tau_vals Numeric vector of filtration parameter values.
#'
#' @return Logical vector indicating consistency at each `tau` value, or `NULL`
#'   if an error occurs.
#' @examples
#' euler <- list(mean = c(1, 0))
#' betti <- list(list(mean = c(1, 0)))
#' check_betti_euler_consistency(euler, betti, dimensions = 0, tau_vals = c(0, 1))
#' @keywords internal
check_betti_euler_consistency <- function(euler_curve, betti_curves, dimensions, tau_vals) {
  log_message("Starting Betti-Euler consistency check.")
  if (!is.list(betti_curves) ||
      any(sapply(betti_curves, function(x)!is.list(x) || !"mean" %in% names(x)))) {
    stop("Invalid Betti curves structure. Ensure each dimension has a 'mean' component.")
  }
  tryCatch({
    consistency_results <- lapply(seq_along(tau_vals), function(tau_idx) {
      tau <- tau_vals[tau_idx]
      euler_at_tau <- euler_curve$mean[tau_idx]
      betti_sum_at_tau <- sum(sapply(seq_along(dimensions), function(dim_idx) {
        if (!is.null(betti_curves[[dim_idx]]$mean))
          betti_curves[[dim_idx]]$mean[tau_idx] * (-1) ^ dimensions[dim_idx]
        else 0
      }))
      abs(euler_at_tau - betti_sum_at_tau) < 1e-6
    })
    log_message("Completed Betti-Euler consistency check.")
    consistency_results
  }, error = function(e) {
    log_message(paste("Error during Betti-Euler consistency check:", e$message))
    return(NULL)
  })
}

compute_null_stats <- function(pd_list, metadata_list, tau_vals, dimensions,
                               curve_type = "Euler", dataset_name, base_sigma,
                               grid_points, tau_max, results_folder, verbose = TRUE,
                               num_cores = 8) {
  integrated_diff <- function(curve_diff, tau_vals) {
    delta_tau <- diff(tau_vals)
    integrated_sq <- sum(delta_tau * (curve_diff[-length(curve_diff)] ^ 2 +
                                        curve_diff[-1] ^ 2) / 2, na.rm = TRUE)
    sqrt(integrated_sq)
  }
  log_message <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", Sys.time(), msg))
  }
  combined_meta <- do.call(rbind, metadata_list)
  bs_cols <- grep("^random_group_bootstrap_", colnames(combined_meta),
                  value = TRUE, ignore.case = TRUE)
  results <- parallel::mclapply(bs_cols, function(col) {
    idx1 <- which(combined_meta[[col]] == 1)
    idx2 <- which(combined_meta[[col]] == 2)
    if (length(idx1) == 0 || length(idx2) == 0) return(NULL)
    pd_grp1 <- pd_list[idx1]
    pd_grp2 <- pd_list[idx2]
    if (curve_type == "Euler") {
      curve1 <- compute_euler_curve(pd_grp1, dimensions, "group1", dataset_name,
                                    grid_points, tau_max, base_sigma, n_bootstrap = 50,
                                    results_folder = results_folder)
      curve2 <- compute_euler_curve(pd_grp2, dimensions, "group2", dataset_name,
                                    grid_points, tau_max, base_sigma, n_bootstrap = 50,
                                    results_folder = results_folder)
      diff_val <- integrated_diff(curve1$mean - curve2$mean, tau_vals)
      denom <- sd(c(curve1$mean, curve2$mean), na.rm = TRUE)
      eff_size <- if (denom == 0) 0 else diff_val / denom
      ind1 <- do.call(cbind, curve1$individual_euler_curves)
      ind2 <- do.call(cbind, curve2$individual_euler_curves)
      ks_stats <- sapply(1:nrow(ind1), function(i) {
        res <- tryCatch(ks.test(ind1[i,], ind2[i,]), error = function(e) NULL)
        if (is.null(res)) NA else as.numeric(res$statistic)
      })
      return(list(eff_size = eff_size, ks_stat = max(ks_stats, na.rm = TRUE)))
    } else {
      return(NULL)
    }
  }, mc.cores = num_cores)
  results <- results[!sapply(results, is.null)]
  null_effects <- sapply(results, function(x) x$eff_size)
  null_ks <- sapply(results, function(x) x$ks_stat)
  list(
    null_effect = if (length(null_effects) > 0) list(
      mean = mean(null_effects, na.rm = TRUE),
      sd = sd(null_effects, na.rm = TRUE),
      quantiles = quantile(null_effects, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    ) else NULL,
    null_ks = if (length(null_ks) > 0) list(
      mean = mean(null_ks, na.rm = TRUE),
      sd = sd(null_ks, na.rm = TRUE),
      quantiles = quantile(null_ks, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    ) else NULL
  )
}

#' Generate a hash for curve parameters
#'
#' Creates a SHA-256 hash based on the persistence diagram and curve parameters.
#'
#' @param pd Persistence diagram as a matrix or data frame.
#' @param dimension Integer dimension for which the curve is computed.
#' @param dataset_name Name of the dataset as a string.
#' @param base_sigma Minimum standard deviation for the Gaussian kernel.
#' @param grid_points Number of evaluation points along the filtration scale.
#' @param tau_max Maximum filtration value.
#'
#' @return Character string representing the hash.
#' @examples
#' generate_hash(matrix(c(0, 0.1, 0.4), ncol = 3), 0, "dataset", 0.05, 10, 1)
#' @keywords internal
generate_hash <- function(pd, dimension, dataset_name, base_sigma, grid_points, tau_max) {
  digest(list(pd, dimension, dataset_name, base_sigma, grid_points, tau_max), algo = "sha256")
}

#' Load a cached curve from disk
#'
#' Retrieves a previously cached curve using the dataset name and parameter hash
#' to locate the appropriate file.
#'
#' @param dataset_name Name of the dataset as a string.
#' @param hash Hash string identifying the curve parameters.
#' @param results_folder Base directory where results are stored.
#'
#' @return The cached curve if available, otherwise `NULL`.
#' @examples
#' \dontrun{
#' load_curve_cache("dataset", "abc123", tempdir())
#' }
#' @keywords internal
load_curve_cache <- function(dataset_name, hash, results_folder) {
  cache_file <- generate_cache_filename(dataset_name, hash, results_folder)
  if (file.exists(cache_file)) {
    log_message(paste("Loaded cached curve from:", cache_file))
    return(readRDS(cache_file))
  }
  NULL
}

#' Save a curve to the cache
#'
#' Stores a computed curve to disk so that subsequent requests with the same
#' parameters can reuse the cached result.
#'
#' @param curve Numeric vector representing the curve to cache.
#' @param dataset_name Name of the dataset as a string.
#' @param hash Hash string identifying the curve parameters.
#' @param results_folder Base directory where results are stored.
#' @return Invisibly returns `NULL` after saving the file.
#' @examples
#' \dontrun{
#' save_curve_cache(1:5, "dataset", "abc123", tempdir())
#' }
#' @keywords internal
save_curve_cache <- function(curve, dataset_name, hash, results_folder) {
  cache_file <- generate_cache_filename(dataset_name, hash, results_folder)
  ensure_directory(dirname(cache_file))
  saveRDS(curve, cache_file)
  log_message(paste("Saved cached curve to:", cache_file))
}

#' Retrieve or compute a Betti curve
#'
#' Attempts to load a cached Betti curve; if not available, computes the curve
#' and stores it for future use.
#'
#' @param pd Persistence diagram as a matrix or data frame.
#' @param dimension Integer dimension for which to compute the curve.
#' @param dataset_name Name of the dataset as a string.
#' @param base_sigma Minimum standard deviation for the Gaussian kernel.
#' @param grid_points Number of evaluation points along the filtration scale.
#' @param tau_max Maximum filtration value.
#' @param results_folder Base directory where results are stored.
#' @return Normalized Betti curve as a numeric vector.
#' @examples
#' \dontrun{
#' pd <- matrix(c(0, 0.1, 0.4), ncol = 3)
#' get_or_compute_curve(pd, 0, "dataset", 0.05, 10, 1, tempdir())
#' }
#' @keywords internal
get_or_compute_curve <- function(pd, dimension, dataset_name, base_sigma, grid_points, tau_max, results_folder) {
  hash <- generate_hash(pd, dimension, dataset_name, base_sigma, grid_points, tau_max)
  cached_curve <- load_curve_cache(dataset_name, hash, results_folder)
  if (!is.null(cached_curve)) return(cached_curve)
  log_message("Cache miss for curve. Computing curve.")
  curve <- compute_betti_curve(pd, dimension, grid_points, tau_max, base_sigma)
  norm_curve <- normalize_curve(curve, grid_points)
  save_curve_cache(norm_curve, dataset_name, hash, results_folder)
  norm_curve
}

#' Bootstrap Betti curves for a subset of persistence diagrams
#'
#' Generates bootstrap samples of Betti curves by resampling persistence
#' diagrams with replacement.
#'
#' @param pd_subset List of persistence diagrams.
#' @param dimension Integer dimension for which to compute Betti curves.
#' @param group_name Character label for the group being processed.
#' @param dataset_name Name of the dataset as a string.
#' @param grid_points Number of evaluation points along the filtration scale.
#' @param tau_max Maximum filtration value.
#' @param base_sigma Minimum standard deviation for the Gaussian kernel.
#' @param n_bootstrap Number of bootstrap samples to draw.
#' @param seed Random seed for reproducibility.
#' @param results_folder Base directory where results are stored.
#' @param num_cores Number of cores for parallel computation.
#'
#' @return List containing the mean curve, lower and upper confidence bands, and
#'   individual persistence diagram curves.
#' @examples
#' \dontrun{
#' pd <- list(matrix(c(0, 0.1, 0.4), ncol = 3))
#' bootstrap_curve(pd, 0, "grp", "dataset", 10, 1, 0.05, 5, results_folder = tempdir())
#' }
#' @keywords internal
bootstrap_curve <- function(pd_subset, dimension, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, seed = 42, results_folder, num_cores = 8) {
  log_message(paste("Starting bootstrapping for dimension", dimension, "with", n_bootstrap, "samples."))
  individual_pd_curves <- parallel::mclapply(seq_along(pd_subset), function(i) {
    pd <- pd_subset[[i]]
    if (is.matrix(pd) && nrow(pd) > 0) {
      curv <- get_or_compute_curve(pd, dimension, dataset_name, base_sigma, grid_points, tau_max, results_folder)
      normalize_curve(curv, grid_points)
    } else {
      rep(0, grid_points)
    }
  }, mc.cores = num_cores)
  names(individual_pd_curves) <- paste0("pd_", seq_along(pd_subset))
  bootstrap_curves <- parallel::mclapply(1:n_bootstrap, function(rep) {
    set.seed(seed + rep)
    sampled <- pd_subset[sample(seq_along(pd_subset), replace = TRUE)]
    agg <- numeric(grid_points)
    for (pd in sampled) {
      if (is.matrix(pd) && nrow(pd) > 0) {
        agg <- agg + get_or_compute_curve(pd, dimension, dataset_name, base_sigma, grid_points, tau_max, results_folder)
      }
    }
    normalize_curve(agg, grid_points)
  }, mc.cores = num_cores)
  valid <- bootstrap_curves[!sapply(bootstrap_curves, is.null)]
  if (length(valid) == 0) {
    return(list(mean = rep(0, grid_points), lower = rep(0, grid_points), upper = rep(0, grid_points),
                individual_pd_curves = individual_pd_curves))
  }
  list(
    mean = rowMeans(do.call(cbind, valid), na.rm = TRUE),
    lower = apply(do.call(cbind, valid), 1, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(do.call(cbind, valid), 1, quantile, probs = 0.975, na.rm = TRUE),
    individual_pd_curves = individual_pd_curves
  )
}

#' Compute the Euler curve from persistence diagrams
#'
#' Aggregates Betti curves across dimensions with alternating signs to obtain
#' the Euler characteristic curve for a group of persistence diagrams.
#'
#' @param pd_subset List of persistence diagrams.
#' @param dimensions Integer vector of homology dimensions to include.
#' @param group_name Character label for the group being processed.
#' @param dataset_name Name of the dataset as a string.
#' @param grid_points Number of evaluation points along the filtration scale.
#' @param tau_max Maximum filtration value.
#' @param base_sigma Minimum standard deviation for the Gaussian kernel.
#' @param n_bootstrap Number of bootstrap samples to draw.
#' @param results_folder Base directory where results are stored.
#'
#' @return List containing the mean Euler curve, lower and upper confidence
#'   bands, and individual Euler curves.
#' @examples
#' \dontrun{
#' pd <- list(matrix(c(0, 0.1, 0.4), ncol = 3))
#' compute_euler_curve(pd, 0:1, "grp", "dataset", 10, 1, 0.05, 5, tempdir())
#' }
#' @keywords internal
compute_euler_curve <- function(pd_subset, dimensions, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, results_folder) {
  log_message("Computing Euler curve with bootstrapping.")
  individual_euler <- lapply(seq_along(pd_subset), function(i) {
    pd <- pd_subset[[i]]
    if (is.matrix(pd) && nrow(pd) > 0) {
      euler <- rep(0, grid_points)
      for (d in dimensions) {
        curv <- get_or_compute_curve(pd, d, dataset_name, base_sigma, grid_points, tau_max, results_folder)
        euler <- euler + ((-1) ^ d) * curv
      }
      normalize_curve(euler, grid_points)
    } else {
      rep(0, grid_points)
    }
  })
  names(individual_euler) <- paste0("pd_", seq_along(pd_subset))
  bootstrapped <- lapply(dimensions, function(d) {
    bootstrap_curve(pd_subset, d, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, seed = 42, results_folder = results_folder)
  })
  euler_mean <- rep(0, grid_points)
  euler_lower <- rep(0, grid_points)
  euler_upper <- rep(0, grid_points)
  for (i in seq_along(dimensions)) {
    euler_mean <- euler_mean + ((-1) ^ dimensions[i]) * bootstrapped[[i]]$mean
    euler_lower <- euler_lower + ((-1) ^ dimensions[i]) * bootstrapped[[i]]$lower
    euler_upper <- euler_upper + ((-1) ^ dimensions[i]) * bootstrapped[[i]]$upper
  }
  list(mean = euler_mean, lower = euler_lower, upper = euler_upper, individual_euler_curves = individual_euler)
}

#' Compute integrated difference between two curves
#'
#' Approximates the L2 distance between consecutive points of a difference curve
#' using the trapezoidal rule.
#'
#' @param curve_diff Numeric vector of differences between two curves.
#' @param tau_vals Numeric vector of corresponding filtration values.
#'
#' @return Single numeric value giving the integrated difference.
#' @examples
#' integrated_diff(c(0, 1, 0), c(0, 0.5, 1))
#' @keywords internal
integrated_diff <- function(curve_diff, tau_vals) {
  if (length(curve_diff) < 2) return(0)
  delta_tau <- diff(tau_vals)
  sqrt(sum(delta_tau * (curve_diff[-length(curve_diff)] ^ 2 + curve_diff[-1] ^ 2) / 2, na.rm = TRUE))
}

#' Summarize a persistence landscape
#'
#' Computes the area under the curve of an aggregated persistence landscape.
#'
#' @param landscape_obj Persistence landscape object with components `dim0` and
#'   `dim1`.
#' @param grid Numeric vector of evaluation points.
#'
#' @return Numeric value representing the area under the landscape curve.
#' @examples
#' \dontrun{
#' land <- list(dim0 = matrix(1, 10, 1), dim1 = matrix(0, 10, 1))
#' summary_landscape(land, grid = seq(0, 1, length.out = 10))
#' }
#' @keywords internal
summary_landscape <- function(landscape_obj, grid = seq(0, 1, length.out = grid_points)) {
  curve <- compute_landscape_curve(landscape_obj, grid = grid)
  delta <- diff(grid)
  auc <- sum(((curve[-length(curve)] + curve[-1]) / 2) * delta)
  auc
}

bootstrap_null_stats_landscape <- function(landscape_list, n_bootstrap = 50, grid_points = 500) {
  if (!exists("summary_landscape")) {
    summary_landscape <- function(landscape_obj, grid) {
      agg <- compute_landscape_curve(landscape_obj, grid = grid)
      delta <- diff(grid)
      sum(((agg[-length(agg)] + agg[-1]) / 2) * delta)
    }
  }
  grid <- seq(0, 1, length.out = grid_points)
  summary_values <- sapply(landscape_list, function(land) {
    summary_landscape(land, grid = grid)
  })
  boot_effects <- numeric(n_bootstrap)
  boot_ks <- numeric(n_bootstrap)
  N <- length(summary_values)
  for (rep in 1:n_bootstrap) {
    set.seed(100 + rep)
    indices <- sample(N, replace = TRUE)
    half <- floor(length(indices) / 2)
    group1_vals <- summary_values[indices[1:half]]
    group2_vals <- summary_values[indices[(half + 1):length(indices)]]
    diff_val <- abs(mean(group1_vals) - mean(group2_vals))
    sd_pooled <- sd(c(group1_vals, group2_vals), na.rm = TRUE)
    boot_effects[rep] <- ifelse(sd_pooled == 0, 0, diff_val / sd_pooled)
    if (sum(group1_vals) == 0 || sum(group2_vals) == 0) {
      boot_ks[rep] <- NA
    } else {
      cdf1 <- cumsum(group1_vals / sum(group1_vals))
      cdf2 <- cumsum(group2_vals / sum(group2_vals))
      boot_ks[rep] <- max(abs(cdf1 - cdf2), na.rm = TRUE)
    }
  }
  list(
    null_effect = list(
      mean = mean(boot_effects, na.rm = TRUE),
      sd = sd(boot_effects, na.rm = TRUE),
      quantiles = quantile(boot_effects, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    ),
    null_ks = list(
      mean = mean(boot_ks, na.rm = TRUE),
      sd = sd(boot_ks, na.rm = TRUE),
      quantiles = quantile(boot_ks, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    )
  )
}

compute_null_stats_individual_landscapes <- function(landscape_list, n_bootstrap = 50, grid_points = 500) {
  N <- length(landscape_list)
  if (N < 2) {
    return(NULL)
  }
  summary_values <- sapply(landscape_list, function(land) {
    summary_landscape(land, grid = seq(0, 1, length.out = grid_points))
  })
  boot_effects <- numeric(n_bootstrap)
  boot_ks <- numeric(n_bootstrap)
  for (rep in 1:n_bootstrap) {
    set.seed(100 + rep)
    indices <- sample(N, replace = TRUE)
    half <- floor(length(indices) / 2)
    group1_vals <- summary_values[indices[1:half]]
    group2_vals <- summary_values[indices[(half + 1):length(indices)]]
    diff_val <- abs(mean(group1_vals) - mean(group2_vals))
    sd_pooled <- sd(c(group1_vals, group2_vals), na.rm = TRUE)
    effect <- ifelse(sd_pooled == 0, 0, diff_val / sd_pooled)
    boot_effects[rep] <- effect
    ks_result <- tryCatch(ks.test(group1_vals, group2_vals), error = function(e) NULL)
    ks_stat <- if (!is.null(ks_result)) as.numeric(ks_result$statistic) else NA
    boot_ks[rep] <- ks_stat
  }
  list(
    null_effect = list(
      mean = mean(boot_effects, na.rm = TRUE),
      sd = sd(boot_effects, na.rm = TRUE),
      lower = quantile(boot_effects, probs = 0.025, na.rm = TRUE),
      median = quantile(boot_effects, probs = 0.5, na.rm = TRUE),
      upper = quantile(boot_effects, probs = 0.975, na.rm = TRUE)
    ),
    null_ks = list(
      mean = mean(boot_ks, na.rm = TRUE),
      sd = sd(boot_ks, na.rm = TRUE),
      lower = quantile(boot_ks, probs = 0.025, na.rm = TRUE),
      median = quantile(boot_ks, probs = 0.5, na.rm = TRUE),
      upper = quantile(boot_ks, probs = 0.975, na.rm = TRUE)
    )
  )
}

#' Compute an aggregated persistence landscape curve
#'
#' Combines dimension 0 and dimension 1 landscape components into a single curve
#' and optionally interpolates it to a specified grid.
#'
#' @param landscape_obj Persistence landscape object with `dim0` and `dim1`
#'   components, each either a vector or matrix.
#' @param grid Numeric vector of evaluation points.
#'
#' @return Numeric vector representing the aggregated landscape curve.
#' @examples
#' \dontrun{
#' land <- list(dim0 = matrix(1, 5, 1), dim1 = matrix(0, 5, 1))
#' compute_landscape_curve(land, grid = seq(0, 1, length.out = 5))
#' }
#' @keywords internal
compute_landscape_curve <- function(landscape_obj, grid = seq(0, 1, length.out = grid_points)) {
  if (is.null(landscape_obj)) return(rep(0, length(grid)))
  curve0 <- if (is.matrix(landscape_obj$dim0)) {
    rowMeans(landscape_obj$dim0, na.rm = TRUE)
  } else {
    landscape_obj$dim0
  }
  curve1 <- if (is.matrix(landscape_obj$dim1)) {
    rowMeans(landscape_obj$dim1, na.rm = TRUE)
  } else {
    landscape_obj$dim1
  }
  aggregated_curve <- sqrt(curve0 ^ 2 + curve1 ^ 2)
  if (length(aggregated_curve) != length(grid)) {
    aggregated_curve <- approx(
      x = seq(min(grid), max(grid), length.out = length(aggregated_curve)),
      y = aggregated_curve,
      xout = grid,
      rule = 2
    )$y
  }
  aggregated_curve
}

#' Compute an aggregated landscape curve for a group
#'
#' Calculates the mean landscape curve across a list of persistence landscapes.
#'
#' @param landscape_list_group List of persistence landscape objects.
#' @param grid Numeric vector of evaluation points.
#'
#' @return Numeric vector representing the averaged landscape curve.
#' @examples
#' \dontrun{
#' lands <- list(list(dim0 = matrix(1,5,1), dim1 = matrix(0,5,1)))
#' compute_aggregated_landscape_curve(lands, grid = seq(0,1,length.out=5))
#' }
#' @keywords internal
compute_aggregated_landscape_curve <- function(landscape_list_group, grid = seq(0, 1, length.out = grid_points)) {
  curves <- sapply(landscape_list_group, function(land) {
    compute_landscape_curve(land, grid = grid)
  })
  if (is.null(dim(curves))) {
    curves
  } else {
    rowMeans(curves, na.rm = TRUE)
  }
}

compute_null_stats_aggregated_landscapes <- function(pd_list, landscape_list, tau_vals, n_bootstrap = 50, num_cores = 8, grid_points = 500) {
  integrated_diff_local <- function(curve_diff, grid) {
    delta <- diff(grid)
    sqrt(sum(delta * (curve_diff[-length(curve_diff)] ^ 2 + curve_diff[-1] ^ 2) / 2, na.rm = TRUE))
  }
  N <- length(pd_list)
  if (N < 2) return(NULL)
  boot_effects <- numeric(n_bootstrap)
  boot_ks <- numeric(n_bootstrap)
  grid <- seq(0, 1, length.out = grid_points)
  for (rep in 1:n_bootstrap) {
    set.seed(100 + rep)
    indices <- sample(N)
    half <- floor(N / 2)
    group1_ids <- indices[1:half]
    group2_ids <- indices[(half + 1):N]
    agg1 <- compute_aggregated_landscape_curve(landscape_list[group1_ids], grid = grid)
    agg2 <- compute_aggregated_landscape_curve(landscape_list[group2_ids], grid = grid)
    diff_val <- integrated_diff_local(agg1 - agg2, grid)
    sd_pooled <- sd(c(agg1, agg2), na.rm = TRUE)
    effect <- ifelse(sd_pooled == 0, 0, diff_val / sd_pooled)
    boot_effects[rep] <- effect
    if (sum(agg1) == 0 || sum(agg2) == 0) {
      ks_stat <- NA
    } else {
      cdf1 <- cumsum(agg1 / sum(agg1))
      cdf2 <- cumsum(agg2 / sum(agg2))
      ks_stat <- max(abs(cdf1 - cdf2), na.rm = TRUE)
    }
    boot_ks[rep] <- ks_stat
  }
  list(
    null_effect = list(
      mean = mean(boot_effects, na.rm = TRUE),
      sd = sd(boot_effects, na.rm = TRUE),
      quantiles = quantile(boot_effects, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    ),
    null_ks = list(
      mean = mean(boot_ks, na.rm = TRUE),
      sd = sd(boot_ks, na.rm = TRUE),
      quantiles = quantile(boot_ks, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    )
  )
}

#' Bootstrap an aggregated landscape curve
#'
#' Creates confidence bands for an aggregated landscape curve via bootstrap
#' resampling of individual curves.
#'
#' @param individual_curves List of individual landscape curves (numeric
#'   vectors).
#' @param grid_points Number of evaluation points along each curve.
#' @param n_bootstrap Number of bootstrap samples to draw.
#'
#' @return List containing the mean curve, lower and upper confidence bands, and
#'   the bootstrap samples.
#' @examples
#' \dontrun{
#' curves <- list(runif(10), runif(10))
#' bootstrap_aggregated_landscape_curve(curves, grid_points = 10, n_bootstrap = 50)
#' }
#' @keywords internal
bootstrap_aggregated_landscape_curve <- function(individual_curves, grid_points, n_bootstrap = 100) {
  boot_samples <- replicate(n_bootstrap, {
    sampled <- sample(individual_curves, replace = TRUE)
    rowMeans(do.call(cbind, sampled), na.rm = TRUE)
  })
  mean_curve <- rowMeans(boot_samples, na.rm = TRUE)
  lower_curve <- apply(boot_samples, 1, quantile, probs = 0.025, na.rm = TRUE)
  upper_curve <- apply(boot_samples, 1, quantile, probs = 0.975, na.rm = TRUE)
  log_message("Checking bootstrap results...")
  log_message(paste("Are the mean and upper curves identical?", all.equal(mean_curve, upper_curve)))
  list(
    mean = mean_curve,
    lower = lower_curve,
    upper = upper_curve,
    bootstrap_curves = boot_samples
  )
}

# --- UNIFIED SCALED PLOTTING FUNCTION ---
analyze_and_plot_curves <- function(curve1, curve2, grp1, grp2, curve_type, tau_vals, tau_max,
                                    alpha, effect_threshold, num_permutations,
                                    null_effect_stats = NULL, null_ks_stats = NULL,
                                    verbose = TRUE) {
  log_message <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", Sys.time(), msg))
  }
  
  integrated_diff_local <- function(curve_diff, tau_vals) {
    if (length(curve_diff) < 2) {
      log_message("curve_diff has less than 2 points; returning 0 for integrated difference.")
      return(0)
    }
    if (length(curve_diff) != length(tau_vals)) {
      log_message("Interpolating curve_diff to match tau_vals length.")
      curve_diff <- approx(
        x = seq(0, 1, length.out = length(curve_diff)),
        y = curve_diff,
        xout = seq(0, 1, length.out = length(tau_vals)),
        rule = 2)$y
    }
    delta_tau <- diff(tau_vals)
    sqrt(sum(delta_tau * (curve_diff[-length(curve_diff)] ^ 2 + curve_diff[-1] ^ 2) / 2, na.rm = TRUE))
  }
  
  build_plot_data <- function(curve, grp) {
    data.frame(Tau = tau_vals, Mean = curve$mean,
               Lower = curve$lower, Upper = curve$upper, Group = grp)
  }
  
  pd1 <- build_plot_data(curve1, grp1)
  pd2 <- build_plot_data(curve2, grp2)
  plot_data <- rbind(pd1, pd2)
  
  if (grepl("Landscape", curve_type, ignore.case = TRUE)) {
    if (!is.null(curve1$individual_landscape_curves) && !is.null(curve2$individual_landscape_curves)) {
      inds1 <- do.call(cbind, curve1$individual_landscape_curves)
      inds2 <- do.call(cbind, curve2$individual_landscape_curves)
    } else {
      inds1 <- inds2 <- NULL
    }
  } else if (grepl("Betti", curve_type, ignore.case = TRUE)) {
    if (!is.null(curve1$individual_pd_curves) && !is.null(curve2$individual_pd_curves)) {
      inds1 <- do.call(cbind, curve1$individual_pd_curves)
      inds2 <- do.call(cbind, curve2$individual_pd_curves)
    } else {
      inds1 <- inds2 <- NULL
    }
  } else if (grepl("Euler", curve_type, ignore.case = TRUE)) {
    if (!is.null(curve1$individual_euler_curves) && !is.null(curve2$individual_euler_curves)) {
      inds1 <- do.call(cbind, curve1$individual_euler_curves)
      inds2 <- do.call(cbind, curve2$individual_euler_curves)
    } else {
      inds1 <- inds2 <- NULL
    }
  } else {
    inds1 <- inds2 <- NULL
  }
  
  if (!is.null(inds1) && !is.null(inds2) && is.matrix(inds1) && is.matrix(inds2) &&
      ncol(inds1) > 0 && ncol(inds2) > 0) {
    m1 <- rowMeans(inds1, na.rm = TRUE)
    m2 <- rowMeans(inds2, na.rm = TRUE)
  } else {
    log_message("Falling back to using the aggregated 'mean' fields for stats.")
    m1 <- curve1$mean
    m2 <- curve2$mean
  }
  
  if (length(m1) < 2) m1 <- rep(0, length(tau_vals))
  if (length(m2) < 2) m2 <- rep(0, length(tau_vals))
  
  if (length(m1) != length(tau_vals)) {
    m1 <- approx(x = seq(0, 1, length.out = length(m1)), y = m1, xout = seq(0, 1, length.out = length(tau_vals)), rule = 2)$y
  }
  if (length(m2) != length(tau_vals)) {
    m2 <- approx(x = seq(0, 1, length.out = length(m2)), y = m2, xout = seq(0, 1, length.out = length(tau_vals)), rule = 2)$y
  }

  tau_rescaled <- tau_vals / max(tau_vals)
  wass <- tryCatch({
    a <- m1 / sum(m1)
    b <- m2 / sum(m2)
    F_a <- cumsum(a)
    F_b <- cumsum(b)
    dx <- diff(tau_rescaled)
    diff_cdf <- abs(F_a - F_b)
    sum(0.5 * (diff_cdf[-length(diff_cdf)] + diff_cdf[-1]) * dx)
  }, error = function(e) {
    log_message(paste("Wasserstein distance failed:", e$message))
    NA
  })
  
  obs_diff <- integrated_diff_local(m1 - m2, tau_vals)
  
  if (!is.null(inds1) && !is.null(inds2) && is.matrix(inds1) && is.matrix(inds2) &&
      ncol(inds1) > 0 && ncol(inds2) > 0) {
    all_inds <- cbind(inds1, inds2)
    n1 <- ncol(inds1)
    n_total <- ncol(all_inds)
    perm_diffs <- replicate(num_permutations, {
      perm <- sample(n_total)
      m1p <- rowMeans(all_inds[, perm[1:n1], drop = FALSE], na.rm = TRUE)
      m2p <- rowMeans(all_inds[, perm[(n1 + 1):n_total], drop = FALSE], na.rm = TRUE)
      integrated_diff_local(m1p - m2p, tau_vals)
    })
    perm_p <- mean(perm_diffs >= obs_diff, na.rm = TRUE)
    ks_res <- lapply(1:nrow(inds1), function(i) {
      res <- tryCatch(ks.test(inds1[i,], inds2[i,]), error = function(e) NULL)
      if (is.null(res)) return(c(stat = NA, p = NA))
      c(stat = as.numeric(res$statistic), p = as.numeric(res$p.value))
    })
    ks_stats <- sapply(ks_res, function(x) x["stat"])
    ks_pvals <- sapply(ks_res, function(x) x["p"])
    ks_combined_stat <- max(ks_stats, na.rm = TRUE)
    ks_adj_pval <- min(p.adjust(ks_pvals, method = "BH"), na.rm = TRUE)
  } else {
    perm_p <- NA
    ks_combined_stat <- NA
    ks_adj_pval <- NA
  }
  
  eff_size <- ifelse(sd(c(m1, m2), na.rm = TRUE) == 0, 0, obs_diff / sd(c(m1, m2), na.rm = TRUE))
  
  if (grepl("Euler", curve_type, ignore.case = TRUE) || grepl("Landscape", curve_type, ignore.case = TRUE)) {
    extract_null_values <- function(null_stats) {
      if (is.null(null_stats)) return(list(median = NA, lower = NA, upper = NA))
      if (!is.null(null_stats$median)) {
        list(median = null_stats$median, lower = null_stats$lower, upper = null_stats$upper)
      } else if (!is.null(null_stats$quantiles)) {
        list(median = as.numeric(null_stats$quantiles["50%"]),
             lower = as.numeric(null_stats$quantiles["2.5%"]),
             upper = as.numeric(null_stats$quantiles["97.5%"]))
      } else {
        list(median = NA, lower = NA, upper = NA)
      }
    }
    null_effect_vals <- extract_null_values(null_effect_stats)
    null_ks_vals <- extract_null_values(null_ks_stats)
    annot_null <- paste("\nNull median effect size:",
                        ifelse(is.na(null_effect_vals$median), "NA", signif(null_effect_vals$median, 3)),
                        "\nNull effect 95% CI: [",
                        ifelse(is.na(null_effect_vals$lower), "NA", signif(null_effect_vals$lower, 3)),
                        ",",
                        ifelse(is.na(null_effect_vals$upper), "NA", signif(null_effect_vals$upper, 3)),
                        "]",
                        "\nNull median KS:",
                        ifelse(is.na(null_ks_vals$median), "NA", signif(null_ks_vals$median, 3)),
                        "\nNull KS 95% CI: [",
                        ifelse(is.na(null_ks_vals$lower), "NA", signif(null_ks_vals$lower, 3)),
                        ",",
                        ifelse(is.na(null_ks_vals$upper), "NA", signif(null_ks_vals$upper, 3)),
                        "]")
  } else {
    annot_null <- ""
  }
  
  sig <- ifelse(!is.na(perm_p) && perm_p < alpha, "sig", "notsig")
  annot <- paste("Wasserstein:", signif(wass, 3),
                 "\nPerm p:", ifelse(is.na(perm_p), "NA", signif(perm_p, 3)),
                 "\nEffect size:", signif(eff_size, 3),
                 annot_null,
                 "\nMax KS stat:", ifelse(is.na(ks_combined_stat), "NA", signif(ks_combined_stat, 3)),
                 "\nSignificance:", sig)
  
  p <- ggplot(plot_data, aes(x = Tau, y = Mean, color = Group, group = Group)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group), alpha = 0.2) +
    labs(
      title = paste(curve_type, "Comparison:", grp1, "vs.", grp2),
      x = "Normalized Filtration Scale (Tau)",
      y = paste(curve_type, "Value")
    ) +
    # scale_x_continuous(labels = scales::percent_format(scale = 100, accuracy = 1), limits = c(0, 1)) +
    annotate("text", x = 0.6, y = Inf, label = annot, hjust = 0, vjust = 1.2, size = 4) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
  scale_x_log10(
    "Normalized Filtration Scale (Tau, log scale)",
    labels = function(breaks) scales::percent(breaks / tau_max, accuracy = 0.01)
  )

  list(stats = list(wasserstein = wass,
                    perm_p = perm_p,
                    effect_size = eff_size,
                    ks_stat = ks_combined_stat,
                    ks_adj_p = ks_adj_pval), 
       plot = p)
}

# --- MODIFICATION: Added tau_max to the function definition ---
perform_pairwise_comparisons <- function(curves, euler_curves, groups, tau_vals, tau_max, dimensions = c(0, 1),
                                         dataset_name = NULL, group_by_col = NULL, plot_output_dir,
                                         alpha = 0.05, effect_threshold = 0.5, num_permutations = 1000,
                                         null_effect_stats = NULL, null_ks_stats = NULL) {
  log_message("Performing pairwise comparisons.")
  sanitize_group <- function(name) gsub(" ", "_", name)
  valid_groups <- groups[!is.na(groups) & trimws(groups) != ""]
  if (length(valid_groups) < 1) {
    log_message("No valid groups.")
    return(NULL)
  }
  res <- list()
  for (i in seq_along(valid_groups)) {
    for (j in i:length(valid_groups)) {
      grp1 <- valid_groups[i]
      grp2 <- valid_groups[j]
      log_message(paste("Comparing", grp1, "vs", grp2))
      comp <- list()
      plots <- list()
      if (is.null(dimensions)) {
        c1 <- tryCatch(curves[[grp1]], error = function(e) NULL)
        c2 <- tryCatch(curves[[grp2]], error = function(e) NULL)
        if (!is.null(c1) && is.atomic(c1)) c1 <- list(mean = c1, lower = c1, upper = c1)
        if (!is.null(c2) && is.atomic(c2)) c2 <- list(mean = c2, lower = c2, upper = c2)
        # --- MODIFICATION: Pass tau_max to the plotting function ---
        ar <- analyze_and_plot_curves(c1, c2, grp1, grp2, "Landscape", tau_vals, tau_max,
                                      alpha, effect_threshold, num_permutations,
                                      null_effect_stats, null_ks_stats)
        if (!is.null(ar$stats)) comp[["landscape"]] <- ar$stats
        if (!is.null(ar$plot)) plots <- c(plots, list(ar$plot))
      } else {
        for (d in dimensions) {
          c1 <- tryCatch(curves[[grp1]][[d + 1]], error = function(e) NULL)
          c2 <- tryCatch(curves[[grp2]][[d + 1]], error = function(e) NULL)
          if (is.null(c1) || is.null(c2)) next
          # --- MODIFICATION: Pass tau_max to the plotting function ---
          ar <- analyze_and_plot_curves(c1, c2, grp1, grp2, paste("Betti Dim", d), tau_vals, tau_max,
                                        alpha, effect_threshold, num_permutations)
          if (!is.null(ar$stats)) comp[[paste0("dimension_", d)]] <- ar$stats
          if (!is.null(ar$plot)) plots <- c(plots, list(ar$plot))
        }
        e1 <- tryCatch(euler_curves[[grp1]], error = function(e) NULL)
        e2 <- tryCatch(euler_curves[[grp2]], error = function(e) NULL)
        if (!is.null(e1) && !is.null(e2)) {
          # --- MODIFICATION: Pass tau_max to the plotting function ---
          er <- analyze_and_plot_curves(e1, e2, grp1, grp2, "Euler", tau_vals, tau_max,
                                        alpha, effect_threshold, num_permutations,
                                        null_effect_stats, null_ks_stats)
          comp[["euler"]] <- er$stats
          plots <- c(plots, list(er$plot))
        }
      }
      if (length(plots) > 0) {
        for (p_idx in seq_along(plots)) {
          p <- plots[[p_idx]]
          plot_title <- p$labels$title
          curve_type_sanitized <- gsub("[[:space:]]+", "_", plot_title) 
          curve_type_sanitized <- gsub("[^A-Za-z0-9_vs-]", "", curve_type_sanitized) 
          png_file_name <- paste0(dataset_name, "_", group_by_col, "_", curve_type_sanitized, ".png")
          png_file_path <- file.path(plot_output_dir, group_by_col, png_file_name)
          ggsave(
            filename = png_file_path,
            plot = p,
            width = 12,
            height = 10,
            dpi = 300
          )
          log_message(paste("Saved individual PNG plot to", png_file_path))
        }
        if (is.null(dimensions)) {
          file_name_pdf <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, "_",
                                  sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.pdf")
          file_name_svg <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, "_",
                                  sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.svg")
        } else {
          file_name_pdf <- paste0(dataset_name, "_", group_by_col, "_",
                                  sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.pdf")
          file_name_svg <- paste0(dataset_name, "_", group_by_col, "_",
                                  sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.svg")
        }
        plot_file_pdf <- file.path(plot_output_dir, group_by_col, file_name_pdf)
        ensure_directory(dirname(plot_file_pdf))
        pdf(plot_file_pdf, width = 12, height = 10 * length(plots))
        gridExtra::grid.arrange(grobs = plots, ncol = 1)
        dev.off()
        log_message(paste("Saved plot to", plot_file_pdf))
        plot_file_svg <- file.path(plot_output_dir, group_by_col, file_name_svg)
        ensure_directory(dirname(plot_file_svg))
        svglite::svglite(plot_file_svg, width = 12, height = 10 * length(plots))
        gridExtra::grid.arrange(grobs = plots, ncol = 1)
        dev.off()
        log_message(paste("Saved plot to", plot_file_svg))
      }
      res[[paste0(grp1, "_vs_", grp2)]] <- comp
    }
  }
  res
}

# --- UNIFIED SCALED PLOTTING FUNCTION ---
plot_aggregated_curves <- function(group_curves, groups, dimensions, tau_vals, tau_max, dataset_name, group_by_col, euler_curves = NULL, results_folder) {
  log_message("Generating aggregated plots.")
  
  if (is.null(dimensions)) {
    out_dir <- file.path(results_folder, "plots", "landscape_curves", "aggregated_curves")
    curve_label <- "Landscape"
  } else {
    out_dir <- file.path(results_folder, "plots", "betti_plots", "aggregated_curves")
    curve_label <- NULL
  }
  ensure_directory(out_dir)
  
  gen_plot <- function(data, curve_type, dim = NULL) {
    data$Group <- as.factor(data$Group)
    
    ggplot(data, aes(x = Tau, y = Mean, color = Group, group = Group)) +
      geom_line(linewidth = 1.2, alpha = 0.6) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group), alpha = 0.2) +
      labs(
        title = paste(curve_type, "Curve", if (!is.null(dim)) paste("(Dim", dim, ")") else ""),
        x = "Normalized Filtration Scale (Tau)", 
        y = paste(curve_type, "Value")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      ) +
      # scale_x_continuous(labels = scales::percent_format(scale = 100, accuracy = 1), limits = c(0, 1)) +
      # coord_trans(x = "log10") + # Add this line for the log scale
      scale_x_log10(
      "Normalized Filtration Scale (Tau, log scale)",
      labels = function(breaks) scales::percent(breaks / tau_max, accuracy = 0.01)
    ) +
      scale_color_viridis_d() + 
      scale_fill_viridis_d()    
  }
  
  plots <- list()
  if (is.null(dimensions)) {
    aggregated_data_list <- list()
    for (group in groups) {
      curve <- tryCatch(group_curves[[group]], error = function(e) NULL)
      if (!is.null(curve)) {
        if (is.atomic(curve)) curve <- list(mean = curve, lower = curve, upper = curve)
        df <- data.frame(Tau = tau_vals, Mean = curve$mean, Lower = curve$lower, Upper = curve$upper, Group = group)
        aggregated_data_list[[group]] <- df
      }
    }
    if (length(aggregated_data_list) > 0) {
      all_data <- do.call(rbind, aggregated_data_list)
      plots <- list(gen_plot(all_data, "Landscape"))
    }
  } else {
    betti_data_list <- list("0" = data.frame(), "1" = data.frame())
    euler_data <- data.frame()
    for (group in groups) {
      for (i in 1:2) {
        curve <- tryCatch(group_curves[[group]][[i]], error = function(e) NULL)
        if (!is.null(curve)) {
          dim_label <- as.character(i - 1)
          df <- data.frame(Tau = tau_vals, Mean = curve$mean, Lower = curve$lower, Upper = curve$upper, Group = group)
          betti_data_list[[dim_label]] <- rbind(betti_data_list[[dim_label]], df)
        }
      }
      ec <- tryCatch(euler_curves[[group]], error = function(e) NULL)
      if (!is.null(ec)) {
        df <- data.frame(Tau = tau_vals, Mean = ec$mean, Lower = ec$lower, Upper = ec$upper, Group = group)
        euler_data <- rbind(euler_data, df)
      }
    }
    for (dim_label in names(betti_data_list)) {
      if (nrow(betti_data_list[[dim_label]]) > 0) {
        plots <- c(plots, list(gen_plot(betti_data_list[[dim_label]], "Betti", dim_label)))
      }
    }
    if (nrow(euler_data) > 0) {
      plots <- c(plots, list(gen_plot(euler_data, "Euler")))
    }
  }
  
  if (length(plots) > 0) {
    for (p_idx in seq_along(plots)) {
      p <- plots[[p_idx]]
      plot_title <- p$labels$title
      curve_type_sanitized <- gsub("[[:space:]]+", "_", plot_title)
      curve_type_sanitized <- gsub("[^A-Za-z0-9_]", "", curve_type_sanitized)
      png_file_name <- paste0(dataset_name, "_aggregated_", group_by_col, "_", curve_type_sanitized, ".png")
      png_file_path <- file.path(out_dir, png_file_name)
      ggsave(
        filename = png_file_path,
        plot = p,
        width = 12,
        height = 10,
        dpi = 300
      )
      log_message(paste("Saved individual aggregated PNG plot to", png_file_path))
    }
    if (is.null(dimensions)) {
      file_name_pdf <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, ".pdf")
      file_name_svg <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, ".svg")
    } else {
      file_name_pdf <- paste0(dataset_name, "_betti_and_euler_curves_by_", group_by_col, ".pdf")
      file_name_svg <- paste0(dataset_name, "_betti_and_euler_curves_by_", group_by_col, ".svg")
    }
    plot_file_pdf <- file.path(out_dir, file_name_pdf)
    ensure_directory(dirname(plot_file_pdf))
    pdf(plot_file_pdf, width = 12, height = 10 * length(plots))
    gridExtra::grid.arrange(grobs = plots, ncol = 1)
    dev.off()
    log_message(paste("Saved aggregated plot to", plot_file_pdf))
    plot_file_svg <- file.path(out_dir, file_name_svg)
    ensure_directory(dirname(plot_file_svg))
    svglite::svglite(plot_file_svg, width = 12, height = 10 * length(plots))
    gridExtra::grid.arrange(grobs = plots, ncol = 1)
    dev.off()
    log_message(paste("Saved aggregated plot to", plot_file_svg))
  } else {
    log_message("No plots to save.")
  }
}

compute_null_stats_cross_iterations <- function(pd_list, metadata_list, tau_vals,
                                                dimensions = c(0, 1),
                                                dataset_name, base_sigma, grid_points, tau_max,
                                                results_folder, verbose = TRUE, num_cores = 16,
                                                landscape_list = NULL) {
  integrated_diff <- function(curve_diff, tau_vals) {
    delta_tau <- diff(tau_vals)
    integrated_sq <- sum(delta_tau * (curve_diff[-length(curve_diff)] ^ 2 + curve_diff[-1] ^ 2) / 2)
    sqrt(integrated_sq)
  }
  log_message <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", Sys.time(), msg))
  }
  combined_meta <- do.call(rbind, metadata_list)
  bs_cols <- grep("^random_group_bootstrap_", colnames(combined_meta), value = TRUE, ignore.case = TRUE)
  results_list <- parallel::mclapply(bs_cols, function(col) {
    group1_idx <- which(combined_meta[[col]] == 1)
    group2_idx <- which(combined_meta[[col]] == 2)
    if (length(group1_idx) == 0 || length(group2_idx) == 0) return(NULL)
    pd_group1 <- pd_list[group1_idx]
    pd_group2 <- pd_list[group2_idx]
    curve1 <- compute_euler_curve(
      pd_subset = pd_group1, dimensions = dimensions, group_name = "group1",
      dataset_name = dataset_name, grid_points = grid_points,
      tau_max = tau_max, base_sigma = base_sigma, n_bootstrap = 50, results_folder = results_folder
    )
    curve2 <- compute_euler_curve(
      pd_subset = pd_group2, dimensions = dimensions, group_name = "group2",
      dataset_name = dataset_name, grid_points = grid_points,
      tau_max = tau_max, base_sigma = base_sigma, n_bootstrap = 50, results_folder = results_folder
    )
    mean_curve1 <- curve1$mean
    mean_curve2 <- curve2$mean
    diff_val <- integrated_diff(mean_curve1 - mean_curve2, tau_vals)
    denom <- sd(c(mean_curve1, mean_curve2), na.rm = TRUE)
    effect_size <- if (denom == 0) 0 else diff_val / denom
    individual_curves1 <- do.call(cbind, curve1$individual_euler_curves)
    individual_curves2 <- do.call(cbind, curve2$individual_euler_curves)
    n_tau <- nrow(individual_curves1)
    ks_stats <- numeric(n_tau)
    for (i in seq_len(n_tau)) {
      test_result <- tryCatch(ks.test(individual_curves1[i,], individual_curves2[i,]), error = function(e) NULL)
      ks_stats[i] <- ifelse(is.null(test_result), NA, as.numeric(test_result$statistic))
    }
    ks_stat <- max(ks_stats, na.rm = TRUE)
    list(effect_size = effect_size, ks_stat = ks_stat)
  }, mc.cores = num_cores)
  results_list <- results_list[!sapply(results_list, is.null)]
  null_effects <- sapply(results_list, function(x) x$effect_size)
  null_ks <- sapply(results_list, function(x) x$ks_stat)
  betti_euler_nulls <- list(
    null_effect = if (length(null_effects) > 0) {
      list(mean = mean(null_effects), sd = sd(null_effects),
           quantiles = quantile(null_effects, probs = c(0.025, 0.5, 0.975)))
    } else { NULL },
    null_ks = if (length(null_ks) > 0) {
      list(mean = mean(null_ks), sd = sd(null_ks),
           quantiles = quantile(null_ks, probs = c(0.025, 0.5, 0.975)))
    } else { NULL }
  )
  landscape_nulls <- NULL
  if (!is.null(landscape_list)) {
    landscape_nulls <- bootstrap_null_stats_landscape(landscape_list, n_bootstrap = 50, grid_points = grid_points)
  }
  list(betti_euler = betti_euler_nulls, landscape = landscape_nulls)
}

compute_and_compare_betti_iterations <- function(
  pd_list, metadata_list, group_by_col = NULL,
  dimensions = c(0, 1), grid_points = 500, num_permutations = 10000,
  bootstrap_samples = 100, dataset_name = NULL, comparison_type = "group",
  base_sigma = 1, num_cores = 20, results_folder = "results", verbose = TRUE,
  landscape_list = NULL
) {
  log_message <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", Sys.time(), msg))
  }
  generate_cache_filename <- function(dataset_name, hash, results_folder) {
    path <- file.path(results_folder, "plots", "betti_plots", "betti_cache", dataset_name, paste0("cache_", hash, ".rds"))
    ensure_directory(dirname(path))
    path
  }
  compute_betti_curve <- function(pd, dimension, grid_points, tau_max, base_sigma) {
    pd <- as.data.frame(pd)
    valid_pd <- pd[pd[, 1] == dimension & pd[, 2] < pd[, 3],]
    if (nrow(valid_pd) == 0) return(rep(0, grid_points))
    tau_grid <- seq(0, tau_max, length.out = grid_points)
    gaussian <- function(x, mu, sigma) {
      exp(-((x - mu) ^ 2) / (2 * sigma ^ 2)) / (sqrt(2 * pi) * sigma)
    }
    sapply(tau_grid, function(tau) {
      sum(sapply(1:nrow(valid_pd), function(i) {
        b <- valid_pd[i, 2]
        d <- valid_pd[i, 3]
        persistence <- d - b
        interval_sigma <- max(base_sigma, persistence / 2)
        persistence * gaussian(tau, (b + d) / 2, interval_sigma)
      }))
    })
  }
  generate_hash <- function(pd, dimension, dataset_name, base_sigma, grid_points, tau_max) {
    digest(list(pd, dimension, dataset_name, base_sigma, grid_points, tau_max), algo = "sha256")
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
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(curve, cache_file)
    log_message(paste("Saved cached curve to:", cache_file))
  }
  get_or_compute_curve <- function(pd, dimensions, dataset_name, base_sigma, grid_points, tau_max, results_folder) {
    hash <- generate_hash(pd, dimensions, dataset_name, base_sigma, grid_points, tau_max)
    cached_curve <- load_curve_cache(dataset_name, hash, results_folder)
    if (!is.null(cached_curve)) return(cached_curve)
    log_message("Cache miss for curve. Computing curve.")
    curve <- compute_betti_curve(pd, dimensions[1], grid_points, tau_max, base_sigma)
    save_curve_cache(curve, dataset_name, hash, results_folder)
    curve
  }
  bootstrap_curve <- function(pd_subset, dimension, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, seed = 42, results_folder) {
    log_message(paste("Starting bootstrapping for Betti curve in dimension", dimension, "with", n_bootstrap, "samples."))
    individual_pd_curves <- lapply(seq_along(pd_subset), function(pd_idx) {
      pd <- pd_subset[[pd_idx]]
      if (is.matrix(pd) && nrow(pd) > 0) {
        curve <- get_or_compute_curve(pd, dimension, dataset_name, base_sigma, grid_points, tau_max, results_folder)
        return(curve / sum(curve))
      } else {
        return(rep(0, grid_points))
      }
    })
    names(individual_pd_curves) <- paste0("pd_", seq_along(pd_subset))
    bootstrap_curves <- lapply(1:n_bootstrap, function(rep) {
      set.seed(seed + rep)
      sampled_pds <- pd_subset[sample(seq_along(pd_subset), replace = TRUE)]
      aggregated_curve <- numeric(grid_points)
      for (pd in sampled_pds) {
        if (is.matrix(pd) && nrow(pd) > 0) {
          curve <- get_or_compute_curve(pd, dimension, dataset_name, base_sigma, grid_points, tau_max, results_folder)
          aggregated_curve <- aggregated_curve + curve
        }
      }
      aggregated_curve / sum(aggregated_curve)
    })
    valid_curves <- bootstrap_curves[!sapply(bootstrap_curves, is.null)]
    if (length(valid_curves) == 0) {
      return(list(mean = rep(0, grid_points), lower = rep(0, grid_points),
                  upper = rep(0, grid_points), individual_pd_curves = individual_pd_curves))
    }
    aggregated_results <- list(
      mean = rowMeans(do.call(cbind, valid_curves), na.rm = TRUE),
      lower = apply(do.call(cbind, valid_curves), 1, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(do.call(cbind, valid_curves), 1, quantile, probs = 0.975, na.rm = TRUE),
      individual_pd_curves = individual_pd_curves
    )
    log_message(paste("Completed bootstrapping for Betti curve in dimension", dimension))
    aggregated_results
  }
  tau_max <- max(unlist(lapply(pd_list, function(pd) max(pd[, 3], na.rm = TRUE))), na.rm = TRUE)
  tau_vals <- seq(0.0001, tau_max, length.out = grid_points)
  group_metadata_mapping <- do.call(rbind, lapply(seq_along(metadata_list), function(i) {
    metadata <- metadata_list[[i]]
    pd_names <- names(pd_list)[names(pd_list) %in% metadata$orig.ident.Iter]
    metadata[metadata$orig.ident.Iter %in% pd_names,]
  }))
  groups <- unique(group_metadata_mapping[[group_by_col]])
  log_message(paste("Groups identified:", paste(groups, collapse = ", ")))
  group_results <- parallel::mclapply(as.character(groups), function(group) {
    log_message(paste("Processing group:", group))
    group_pd <- pd_list[group_metadata_mapping$orig.ident.Iter[group_metadata_mapping[[group_by_col]] == group]]
    if (length(group_pd) == 0) {
      log_message(paste("No PDs for group:", group))
      return(NULL)
    }
    tryCatch({
      bootstrap_curves <- lapply(dimensions, function(dim) {
        log_message(paste("Bootstrapping Betti for dimension", dim, "in group", group))
        bootstrap_curve(group_pd, dim, group, dataset_name, grid_points, tau_max, base_sigma, bootstrap_samples, results_folder = results_folder)
      })
      bootstrap_curves <- list(betti = bootstrap_curves)
      log_message(paste("Computing Euler curve for group:", group))
      ec <- compute_euler_curve(group_pd, dimensions, group, dataset_name, grid_points, tau_max, base_sigma, bootstrap_samples, results_folder)
      ec <- list(euler = ec)
      consistency <- check_betti_euler_consistency(ec$euler, bootstrap_curves$betti, dimensions, tau_vals)
      list(group = group, betti_curves = setNames(bootstrap_curves, group), euler_curve = setNames(ec, group), consistency = consistency)
    }, error = function(e) {
      log_message(paste("Error processing group:", group, "-", e$message))
      NULL
    })
  }, mc.cores = num_cores)
  group_results <- Filter(Negate(is.null), group_results)
  betti_curves_by_group <- setNames(lapply(group_results, function(res) {
    curves <- res$betti_curves[[res$group]]
    setNames(curves, paste0("dimension_", seq_along(curves)))
  }), vapply(group_results, `[[`, "", "group"))
  euler_curves_by_group <- setNames(lapply(group_results, function(res) {
    res$euler_curve[[res$group]]
  }), vapply(group_results, `[[`, "", "group"))
  aggregated_landscape_curves_by_group <- NULL
  if (!is.null(landscape_list)) {
    aggregated_landscape_curves_by_group <- setNames(lapply(as.character(groups), function(group) {
      group_ids <- group_metadata_mapping$orig.ident.Iter[group_metadata_mapping[[group_by_col]] == group]
      grid <- seq(0, 1, length.out = grid_points) # Use a 0-1 grid for landscapes
      if (length(group_ids) == 0) {
        return(list(
          mean = rep(0, grid_points), lower = rep(0, grid_points),
          upper = rep(0, grid_points), individual_landscape_curves = list()
        ))
      }
      individual_curves <- lapply(landscape_list[group_ids], function(land) {
        compute_landscape_curve(land, grid = grid)
      })
      bootstrapped_agg <- bootstrap_aggregated_landscape_curve(individual_curves, grid_points, n_bootstrap = bootstrap_samples)
      list(
        mean = bootstrapped_agg$mean, lower = bootstrapped_agg$lower,
        upper = bootstrapped_agg$upper, individual_landscape_curves = individual_curves,
        bootstrap_aggregates = bootstrapped_agg$bootstrap_curves
      )
    }), groups)
  }
  bootstrap_nulls <- compute_null_stats_cross_iterations(
    pd_list = pd_list, metadata_list = metadata_list, tau_vals = tau_vals,
    dimensions = dimensions, dataset_name = dataset_name,
    base_sigma = base_sigma, grid_points = grid_points, tau_max = tau_max, results_folder = results_folder
  )
  # --- MODIFICATION: Updated function call ---
  pairwise_results <- perform_pairwise_comparisons(
    curves = betti_curves_by_group, 
    euler_curves = euler_curves_by_group,
    groups = groups, 
    tau_vals = tau_vals, 
    tau_max = tau_max, # Pass tau_max
    dimensions = dimensions,
    dataset_name = dataset_name, 
    group_by_col = group_by_col,
    plot_output_dir = file.path(results_folder, "plots", "betti_plots", "cross_iteration_comparisons", dataset_name, "pairwise_comparisons"),
    null_effect_stats = bootstrap_nulls$betti_euler$null_effect,
    null_ks_stats = bootstrap_nulls$betti_euler$null_ks,
    alpha = 0.05, 
    effect_threshold = 0.5, 
    num_permutations = 1000
  )
  landscape_pairwise_results <- NULL
  if (!is.null(aggregated_landscape_curves_by_group)) {
    land_tau <- seq(0, 1, length.out = grid_points) # Landscape tau is 0-1
    landscape_nulls <- bootstrap_null_stats_landscape(
      landscape_list,
      n_bootstrap = 50,
      grid_points = grid_points
    )
    # --- MODIFICATION: Updated function call ---
    landscape_pairwise_results <- perform_pairwise_comparisons(
      curves = aggregated_landscape_curves_by_group,
      euler_curves = NULL,
      groups = groups,
      tau_vals = land_tau,
      tau_max = 1.0, # Pass tau_max for landscapes (which is 1.0)
      dimensions = NULL,
      dataset_name = dataset_name,
      group_by_col = group_by_col,
      plot_output_dir = file.path(
        results_folder,
        "plots", "landscape_curves", dataset_name,
        "pairwise_comparisons"
      ),
      null_effect_stats = landscape_nulls$null_effect,
      null_ks_stats = landscape_nulls$null_ks,
      alpha = 0.05,
      effect_threshold = 0.5,
      num_permutations = 1000
    )
  }
  log_message("Generating aggregated plots for Betti/Euler curves.")
  # --- MODIFICATION: Updated function call ---
  plot_aggregated_curves(
    group_curves = betti_curves_by_group, 
    groups = groups, 
    dimensions = dimensions,
    tau_vals = tau_vals, 
    tau_max = tau_max, # Pass tau_max
    dataset_name = dataset_name, 
    group_by_col = group_by_col,
    euler_curves = euler_curves_by_group, 
    results_folder = results_folder
  )
  if (!is.null(aggregated_landscape_curves_by_group)) {
    log_message("Generating aggregated landscape plots.")
    # --- MODIFICATION: Updated function call ---
    plot_aggregated_curves(
      group_curves = aggregated_landscape_curves_by_group, 
      groups = groups, 
      dimensions = NULL,
      tau_vals = seq(0, 1, length.out = grid_points), 
      tau_max = 1.0, # Pass tau_max for landscapes
      dataset_name = paste0(dataset_name, "_landscape"),
      group_by_col = group_by_col, 
      euler_curves = NULL, 
      results_folder = results_folder
    )
  }
  list(
    group_results = group_results,
    pairwise_results = pairwise_results,
    aggregated_landscape_curves = aggregated_landscape_curves_by_group,
    landscape_pairwise_results = landscape_pairwise_results,
    bootstrap_nulls = bootstrap_nulls
  )
}

cross_iteration_comparison_with_betti <- function(data_iterations,
                                                  group_by_col = "Tissue",
                                                  orig_ident_columns = c("orig.ident"),
                                                  tau_vals = seq(0, 1, length.out = 500),
                                                  dimensions = c(0, 1), grid_points = 500,
                                                  num_permutations = 10000, bootstrap_samples = 100,
                                                  base_sigma = 1,
                                                  output_folder = file.path(results_folder, "cross_iteration_comparisons"),
                                                  num_cores = 16,
                                                  verbose = TRUE,
                                                  metadata = NULL) {
    if (!all(colnames(seurat_obj) == reference_cell_names)) {
      if (ncol(seurat_obj) == length(reference_cell_names)) {
        colnames(seurat_obj) <- reference_cell_names
        data_iterations[[i]]$seurat_obj <- seurat_obj
        log_message(paste0("Iteration ", data_iterations[[i]]$name, " renamed to match reference cell names."))
      } else {
        log_message(paste0("Iteration ", data_iterations[[i]]$name, " has a different number of cells; skipping renaming."))
      }
    }
  process_group <- function(group_value) {
    log_message(paste("Processing", group_by_col, ":", group_value))
    combined_pd_list <- list()
    metadata_list <- list()
    combined_landscape_list <- list()
    iter_results <- parallel::mclapply(data_iterations, function(iter) {
      pd_list <- readRDS(iter$pd_list)
      if (!all(startsWith(names(pd_list), paste0(iter$name, "_")))) {
        valid_names <- paste0(iter$name, "_", names(pd_list))
        names(pd_list) <- valid_names
      } else {
        valid_names <- names(pd_list)
      }
      if (!is.null(iter$landscape_list)) {
        landscape_list <- readRDS(iter$landscape_list)
        if (!all(startsWith(names(landscape_list), paste0(iter$name, "_")))) {
          valid_landscape_names <- paste0(iter$name, "_", names(landscape_list))
          names(landscape_list) <- valid_landscape_names
        } else {
          valid_landscape_names <- names(landscape_list)
        }
      } else {
        landscape_list <- NULL
      }
      metadata_local <- iter$seurat_obj@meta.data
      metadata_local$Iteration <- iter$name
      if (!("orig.ident" %in% names(metadata_local)) || !identical(orig_ident_columns, "orig.ident")) {
        metadata_local$orig.ident <- do.call(paste, c(metadata_local[orig_ident_columns], sep = "_"))
      }
      metadata_local$orig.ident.Iter <- paste(iter$name, metadata_local$orig.ident, sep = "_")
      metadata_local$orig.ident.Iter <- as.character(metadata_local$orig.ident.Iter)
      group_metadata <- metadata_local[metadata_local[[group_by_col]] == group_value,]
      group_pd_list <- pd_list[valid_names %in% group_metadata$orig.ident.Iter]
      if (length(group_pd_list) == 0) {
        log_message(paste("No persistence diagrams for", group_by_col, ":", group_value, "in iteration:", iter$name))
        return(NULL)
      }
      pd_names_filtered <- valid_names[valid_names %in% group_metadata$orig.ident.Iter]
      names(group_pd_list) <- pd_names_filtered
      aligned_metadata <- lapply(pd_names_filtered, function(nm) {
        rows <- group_metadata[group_metadata$orig.ident.Iter == nm,]
        if (nrow(rows) > 0) return(rows[1,, drop = FALSE]) else return(NULL)
      })
      valid_idx <- which(!sapply(aligned_metadata, is.null))
      if (length(valid_idx) != length(pd_names_filtered)) {
        pd_names_filtered <- pd_names_filtered[valid_idx]
        group_pd_list <- group_pd_list[valid_idx]
        aligned_metadata <- aligned_metadata[valid_idx]
      }
      if (!is.null(landscape_list)) {
        group_landscape_list <- landscape_list[valid_landscape_names %in% group_metadata$orig.ident.Iter]
      } else {
        group_landscape_list <- NULL
      }
      list(
        pd_list = group_pd_list,
        metadata = setNames(aligned_metadata, pd_names_filtered),
        landscape_list = group_landscape_list
      )
    }, mc.cores = num_cores)
    iter_results <- Filter(Negate(is.null), iter_results)
    for (res in iter_results) {
      combined_pd_list <- c(combined_pd_list, res$pd_list)
      metadata_list <- c(metadata_list, res$metadata)
      if (!is.null(res$landscape_list)) {
        combined_landscape_list <- c(combined_landscape_list, res$landscape_list)
      }
    }
    log_message(paste("Combined pd_list entries:", length(combined_pd_list)))
    log_message(paste("Combined landscape entries:", length(combined_landscape_list)))
    if (length(metadata_list) < 2) {
      log_message(paste("Skipping", group_by_col, group_value, "- fewer than two iterations available."))
      return(NULL)
    }
    common_cols <- Reduce(intersect, lapply(metadata_list, colnames))
    metadata_list <- lapply(metadata_list, function(df) df[, common_cols, drop = FALSE])
    if (length(combined_pd_list) != length(metadata_list))
      stop("The length of pd_list and metadata_list must be the same after alignment.")
    results <- compute_and_compare_betti_iterations(
      pd_list = combined_pd_list,
      metadata_list = metadata_list,
      group_by_col = "Iteration",
      dimensions = dimensions,
      grid_points = grid_points,
      num_permutations = num_permutations,
      bootstrap_samples = bootstrap_samples,
      dataset_name = paste(group_by_col, group_value, sep = "_"),
      base_sigma = base_sigma,
      num_cores = num_cores,
      results_folder = output_folder,
      landscape_list = if (length(combined_landscape_list) > 0) combined_landscape_list else NULL
    )
    output_file <- file.path(output_folder, paste0("comparison_", group_by_col, "_", group_value, ".rds"))
    saveRDS(results, output_file)
    log_message(paste("Completed comparison for", group_by_col, ":", group_value, "- results saved to:", output_file))
    results
  }
  results <- parallel::mclapply(group_values, process_group, mc.cores = num_cores)
  log_message("Flattening all pairwise curve-comparison stats into one data.frame")
  metrics <- c("wasserstein", "perm_p", "effect_size", "ks_stat", "ks_adj_p")
  extract_row <- function(st, grp, pair, comp) {
    vals <- setNames(lapply(metrics, function(m) if (!is.null(st[[m]])) st[[m]] else NA_real_), metrics)
    data.frame(
      Group = grp,
      Pair = pair,
      Comparison = comp,
      as.data.frame(vals, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
  }
  betti_df <- do.call(rbind, lapply(seq_along(results), function(i) {
    grp <- group_values[i]
    prs <- results[[i]]$pairwise_results
    if (is.null(prs)) return(NULL)
    do.call(rbind, lapply(names(prs), function(pair) {
      comps <- prs[[pair]]
      do.call(rbind, lapply(names(comps), function(comp) {
        extract_row(comps[[comp]], grp, pair, comp)
      }))
    }))
  }))
  land_df <- do.call(rbind, lapply(seq_along(results), function(i) {
    grp <- group_values[i]
    lprs <- results[[i]]$landscape_pairwise_results
    if (is.null(lprs)) return(NULL)
    do.call(rbind, lapply(names(lprs), function(pair) {
      stats <- lprs[[pair]][["landscape"]]
      extract_row(
        stats,
        grp, pair, "landscape"
      )
    }))
  }))
  all_pairs <- rbind(betti_df, land_df)
  csv_out <- file.path(output_folder, paste(group_by_col, "all_cross_iteration_pairwise_stats.csv"))
  write.csv(all_pairs, csv_out, row.names = FALSE)
  log_message(paste("Wrote all curve-comparison stats to:", csv_out))
  log_message("All group-specific comparisons completed.")
}

#' @title Run cross-iteration analysis
#' @description
#' A lightweight wrapper around `cross_iteration_comparison_with_betti()`
#' that performs the comparison using default settings and writes the
#' results to `results_folder`. The function is typically invoked from
#' `run_modular_analysis()`.
#'
#' @param data_iterations List of iteration objects created by the
#'   preprocessing pipeline.
#' @param results_folder Directory where output files should be written.
#' @param group_by_col Metadata column used to group cells when comparing
#'   Betti curves across iterations.
#' @param ... Additional parameters passed to
#'   `cross_iteration_comparison_with_betti()`.
#'
#' @return The result object returned by
#'   `cross_iteration_comparison_with_betti()`.
#'
#' @examples
#' \dontrun{
#' run_cross_iteration(data_iterations)
#' }
#' @export
run_cross_iteration <- function(data_iterations,
                                results_folder = "results",
                                group_by_col = "Tissue",
                                ...) {
  if (!dir.exists(results_folder)) {
    dir.create(results_folder, recursive = TRUE)
  }
  cross_iteration_comparison_with_betti(
    data_iterations = data_iterations,
    group_by_col = group_by_col,
    output_folder = file.path(results_folder, "cross_iteration_comparisons"),
    ...
  )
}

