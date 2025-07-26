# Make sure you have cowplot and other required packages installed

# Helper function for logging
log_message <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

# ----------------------------------------------------
# Betti/Euler/Landscape Curve Comparison Function
# ----------------------------------------------------

# ----------------------------------------------------
# Betti/Euler/Landscape Curve Comparison Function
# ----------------------------------------------------
compute_and_compare_betti_curves <- function(pd_list, landscape_list, seurat_objects = NULL, group_by_col = NULL,
                                             dimensions = c(0, 1), grid_points = 500, num_permutations = 10000,
                                             bootstrap_samples = 100, dataset_name = NULL, comparison_type = "group",
                                             base_sigma = 1, num_cores = 16, results_folder = 'results') {

  # Simple logging function.
  log_message <- function(msg) {
    cat(sprintf("[%s] %s\n", Sys.time(), msg))
  }

  log_debug_to_file <- function(msg, file_path) {
    formatted_msg <- sprintf("[%s] %s\n", Sys.time(), msg)
    cat(formatted_msg, file = file_path, append = TRUE)
  }
  debug_log_file <- file.path(results_folder, "landscape_debug_run.txt")

  # --- Helper: Directory Creation ---
  ensure_directory <- function(path) {
    if (length(path) == 0 || is.na(path) || path == "") {
      stop("Invalid directory path provided.")
    }
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }

  # --- Helper: Normalization with Division-by-Zero Check ---
  normalize_curve <- function(curve, grid_points) {
    tot <- sum(curve)
    if (tot == 0) rep(0, grid_points) else curve / tot
  }

  # --- Curve Caching ---
  generate_cache_filename <- function(dataset_name, hash, results_folder) {
    path <- file.path(results_folder, "plots", "betti_plots", "betti_cache", dataset_name, paste0("cache_", hash, ".rds"))
    ensure_directory(dirname(path))
    path
  }

  # --- Compute Betti Curve ---
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

  # --- UNIFIED NULL STATISTICS (for Betti/Euler) ---
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
      if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
    }
    if (is.null(metadata_list) || length(metadata_list) == 0 || !all(sapply(metadata_list, is.data.frame))) {
        log_message("Warning: metadata_list is invalid or empty. Cannot compute metadata-driven null stats. Returning NULL.")
        return(NULL)
    }
    combined_meta <- do.call(rbind, metadata_list)
    bs_cols <- grep("^random_group_bootstrap_", colnames(combined_meta),
                    value = TRUE, ignore.case = TRUE)
    if (length(bs_cols) == 0) {
        log_message("Warning: No 'random_group_bootstrap_' columns found in metadata. Cannot compute null stats. Returning NULL.")
        return(NULL)
    }
    
    results <- mclapply(bs_cols, function(col) {
      if (!"orig.ident" %in% colnames(combined_meta)) {
          log_message("Error: 'orig.ident' not found in combined metadata. Cannot match PDs to metadata for null stats.")
          return(NULL)
      }
      meta_idents <- combined_meta[["orig.ident"]]
      pd_list_names <- names(pd_list)
      
      pd_indices_in_meta <- match(pd_list_names, meta_idents)
      pd_indices_in_meta <- pd_indices_in_meta[!is.na(pd_indices_in_meta)]
      
      if(length(pd_indices_in_meta) == 0) {
        log_message("Warning: No overlap between pd_list names and metadata 'orig.ident'.")
        return(NULL)
      }
      
      random_groups_for_pds <- combined_meta[[col]][pd_indices_in_meta]

      idx1_pd_names <- pd_list_names[which(random_groups_for_pds == 1)]
      idx2_pd_names <- pd_list_names[which(random_groups_for_pds == 2)]
      
      if (length(idx1_pd_names) == 0 || length(idx2_pd_names) == 0) return(NULL)
      
      pd_grp1 <- pd_list[idx1_pd_names]
      pd_grp2 <- pd_list[idx2_pd_names]
      
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
    ensure_directory(dirname(cache_file))
    saveRDS(curve, cache_file)
    log_message(paste("Saved cached curve to:", cache_file))
  }

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

  bootstrap_curve <- function(pd_subset, dimension, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, seed = 42, results_folder, num_cores = 8) {
    log_message(paste("Starting bootstrapping for dimension", dimension, "with", n_bootstrap, "samples."))
    individual_pd_curves <- mclapply(seq_along(pd_subset), function(i) {
      pd <- pd_subset[[i]]
      if (is.matrix(pd) && nrow(pd) > 0) {
        curv <- get_or_compute_curve(pd, dimension, dataset_name, base_sigma, grid_points, tau_max, results_folder)
        normalize_curve(curv, grid_points)
      } else {
        rep(0, grid_points)
      }
    }, mc.cores = num_cores)
    names(individual_pd_curves) <- paste0("pd_", seq_along(pd_subset))
    bootstrap_curves <- mclapply(1:n_bootstrap, function(rep) {
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

  integrated_diff <- function(curve_diff, tau_vals) {
    if (length(curve_diff) < 2) return(0)
    delta_tau <- diff(tau_vals)
    sqrt(sum(delta_tau * (curve_diff[-length(curve_diff)] ^ 2 + curve_diff[-1] ^ 2) / 2, na.rm = TRUE))
  }

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
        x = seq(0, 1, length.out = length(aggregated_curve)),
        y = aggregated_curve,
        xout = grid,
        rule = 2
      )$y
    }
    aggregated_curve
  }
  
  summary_landscape <- function(landscape_obj, grid = seq(0, 1, length.out = grid_points)) {
    curve <- compute_landscape_curve(landscape_obj, grid = grid)
    delta <- diff(grid)
    auc <- sum(((curve[-length(curve)] + curve[-1]) / 2) * delta)
    auc
  }

  compute_null_stats_landscapes <- function(landscape_list, n_bootstrap = 50, grid_points = 500) {
    if (length(landscape_list) < 2) {
        log_message("Warning: Fewer than 2 landscapes provided. Cannot compute null stats.")
        return(NULL)
    }
    summary_values <- sapply(landscape_list, function(land) {
        summary_landscape(land, grid = seq(0, 1, length.out = grid_points))
    })
    boot_effects <- numeric(n_bootstrap)
    boot_ks <- numeric(n_bootstrap)
    N <- length(summary_values)
    for (rep in 1:n_bootstrap) {
        set.seed(100 + rep)
        indices <- sample(N, replace = TRUE)
        half <- floor(length(indices) / 2)
        if (half == 0 || half == N) next 
        
        group1_vals <- summary_values[indices[1:half]]
        group2_vals <- summary_values[indices[(half + 1):length(indices)]]
        
        diff_val <- abs(mean(group1_vals) - mean(group2_vals))
        sd_pooled <- sd(c(group1_vals, group2_vals), na.rm = TRUE)
        boot_effects[rep] <- ifelse(sd_pooled == 0, 0, diff_val / sd_pooled)
        
        ks_result <- tryCatch(ks.test(group1_vals, group2_vals), error = function(e) NULL)
        boot_ks[rep] <- if (!is.null(ks_result)) as.numeric(ks_result$statistic) else NA
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

  bootstrap_aggregated_landscape_curve <- function(individual_curves, grid_points, n_bootstrap = 100) {
    boot_samples <- replicate(n_bootstrap, {
      sampled <- sample(individual_curves, replace = TRUE)
      rowMeans(do.call(cbind, sampled), na.rm = TRUE)
    })
    mean_curve <- rowMeans(boot_samples, na.rm = TRUE)
    lower_curve <- apply(boot_samples, 1, quantile, probs = 0.025, na.rm = TRUE)
    upper_curve <- apply(boot_samples, 1, quantile, probs = 0.975, na.rm = TRUE)
    log_debug_to_file("Checking bootstrap results...", debug_log_file)
    log_debug_to_file(paste("Are the mean and upper curves identical?", all.equal(mean_curve, upper_curve)), debug_log_file)
    list(
      mean = mean_curve,
      lower = lower_curve,
      upper = upper_curve,
      bootstrap_curves = boot_samples
    )
  }

  # --- UPDATED FUNCTION (from Canvas) ---
  analyze_and_plot_curves <- function(curve1, curve2, grp1, grp2, curve_type, tau_vals, tau_max,
                                    alpha, effect_threshold, num_permutations,
                                    null_effect_stats = NULL, null_ks_stats = NULL,
                                    verbose = TRUE) {
    log_message <- function(msg) {
      if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
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
    
    p <- ggplot(plot_data, aes(x = Tau / tau_max, y = Mean, color = Group, group = Group)) +
      geom_line(linewidth = 1.2) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group), alpha = 0.2) +
      labs(title = paste(curve_type, "Comparison:", grp1, "vs.", grp2),
           x = "Normalized Filtration Scale (Tau)",
           y = paste(curve_type, "Value")) +
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
        labels = function(breaks) scales::percent(breaks / tau_max, accuracy = 0.0000001)
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
            ensure_directory(dirname(png_file_path))
            ggsave(filename = png_file_path, plot = p, width = 12, height = 10, dpi = 300)
            log_message(paste("Saved individual PNG plot to", png_file_path))
          }
          if (is.null(dimensions)) {
            file_name_pdf <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, "_", sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.pdf")
            file_name_svg <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, "_", sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.svg")
          } else {
            file_name_pdf <- paste0(dataset_name, "_", group_by_col, "_", sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.pdf")
            file_name_svg <- paste0(dataset_name, "_", group_by_col, "_", sanitize_group(grp1), "_vs_", sanitize_group(grp2), "_comparison.svg")
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

  # --- UPDATED FUNCTION (from Canvas) ---
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
      
      ggplot(data, aes(x = Tau / tau_max, y = Mean, color = Group, group = Group)) +
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
        scale_x_log10(
          "Normalized Filtration Scale (Tau, log scale)",
          labels = function(breaks) scales::percent(breaks / tau_max, accuracy = 0.0000001)
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
        ensure_directory(dirname(png_file_path))
        ggsave(filename = png_file_path, plot = p, width = 12, height = 10, dpi = 300)
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

  # --- Main Execution Logic ---
  log_message("Starting computation.")
  tau_max <- max(unlist(lapply(pd_list, function(pd) if(is.matrix(pd) && ncol(pd) >= 3) max(pd[, 3], na.rm = TRUE) else -Inf)), na.rm = TRUE)
  if (is.infinite(tau_max)) {
      log_message("Warning: Could not determine a valid tau_max. Defaulting to 1.")
      tau_max <- 1.0
  }
  tau_vals <- seq(0.0001, tau_max, length.out = grid_points)
  
  if (is.null(seurat_objects) || length(seurat_objects) == 0) stop("Seurat objects list is required for metadata.")
  
  group_metadata_mapping <- do.call(rbind, lapply(seurat_objects, function(obj) {
    meta <- obj@meta.data
    unique(meta[, c("orig.ident", group_by_col)])
  }))
  if (nrow(group_metadata_mapping) == 0) stop("No valid metadata mapping found.")
  groups <- unique(group_metadata_mapping[[group_by_col]])
  log_message(paste("Groups identified:", paste(groups, collapse = ", ")))
  
  metadata_list <- lapply(seurat_objects, function(obj) obj@meta.data)
  
  log_message("Computing null statistics for Betti/Euler curves...")
  bootstrap_nulls <- compute_null_stats(pd_list, metadata_list, tau_vals, dimensions, "Euler",
                                        dataset_name, base_sigma, grid_points, tau_max, results_folder, verbose = TRUE)
  log_message("Completed null statistics for Betti/Euler curves.")
  
  log_message("Computing null statistics for Landscape curves...")
  landscape_nulls <- compute_null_stats_landscapes(landscape_list, n_bootstrap = 50, grid_points = grid_points)
  log_message("Completed null statistics for Landscape curves.")

  log_message("Aggregating and analyzing Landscape curves...")
  aggregated_landscape_curves_by_group <- setNames(lapply(as.character(groups), function(group) {
    group_ids <- group_metadata_mapping$orig.ident[group_metadata_mapping[[group_by_col]] == group]
    grid <- seq(0, 1, length.out = grid_points)
    if (length(group_ids) == 0) {
      return(list(
        mean = rep(0, grid_points), lower = rep(0, grid_points), upper = rep(0, grid_points),
        individual_landscape_curves = list()
      ))
    }
    individual_curves <- lapply(landscape_list[group_ids], function(land) {
      compute_landscape_curve(land, grid = grid)
    })
    bootstrapped_agg <- bootstrap_aggregated_landscape_curve(individual_curves, grid_points, n_bootstrap = bootstrap_samples)
    list(
      mean = bootstrapped_agg$mean, lower = bootstrapped_agg$lower, upper = bootstrapped_agg$upper,
      individual_landscape_curves = individual_curves,
      bootstrap_aggregates = bootstrapped_agg$bootstrap_curves
    )
  }), groups)

  # --- MODIFICATION: Updated function call ---
  landscape_pairwise_results <- perform_pairwise_comparisons(aggregated_landscape_curves_by_group,
                                                             NULL,
                                                             groups, 
                                                             seq(0, 1, length.out=grid_points), 
                                                             tau_max = 1.0, # Landscapes are on a 0-1 scale
                                                             dimensions = NULL,
                                                             dataset_name = dataset_name,
                                                             group_by_col = group_by_col,
                                                             plot_output_dir = file.path(results_folder, "plots", "landscape_curves", dataset_name, "pairwise_comparisons"),
                                                             null_effect_stats = landscape_nulls$null_effect,
                                                             null_ks_stats = landscape_nulls$null_ks)
  
  # --- MODIFICATION: Updated function call ---
  plot_aggregated_curves(aggregated_landscape_curves_by_group, 
                         groups, 
                         dimensions = NULL, 
                         seq(0, 1, length.out=grid_points),
                         tau_max = 1.0, # Landscapes are on a 0-1 scale
                         dataset_name, 
                         group_by_col, 
                         euler_curves = NULL, 
                         results_folder)
  
  log_message("Aggregating and analyzing Betti/Euler curves...")
  group_results <- mclapply(as.character(groups), function(group) {
    log_message(paste("Processing group:", group))
    group_pd_ids <- group_metadata_mapping$orig.ident[group_metadata_mapping[[group_by_col]] == group]
    group_pd <- pd_list[names(pd_list) %in% group_pd_ids]
    if (length(group_pd) == 0) {
      log_message(paste("No PDs for group:", group))
      return(NULL)
    }
    tryCatch({
      bc <- lapply(dimensions, function(d) {
        log_message(paste("Bootstrapping Betti for dim", d, "in group", group))
        bootstrap_curve(group_pd, d, group, dataset_name, grid_points, tau_max, base_sigma, bootstrap_samples, results_folder = results_folder)
      })
      bc <- list(betti = bc)
      ec <- compute_euler_curve(group_pd, dimensions, group, dataset_name, grid_points, tau_max, base_sigma, bootstrap_samples, results_folder)
      ec <- list(euler = ec)
      cons <- check_betti_euler_consistency(ec$euler, bc$betti, dimensions, tau_vals)
      list(group = group, betti_curves = setNames(bc, group), euler_curve = setNames(ec, group), consistency = cons)
    }, error = function(e) {
      log_message(paste("Error in group", group, ":", e$message))
      NULL
    })
  }, mc.cores = num_cores)
  
  group_results <- Filter(Negate(is.null), group_results)
  if (length(group_results) == 0) { log_message("No results computed for Betti/Euler."); return(NULL) }
  
  betti_curves_by_group <- setNames(lapply(group_results, function(res) {
    res$betti_curves[[res$group]]
  }), vapply(group_results, `[[`, "", "group"))
  
  euler_curves_by_group <- setNames(lapply(group_results, function(res) {
    res$euler_curve[[res$group]]
  }), vapply(group_results, `[[`, "", "group"))
  
  # --- MODIFICATION: Updated function call ---
  pairwise_results <- perform_pairwise_comparisons(betti_curves_by_group, 
                                                   euler_curves_by_group, 
                                                   groups, 
                                                   tau_vals, 
                                                   tau_max, # Pass the calculated tau_max
                                                   dimensions,
                                                   dataset_name, 
                                                   group_by_col,
                                                   plot_output_dir = file.path(results_folder, "plots", "betti_plots", dataset_name, "pairwise_comparisons"),
                                                   null_effect_stats = bootstrap_nulls$null_effect, 
                                                   null_ks_stats = bootstrap_nulls$null_ks)
                                                   
  # --- MODIFICATION: Updated function call ---
  plot_aggregated_curves(betti_curves_by_group, 
                         groups, 
                         dimensions, 
                         tau_vals, 
                         tau_max, # Pass the calculated tau_max
                         dataset_name, 
                         group_by_col, 
                         euler_curves_by_group, 
                         results_folder)
  
  list(
    group_results = group_results,
    pairwise_results = pairwise_results,
    aggregated_landscape_curves = aggregated_landscape_curves_by_group,
    landscape_pairwise_results = landscape_pairwise_results
  )
}
