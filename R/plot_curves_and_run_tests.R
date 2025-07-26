#orig_ident_columns = c("orig.ident")
# orig_ident_columns <- c("SRA Number", "SRS Number")
# ---------------------------
# Heatmap Generation Function
# ---------------------------
generate_heatmaps_for_all_statistics <- function(
  results_list, dataset_name, plots_folder,
  null_summary = list(euler_effect_size = NULL, euler_ks = NULL,
                      landscape_effect_size = NULL, landscape_ks = NULL),
  plot_output_dir = file.path(plots_folder, "betti_plots", dataset_name, "statistical_comparison_heatmaps")
) {
  log_message("Generating heatmaps for all statistics in pairwise comparison results.")

  ensure_directory <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }
  ensure_directory(plot_output_dir)

  all_combined_csvs <- list() # For final combined CSV across all comparisons and dimensions

  # Build null summary strings
  null_summary_str_euler_effect <- if (!is.null(null_summary$euler_effect_size)) {
    paste("Null Eff Mean:", signif(null_summary$euler_effect_size$mean, 3),
          "SD:", signif(null_summary$euler_effect_size$sd, 3),
          "Quantiles:", paste(signif(null_summary$euler_effect_size$quantiles, 3), collapse = ", "))
  } else { "" }

  null_summary_str_euler_ks <- if (!is.null(null_summary$euler_ks)) {
    paste("Null KS Mean:", signif(null_summary$euler_ks$mean, 3),
          "SD:", signif(null_summary$euler_ks$sd, 3),
          "Quantiles:", paste(signif(null_summary$euler_ks$quantiles, 3), collapse = ", "))
  } else { "" }

  null_summary_str_landscape_effect <- if (!is.null(null_summary$landscape_effect_size)) {
    paste("Null Eff Mean:", signif(null_summary$landscape_effect_size$mean, 3),
          "SD:", signif(null_summary$landscape_effect_size$sd, 3),
          "Quantiles:", paste(signif(null_summary$landscape_effect_size$quantiles, 3), collapse = ", "))
  } else { "" }

  null_summary_str_landscape_ks <- if (!is.null(null_summary$landscape_ks)) {
    paste("Null KS Mean:", signif(null_summary$landscape_ks$mean, 3),
          "SD:", signif(null_summary$landscape_ks$sd, 3),
          "Quantiles:", paste(signif(null_summary$landscape_ks$quantiles, 3), collapse = ", "))
  } else { "" }

  for (comparison_type in names(results_list)) {
    log_message(paste("Processing heatmaps for", comparison_type))

    if (is.null(results_list[[comparison_type]]) || !is.list(results_list[[comparison_type]])) {
      log_message(paste("Invalid or missing pairwise results for", comparison_type, "- Skipping heatmaps."))
      next
    }

    heatmap_list <- list()

    for (dimension in c("dimension_0", "dimension_1", "euler", "landscape")) {
      log_message(paste("Processing dimension:", dimension, "for", comparison_type))

      pairwise_stats_df <- list() # Collect stats for wide CSV

      current_results <- if (tolower(dimension) == "landscape") {
        results_list[[comparison_type]]$landscape_pairwise_results
      } else {
        results_list[[comparison_type]]$pairwise_results
      }

      if (is.null(current_results)) {
        log_message(paste("No pairwise results available for", dimension, "in", comparison_type))
        next
      }

      available_statistics <- if (tolower(dimension) == "landscape") {
        unique(unlist(lapply(current_results, function(pair) {
          if (!is.null(pair[["landscape"]])) names(pair[["landscape"]])
        })))
      } else {
        unique(unlist(lapply(current_results, function(pair) {
          if (!is.null(pair[[dimension]])) names(pair[[dimension]])
        })))
      }

      for (stat in available_statistics) {
        log_message(paste("Generating heatmap for statistic:", stat, "in", dimension, "for", comparison_type))

        fill_label <- stat
        subtitle_text <- NULL
        if (tolower(dimension) %in% c("euler", "landscape")) {
          if (tolower(stat) == "effect_size") {
            subtitle_text <- if (dimension == "euler") null_summary_str_euler_effect else null_summary_str_landscape_effect
          } else if (tolower(stat) %in% c("ks_result", "ks_stat")) {
            fill_label <- "-log10(p-value)"
            subtitle_text <- if (dimension == "euler") null_summary_str_euler_ks else null_summary_str_landscape_ks
          } else if (tolower(stat) == "perm_p") {
            fill_label <- "-log10(perm p-value)"
            subtitle_text <- "Fill = -log10(perm p-value)."
          }
        } else {
          if (tolower(stat) %in% c("ks_result", "ks_stat")) {
            fill_label <- "-log10(p-value)"
            subtitle_text <- "Fill = -log10(p-value); Block label = KS statistic."
          } else if (tolower(stat) == "perm_p") {
            fill_label <- "-log10(perm p-value)"
            subtitle_text <- "Fill = -log10(perm p-value)."
          }
        }

        pairwise_data <- tryCatch({
          do.call(rbind, lapply(names(current_results), function(pair) {
            result <- current_results[[pair]]
            group_names <- strsplit(pair, "_vs_")[[1]]
            value <- NULL
            label <- NA

            if (tolower(dimension) == "landscape") {
              value <- result[["landscape"]][[stat]]
            } else {
              value <- result[[dimension]][[stat]]
            }

            if (tolower(stat) %in% c("ks_result", "ks_stat")) {
              if (is.list(value)) {
                p_val <- value[["combined_p"]];
                if (p_val == 0) p_val <- 1e-16
                value <- -log10(p_val)
                label <- round(as.numeric(value[["combined_stat"]])[1], 3)
              } else {
                value <- -log10(ifelse(value == 0, 1e-16, value))
              }
            } else if (tolower(stat) == "perm_p") {
              if (is.list(value) && !is.null(value[["perm_p"]])) {
                value <- value[["perm_p"]]
              }
              value <- -log10(ifelse(value == 0, 1e-16, value))
            } else if (inherits(value, "htest")) {
              value <- value$p.value
            }

            data.frame(Group1 = group_names[1], Group2 = group_names[2], Value = as.numeric(value), Label = as.numeric(label))
          }))
        }, error = function(e) {
          log_message(paste("Error collecting data for", stat, "in", dimension, "for", comparison_type, ":", conditionMessage(e)))
          NULL
        })

        if (is.null(pairwise_data) || nrow(pairwise_data) == 0) {
          log_message(paste("No valid data for statistic", stat, "in", dimension, "for", comparison_type, "- Skipping heatmap."))
          next
        }

        pairwise_stat <- pairwise_data[, c("Group1", "Group2", "Value")]
        colnames(pairwise_stat)[3] <- stat
        pairwise_stats_df[[stat]] <- pairwise_stat

        groups <- unique(c(pairwise_data$Group1, pairwise_data$Group2))
        heatmap_matrix <- matrix(NA, length(groups), length(groups), dimnames = list(groups, groups))
        label_matrix <- matrix(NA, length(groups), length(groups), dimnames = list(groups, groups))
        for (i in seq_len(nrow(pairwise_data))) {
          g1 <- pairwise_data$Group1[i];
          g2 <- pairwise_data$Group2[i]
          heatmap_matrix[g1, g2] <- pairwise_data$Value[i]
          heatmap_matrix[g2, g1] <- pairwise_data$Value[i]
          if (!is.na(pairwise_data$Label[i])) {
            label_matrix[g1, g2] <- pairwise_data$Label[i]
            label_matrix[g2, g1] <- pairwise_data$Label[i]
          }
        }

        heatmap_data <- as.data.frame(as.table(heatmap_matrix))
        colnames(heatmap_data) <- c("Group1", "Group2", "Value")

        if (tolower(stat) %in% c("ks_result", "ks_stat", "perm_p")) {
          label_data <- as.data.frame(as.table(label_matrix))
          colnames(label_data) <- c("Group1", "Group2", "Label")
          heatmap_data <- merge(heatmap_data, label_data, by = c("Group1", "Group2"), all.x = TRUE)
        }

        fill_scale <- scale_fill_gradient(low = "white", high = "blue", na.value = "grey50")

        heatmap_plot <- tryCatch({
          p <- ggplot(heatmap_data, aes(x = Group1, y = Group2, fill = Value)) +
            geom_tile(color = "white") +
            fill_scale +
            labs(
              title = paste("Pairwise Heatmap for", comparison_type, "-", dimension, "-", stat),
              subtitle = subtitle_text,
              x = "Group 1", y = "Group 2", fill = fill_label
            ) +
            theme_minimal(base_size = 14) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
              axis.text.y = element_text(size = 12),
              plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14, face = "bold")
            )
          
          if ("Label" %in% names(heatmap_data)) {
            p <- p + geom_text(aes(label = Label), color = "black", size = 4.5, na.rm = TRUE)
          }
          p
        }, error = function(e) {
          log_message(paste("Error generating plot for", stat, "in", dimension, "for", comparison_type, ":", conditionMessage(e)))
          NULL
        })

        if (!is.null(heatmap_plot)) {
          heatmap_list[[paste(dimension, stat, sep = "_")]] <- heatmap_plot
        }
      }

      if (length(pairwise_stats_df) > 0) {
        csv_output_dir <- file.path(plot_output_dir, "statistical_csvs", comparison_type)
        ensure_directory(csv_output_dir)
        wide_df <- Reduce(function(x, y) merge(x, y, by = c("Group1", "Group2"), all = TRUE), pairwise_stats_df)
        csv_file <- file.path(csv_output_dir, paste0(dataset_name, "_", comparison_type, "_", dimension, "_all_statistics.csv"))
        write.csv(wide_df, csv_file, row.names = FALSE)
        log_message(paste("Saved wide-format CSV for", dimension, "in", comparison_type, "to", csv_file))

        wide_df$ComparisonType <- comparison_type
        wide_df$Dimension <- dimension
        all_combined_csvs[[paste(comparison_type, dimension, sep = "_")]] <- wide_df
      }
    }

    if (length(heatmap_list) > 0) {
      combined <- tryCatch({
        gridExtra::marrangeGrob(grobs = heatmap_list, nrow = 2, ncol = 2)
      }, error = function(e) {
        log_message(paste("Error combining heatmaps for", comparison_type, ":", conditionMessage(e)))
        NULL
      })
      if (!is.null(combined)) {
        combined_file <- file.path(plot_output_dir, paste0(dataset_name, "_", comparison_type, "_combined_heatmaps.pdf"))
        ggsave(combined_file, combined, width = 16, height = 12)
        log_message(paste("Saved combined heatmaps for", comparison_type, "to", combined_file))
      }
    }

    log_message(paste("Completed heatmap generation for", comparison_type))
  }

  if (length(all_combined_csvs) > 0) {
    all_stats_df <- do.call(rbind, all_combined_csvs)
    master_file <- file.path(plot_output_dir, paste0(dataset_name, "_all_pairwise_statistics_combined.csv"))
    write.csv(all_stats_df, master_file, row.names = FALSE)
    log_message(paste("Saved final combined CSV for all stats to", master_file))
  }

  log_message("Completed heatmap generation for all statistics.")
  return(all_stats_df)
}

# Functions are available once the package is loaded

# Prepare a list to hold each iteration’s results
meta_list <- list()

for (iter in data_iterations) {
  dataset_name <- iter$name
  pd_list <- readRDS(iter$pd_list)
  landscape_list <- readRDS(iter$landscape_list)
  seurat_obj <- iter$seurat_obj

  results_list <- list()

  # Main comparisons
  if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
    results_list$tissue_comparison <- tryCatch(
        compute_and_compare_betti_curves(
          pd_list = pd_list,
          landscape_list = landscape_list,
          seurat_obj = list(seurat_obj),
          group_by_col = "Tissue",
          grid_points = 500,
          num_permutations = 1000,
          dataset_name = dataset_name,
          comparison_type = "group",
          results_folder = results_folder
        ),
        error = function(e) {
          cat("Error in Tissue-based comparison:", conditionMessage(e), "\n")
          NULL
        }
      )
  } else {
    cat("Warning: 'Tissue' column not found in seurat_obj metadata. Skipping tissue-based comparison.\n")
  }

  if ("SRA" %in% colnames(seurat_obj@meta.data)) {
    results_list$sra_comparison <- tryCatch(
        compute_and_compare_betti_curves(
          pd_list = pd_list,
          landscape_list = landscape_list,
          seurat_obj = list(seurat_obj),
          group_by_col = "SRA",
          grid_points = 500,
          num_permutations = 1000,
          dataset_name = dataset_name,
          comparison_type = "group",
          results_folder = results_folder
        ),
        error = function(e) {
          cat("Error in SRA-based comparison:", conditionMessage(e), "\n")
          NULL
        }
      )
  } else {
    cat("Warning: 'SRA' column not found in seurat_obj metadata. Skipping SRA-based comparison.\n")
  }

  if ("Approach" %in% colnames(seurat_obj@meta.data)) {
    results_list$approach_comparison <- tryCatch(
        compute_and_compare_betti_curves(
          pd_list = pd_list,
          landscape_list = landscape_list,
          seurat_obj = list(seurat_obj),
          group_by_col = "Approach",
          grid_points = 500,
          num_permutations = 1000,
          dataset_name = dataset_name,
          comparison_type = "group",
          results_folder = results_folder
        ),
        error = function(e) {
          cat("Error in Approach-based comparison:", conditionMessage(e), "\n")
          NULL
        }
      )
  } else {
    cat("Warning: 'Approach' column not found in seurat_obj metadata. Skipping Approach-based comparison.\n")
  }

  if ("sample" %in% colnames(seurat_obj@meta.data)) {
    results_list$sample_comparison <- tryCatch(
        compute_and_compare_betti_curves(
          pd_list = pd_list,
          landscape_list = landscape_list,
          seurat_obj = list(seurat_obj),
          group_by_col = "sample",
          grid_points = 500,
          num_permutations = 1000,
          dataset_name = dataset_name,
          comparison_type = "group",
          results_folder = results_folder
        ),
        error = function(e) {
          cat("Error in sample-based comparison:", conditionMessage(e), "\n")
          NULL
        }
      )
  } else {
    cat("Warning: 'sample' column not found in seurat_obj metadata. Skipping sample-based comparison.\n")
  }


  # Loop over bootstrapped random group columns
  random_group_cols <- grep("^Random_Group", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
  if (length(random_group_cols) > 0) {
    for (rg in random_group_cols[[1]]) {
      results_list[[paste0("random_group_comparison_", rg)]] <- tryCatch(
          compute_and_compare_betti_curves(
            pd_list = pd_list,
            landscape_list = landscape_list,
            seurat_obj = list(seurat_obj),
            group_by_col = rg,
            grid_points = 500,
            num_permutations = 1000,
            dataset_name = dataset_name,
            comparison_type = "group",
            results_folder = results_folder
          ),
          error = function(e) {
            cat("Error in", rg, "based comparison:", conditionMessage(e), "\n")
            NULL
          }
        )
    }
  } else {
    cat("Warning: No 'Random_Group' column found in seurat_obj metadata. Skipping Random_Group-based comparison.\n")
  }

  # cluster_cols <- grep("hierarchical", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
  # if (length(cluster_cols) > 0) {
  #   for (cc in cluster_cols) {
  #     results_list[[paste0("cluster_comparison_", cc)]] <- tryCatch(
  #         compute_and_compare_betti_curves(
  #           pd_list = pd_list,
  #           landscape_list = landscape_list,
  #           seurat_obj = list(seurat_obj),
  #           group_by_col = cc,
  #           grid_points = 500,
  #           num_permutations = 10000,
  #           dataset_name = dataset_name,
  #           comparison_type = "group",
  #           results_folder = results_folder
  #         ),
  #         error = function(e) {
  #           cat("Error in", rg, "based comparison:", conditionMessage(e), "\n")
  #           NULL
  #         }
  #       )
  #   }
  # } else {
  #   cat("Warning: No 'cluster' column found in seurat_obj metadata. Skipping cluster-based comparison.\n")
  # }

  # cluster_cols <- grep("kmeans", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
  # if (length(cluster_cols) > 0) {
  #   for (cc in cluster_cols) {
  #     results_list[[paste0("cluster_comparison_", cc)]] <- tryCatch(
  #         compute_and_compare_betti_curves(
  #           pd_list = pd_list,
  #           landscape_list = landscape_list,
  #           seurat_obj = list(seurat_obj),
  #           group_by_col = cc,
  #           grid_points = 500,
  #           num_permutations = 10000,
  #           dataset_name = dataset_name,
  #           comparison_type = "group",
  #           results_folder = results_folder
  #         ),
  #         error = function(e) {
  #           cat("Error in", rg, "based comparison:", conditionMessage(e), "\n")
  #           NULL
  #         }
  #       )
  #   }
  # } else {
  #   cat("Warning: No 'cluster' column found in seurat_obj metadata. Skipping cluster-based comparison.\n")
  # }

  # --- Aggregate bootstrap statistics from random group comparisons using Euler metrics ---
  random_group_keys <- grep("^random_group_comparison", names(results_list), value = TRUE, ignore.case = TRUE)
  euler_null_stats_effect <- c()
  euler_null_stats_ks <- c()
  landscape_null_stats_effect <- c()
  landscape_null_stats_ks <- c()

  for (key in random_group_keys) {
    rg_result <- results_list[[key]]$pairwise_results
    rg_result_landscape <- results_list[[key]]$landscape_pairwise_results
    if (!is.null(rg_result)) {
      # Loop over each pairwise comparison in this random group result
      for (comp_name in names(rg_result)) {
        comp <- rg_result[[comp_name]]
        comp_landscape <- rg_result_landscape[[comp_name]]
        # Use the Euler curve result since it captures all dimensions
        if (is.list(comp) && !is.null(comp$euler)) {
          # Aggregate effect size if available
          if (!is.null(comp$euler$effect_size)) {
            euler_null_stats_effect <- c(euler_null_stats_effect, comp$euler$effect_size)
          }
          # Aggregate Kolmogorov–Smirnov statistic if available
          if (!is.null(comp$euler)) {
            euler_null_stats_ks <- c(euler_null_stats_ks, comp$euler$ks_stat)
            # euler_null_stats_ks <- c(euler_null_stats_ks, comp$euler$ks_result$combined_stat)

          }
        }

        # Check if the landscape comparisons are stored under "dimension_landscape"
        if (!is.null(comp_landscape$landscape)) {
          # Extract effect size if available
          if (!is.null(comp_landscape$landscape$effect_size)) {
            landscape_null_stats_effect <- c(landscape_null_stats_effect, comp_landscape$landscape$effect_size)
          }
          # Extract KS statistic if available
          if (!is.null(comp_landscape$landscape$ks_stat)) {
            landscape_null_stats_ks <- c(landscape_null_stats_ks, comp_landscape$landscape$ks_stat)
          }
        }

      }
    }
  }

  #  # Remove NA values for each metric
  #   euler_null_stats_effect <- euler_null_stats_effect[!is.na(euler_null_stats_effect)]
  #   euler_null_stats_ks <- euler_null_stats_ks[!is.na(euler_null_stats_ks)]
  # Remove NA values
  landscape_null_stats_effect <- landscape_null_stats_effect[!is.na(landscape_null_stats_effect)]
  landscape_null_stats_ks <- landscape_null_stats_ks[!is.na(landscape_null_stats_ks)]


  # Compute summary statistics for effect size
  if (length(euler_null_stats_effect) > 0) {
    euler_null_stats_effect <- mean(euler_null_stats_effect)
    euler_null_sd_effect <- sd(euler_null_stats_effect)
    euler_null_quantiles_effect <- quantile(euler_null_stats_effect, probs = c(0.025, 0.5, 0.975))
    cat("Random Group Null Distribution (Effect Size) Mean:", euler_null_stats_effect, "\n")
    cat("Random Group Null Distribution (Effect Size) SD:", euler_null_sd_effect, "\n")
    cat("Random Group Null Distribution (Effect Size) Quantiles:", euler_null_quantiles_effect, "\n")

    euler_random_group_null_effect <- list(mean = euler_null_stats_effect, sd = euler_null_sd_effect, quantiles = euler_null_quantiles_effect)
  } else {
    cat("No valid random group effect size statistics were found to build a null distribution.\n")
    euler_random_group_null_effect <- NULL
  }

  # Compute summary statistics for KS statistic
  if (length(euler_null_stats_ks) > 0) {
    euler_null_mean_ks <- mean(euler_null_stats_ks)
    euler_null_sd_ks <- sd(euler_null_stats_ks)
    euler_null_quantiles_ks <- quantile(euler_null_stats_ks, probs = c(0.025, 0.5, 0.975))
    cat("Random Group Null Distribution (KS) Mean:", euler_null_mean_ks, "\n")
    cat("Random Group Null Distribution (KS) SD:", euler_null_sd_ks, "\n")
    cat("Random Group Null Distribution (KS) Quantiles:", euler_null_quantiles_ks, "\n")

    euler_random_group_null_ks <- list(mean = euler_null_mean_ks, sd = euler_null_sd_ks, quantiles = euler_null_quantiles_ks)
  } else {
    cat("No valid random group KS statistics were found to build a null distribution.\n")
    euler_random_group_null_ks <- NULL
  }


  # Compute summary statistics for the landscape effect size
  if (length(landscape_null_stats_effect) > 0) {
    landscape_null_mean_effect <- mean(landscape_null_stats_effect)
    landscape_null_sd_effect <- sd(landscape_null_stats_effect)
    landscape_null_quantiles_effect <- quantile(landscape_null_stats_effect, probs = c(0.025, 0.5, 0.975))
    cat("Random Group Landscape Null Distribution (Effect Size) Mean:", landscape_null_mean_effect, "\n")
    cat("Random Group Landscape Null Distribution (Effect Size) SD:", landscape_null_sd_effect, "\n")
    cat("Random Group Landscape Null Distribution (Effect Size) Quantiles:", landscape_null_quantiles_effect, "\n")

    landscape_random_group_null_effect <- list(mean = landscape_null_mean_effect,
                                                  sd = landscape_null_sd_effect,
                                                  quantiles = landscape_null_quantiles_effect)
  } else {
    cat("No valid random group landscape effect size statistics were found.\n")
    landscape_random_group_null_effect <- NULL
  }

  # Compute summary statistics for the landscape KS statistic
  if (length(landscape_null_stats_ks) > 0) {
    landscape_null_mean_ks <- mean(landscape_null_stats_ks)
    landscape_null_sd_ks <- sd(landscape_null_stats_ks)
    landscape_null_quantiles_ks <- quantile(landscape_null_stats_ks, probs = c(0.025, 0.5, 0.975))
    cat("Random Group Landscape Null Distribution (KS) Mean:", landscape_null_mean_ks, "\n")
    cat("Random Group Landscape Null Distribution (KS) SD:", landscape_null_sd_ks, "\n")
    cat("Random Group Landscape Null Distribution (KS) Quantiles:", landscape_null_quantiles_ks, "\n")

    landscape_random_group_null_ks <- list(mean = landscape_null_mean_ks,
                                              sd = landscape_null_sd_ks,
                                              quantiles = landscape_null_quantiles_ks)
  } else {
    cat("No valid random group landscape KS statistics were found.\n")
    landscape_random_group_null_ks <- NULL
  }

  # Combine the null stats from Euler (if already computed) with the landscape null stats.
  # Here we assume you have previously computed euler_random_group_null_effect and euler_random_group_null_ks.
  overall_null_summary <- list(
      euler_effect_size = euler_random_group_null_effect, # from your Euler extraction
      euler_ks = euler_random_group_null_ks, # from your Euler extraction
      landscape_effect_size = landscape_random_group_null_effect,
      landscape_ks = landscape_random_group_null_ks
    )


  # Generate heatmaps for all statistics, passing the null summaries as a list
  tryCatch(
      df_iter <- generate_heatmaps_for_all_statistics(
        results_list = results_list,
        dataset_name = dataset_name,
        plots_folder = plots_folder,
        null_summary = overall_null_summary
      ),
      error = function(e) {
        cat("Error in generating heatmaps:", conditionMessage(e), "\n")
      }
    )

  df_iter <- df_iter %>%
       mutate(Dataset = dataset_name)

  meta_list[[dataset_name]] <- df_iter

}

# bind all rows together
meta_master_df <- bind_rows(meta_list)

# write to disk
write.csv(
  meta_master_df,
  file.path(plots_folder, "betti_plots", "ALL_datasets_all_pairwise_statistics_meta_master.csv"),
  row.names = FALSE
)