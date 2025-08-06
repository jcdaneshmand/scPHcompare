# ---------------------------------------------------------------
# Heatmap Generation Helper
# ---------------------------------------------------------------
generate_heatmaps_for_all_statistics <- function(
  results_list, dataset_name, plots_folder,
  null_summary = list(euler_effect_size = NULL, euler_ks = NULL,
                      landscape_effect_size = NULL, landscape_ks = NULL),
  plot_output_dir = file.path(plots_folder, "betti_plots", dataset_name, "statistical_comparison_heatmaps"),
  verbose = TRUE
) {
  log_message <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", Sys.time(), msg))
  }
  log_message("Generating heatmaps for all statistics in pairwise comparison results.")

  ensure_directory <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }
  ensure_directory(plot_output_dir)

  all_combined_csvs <- list()

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

      pairwise_stats_df <- list()

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
                p_val <- value[["combined_p"]]
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
          g1 <- pairwise_data$Group1[i]
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

        fill_scale <- ggplot2::scale_fill_gradient(low = "white", high = "blue", na.value = "grey50")

        heatmap_plot <- tryCatch({
          p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = Group1, y = Group2, fill = Value)) +
            ggplot2::geom_tile(color = "white") +
            fill_scale +
            ggplot2::labs(
              title = paste("Pairwise Heatmap for", comparison_type, "-", dimension, "-", stat),
              subtitle = subtitle_text,
              x = "Group 1", y = "Group 2", fill = fill_label
            ) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
              axis.text.y = ggplot2::element_text(size = 12),
              plot.title = ggplot2::element_text(size = 14, face = "bold"),
              plot.subtitle = ggplot2::element_text(size = 12),
              legend.text = ggplot2::element_text(size = 12),
              legend.title = ggplot2::element_text(size = 14, face = "bold")
            )

          if ("Label" %in% names(heatmap_data)) {
            p <- p + ggplot2::geom_text(ggplot2::aes(label = Label), color = "black", size = 4.5, na.rm = TRUE)
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
        utils::write.csv(wide_df, csv_file, row.names = FALSE)
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
        ggplot2::ggsave(combined_file, combined, width = 16, height = 12)
        log_message(paste("Saved combined heatmaps for", comparison_type, "to", combined_file))
      }
    }

    log_message(paste("Completed heatmap generation for", comparison_type))
  }

  if (length(all_combined_csvs) > 0) {
    all_stats_df <- do.call(rbind, all_combined_csvs)
    master_file <- file.path(plot_output_dir, paste0(dataset_name, "_all_pairwise_statistics_combined.csv"))
    utils::write.csv(all_stats_df, master_file, row.names = FALSE)
    log_message(paste("Saved final combined CSV for all stats to", master_file))
  }

  log_message("Completed heatmap generation for all statistics.")
  return(all_stats_df)
}


#' Run post-processing modules
#'
#' This helper runs the optional post-processing analyses on the
#' output of `process_datasets_PH()`. Depending on the supplied
#' arguments, clustering comparisons, Betti curve analyses and a
#' cross-iteration comparison can be performed.
#'
#' @param ph_results The list returned by `process_datasets_PH()`.
#' @param results_dir Directory where results should be written.
#' @param run_cluster Logical, run clustering comparisons if `TRUE`.
#' @param run_betti Logical, compute Betti curves if `TRUE`.
#' @param run_cross_iteration Logical, run cross-iteration analysis.
#' @param SRA_col Metadata column containing sample identifiers.
#' @param Tissue_col Metadata column with tissue labels.
#' @param Approach_col Metadata column describing experimental approach.
#' @param ... Additional parameters passed to the downstream functions.
#'
#' @return A list containing any results generated by the chosen modules.
#'
#' @examples
#' \dontrun{
#' ph_results <- process_datasets_PH(read.csv("metadata.csv"))
#' run_modular_analysis(ph_results, run_cluster = TRUE)
#' }
#' @export
run_modular_analysis <- function(ph_results,
                                 results_dir = "results",
                                 run_cluster = FALSE,
                                 run_betti = FALSE,
                                 run_cross_iteration = TRUE,
                                 SRA_col = "SRA",
                                 Tissue_col = "Tissue",
                                 Approach_col = "Approach",
                                 ...) {
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  # Required modules are loaded with the package

  results <- list()
  data_iterations <- ph_results$data_iterations

  if (run_cluster && exists("run_cluster_comparison") && !is.null(data_iterations)) {
    results$cluster <- try(run_cluster_comparison(data_iterations,
                                                 results_folder = results_dir,
                                                 SRA_col = SRA_col,
                                                 ...),
                           silent = TRUE)
  }

  if (run_betti && !is.null(data_iterations)) {
    betti_results <- list()
    meta_list <- list()
    for (iter in data_iterations) {
      pd_list <- readRDS(iter$pd_list)
      landscape_list <- if (!is.null(iter$landscape_list)) readRDS(iter$landscape_list) else NULL

      iter_results <- list()

      if (Tissue_col %in% colnames(iter$seurat_obj@meta.data)) {
        iter_results$Tissue <- try(
          compute_and_compare_betti_curves(
            pd_list = pd_list,
            landscape_list = landscape_list,
            seurat_objects = list(iter$seurat_obj),
            group_by_col = Tissue_col,
            dataset_name = iter$name,
            results_folder = results_dir,
            ...
          ),
          silent = TRUE
        )
      }

      if (SRA_col %in% colnames(iter$seurat_obj@meta.data)) {
        iter_results$SRA <- try(
          compute_and_compare_betti_curves(
            pd_list = pd_list,
            landscape_list = landscape_list,
            seurat_objects = list(iter$seurat_obj),
            group_by_col = SRA_col,
            dataset_name = iter$name,
            results_folder = results_dir,
            ...
          ),
          silent = TRUE
        )
      }

      if (Approach_col %in% colnames(iter$seurat_obj@meta.data)) {
        iter_results$Approach <- try(
          compute_and_compare_betti_curves(
            pd_list = pd_list,
            landscape_list = landscape_list,
            seurat_objects = list(iter$seurat_obj),
            group_by_col = Approach_col,
            dataset_name = iter$name,
            results_folder = results_dir,
            ...
          ),
          silent = TRUE
        )
      }

      # if ("sample" %in% colnames(iter$seurat_obj@meta.data)) {
      #   iter_results$sample <- try(
      #     compute_and_compare_betti_curves(
      #       pd_list = pd_list,
      #       landscape_list = landscape_list,
      #       seurat_objects = list(iter$seurat_obj),
      #       group_by_col = "sample",
      #       dataset_name = iter$name,
      #       results_folder = results_dir,
      #       ...
      #     ),
      #     silent = TRUE
      #   )
      # }

      # random_group_cols <- grep("^Random_Group", colnames(iter$seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
      # if (length(random_group_cols) > 0) {
      #   for (rg in random_group_cols) {
      #     iter_results[[paste0("random_group_comparison_", rg)]] <- try(
      #       compute_and_compare_betti_curves(
      #         pd_list = pd_list,
      #         landscape_list = landscape_list,
      #         seurat_objects = list(iter$seurat_obj),
      #         group_by_col = rg,
      #         dataset_name = iter$name,
      #         results_folder = results_dir,
      #         ...
      #       ),
      #       silent = TRUE
      #     )
      #   }
      # }

      random_group_keys <- grep("^random_group_comparison", names(iter_results), value = TRUE, ignore.case = TRUE)
      euler_effect <- c()
      euler_ks <- c()
      landscape_effect <- c()
      landscape_ks <- c()
      for (key in random_group_keys) {
        rg_result <- iter_results[[key]]$pairwise_results
        rg_landscape <- iter_results[[key]]$landscape_pairwise_results
        if (!is.null(rg_result)) {
          for (comp_name in names(rg_result)) {
            comp <- rg_result[[comp_name]]
            comp_landscape <- rg_landscape[[comp_name]]
            if (is.list(comp) && !is.null(comp$euler)) {
              if (!is.null(comp$euler$effect_size)) euler_effect <- c(euler_effect, comp$euler$effect_size)
              if (!is.null(comp$euler$ks_stat)) euler_ks <- c(euler_ks, comp$euler$ks_stat)
            }
            if (!is.null(comp_landscape$landscape)) {
              if (!is.null(comp_landscape$landscape$effect_size)) landscape_effect <- c(landscape_effect, comp_landscape$landscape$effect_size)
              if (!is.null(comp_landscape$landscape$ks_stat)) landscape_ks <- c(landscape_ks, comp_landscape$landscape$ks_stat)
            }
          }
        }
      }

      euler_effect <- euler_effect[!is.na(euler_effect)]
      euler_ks <- euler_ks[!is.na(euler_ks)]
      landscape_effect <- landscape_effect[!is.na(landscape_effect)]
      landscape_ks <- landscape_ks[!is.na(landscape_ks)]

      euler_random_group_null_effect <- if (length(euler_effect) > 0) list(mean = mean(euler_effect), sd = sd(euler_effect), quantiles = stats::quantile(euler_effect, probs = c(0.025, 0.5, 0.975))) else NULL
      euler_random_group_null_ks <- if (length(euler_ks) > 0) list(mean = mean(euler_ks), sd = sd(euler_ks), quantiles = stats::quantile(euler_ks, probs = c(0.025, 0.5, 0.975))) else NULL
      landscape_random_group_null_effect <- if (length(landscape_effect) > 0) list(mean = mean(landscape_effect), sd = sd(landscape_effect), quantiles = stats::quantile(landscape_effect, probs = c(0.025, 0.5, 0.975))) else NULL
      landscape_random_group_null_ks <- if (length(landscape_ks) > 0) list(mean = mean(landscape_ks), sd = sd(landscape_ks), quantiles = stats::quantile(landscape_ks, probs = c(0.025, 0.5, 0.975))) else NULL

      overall_null_summary <- list(
        euler_effect_size = euler_random_group_null_effect,
        euler_ks = euler_random_group_null_ks,
        landscape_effect_size = landscape_random_group_null_effect,
        landscape_ks = landscape_random_group_null_ks
      )

      df_iter <- try(
        generate_heatmaps_for_all_statistics(
          results_list = iter_results,
          dataset_name = iter$name,
          plots_folder = results_dir,
          null_summary = overall_null_summary
        ),
        silent = TRUE
      )

      if (!inherits(df_iter, "try-error") && !is.null(df_iter)) {
        df_iter <- dplyr::mutate(df_iter, Dataset = iter$name)
        meta_list[[iter$name]] <- df_iter
      }

      betti_results[[iter$name]] <- iter_results
    }
    results$betti <- betti_results

    if (length(meta_list) > 0) {
      meta_master_df <- dplyr::bind_rows(meta_list)
      utils::write.csv(
        meta_master_df,
        file.path(results_dir, "betti_plots", "ALL_datasets_all_pairwise_statistics_meta_master.csv"),
        row.names = FALSE
      )
      results$betti_meta <- meta_master_df
    }
  }

  if (run_cross_iteration && exists("run_cross_iteration") && !is.null(data_iterations)) {
    results$cross_iteration <- try(run_cross_iteration(data_iterations,
                                                      results_folder = results_dir,
                                                      ...),
                                   silent = TRUE)
  }

  invisible(results)
}

#' Full post-processing pipeline
#'
#' Performs matrix calculations, clustering and all downstream analyses
#' for the persistent homology results returned by
#' `process_datasets_PH()`.
#'
#' @param ph_results Result object from `process_datasets_PH()`.
#' @param results_dir Directory where output should be saved.
#' @param num_cores Number of cores used when calculating matrices.
#' @param run_standard_seurat_clustering Logical, run Seurat clustering.
#' @param run_kmeans_clustering Logical, run k-means clustering.
#' @param run_hierarchical_ph_clustering Logical, perform hierarchical
#'   clustering on PH distances.
#' @param run_spectral_clustering Logical, run spectral clustering.
#' @param run_visualizations Logical, generate plots for each iteration.
#' @param run_sample_level_heatmap Logical, create sample-level heatmaps.
#' @param run_cluster Logical, run clustering comparisons if `TRUE`.
#' @param run_betti Logical, compute Betti curves if `TRUE`.
#' @param run_cross_iteration Logical, run cross-iteration analysis.
#' @param metadata_path Optional path to metadata CSV for plotting.
#' @param SRA_col Metadata column with sample identifiers.
#' @param Tissue_col Metadata column containing tissue labels.
#' @param Approach_col Metadata column describing experimental approach.
#' @param ... Additional arguments passed to helper functions.
#'
#' @return Invisibly returns a list with generated results.
#'
#' @examples
#' \dontrun{
#' ph_results <- process_datasets_PH(read.csv("metadata.csv"))
#' run_postprocessing_pipeline(ph_results, num_cores = 2)
#' }
#' @export
run_postprocessing_pipeline <- function(ph_results,
                                        results_dir = "results",
                                        num_cores = parallel::detectCores(),
                                        run_standard_seurat_clustering = TRUE,
                                        run_kmeans_clustering = TRUE,
                                        run_hierarchical_ph_clustering = TRUE,
                                        run_spectral_clustering = TRUE,
                                        run_visualizations = TRUE,
                                        run_sample_level_heatmap = TRUE,
                                        run_cluster = TRUE,
                                        run_betti = TRUE,
                                        run_cross_iteration = TRUE,
                                        metadata_path = "./data/VastlyDifferentTissues/metadata.csv",
                                        SRA_col = "orig.ident",
                                        Tissue_col = "Tissue",
                                        Approach_col = "Approach",
                                        ...) {
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  metadata <- if (!is.null(metadata_path) && file.exists(metadata_path)) {
    readr::read_csv(metadata_path)
  } else {
    NULL
  }

  data_iterations <- ph_results$data_iterations
  if (is.null(data_iterations) || length(data_iterations) == 0) {
    stop("'ph_results' must contain data_iterations")
  }

  for (i in seq_along(data_iterations)) {
    iter <- data_iterations[[i]]
    log_message(paste("Processing", iter$name))

    # compute matrices if necessary
    process_iteration_calculate_matrices(iteration = iter,
                                         num_cores = num_cores,
                                         log_message = log_message)

    seurat_obj <- iter$seurat_obj
    assay <- iter$assay
    variable_features_path <- if ("variable_features" %in% names(iter)) iter$variable_features else NULL

    # Assign balanced random groups for permutation-based analyses
    seurat_obj <- assignRandomGroup(seurat_obj, k = 5, new_col_name = "Random_Group")

    if (run_standard_seurat_clustering) {
      seurat_obj <- perform_standard_seurat_clustering(seurat_obj, assay, variable_features_path)
      seurat_obj@meta.data[[paste0("seurat_cluster_", tolower(iter$name))]] <- Idents(seurat_obj)
      Idents(seurat_obj) <- seurat_obj@meta.data[[paste0("seurat_cluster_", tolower(iter$name))]]
      seurat_obj@meta.data$seurat_clusters <- NULL
    }

    bdm <- if (!is.null(iter$bdm_matrix) && file.exists(iter$bdm_matrix)) as.matrix(readRDS(iter$bdm_matrix)) else NULL
    sdm <- if (!is.null(iter$sdm_matrix) && file.exists(iter$sdm_matrix)) as.matrix(readRDS(iter$sdm_matrix)) else NULL
    landscape <- if (!is.null(iter$landscape_l2_distance_matrix) && file.exists(iter$landscape_l2_distance_matrix)) as.matrix(readRDS(iter$landscape_l2_distance_matrix)) else NULL

    seurat_obj <- apply_all_clustering_methods(
      seurat_obj,
      dataset_name = iter$name,
      assay = assay,
      bdm_matrix = bdm,
      sdm_matrix = sdm,
      landscape_matrix = landscape,
      run_kmeans_clustering = run_kmeans_clustering,
      run_hierarchical_ph_clustering = run_hierarchical_ph_clustering,
      run_spectral_clustering = run_spectral_clustering,
      SRA_col = SRA_col,
      Tissue_col = Tissue_col,
      Approach_col = Approach_col
    )

    seurat_obj <- generate_visualizations_for_iteration(
      seurat_obj = seurat_obj,
      dataset_name = iter$name,
      assay = assay,
      bdm_matrix = bdm,
      sdm_matrix = sdm,
      landscape_matrix = landscape,
      metadata = metadata,
      plots_folder = file.path(results_dir, "plots"),
      run_visualizations = run_visualizations,
      run_sample_level_heatmap = run_sample_level_heatmap,
      SRA_col = SRA_col,
      Tissue_col = Tissue_col,
      Approach_col = Approach_col
    )

    save_path <- file.path(results_dir, "seurat_objects",
                           paste0(tolower(iter$name), "_seurat_object.rds"))
    if (!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)
    saveRDS(seurat_obj, save_path)

    data_iterations[[i]]$seurat_obj <- seurat_obj
  }

  ph_results$data_iterations <- data_iterations

  if (run_cluster || run_betti || run_cross_iteration) {
    run_modular_analysis(ph_results,
                         results_dir = results_dir,
                         run_cluster = run_cluster,
                         run_betti = run_betti,
                         run_cross_iteration = run_cross_iteration,
                         SRA_col = SRA_col,
                         Tissue_col = Tissue_col,
                         Approach_col = Approach_col,
                         ...)
  }

  invisible(ph_results)
}


compute_and_save_distance_matrices <- function(
    pd_list_path, expr_list, bdm_path, sdm_path, dataset_name, num_cores = 32,
    log_message, recompute_bdm = TRUE, recompute_sdm = TRUE) {

  # Function to check if a matrix is valid
  is_valid_matrix <- function(matrix_path, expected_length = NULL) {
    if (!file.exists(matrix_path)) {
      return(FALSE)
    }
    tryCatch({
      matrix_data <- readRDS(matrix_path)
      if (!is.matrix(matrix_data) && !is.data.frame(matrix_data)) {
        return(FALSE)
      }
      if (any(!is.finite(as.numeric(matrix_data)))) {
        return(FALSE)
      }
      if (!is.null(expected_length) && nrow(matrix_data) != expected_length) {
        return(FALSE)
      }
      return(TRUE)
    }, error = function(e) {
      return(FALSE)
    })
  }

  # Load the PD list (needed for computation) and get names from the provided expr_list
  pd_list <- readRDS(pd_list_path)
  expr_list_names <- names(expr_list)
  names(pd_list) <- expr_list_names
  saveRDS(pd_list, file = pd_list_path)

  # Unique progress file for BDM
  bdm_progress_file <- paste0("BDM_progress_", dataset_name, ".rds")

  # Check or load/recompute BDM
  if (is_valid_matrix(bdm_path, expected_length = length(expr_list_names)) && !recompute_bdm) {
    log_message(paste("Loading preexisting and valid BDM for", dataset_name))
    BDM <- readRDS(bdm_path)
    # Update row and column names if they differ or are missing
    if (is.null(rownames(BDM)) || is.null(colnames(BDM)) ||
        !identical(rownames(BDM), expr_list_names) || !identical(colnames(BDM), expr_list_names)) {
      log_message("Updating row and column names for existing BDM.")
      rownames(BDM) <- expr_list_names
      colnames(BDM) <- expr_list_names
      BDM_df <- as.data.frame(BDM)
      saveRDS(BDM_df, file = bdm_path)
      write.csv(BDM_df, file = sub(".rds", ".csv", bdm_path))
    }
  } else {
    log_message(paste("Creating Bottleneck Distance Matrix for", dataset_name))
    BDM <- tryCatch({
      CreateBottleneckDistanceMatrixParallel(
          PD = pd_list,
          max_cores = num_cores,
          log_message = log_message,
          progress_file = bdm_progress_file,
          dataset_name = dataset_name
        )
    },
      error = function(e) {
        log_message(paste("Error in creating BDM for", dataset_name, ":", e$message))
        NULL
      }
    )

    if (!is.null(BDM)) {
      tryCatch({
        BDM_df <- as.data.frame(BDM)
        rownames(BDM_df) <- expr_list_names
        colnames(BDM_df) <- expr_list_names

        saveRDS(BDM_df, file = bdm_path)
        write.csv(BDM_df, file = sub(".rds", ".csv", bdm_path))
        log_message(paste("BDM for", dataset_name, "saved successfully."))
      },
        error = function(e) {
          log_message(paste("Error in saving BDM for", dataset_name, ":", e$message))
        }
      )
    }
  }

  # Check or recompute SDM
  if (is_valid_matrix(sdm_path, expected_length = length(expr_list_names)) && !recompute_sdm) {
    log_message(paste("Loading preexisting and valid SDM for", dataset_name))
    SDM <- readRDS(sdm_path)
    # Update row and column names if they differ or are missing
    if (is.null(rownames(SDM)) || is.null(colnames(SDM)) ||
        !identical(rownames(SDM), expr_list_names) || !identical(colnames(SDM), expr_list_names)) {
      log_message("Updating row and column names for existing SDM.")
      rownames(SDM) <- expr_list_names
      colnames(SDM) <- expr_list_names
      SDM_df <- as.data.frame(SDM)
      saveRDS(SDM_df, file = sdm_path)
      write.csv(SDM_df, file = sub(".rds", ".csv", sdm_path))
    }
  } else {
    log_message(paste("Creating Spectral Distance Matrix for", dataset_name))
    SDM <- tryCatch({
      CreateSpectralDistanceMatrixFromPD(
          pd_list = pd_list,
          bdm_matrix = BDM, # Pass validated BDM
          num_eigen = 50,
          log_message = log_message,
          dataset_name = dataset_name
        )
    },
      error = function(e) {
        log_message(paste("Error in creating SDM for", dataset_name, ":", e$message))
        NULL
      }
    )

    if (!is.null(SDM)) {
      tryCatch({
        SDM_df <- as.data.frame(SDM)
        rownames(SDM_df) <- expr_list_names
        colnames(SDM_df) <- expr_list_names

        saveRDS(SDM_df, file = sdm_path)
        write.csv(SDM_df, file = sub(".rds", ".csv", sdm_path))
        log_message(paste("SDM for", dataset_name, "saved successfully."))
      },
        error = function(e) {
          log_message(paste("Error in saving SDM for", dataset_name, ":", e$message))
        }
      )
    }
  }

  # Final consistency check: load, update (if needed), and resave the matrices using names from expr_list
  final_BDM <- tryCatch(readRDS(bdm_path), error = function(e) NULL)
  if (!is.null(final_BDM)) {
    if (is.null(rownames(final_BDM)) || is.null(colnames(final_BDM)) ||
        !identical(rownames(final_BDM), expr_list_names) || !identical(colnames(final_BDM), expr_list_names)) {
      log_message("Final check: Updating row and column names for BDM.")
      rownames(final_BDM) <- expr_list_names
      colnames(final_BDM) <- expr_list_names
      final_BDM_df <- as.data.frame(final_BDM)
      saveRDS(final_BDM_df, file = bdm_path)
      write.csv(final_BDM_df, file = sub(".rds", ".csv", bdm_path))
    }
  } else {
    log_message("Final check: BDM could not be loaded for final validation.")
  }

  final_SDM <- tryCatch(readRDS(sdm_path), error = function(e) NULL)
  if (!is.null(final_SDM)) {
    if (is.null(rownames(final_SDM)) || is.null(colnames(final_SDM)) ||
        !identical(rownames(final_SDM), expr_list_names) || !identical(colnames(final_SDM), expr_list_names)) {
      log_message("Final check: Updating row and column names for SDM.")
      rownames(final_SDM) <- expr_list_names
      colnames(final_SDM) <- expr_list_names
      final_SDM_df <- as.data.frame(final_SDM)
      saveRDS(final_SDM_df, file = sdm_path)
      write.csv(final_SDM_df, file = sub(".rds", ".csv", sdm_path))
    }
  } else {
    log_message("Final check: SDM could not be loaded for final validation.")
  }
}

# ---------------------------
# New Helper Function: Compute Persistence Landscapes for Dimensions 0 and 1
# ---------------------------
ComputePersistenceLandscapes <- function(pd, grid = seq(0, 1, length.out = 100)) {
  res <- tryCatch({
    landscape0 <- TDA::landscape(Diag = pd, dimension = 0, tseq = grid)
    landscape1 <- TDA::landscape(Diag = pd, dimension = 1, tseq = grid)
    list(dim0 = landscape0, dim1 = landscape1)
  }, error = function(e) {
    log_message(paste("Error computing persistence landscape for a PD:", e$message))
    return(NULL)
  })
  return(res)
}

# ---------------------------
# New Function: Compute and Save Landscape List and Combined Distance Matrix
# ---------------------------
compute_and_save_landscape_matrices <- function(
    pd_list_path,
    landscape_list_path,
    landscape_distance_matrix_path,
    log_message,
    num_cores = 6) {

  # Load the PD list
  pd_list <- readRDS(pd_list_path)

  # Compute persistence landscapes for each PD in parallel
  landscape_list <- mclapply(pd_list, function(pd) {
    ComputePersistenceLandscapes(pd)
  }, mc.cores = num_cores)

  # Define expected grid length (as used in ComputePersistenceLandscapes)
  grid_length <- length(seq(0, 1, length.out = 100))

  # Validate structure of computed landscape_list before saving
  valid_list <- TRUE
  for (i in seq_along(landscape_list)) {
    if (!is.list(landscape_list[[i]]) || is.null(landscape_list[[i]]$dim0) || is.null(landscape_list[[i]]$dim1)) {
      log_message(paste("Validation Error: Landscape element", i, "is missing dim0 and/or dim1 components."))
      valid_list <- FALSE
    } else {
      if (length(landscape_list[[i]]$dim0) != grid_length || length(landscape_list[[i]]$dim1) != grid_length) {
        log_message(paste("Validation Error: Landscape element", i, "does not have the expected grid length of", grid_length))
        valid_list <- FALSE
      }
    }
  }
  if (!valid_list) {
    log_message("Aborting save: Computed landscape list structure is invalid.")
    return(invisible(NULL))
  }

  # Save the computed landscape list
  saveRDS(landscape_list, file = landscape_list_path)
  log_message(paste("Saved persistence landscape list to", landscape_list_path))

  # Validate that the landscape list was saved correctly by re-loading it
  if (file.exists(landscape_list_path)) {
    loaded_landscape_list <- tryCatch({
      readRDS(landscape_list_path)
    }, error = function(e) {
      log_message(paste("Validation Error: Failed to load landscape list from", landscape_list_path, ":", e$message))
      NULL
    })

    if (!is.null(loaded_landscape_list)) {
      if (isTRUE(all.equal(loaded_landscape_list, landscape_list))) {
        log_message("Validation passed: Landscape list saved and loaded correctly.")
      } else {
        log_message("Validation Warning: Reloaded landscape list differs from the in-memory object.")
      }
    }
  } else {
    log_message(paste("Validation Error: Landscape list file", landscape_list_path, "does not exist."))
  }

  # Number of PDs
  n <- length(landscape_list)

  # Compute pairwise distances for each pair of PDs to form the combined distance matrix
  combined_distance_matrix <- matrix(0, n, n)
  rownames(combined_distance_matrix) <- names(pd_list)
  colnames(combined_distance_matrix) <- names(pd_list)

  for (i in 1:n) {
    for (j in i:n) {
      if (any(is.na(landscape_list[[i]]$dim0)) || any(is.na(landscape_list[[j]]$dim0))) {
        combined_dist <- NA
      } else {
        dist_dim0 <- sqrt(sum((landscape_list[[i]]$dim0 - landscape_list[[j]]$dim0) ^ 2))
        dist_dim1 <- sqrt(sum((landscape_list[[i]]$dim1 - landscape_list[[j]]$dim1) ^ 2))
        combined_dist <- sqrt(dist_dim0 ^ 2 + dist_dim1 ^ 2)
      }
      combined_distance_matrix[i, j] <- combined_dist
      combined_distance_matrix[j, i] <- combined_dist
    }
  }

  # Validate structure of combined_distance_matrix before saving
  if (!is.matrix(combined_distance_matrix)) {
    log_message("Validation Error: Combined distance matrix is not a matrix. Aborting save.")
    return(invisible(NULL))
  }
  if (nrow(combined_distance_matrix) != ncol(combined_distance_matrix)) {
    log_message("Validation Error: Combined distance matrix is not square. Aborting save.")
    return(invisible(NULL))
  }
  if (is.null(rownames(combined_distance_matrix)) || is.null(colnames(combined_distance_matrix))) {
    log_message("Validation Error: Combined distance matrix does not have row/column names. Aborting save.")
    return(invisible(NULL))
  }
  if (!isTRUE(all.equal(combined_distance_matrix, t(combined_distance_matrix)))) {
    log_message("Validation Warning: Combined distance matrix is not symmetric.")
  } else {
    log_message("Validation passed: Combined distance matrix is a symmetric square matrix.")
  }

  # Save the combined distance matrix (RDS and CSV)
  saveRDS(combined_distance_matrix, file = landscape_distance_matrix_path)
  write.csv(combined_distance_matrix, file = sub(".Rds", ".csv", landscape_distance_matrix_path))
  log_message(paste("Saved combined landscape distance matrix to", landscape_distance_matrix_path))

  # Validate that the combined distance matrix was saved correctly by re-loading it
  if (file.exists(landscape_distance_matrix_path)) {
    loaded_combined_matrix <- tryCatch({
      readRDS(landscape_distance_matrix_path)
    }, error = function(e) {
      log_message(paste("Validation Error: Failed to load combined distance matrix from",
                        landscape_distance_matrix_path, ":", e$message))
      NULL
    })

    if (!is.null(loaded_combined_matrix)) {
      if (isTRUE(all.equal(loaded_combined_matrix, combined_distance_matrix))) {
        log_message("Validation passed: Combined distance matrix saved and loaded correctly.")
      } else {
        log_message("Validation Warning: Reloaded combined distance matrix differs from the in-memory object.")
      }
    }
  } else {
    log_message(paste("Validation Error: Combined distance matrix file", landscape_distance_matrix_path, "does not exist."))
  }

  # Additional file check: Log file size for troubleshooting
  if (file.exists(landscape_distance_matrix_path)) {
    file_size <- file.info(landscape_distance_matrix_path)$size
    log_message(paste("File", landscape_distance_matrix_path, "exists with size:", file_size, "bytes."))
  }
}

# ---------------------------
# Updated Per-Iteration Function: Integrate All Computations
# ---------------------------
process_iteration_calculate_matrices <- function(iteration, num_cores, log_message) {
  log_message(paste("Processing dataset:", iteration$name))

  # Compute Bottleneck and Spectral Distance Matrices
  compute_and_save_distance_matrices(
    pd_list_path = iteration$pd_list,
    expr_list = iteration$expr_list,
    bdm_path = iteration$bdm_matrix,
    sdm_path = iteration$sdm_matrix,
    dataset_name = iteration$name,
    num_cores = num_cores,
    log_message = log_message
  )

  # Compute Persistence Landscapes and Combined Distance Matrix for Dimensions 0 and 1
  log_message(paste("Computing persistence landscapes and combined distance matrix for", iteration$name))
  compute_and_save_landscape_matrices(
    pd_list_path = iteration$pd_list,
    landscape_list_path = iteration$landscape_list,
    landscape_distance_matrix_path = iteration$landscape_l2_distance_matrix,
    log_message = log_message,
    num_cores = num_cores
  )
}

# Validate that an object is a matrix (or data frame) with correct dimensions and names
validate_matrix <- function(mat, expected_dim = NULL, expected_names = NULL, verbose = TRUE) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    if (verbose) message("Validation failed: Object is not a matrix or data.frame.")
    return(FALSE)
  }

  dims <- dim(mat)
  if (verbose) message("Matrix dimensions: ", paste(dims, collapse = " x "))
  if (!is.null(expected_dim) && !all(dims == expected_dim)) {
    if (verbose) message("Dimension mismatch. Expected: ", paste(expected_dim, collapse = " x "))
    return(FALSE)
  }

  if (!is.null(expected_names)) {
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      if (verbose) message("Row or column names are missing.")
      return(FALSE)
    }
    if (!all(rownames(mat) == expected_names)) {
      if (verbose) message("Row names do not match expected names.")
      return(FALSE)
    }
    if (!all(colnames(mat) == expected_names)) {
      if (verbose) message("Column names do not match expected names.")
      return(FALSE)
    }
  }

  if (!isTRUE(all.equal(mat, t(mat)))) {
    if (verbose) message("Warning: Matrix is not symmetric.")
  } else if (verbose) {
    message("Matrix is symmetric.")
  }

  TRUE
}

# Validate that the landscape list is structured properly.
validate_landscape_list <- function(landscape_list, expected_length = NULL, grid_length = 100, verbose = TRUE) {
  if (!is.list(landscape_list)) {
    if (verbose) message("Validation failed: Landscape object is not a list.")
    return(FALSE)
  }

  n <- length(landscape_list)
  if (verbose) message("Length of landscape list: ", n)
  if (!is.null(expected_length) && n != expected_length) {
    if (verbose) message("Landscape list length does not match expected value: ", expected_length)
    return(FALSE)
  }

  valid <- TRUE
  for (i in seq_along(landscape_list)) {
    elem <- landscape_list[[i]]
    if (!is.list(elem) || is.null(elem$dim0) || is.null(elem$dim1)) {
      if (verbose) message("Element ", i, " is not properly structured (missing dim0 or dim1).")
      valid <- FALSE
    } else {
      len0 <- length(elem$dim0)
      len1 <- length(elem$dim1)
      if (verbose) message("Element ", i, " has dim0 length ", len0, " and dim1 length ", len1)
      if (len0 != grid_length || len1 != grid_length) {
        if (verbose) message("Element ", i, " does not have the expected grid length of ", grid_length)
        valid <- FALSE
      }
    }
  }
  valid
}


# ---------------------------
# Functions for Modular Analysis
# ---------------------------


# Function: Perform standard Seurat clustering, handling various assays and precomputed variable features
perform_standard_seurat_clustering <- function(seurat_obj,
                                               assay = "integrated",
                                               resolution = 0.5,
                                               variable_features_path = NULL,
                                               verbose = TRUE) {
  # Set active assay
  DefaultAssay(seurat_obj) <- assay

  # Logging helper function
  log_message <- function(message) {
    if (verbose) message(sprintf("%s - %s", Sys.time(), message))
  }

  log_message(paste("Starting clustering for assay:", assay))

  # Assay-specific handling
  if (assay == "RNA") {
    log_message("Assay type: RNA. Running FindVariableFeatures, NormalizeData, and ScaleData.")
    seurat_obj <- FindVariableFeatures(seurat_obj, assay = assay)
    seurat_obj <- NormalizeData(seurat_obj, assay = assay, verbose = TRUE)
    seurat_obj <- ScaleData(seurat_obj, assay = assay, verbose = TRUE)
  } else if (assay %in% c("SCT", "SCT_Ind")) {
    log_message("Assay type: SCT or SCT_Ind.")
    if (length(VariableFeatures(seurat_obj)) < 2) {
      log_message("Insufficient variable features found. Running FindVariableFeatures.")
      seurat_obj <- FindVariableFeatures(seurat_obj, assay = assay)
    }
    seurat_obj <- ScaleData(seurat_obj, assay = assay, verbose = TRUE)
  } else if (assay == "integrated") {
    log_message("Assay type: Integrated. Checking for SCT models and handling variable features.")

    # For Seurat v5 integrated assays, an SCTModel.list is not expected.
    if ("SCTModel.list" %in% slotNames(seurat_obj@assays[[assay]])) {
      sct_models <- seurat_obj@assays[[assay]]@SCTModel.list
      if (length(sct_models) > 1) {
        log_message("Multiple SCT models found. Merging models.")
        feature_data_frames <- list()
        all_columns <- unique(unlist(lapply(sct_models, function(model) colnames(model@feature.attributes))))
        for (model in sct_models) {
          feature_data <- model@feature.attributes
          missing_columns <- setdiff(all_columns, colnames(feature_data))
          feature_data[missing_columns] <- NA
          feature_data_frames <- append(feature_data_frames, list(feature_data))
        }
        consolidated_features <- do.call(rbind, feature_data_frames)
        consolidated_features <- consolidated_features[!duplicated(rownames(consolidated_features)),]
        consolidated_model <- sct_models[[1]]
        consolidated_model@feature.attributes <- consolidated_features
        seurat_obj@assays[[assay]]@SCTModel.list <- list(consolidated_model)
        log_message("Successfully merged SCT models into a single consolidated model.")
      } else {
        log_message("Only one SCT model is present; no merging is required.")
      }
    } else {
      log_message("No SCTModel.list slot found in the integrated assay (expected for Seurat v5). Skipping model merging.")
    }

    # If a variable features file is provided, load and set it.
    if (!is.null(variable_features_path)) {
      variable_features <- readRDS(variable_features_path)
      if (!is.null(variable_features) && is.character(variable_features)) {
        compatible_features <- intersect(variable_features, rownames(seurat_obj[["integrated"]]))
        seurat_obj@assays[["integrated"]]@var.features <- compatible_features
        log_message("Variable features successfully set on the integrated data.")
      } else {
        log_message("The variable features file did not load correctly or is not a character vector.")
      }
    } else {
      log_message("Variable features path is NULL. Using all features for PCA.")
    }

    # Strategy for feature selection for PCA:
    # 1. If scale.data exists and has proper rownames, use that.
    # 2. Otherwise, run ScaleData to create scale.data.
    # 3. If after running ScaleData scale.data is still not proper, fall back to the data slot.
    if (is.null(seurat_obj@assays[[assay]]@scale.data) ||
      nrow(seurat_obj@assays[[assay]]@scale.data) == 0 ||
      is.null(rownames(seurat_obj@assays[[assay]]@scale.data))) {
      log_message("scale.data is missing or lacks proper rownames. Running ScaleData to generate scale.data.")
      seurat_obj <- ScaleData(seurat_obj, assay = assay, verbose = TRUE)
    }

    if (!is.null(seurat_obj@assays[[assay]]@scale.data) &&
      nrow(seurat_obj@assays[[assay]]@scale.data) > 0 &&
      !is.null(rownames(seurat_obj@assays[[assay]]@scale.data))) {
      all_features <- rownames(seurat_obj@assays[[assay]]@scale.data)
      log_message("Using features from scale.data.")
    } else if (!is.null(seurat_obj@assays[[assay]]@data) &&
      nrow(seurat_obj@assays[[assay]]@data) > 0 &&
      !is.null(rownames(seurat_obj@assays[[assay]]@data))) {
      all_features <- rownames(seurat_obj@assays[[assay]]@data)
      log_message("scale.data still not available; falling back to using features from the data slot.")
    } else {
      stop("No features found in the integrated assay even after attempting ScaleData and fallback to the data slot. Cannot proceed with PCA.")
    }
  } else {
    stop("Unsupported assay type.")
  }

  # Run PCA if not already computed
  if (!"pca" %in% names(seurat_obj@reductions)) {
    log_message("PCA not found, running RunPCA.")
    if (assay == "integrated" && is.null(variable_features_path)) {
      seurat_obj <- RunPCA(seurat_obj, assay = assay, features = all_features, npcs = 50, verbose = TRUE)
    } else {
      seurat_obj <- RunPCA(seurat_obj, assay = assay, features = VariableFeatures(seurat_obj), npcs = 50, verbose = TRUE)
    }
  } else {
    log_message("PCA already found, skipping RunPCA step.")
  }

  # Proceed with clustering
  log_message("Finding neighbors and clusters.")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

  log_message("Clustering complete.")
  return(seurat_obj)
}

# Function: K-means Clustering
perform_kmeans_clustering <- function(seurat_obj, assay, dims = 1:50, k = 5, verbose = TRUE) {
  log_message <- function(message) {
    if (verbose) message(sprintf("%s - %s", Sys.time(), message))
  }

  log_message(paste("Starting K-means clustering on assay:", assay, "with", k, "clusters"))

  DefaultAssay(seurat_obj) <- assay
  log_message("Default assay set for Seurat object.")

  # Extract PCA embeddings for specified dimensions
  if ("pca" %in% names(seurat_obj@reductions)) {
    clustering_data <- Embeddings(seurat_obj, reduction = "pca")[, dims]
    log_message("PCA embeddings retrieved for clustering.")
  } else {
    stop("PCA not found in the Seurat object. Run PCA before K-means clustering.")
  }

  # Perform K-means clustering
  log_message("Running K-means clustering.")
  kmeans_result <- kmeans(clustering_data, centers = k, nstart = 25, iter.max = 100)
  log_message("K-means clustering completed.")

  # Assign cluster labels to Seurat object
  seurat_obj$kmeans_cluster <- factor(kmeans_result$cluster)
  log_message("K-means cluster labels assigned to Seurat object.")

  return(seurat_obj)
}

# Function: Hierarchical Clustering using BDM matrix
perform_hierarchical_clustering_ph <- function(bdm_matrix, k, verbose = TRUE) {
  log_message <- function(message) {
    if (verbose) message(sprintf("%s - %s", Sys.time(), message))
  }

  log_message("Starting hierarchical clustering using BDM matrix.")
  if (!is.matrix(bdm_matrix)) {
    stop("BDM input must be a matrix.")
  }

  # Convert BDM matrix to distance object and perform hierarchical clustering
  log_message("Calculating distance and performing hierarchical clustering.")
  hc <- hclust(as.dist(bdm_matrix), method = "ward.D2")
  log_message("Hierarchical clustering completed.")

  # Cut dendrogram to form clusters
  clusters_ph <- cutree(hc, k = k)
  log_message(paste("Clusters formed using cutree with k =", k))

  # Return both the clusters and the hierarchical tree
  return(list(clusters = clusters_ph, tree = hc))
}

# Function: Assign PH Clusters to Seurat Object
assign_ph_clusters <- function(seurat_obj, clusters_ph, new_cluster_col,
                               SRA_col = "orig.ident", verbose = TRUE) {
  log_message <- function(message) {
    if (verbose) message(sprintf("%s - %s", Sys.time(), message))
  }

  log_message(paste("Assigning PH clusters to Seurat object with column name:", new_cluster_col))

  # if ("SRA" %in% colnames(seurat_obj@meta.data)) {
  #   if (length(unique(seurat_obj@meta.data$SRA)) > 1) {
  #     seurat_sample_ids <- seurat_obj@meta.data$SRA
  #   } else if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
  #     seurat_sample_ids <- seurat_obj@meta.data$orig.ident
  #   } else if (any(grepl("__", rownames(seurat_obj@meta.data)))) {
  #     seurat_sample_ids <- sapply(strsplit(rownames(seurat_obj@meta.data), "__"), `[`, 1)
  #   } else {
  #     stop("SRA column exists but contains a single unique value and neither orig.ident nor rownames are informative.")
  #   }
  # } else
  if (SRA_col %in% colnames(seurat_obj@meta.data)) {
    seurat_sample_ids <- seurat_obj@meta.data[[SRA_col]]
  } else if (any(grepl("__", rownames(seurat_obj@meta.data)))) {
    seurat_sample_ids <- sapply(strsplit(rownames(seurat_obj@meta.data), "__"), `[`, 1)
  } else {
    stop("No suitable sample identifier found in seurat_obj@meta.data.")
  }

  log_message("Sample IDs extracted from Seurat object:")
  if (verbose) log_message(paste(head(seurat_sample_ids), collapse = ", "))

  # Identify common sample-level identifiers between Seurat object and clusters_ph
  common_sample_ids <- intersect(seurat_sample_ids, names(clusters_ph))
  log_message(paste("Found", length(common_sample_ids), "common sample identifiers for cluster assignment."))

  # Ensure that there are common sample identifiers; if not, issue a warning and exit
  if (length(common_sample_ids) == 0) {
    warning("No common sample identifiers found between Seurat object and clusters_ph.")
    return(seurat_obj)
  }

  # Initialize the new cluster column with NA values if it doesn’t exist
  if (!new_cluster_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data[[new_cluster_col]] <- NA
  }

  # Assign PH clusters to cells based on the sample-level identifier
  for (sample_id in common_sample_ids) {
    cluster_value <- clusters_ph[sample_id]
    if (!is.na(cluster_value)) {
      # Assign the cluster value from clusters_ph to all cells of this sample in the new column
      seurat_obj@meta.data[[new_cluster_col]][seurat_sample_ids == sample_id] <- cluster_value
    }
  }

  # Convert to factor to align with other Seurat metadata
  seurat_obj@meta.data[[new_cluster_col]] <- factor(seurat_obj@meta.data[[new_cluster_col]], levels = unique(clusters_ph))

  # Log any cells still with NA values in the new cluster column
  missing_cells <- sum(is.na(seurat_obj@meta.data[[new_cluster_col]]))
  if (missing_cells > 0) {
    log_message(paste(missing_cells, "cells could not be assigned a PH cluster."))
  }

  log_message(paste("Assigned PH clusters to column:", new_cluster_col))

  return(seurat_obj)
}


ensure_directory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# ---------------------------
# Additional Clustering Helpers
# ---------------------------

#' Perform spectral clustering on a distance matrix
perform_spectral_clustering <- function(distance_matrix, k) {
  epsilon <- 1e-8
  res <- kernlab::specc(kernlab::as.kernelMatrix(similarity_matrix), centers = k)
  clusters <- res@.Data
  names(clusters) <- rownames(distance_matrix)
  clusters
}

#' Apply all clustering approaches to a Seurat object
apply_all_clustering_methods <- function(seurat_obj, dataset_name, assay,
                                         bdm_matrix = NULL, sdm_matrix = NULL,
                                         landscape_matrix = NULL,
                                         run_kmeans_clustering = TRUE,
                                         run_hierarchical_ph_clustering = TRUE,
                                         run_spectral_clustering = TRUE,
                                         SRA_col = "orig.ident",
                                         Tissue_col = "Tissue",
                                         Approach_col = "Approach") {

  prefix <- tolower(dataset_name)

  k_tissue <- if (Tissue_col %in% colnames(seurat_obj@meta.data))
    length(unique(seurat_obj@meta.data[[Tissue_col]])) else NULL
  k_sra <- if (SRA_col %in% colnames(seurat_obj@meta.data))
    length(unique(seurat_obj@meta.data[[SRA_col]])) else NULL
  k_approach <- if (Approach_col %in% colnames(seurat_obj@meta.data))
    length(unique(seurat_obj@meta.data[[Approach_col]])) else NULL

  ## ---- K-means ---------------------------------------------------------
  if (run_kmeans_clustering) {
    if (!is.null(k_tissue)) {
      seurat_obj <- perform_kmeans_clustering(seurat_obj, assay, k = k_tissue)
      seurat_obj@meta.data[[paste0("kmeans_cluster_", prefix, "_tissue")]] <- seurat_obj$kmeans_cluster
    }
    if (!is.null(k_sra)) {
      seurat_obj <- perform_kmeans_clustering(seurat_obj, assay, k = k_sra)
      seurat_obj@meta.data[[paste0("kmeans_cluster_", prefix, "_sra")]] <- seurat_obj$kmeans_cluster
    }
    if (!is.null(k_approach)) {
      seurat_obj <- perform_kmeans_clustering(seurat_obj, assay, k = k_approach)
      seurat_obj@meta.data[[paste0("kmeans_cluster_", prefix, "_approach")]] <- seurat_obj$kmeans_cluster
    }
    seurat_obj$kmeans_cluster <- NULL
  }

  ## ---- Hierarchical Clustering ---------------------------------------
  if (run_hierarchical_ph_clustering) {
    if (!is.null(bdm_matrix)) {
      if (!is.null(k_tissue)) {
        res <- perform_hierarchical_clustering_ph(bdm_matrix, k_tissue)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_bdm_ph_", prefix, "_tissue"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_sra)) {
        res <- perform_hierarchical_clustering_ph(bdm_matrix, k_sra)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_bdm_ph_", prefix, "_sra"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_approach)) {
        res <- perform_hierarchical_clustering_ph(bdm_matrix, k_approach)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_bdm_ph_", prefix, "_approach"),
          SRA_col = SRA_col)
      }
    }

    if (!is.null(sdm_matrix)) {
      if (!is.null(k_tissue)) {
        res <- perform_hierarchical_clustering_ph(sdm_matrix, k_tissue)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_sdm_ph_", prefix, "_tissue"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_sra)) {
        res <- perform_hierarchical_clustering_ph(sdm_matrix, k_sra)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_sdm_ph_", prefix, "_sra"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_approach)) {
        res <- perform_hierarchical_clustering_ph(sdm_matrix, k_approach)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_sdm_ph_", prefix, "_approach"),
          SRA_col = SRA_col)
      }
    }

    if (!is.null(landscape_matrix)) {
      if (!is.null(k_tissue)) {
        res <- perform_hierarchical_clustering_ph(landscape_matrix, k_tissue)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_landscape_ph_", prefix, "_tissue"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_sra)) {
        res <- perform_hierarchical_clustering_ph(landscape_matrix, k_sra)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_landscape_ph_", prefix, "_sra"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_approach)) {
        res <- perform_hierarchical_clustering_ph(landscape_matrix, k_approach)
        seurat_obj <- assign_ph_clusters(seurat_obj, res$clusters,
          paste0("hierarchical_cluster_landscape_ph_", prefix, "_approach"),
          SRA_col = SRA_col)
      }
    }
  }

  ## ---- Spectral Clustering -------------------------------------------
  if (run_spectral_clustering) {
    if (!is.null(bdm_matrix)) {
      if (!is.null(k_tissue)) {
        cl <- perform_spectral_clustering(bdm_matrix, k_tissue)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_bdm_", prefix, "_tissue"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_sra)) {
        cl <- perform_spectral_clustering(bdm_matrix, k_sra)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_bdm_", prefix, "_sra"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_approach)) {
        cl <- perform_spectral_clustering(bdm_matrix, k_approach)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_bdm_", prefix, "_approach"),
          SRA_col = SRA_col)
      }
    }

    if (!is.null(sdm_matrix)) {
      if (!is.null(k_tissue)) {
        cl <- perform_spectral_clustering(sdm_matrix, k_tissue)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_sdm_", prefix, "_tissue"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_sra)) {
        cl <- perform_spectral_clustering(sdm_matrix, k_sra)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_sdm_", prefix, "_sra"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_approach)) {
        cl <- perform_spectral_clustering(sdm_matrix, k_approach)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_sdm_", prefix, "_approach"),
          SRA_col = SRA_col)
      }
    }

    if (!is.null(landscape_matrix)) {
      if (!is.null(k_tissue)) {
        cl <- perform_spectral_clustering(landscape_matrix, k_tissue)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_landscape_", prefix, "_tissue"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_sra)) {
        cl <- perform_spectral_clustering(landscape_matrix, k_sra)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_landscape_", prefix, "_sra"),
          SRA_col = SRA_col)
      }
      if (!is.null(k_approach)) {
        cl <- perform_spectral_clustering(landscape_matrix, k_approach)
        seurat_obj <- assign_ph_clusters(seurat_obj, cl,
          paste0("spectral_cluster_landscape_", prefix, "_approach"),
          SRA_col = SRA_col)
      }
    }
  }

  seurat_obj
}

generate_heatmaps <- function(dataset_name, metadata, seurat_obj, bdm_matrix, plots_folder,
                              hc_tissue, hc_sra, hc_approach,
                              k_tissue, k_sra, k_approach,
                              SRA_col = "orig.ident",
                              Tissue_col = "Tissue",
                              Approach_col = "Approach") {
  if (!is.matrix(bdm_matrix) || nrow(bdm_matrix) != ncol(bdm_matrix)) {
    stop("BDM must be a square matrix with equal rows and columns.")
  }
  if (is.null(rownames(bdm_matrix)) || is.null(colnames(bdm_matrix))) {
    stop("BDM matrix must have row and column names.")
  }

  metadata2 <- metadata %>% dplyr::mutate(Sample = .data[[SRA_col]])
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::mutate(Sample = .data[[SRA_col]]) %>%
    dplyr::left_join(
      metadata2 %>% dplyr::select(Sample, all_of(c(SRA_col, Tissue_col, Approach_col))),
      by = "Sample",
      suffix = c(".seurat", ".meta")
    ) %>%
    dplyr::mutate(
      !!Tissue_col := coalesce(.data[[paste0(Tissue_col, ".meta")]], .data[[paste0(Tissue_col, ".seurat")]]),
      !!Approach_col := coalesce(.data[[paste0(Approach_col, ".meta")]], .data[[paste0(Approach_col, ".seurat")]])
    ) %>%
    dplyr::select(-matches("\\.meta$|\\.seurat$"))

  sample_level_metadata <- seurat_obj@meta.data %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarize(
      !!SRA_col := dplyr::first(.data[[SRA_col]]),
      !!Tissue_col := dplyr::first(.data[[Tissue_col]]),
      !!Approach_col := dplyr::first(.data[[Approach_col]]),
      across(contains("cluster"), ~ dplyr::first(.)),
      .groups = "drop"
    )
  annot_fields <- c(SRA_col, Tissue_col, Approach_col)

  cluster_cols <- colnames(sample_level_metadata)[
    grepl("cluster", colnames(sample_level_metadata), ignore.case = TRUE) &
      !grepl("random_group", colnames(sample_level_metadata), ignore.case = TRUE)
  ]
  if (length(cluster_cols) == 0) {
    stop("No cluster columns (excluding Random_Group) found in the sample-level metadata.")
  }

  hierarchical_cols <- cluster_cols[grepl("hierarchical", cluster_cols, ignore.case = TRUE)]
  non_hierarchical_cols <- setdiff(cluster_cols, hierarchical_cols)

  sample_annotations <- sample_level_metadata %>%
    dplyr::filter(Sample %in% rownames(bdm_matrix)) %>%
    dplyr::arrange(match(Sample, rownames(bdm_matrix)))

  if (!all(sample_annotations$Sample %in% rownames(bdm_matrix))) {
    stop("Sample identifiers in the metadata do not match the BDM matrix.")
  }
  bdm_matrix <- bdm_matrix[sample_annotations$Sample, sample_annotations$Sample, drop = FALSE]

  annotation_colors <- lapply(sample_annotations, function(col) {
    if (is.factor(col) || is.character(col)) {
      unique_vals <- unique(col)
      setNames(viridis::viridis(length(unique_vals)), unique_vals)
    } else {
      NULL
    }
  })
  annotation_colors <- annotation_colors[!sapply(annotation_colors, is.null)]

  ha <- HeatmapAnnotation(
    df = sample_annotations %>% dplyr::select(all_of(c(annot_fields, cluster_cols))),
    col = annotation_colors,
    annotation_name_side = "left"
  )

  pdf_path <- file.path(plots_folder, paste0("BDM_Heatmaps_", dataset_name, ".pdf"))
  if (file.exists(pdf_path)) {
    file.remove(pdf_path)
  }
  pdf(pdf_path, height = 20, width = 25)

  for (cluster_col in hierarchical_cols) {
    cluster_annotations <- sample_annotations %>% dplyr::pull(!!rlang::sym(cluster_col))
    cluster_colors <- annotation_colors[[cluster_col]]
    if (grepl("tissue", cluster_col, ignore.case = TRUE)) {
      hc <- hc_tissue; k_val <- k_tissue
    } else if (grepl("approach", cluster_col, ignore.case = TRUE)) {
      hc <- hc_approach; k_val <- k_approach
    } else if (grepl("sra", cluster_col, ignore.case = TRUE)) {
      hc <- hc_sra; k_val <- k_sra
    } else {
      hc <- hc_sra; k_val <- NA
    }

    sample_order <- order.dendrogram(as.dendrogram(hc))
    ordered_annotations <- cluster_annotations[sample_order]
    dend <- as.dendrogram(hc)
    colored_dendrogram <- dendextend::set(dend, labels_col = cluster_colors[ordered_annotations])
    unique_colors <- unique(cluster_colors[ordered_annotations])
    colored_dendrogram <- dendextend::color_branches(colored_dendrogram, k = length(unique_colors), col = unique_colors)

    plot(colored_dendrogram,
         main = paste("Dendrogram -", cluster_col, dataset_name),
         xlab = "", ylab = "Height", sub = "")
    text(x = 1, y = 0, labels = paste("Number of Clusters (k) =", k_val), pos = 4, col = "blue", font = 2)

    reordered_bdm_matrix <- bdm_matrix[sample_order, sample_order]

    draw(Heatmap(
      matrix = reordered_bdm_matrix,
      name = paste("Raw Distance -", cluster_col, "(Hierarchical)"),
      top_annotation = ha,
      row_order = sample_order,
      column_order = sample_order,
      cluster_rows = FALSE, cluster_columns = FALSE,
      col = viridis::viridis(100),
      heatmap_legend_param = list(title = "Raw Distance", legend_direction = "horizontal")
    ), annotation_legend_side = "bottom", heatmap_legend_side = "bottom")

    bdm_zscore <- scale(reordered_bdm_matrix, center = TRUE, scale = TRUE)
    draw(Heatmap(
      matrix = bdm_zscore,
      name = paste("Z-score Distance -", cluster_col, "(Hierarchical)"),
      top_annotation = ha,
      row_order = sample_order,
      column_order = sample_order,
      cluster_rows = FALSE, cluster_columns = FALSE,
      col = viridis::viridis(100),
      heatmap_legend_param = list(title = "Z-score Distance", legend_direction = "horizontal")
    ), annotation_legend_side = "bottom", heatmap_legend_side = "bottom")

    draw(Heatmap(
      matrix = bdm_matrix,
      name = paste("Raw Distance -", cluster_col, "(Natural Order)"),
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      col = viridis::viridis(100),
      heatmap_legend_param = list(title = "Raw Distance (Natural)", legend_direction = "horizontal")
    ), annotation_legend_side = "bottom", heatmap_legend_side = "bottom")

    bdm_zscore <- scale(bdm_matrix, center = TRUE, scale = TRUE)
    draw(Heatmap(
      matrix = bdm_zscore,
      name = paste("Z-score Distance -", cluster_col, "(Natural Order)"),
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      col = viridis::viridis(100),
      heatmap_legend_param = list(title = "Z-score Distance (Natural)", legend_direction = "horizontal")
    ), annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
  }

  if (length(non_hierarchical_cols) > 0) {
    rep_col <- non_hierarchical_cols[1]
    draw(Heatmap(
      matrix = bdm_matrix,
      name = paste("Raw Distance -", rep_col, "(Natural Order)"),
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      col = viridis::viridis(100),
      heatmap_legend_param = list(title = "Raw Distance (Natural)", legend_direction = "horizontal")
    ), annotation_legend_side = "bottom", heatmap_legend_side = "bottom")

    bdm_zscore <- scale(bdm_matrix, center = TRUE, scale = TRUE)
    draw(Heatmap(
      matrix = bdm_zscore,
      name = paste("Z-score Distance -", rep_col, "(Natural Order)"),
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      col = viridis::viridis(100),
      heatmap_legend_param = list(title = "Z-score Distance (Natural)", legend_direction = "horizontal")
    ), annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
  }

  dev.off()
  message("All heatmaps and dendrograms saved to: ", pdf_path)
}

#' Generate visualizations and heatmaps for a single iteration
generate_visualizations_for_iteration <- function(seurat_obj, dataset_name, assay,
                                                  bdm_matrix = NULL, sdm_matrix = NULL,
                                                  landscape_matrix = NULL, metadata = NULL,
                                                  plots_folder = "plots",
                                                  run_visualizations = TRUE,
                                                  run_sample_level_heatmap = TRUE,
                                                  SRA_col = "orig.ident",
                                                  Tissue_col = "Tissue",
                                                  Approach_col = "Approach") {
  dataset_lower <- tolower(dataset_name)

  k_tissue <- if (Tissue_col %in% colnames(seurat_obj@meta.data))
    length(unique(seurat_obj@meta.data[[Tissue_col]])) else NULL
  k_sra <- if (SRA_col %in% colnames(seurat_obj@meta.data))
    length(unique(seurat_obj@meta.data[[SRA_col]])) else NULL
  k_approach <- if (Approach_col %in% colnames(seurat_obj@meta.data))
    length(unique(seurat_obj@meta.data[[Approach_col]])) else NULL

  if (run_visualizations) {
    if (!"umap" %in% names(seurat_obj@reductions)) {
      seurat_obj <- RunUMAP(seurat_obj, dims = 1:50, reduction = "pca",
                            assay = assay, verbose = TRUE,
                            umap.method = "uwot", n.neighbors = 200L,
                            min.dist = 0.001, seed.use = 1L)
    }

    random_group_cluster_cols <- grep("random_group", colnames(seurat_obj@meta.data),
                                      ignore.case = TRUE, value = TRUE)

    cluster_types <- c(
      paste0("seurat_cluster_", dataset_lower),
      paste0("kmeans_cluster_", dataset_lower, "_tissue"),
      paste0("kmeans_cluster_", dataset_lower, "_sra"),
      paste0("kmeans_cluster_", dataset_lower, "_approach"),
      paste0("hierarchical_cluster_bdm_ph_", dataset_lower, "_tissue"),
      paste0("hierarchical_cluster_bdm_ph_", dataset_lower, "_sra"),
      paste0("hierarchical_cluster_bdm_ph_", dataset_lower, "_approach"),
      paste0("hierarchical_cluster_sdm_ph_", dataset_lower, "_tissue"),
      paste0("hierarchical_cluster_sdm_ph_", dataset_lower, "_sra"),
      paste0("hierarchical_cluster_sdm_ph_", dataset_lower, "_approach"),
      paste0("hierarchical_cluster_landscape_ph_", dataset_lower, "_tissue"),
      paste0("hierarchical_cluster_landscape_ph_", dataset_lower, "_sra"),
      paste0("hierarchical_cluster_landscape_ph_", dataset_lower, "_approach"),
      paste0("spectral_cluster_bdm_", dataset_lower, "_tissue"),
      paste0("spectral_cluster_bdm_", dataset_lower, "_sra"),
      paste0("spectral_cluster_bdm_", dataset_lower, "_approach"),
      paste0("spectral_cluster_sdm_", dataset_lower, "_tissue"),
      paste0("spectral_cluster_sdm_", dataset_lower, "_sra"),
      paste0("spectral_cluster_sdm_", dataset_lower, "_approach"),
      paste0("spectral_cluster_landscape_", dataset_lower, "_tissue"),
      paste0("spectral_cluster_landscape_", dataset_lower, "_sra"),
      paste0("spectral_cluster_landscape_", dataset_lower, "_approach"),
      random_group_cluster_cols,
      Tissue_col, SRA_col, Approach_col
    )

    filtered_cluster_types <- c()
    boot_base_seen <- c()
    for (clust in cluster_types) {
      if (grepl("random_group_bootstrap", clust, ignore.case = TRUE)) {
        base <- strsplit(clust, "_[Rr]andom_[Gg]roup_bootstrap_")[[1]][1]
        if (!(base %in% boot_base_seen)) {
          filtered_cluster_types <- c(filtered_cluster_types, clust)
          boot_base_seen <- c(boot_base_seen, base)
        }
      } else {
        filtered_cluster_types <- c(filtered_cluster_types, clust)
      }
    }
    cluster_types <- filtered_cluster_types

    pdf(file.path(plots_folder, paste0("UMAP_Plots_", dataset_name, "_All_Clusters.pdf")))
    for (cluster_col in cluster_types) {
      if (cluster_col %in% colnames(seurat_obj@meta.data)) {
        p <- DimPlot(seurat_obj, reduction = "umap", group.by = cluster_col) +
          ggtitle(paste("UMAP -", dataset_name, "(", cluster_col, ")")) +
          xlab("UMAP 1") + ylab("UMAP 2") +
          theme_minimal() +
          theme(legend.title = element_text(size = 12),
                legend.text = element_text(size = 10),
                plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
        grid::grid.draw(p)
      }
    }
    dev.off()
  }

  if (run_sample_level_heatmap && !is.null(metadata) && !is.null(bdm_matrix)) {
    orig_idents_in_seurat <- unique(seurat_obj@meta.data[[SRA_col]])
    metadata_filtered <- metadata %>% dplyr::filter(.data[[SRA_col]] %in% orig_idents_in_seurat)
    ordered_orig_idents <- as.character(metadata_filtered[[SRA_col]])
    dimnames(bdm_matrix) <- list(ordered_orig_idents, ordered_orig_idents)
    if (!is.null(sdm_matrix)) dimnames(sdm_matrix) <- list(ordered_orig_idents, ordered_orig_idents)
    if (!is.null(landscape_matrix)) dimnames(landscape_matrix) <- list(ordered_orig_idents, ordered_orig_idents)

    hc_tissue_bdm <- if (!is.null(k_tissue) && !is.null(bdm_matrix))
      perform_hierarchical_clustering_ph(bdm_matrix, k_tissue)$tree else NULL
    hc_sra_bdm <- if (!is.null(k_sra) && !is.null(bdm_matrix))
      perform_hierarchical_clustering_ph(bdm_matrix, k_sra)$tree else NULL
    hc_approach_bdm <- if (!is.null(k_approach) && !is.null(bdm_matrix))
      perform_hierarchical_clustering_ph(bdm_matrix, k_approach)$tree else NULL

    hc_tissue_sdm <- if (!is.null(k_tissue) && !is.null(sdm_matrix))
      perform_hierarchical_clustering_ph(sdm_matrix, k_tissue)$tree else NULL
    hc_sra_sdm <- if (!is.null(k_sra) && !is.null(sdm_matrix))
      perform_hierarchical_clustering_ph(sdm_matrix, k_sra)$tree else NULL
    hc_approach_sdm <- if (!is.null(k_approach) && !is.null(sdm_matrix))
      perform_hierarchical_clustering_ph(sdm_matrix, k_approach)$tree else NULL

    hc_tissue_landscape <- if (!is.null(k_tissue) && !is.null(landscape_matrix))
      perform_hierarchical_clustering_ph(landscape_matrix, k_tissue)$tree else NULL
    hc_sra_landscape <- if (!is.null(k_sra) && !is.null(landscape_matrix))
      perform_hierarchical_clustering_ph(landscape_matrix, k_sra)$tree else NULL
    hc_approach_landscape <- if (!is.null(k_approach) && !is.null(landscape_matrix))
      perform_hierarchical_clustering_ph(landscape_matrix, k_approach)$tree else NULL

    generate_heatmaps(paste0(dataset_name, "_BDM"), metadata_filtered, seurat_obj,
                      bdm_matrix, plots_folder, hc_tissue_bdm, hc_sra_bdm,
                      hc_approach_bdm, k_tissue, k_sra, k_approach,
                      SRA_col = SRA_col, Tissue_col = Tissue_col,
                      Approach_col = Approach_col)

    if (!is.null(sdm_matrix)) {
      generate_heatmaps(paste0(dataset_name, "_SDM"), metadata_filtered, seurat_obj,
                        sdm_matrix, plots_folder, hc_tissue_sdm, hc_sra_sdm,
                        hc_approach_sdm, k_tissue, k_sra, k_approach,
                        SRA_col = SRA_col, Tissue_col = Tissue_col,
                        Approach_col = Approach_col)
    }

    if (!is.null(landscape_matrix)) {
      generate_heatmaps(paste0(dataset_name, "_Landscape"), metadata_filtered, seurat_obj,
                        landscape_matrix, plots_folder, hc_tissue_landscape,
                        hc_sra_landscape, hc_approach_landscape,
                        k_tissue, k_sra, k_approach,
                        SRA_col = SRA_col, Tissue_col = Tissue_col,
                        Approach_col = Approach_col)
    }
  }

  seurat_obj
}


# Function to plot persistence diagrams, barcodes, and landscapes
plot_persistence <- function(pd, landscape = NULL, grid = NULL, output_file_base, plot_title) {
  tryCatch({
    if (!is.matrix(pd) || ncol(pd) != 3) {
      stop("PD is not a valid matrix with three columns (dimension, birth, death).")
    }

    # Convert persistence diagram to data frame
    pd_df <- data.frame(
      Dimension = pd[, "dimension"],
      Birth = pd[, "birth"],
      Death = pd[, "death"]
    )

    # Determine plot limits
    max_value <- max(c(pd_df$Birth, pd_df$Death), na.rm = TRUE)

    # 1. Persistence Diagram Plot
    pd_plot <- ggplot(pd_df, aes(x = Birth, y = Death, color = as.factor(Dimension))) +
      geom_point(size = 3, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
      scale_x_continuous(limits = c(0, max_value), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, max_value), expand = c(0, 0)) +
      labs(
        title = paste("Persistence Diagram -", plot_title),
        x = "Birth Time",
        y = "Death Time",
        color = "Dimension"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        aspect.ratio = 1
      )

    # Save persistence diagram
    ggsave(paste0(output_file_base, "_diagram.png"), pd_plot, width = 8, height = 6, dpi = 300)

    # 2. Persistence Barcode Plot
    barcode_plot <- ggplot(pd_df, aes(y = as.factor(Dimension), x = Birth, xend = Death, color = as.factor(Dimension))) +
      geom_segment(aes(xend = Death, yend = as.factor(Dimension)), size = 2) +
      scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
      labs(
        title = paste("Persistence Barcode -", plot_title),
        x = "Time",
        y = "Dimension",
        color = "Dimension"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top"
      )

    # Save barcode
    ggsave(paste0(output_file_base, "_barcode.png"), barcode_plot, width = 8, height = 4, dpi = 300)

    # 3. Persistence Landscape Plot (Optional)
    if (!is.null(landscape) && !is.null(grid)) {
      dim0_df <- as.data.frame(landscape$dim0)
      dim0_df$t <- grid
      df0 <- reshape2::melt(dim0_df, id.vars = "t", variable.name = "Level", value.name = "Value")
      df0$Dimension <- "0"

      dim1_df <- as.data.frame(landscape$dim1)
      dim1_df$t <- grid
      df1 <- reshape2::melt(dim1_df, id.vars = "t", variable.name = "Level", value.name = "Value")
      df1$Dimension <- "1"

      df <- rbind(df0, df1)
      df$Level <- as.factor(df$Level)

      landscape_plot <- ggplot(df, aes(x = t, y = Value, color = Level)) +
        geom_line() +
        facet_wrap(~Dimension, scales = "free_y") +
        labs(title = paste("Persistence Landscape -", plot_title),
             x = "t", y = "Landscape Value") +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "top")

      # Save landscape plots
      ggsave(paste0(output_file_base, "_landscape.png"), landscape_plot, width = 8, height = 6, dpi = 300)
      ggsave(paste0(output_file_base, "_landscape.svg"), landscape_plot, width = 8, height = 6)
      ggsave(paste0(output_file_base, "_landscape.pdf"), landscape_plot, width = 8, height = 6)

      log_message(paste("Saved landscape plot to", output_file_base))
    }

  }, error = function(e) {
    log_message(paste("Error in plotting persistence data:", e$message))
  })
}