#' Compare clustering results across methods.
#'
#' This function assumes clustering has already been performed on the
#' provided Seurat objects.  It collects the clustering columns and
#' produces comparative metrics and plots.
run_cluster_comparison <- function(data_iterations, results_folder,
                                   dataset_tag = "dataset",
                                   run_comparative_metrics = TRUE,
                                   include_silhouette = FALSE,
                                   SRA_col = "SRA",
                                   verbose = TRUE) {
  # Ensure required packages are loaded

  # Helper functions included in the package

  log_message <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", Sys.time(), msg))
  }

  normalized_bar_plots_list <- list()
  for (iter in data_iterations) {
    results <- enhanced_cluster_comparison_with_pvals(
      seurat_obj = iter$seurat_obj,
      dataset_name = iter$name,
      plots_folder = file.path(results_folder, "plots/withinIterationClusterComparison"),
      run_comparative_metrics = run_comparative_metrics,
      verbose = verbose,
      SRA_col = SRA_col,
      include_silhouette = include_silhouette
    )

    if (!is.null(results) && !is.null(results$normalized_bar_plot)) {
      normalized_bar_plots_list[[iter$name]] <- results$normalized_bar_plot
      log_message(sprintf("Collected normalized bar plot for %s.", iter$name))
    } else {
      log_message(sprintf(
        "Could not collect normalized bar plot for %s. Results might be NULL or plot creation failed.",
        iter$name
      ))
    }
  }

  num_collected_plots <- length(normalized_bar_plots_list)
  log_message(sprintf(
    "Attempting to combine %d normalized bar plots into a grid with combined legend.",
    num_collected_plots
  ))

  
#==============================================================================
# START: Combined Plot Generation Code
#==============================================================================

if (num_collected_plots > 0) {

  log_message(sprintf("Beginning combined plot generation for %d plots.", num_collected_plots))

  # --- 1. Setup File Paths ---
  plots_folder <- file.path(results_folder, "plots/withinIterationClusterComparison")
  combined_plot_basename <- "Combined_Normalized_Metrics_Aligned_Colored_Publication"
  combined_plot_png <- file.path(plots_folder, paste0(combined_plot_basename, ".png"))
  combined_plot_pdf <- file.path(plots_folder, paste0(combined_plot_basename, ".pdf"))

  # --- 2. Define Color Palettes and Global Settings ---
  bold_colorblind_palette <- c(
    "Purity"  = "#CCBB44",
    "ARI"     = "#4477AA",
    "NMI"     = "#EE6677",
    "VI"      = "#66CCEE",
    "Jaccard" = "#228833"
  )
  metric_order <- names(bold_colorblind_palette)

  ref_colors <- c(
    "Sample"   = "#D95F02",
    "Tissue"   = "#7570B3",
    "Approach" = "#1B9E77"
  )

  # --- 3. Calculate Global Axis Range for Alignment ---
  log_message("Calculating global axis range to align plot origins.")
  all_norm_values <- unlist(lapply(normalized_bar_plots_list, function(p) p$data$Normalized))
  all_norm_values <- all_norm_values[is.finite(all_norm_values)]
  global_min <- min(all_norm_values, na.rm = TRUE)
  global_max <- max(all_norm_values, na.rm = TRUE)
  axis_buffer <- (global_max - global_min) * 0.05
  global_range <- c(global_min - axis_buffer, global_max + axis_buffer)

  # --- 4. Modify Individual Plots ---
  log_message("--- Starting individual plot modification loop. ---")
  plots_for_grid <- vector("list", length = num_collected_plots)

  all_iteration_suffixes <- c("raw", "sct_individual", "sct_whole", "integrated")
  pattern_to_find_suffixes <- paste0("_", all_iteration_suffixes, collapse = "|")
  group_order <- c("Approach", "Tissue", "Sample")

  for (i in seq_len(num_collected_plots)) {
    plot_name <- names(normalized_bar_plots_list)[i]
    suffix <- if (nzchar(dataset_tag)) paste0("_", dataset_tag, "$") else "$"
    clean_plot_name <- gsub("_", " ", gsub(suffix, "", plot_name))
    p <- normalized_bar_plots_list[[i]]

    grouping_vars <- setdiff(names(p$data), c("Metric", "Normalized"))

    plot_data <- p$data %>%
      tidyr::complete(nesting(!!!rlang::syms(grouping_vars)), Metric = metric_order, fill = list(Normalized = 0)) %>%
      ungroup()

    plot_data$MethodRef <- factor(plot_data$MethodRef, levels = unique(plot_data$MethodRef))
    plot_data$Metric <- factor(plot_data$Metric, levels = rev(metric_order))

    axis_label_df <- plot_data %>%
      distinct(MethodRef) %>%
      mutate(
        CleanedMethodRef = gsub("\\n.*", "", MethodRef),
        DisplayGroup = case_when(
          grepl("_approach$", CleanedMethodRef) ~ "Approach",
          grepl("_tissue$", CleanedMethodRef) ~ "Tissue",
          grepl("_sra$", CleanedMethodRef) ~ "Sample"
        )
      ) %>%
      filter(!is.na(DisplayGroup)) %>%
      mutate(DisplayGroup = factor(DisplayGroup, levels = group_order)) %>%
      arrange(DisplayGroup, CleanedMethodRef) %>%
      mutate(
        ShortLabel = gsub("_cluster", "", gsub(paste0("(", pattern_to_find_suffixes, ").*"), "",
                              gsub("_sra$|_tissue$|_approach$", "", CleanedMethodRef))),
        color = ref_colors[as.character(DisplayGroup)],
        HtmlLabel = glue::glue("<b style='color:{color};'>{ShortLabel}</b>")
      )

    final_axis_list <- list()
    if (nrow(axis_label_df) > 0) {
      final_axis_list[[1]] <- axis_label_df[1, ]
      if (nrow(axis_label_df) > 1) {
        for (j in 2:nrow(axis_label_df)) {
          if (axis_label_df$DisplayGroup[j] != axis_label_df$DisplayGroup[j - 1]) {
            spacer_row <- data.frame(MethodRef = paste0("spacer", j), HtmlLabel = " ", stringsAsFactors = FALSE)
            final_axis_list[[length(final_axis_list) + 1]] <- spacer_row
          }
          final_axis_list[[length(final_axis_list) + 1]] <- axis_label_df[j, ]
        }
      }
    }

    final_axis_df <- dplyr::bind_rows(final_axis_list) %>% mutate(y_pos = seq_len(n()))

    divider_data <- final_axis_df %>%
      mutate(is_spacer = grepl("spacer", MethodRef)) %>%
      mutate(prev_is_spacer = dplyr::lag(is_spacer, default = TRUE)) %>%
      filter(!is_spacer & !prev_is_spacer) %>%
      mutate(x_intercept = y_pos - 0.5)

    bracket_data <- final_axis_df %>%
      tidyr::fill(DisplayGroup, color, .direction = "down") %>%
      filter(!is.na(DisplayGroup)) %>%
      group_by(DisplayGroup, color) %>%
      summarise(y_start = min(y_pos) - 0.45, y_end = max(y_pos) + 0.45, y_mid = (min(y_pos) + max(y_pos)) / 2, .groups = "drop")

    x_pos_spine <- global_max + (axis_buffer * 0.8)
    x_pos_ticks <- x_pos_spine - (axis_buffer * 0.4)
    x_pos_label <- x_pos_spine + (axis_buffer * 1.8)

    combined_palette <- c(bold_colorblind_palette, ref_colors)

    final_plot <- ggplot(plot_data, aes(x = MethodRef, y = Normalized, fill = Metric)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
      geom_vline(data = divider_data, aes(xintercept = x_intercept), color = "grey70", linewidth = 0.4) +
      geom_segment(data = bracket_data, aes(x = y_start, xend = y_end, y = x_pos_spine, yend = x_pos_spine, color = DisplayGroup), inherit.aes = FALSE, size = 0.7) +
      geom_segment(data = bracket_data, aes(x = y_start, xend = y_start, y = x_pos_spine, yend = x_pos_ticks, color = DisplayGroup), inherit.aes = FALSE, size = 0.7) +
      geom_segment(data = bracket_data, aes(x = y_end, xend = y_end, y = x_pos_spine, yend = x_pos_ticks, color = DisplayGroup), inherit.aes = FALSE, size = 0.7) +
      geom_text(data = bracket_data, aes(x = y_mid, y = x_pos_label, label = DisplayGroup, color = DisplayGroup), inherit.aes = FALSE, angle = 90, vjust = 0.5, hjust = 0.5, size = 3.5, fontface = "bold.italic") +
      coord_flip(clip = "off") +
      scale_y_continuous(limits = c(global_range[1], x_pos_label + axis_buffer), expand = expansion(mult = c(0.05, 0.1))) +
      scale_x_discrete(limits = final_axis_df$MethodRef, labels = setNames(final_axis_df$HtmlLabel, final_axis_df$MethodRef), drop = FALSE) +
      scale_fill_manual(values = bold_colorblind_palette, name = "Metric", limits = metric_order) +
      scale_color_manual(values = combined_palette, guide = "none") +
      theme_minimal_grid(font_size = 10) +
      theme(
        panel.background    = element_rect(fill = "white", colour = NA),
        plot.background     = element_rect(fill = "white", colour = NA),
        axis.text.y         = element_markdown(size = 7, hjust = 0),
        axis.text.x         = element_text(size = 8),
        axis.title.x        = element_text(size = 10, face = "bold"),
        axis.title.y        = element_blank(),
        axis.ticks.y        = element_blank(),
        plot.margin         = margin(t = 15, r = 25, b = 5, l = 5),
        plot.title          = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor.y  = element_blank()
      ) +
      labs(title = clean_plot_name, x = NULL, y = "Normalized Score (Improvement Over Random)")

    if (i < num_collected_plots) {
      final_plot <- final_plot +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
    }

    plots_for_grid[[i]] <- final_plot
  }

  log_message("\n--- Arranging and saving final plot. ---")
  plots_arranged <- Reduce(`/`, plots_for_grid)
  final_combined_plot <- plots_arranged +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = 'collect') &
    theme(
      legend.position     = 'top',
      legend.direction    = 'horizontal',
      legend.justification = 'left'
    )

  plot_width  <- 8.5
  plot_height <- 11

  log_message(sprintf("Saving final plot object of class '%s' to %s and %s", class(final_combined_plot)[1], combined_plot_png, combined_plot_pdf))

  ggsave(filename  = combined_plot_png, plot = final_combined_plot, width = plot_width, height = plot_height, units = "in", dpi = 300, limitsize = FALSE, type = "cairo")
  ggsave(filename  = combined_plot_pdf, plot = final_combined_plot, width = plot_width, height = plot_height, units = "in", limitsize = FALSE)

  log_message("--- Successfully generated and saved the final plot. ---")

} else {
  log_message("No normalized bar plots were successfully collected to combine.")
}

  invisible(normalized_bar_plots_list)
}
