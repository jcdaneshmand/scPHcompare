# Make sure you have cowplot and other required packages installed
library(cowplot)
library(dplyr)
library(tidyr)
library(ggplot2) # Ensure ggplot2 is loaded
library(ggtext)
library(glue)

source("cluster_comparison_FUNCTION.R")

# List to store the normalized bar plots from each iteration
normalized_bar_plots_list <- list()
for (iter in data_iterations) {
  # to run outside of loop if needed (ie, to adjust plots, debugging)
  results <- enhanced_cluster_comparison_with_pvals(
    seurat_obj = iter$seurat_obj,
    dataset_name = iter$name,
    plots_folder = paste0(results_folder, "/plots/withinIterationClusterComparison"),
    run_comparative_metrics = TRUE,
    verbose = TRUE,
    # SRA_col = "orig.ident",
    SRA_col = "SRA",
    include_silhouette = FALSE
  )

  # Store the normalized bar plot if the function returned results and the plot
  if (!is.null(results) && !is.null(results$normalized_bar_plot)) {
    normalized_bar_plots_list[[iter$name]] <- results$normalized_bar_plot
    log_message(sprintf("Collected normalized bar plot for %s.", iter$name))
  } else {
    log_message(sprintf("Could not collect normalized bar plot for %s. Results might be NULL or plot creation failed.", iter$name))
  }
}

num_collected_plots <- length(normalized_bar_plots_list)
log_message(sprintf("Attempting to combine %d normalized bar plots into a grid with combined legend.", num_collected_plots))

#==============================================================================
# START: Combined Plot Generation Code (with reordered rows by ReferenceCol)
#==============================================================================
if (num_collected_plots > 0) {

  log_message(sprintf("Beginning combined plot generation for %d plots.", num_collected_plots))

  # --- 1. Setup File Paths ---
  plots_folder <- file.path(results_folder, "plots/withinIterationClusterComparison")
  combined_plot_basename <- "Combined_Normalized_Metrics_Aligned_Colored_Publication" ## MODIFICATION: New filename
  combined_plot_png <- file.path(plots_folder, paste0(combined_plot_basename, ".png"))
  combined_plot_pdf <- file.path(plots_folder, paste0(combined_plot_basename, ".pdf"))

  # --- 2. Extract a Shared Legend (with increased font size) ---
  log_message("Extracting a shared legend from the first plot.")
  ## MODIFICATION: Increase legend font size on the plot before extracting the legend
  first_plot_for_legend <- normalized_bar_plots_list[[1]] +
    theme(legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12))
  shared_legend <- tryCatch({ cowplot::get_legend(first_plot_for_legend) }, error = function(e) { NULL })

  # --- 3. Calculate Global Axis Range for Alignment ---
  log_message("Calculating global axis range to align plot origins.")
  all_norm_values <- unlist(lapply(normalized_bar_plots_list, function(p) p$data$Normalized))
  all_norm_values <- all_norm_values[is.finite(all_norm_values)]
  global_min <- min(all_norm_values, na.rm = TRUE)
  global_max <- max(all_norm_values, na.rm = TRUE)
  axis_buffer <- (global_max - global_min) * 0.05
  global_range <- c(global_min - axis_buffer, global_max + axis_buffer)

  # --- 4. Modify Individual Plots: reorder rows, color labels, uniform axes ---
  log_message("Modifying individual plots: reordering rows by ReferenceCol, coloring labels, and setting uniform axes.")
  plots_for_grid <- vector("list", length = num_collected_plots)

  for (i in seq_len(num_collected_plots)) {
    plot_name <- names(normalized_bar_plots_list)[i]
    p <- normalized_bar_plots_list[[i]]
    plot_data <- p$data

    # ---- 4a) Reorder rows so same ReferenceCol are grouped ----
    label_color_lookup <- plot_data %>% distinct(MethodRef, ReferenceCol)
    group_order <- c("SRA", "orig.ident", "Tissue", "Approach", "default")
    ordered_levels <- label_color_lookup %>%
      mutate(ReferenceCol = factor(ReferenceCol, levels = group_order)) %>%
      arrange(ReferenceCol, MethodRef) %>%
      pull(MethodRef)
    p$data$MethodRef <- factor(p$data$MethodRef, levels = ordered_levels)

    # ---- 4b) Color-code Y-axis Labels by ReferenceCol ----
    ref_colors <- c(
      "SRA"        = "#D55E00",
      "orig.ident" = "#D55E00",
      "Tissue"     = "#0072B2",
      "Approach"   = "#009E73",
      "default"    = "black"
    )
    label_color_lookup <- p$data %>%
      distinct(MethodRef, ReferenceCol) %>%
      mutate(color = ref_colors[match(ReferenceCol, names(ref_colors))])
    original_levels <- levels(p$data$MethodRef)
    new_levels <- original_levels
    for (j in seq_len(nrow(label_color_lookup))) {
      idx <- which(original_levels == label_color_lookup$MethodRef[j])
      if (length(idx) > 0) { ## Safety check
        new_levels[idx] <- glue(
          "<b style='color:{label_color_lookup$color[j]};'>{original_levels[idx]}</b>"
        )
      }
    }
    levels(p$data$MethodRef) <- new_levels

    # ---- 4c) Apply coord_flip, uniform axis, theming ----
    ## MODIFICATION: Increased font sizes and applied coord_flip()
    p <- p +
      coord_flip() +
      scale_y_continuous(limits = global_range) +
      theme(
        legend.position = "none",
        plot.title      = element_text(face = "bold", size = 16), # Increased size
        axis.title.y    = element_blank(),
        axis.text.y     = element_markdown(size = 10), # Increased size
        axis.title.x    = element_text(size = 12), # Increased size
        axis.text.x     = element_text(size = 10)  # Increased size
      ) +
      labs(title = plot_name, y = "Normalized Value") # Ensure y-axis title is set

    # Remove X-axis title and text on all but the bottom-most plot
    if (i < num_collected_plots) {
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank() # Also remove ticks
      )
    }

    plots_for_grid[[i]] <- p
  }

  # --- 5. Arrange and Save the Final Plot ---
  log_message("Arranging plots into a final grid.")
  ## MODIFICATION: Added labels = "AUTO" and increased label_size
  plots_arranged <- cowplot::plot_grid(
    plotlist = plots_for_grid,
    ncol     = 1,
    align    = "v",
    axis     = "l",
    labels   = "AUTO",
    label_size = 18
  )

  if (!is.null(shared_legend)) {
    log_message("Combining plot stack with the shared legend.")
    ## MODIFICATION: Changed rel_widths to give legend more space
    final_combined_plot <- cowplot::plot_grid(
      plots_arranged, shared_legend,
      ncol       = 2,
      rel_widths = c(4, 1) # Changed from c(5, 1)
    )
  } else {
    final_combined_plot <- plots_arranged
  }

  ## MODIFICATION: Adjusted plot dimensions for better fit
  plot_width           <- 16 # Increased width
  plot_height_per_plot <- 6  # Increased height per plot
  total_height         <- min(num_collected_plots * plot_height_per_plot, 40)

  log_message(sprintf(
    "Saving final aligned & colored plot to %s and %s",
    combined_plot_png, combined_plot_pdf
  ))
  ggsave(
    filename  = combined_plot_png,
    plot      = final_combined_plot,
    width     = plot_width,
    height    = total_height,
    dpi       = 300,
    limitsize = FALSE
  )
  ggsave(
    filename  = combined_plot_pdf,
    plot      = final_combined_plot,
    width     = plot_width,
    height    = total_height,
    limitsize = FALSE
  )

  log_message("Successfully generated and saved the final plot.")

} else {
  log_message("No normalized bar plots were successfully collected to combine.")
}
#==============================================================================
# END: Combined Plot Generation Code
#==============================================================================
