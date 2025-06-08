# Modular_Analysis_and_PostProcessing.R
# Focus: Runs clustering, visualization, and statistical comparisons using pre-computed inputs.

# ---------------------------
# Libraries (Keep all necessary ones)
# ---------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
# ... (include all libraries needed for the functions below)
library(mclust)
library(aricode)
library(RColorBrewer)
library(gridExtra)
library(cluster)
library(parallel)
library(foreach)
library(doParallel)
library(reshape2)
library(TDA)
library(TDAstats)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dendextend)
library(digest)
library(transport)
library(kernlab) # Added for spectral clustering if used
library(svglite) # Added for saving plots

# ---------------------------
# Source Helper Functions (Optional)
# ---------------------------
# If you prefer to keep analysis functions separate:
# source("PH_Functions.R") 
# source("analysis_comparison_functions.R") # Example for Betti/Cluster compare funcs
# Otherwise, keep the function definitions below.

# ---------------------------
# Configuration & Options (Keep relevant ones)
# ---------------------------
results_folder <- "results_march25" # Or load from a config file/environment
plots_folder <- file.path(results_folder, "plots")
# ... (Define other necessary output subdirectories)

# Ensure directories exist (Can be a separate setup function)
setup_directories <- function(results_folder) {
  plots_folder <- file.path(results_folder, "plots")
  plot_subfolders <- list(
    cluster_size = file.path(plots_folder, "Cluster_Size_Distribution"),
    betti_curves = file.path(plots_folder, "Betti_Curves"),
    persistence_diagrams = file.path(plots_folder, "Persistence_Diagrams"),
    betti_plots_comparisons = file.path(plots_folder, "betti_plots"), # For compute_and_compare... output
    within_iteration_comparison = file.path(plots_folder, "withinIterationClusterComparison"), # For enhanced_cluster_comparison... output
    heatmaps = file.path(plots_folder, "heatmaps"), # For generate_heatmaps output
    umap_plots = file.path(plots_folder, "umap_plots") # Define if needed
  )
  metric_subfolders <- list(
      ari_nmi = file.path(results_folder, "metrics", "ARI_NMI"),
      purity = file.path(results_folder, "metrics", "Purity"),
      silhouette = file.path(results_folder, "metrics", "Silhouette"),
      chi_square = file.path(results_folder, "metrics", "Chi_Square"),
      betti_stats = file.path(results_folder, "metrics", "Betti_Statistics")
  )

  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(results_folder, "seurat_objects"), recursive = TRUE, showWarnings = FALSE) # For saving Seurat objects
  lapply(plot_subfolders, dir.create, recursive = TRUE, showWarnings = FALSE)
  lapply(metric_subfolders, dir.create, recursive = TRUE, showWarnings = FALSE)

  # Return paths needed later
  return(list(
    results_folder = results_folder,
    plots_folder = plots_folder,
    plot_folders = plot_subfolders, # Pass this to functions needing specific plot dirs
    metric_folders = metric_subfolders # Pass this if needed
  ))
}

config_paths <- setup_directories(results_folder)
plots_folder <- config_paths$plots_folder # Main plots folder accessible
plot_folders <- config_paths$plot_folders # List of specific plot folders

log_message <- function(message, task = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  prefix <- if (!is.null(task)) sprintf("[%s] Task: %s -", timestamp, task) else sprintf("[%s]", timestamp)
  cat(paste(prefix, message, "\n"))
}

# --- Run Options ---
run_standard_seurat_clustering <- TRUE
run_kmeans_clustering <- TRUE
run_hierarchical_ph_clustering <- TRUE
run_spectral_clustering <- TRUE
run_visualizations <- TRUE # UMAP plots
run_sample_level_heatmap <- TRUE
run_comparative_statistical_analysis_ph <- TRUE # Betti comparisons
run_cluster_comparison_metrics <- TRUE # ARI, NMI etc. comparisons
save_pd_plots <- TRUE # Persistence diagrams/landscapes from input

num_cores <- 12 # Or load from config

# ---------------------------
# Analysis Function Definitions (Keep or Source)
# ---------------------------
# Keep all the function definitions from your script here, *EXCEPT* for:
# - compute_and_save_distance_matrices()
# - ComputePersistenceLandscapes() 
# - compute_and_save_landscape_matrices()
# - process_iteration_calculate_matrices()
# - validate_matrix()
# - validate_landscape_list()
# - updateApproach()
# - assignRandomGroup() 

# Make sure the remaining functions are defined here or sourced from e.g., PH_Functions.R
# Example: Ensure these are defined:
# perform_standard_seurat_clustering() 
# perform_kmeans_clustering()
# perform_hierarchical_clustering_ph()
# assign_ph_clusters()
# perform_spectral_clustering() # Definition was included in the old script
# generate_heatmaps() # Definition was included
# compute_and_compare_betti_curves() # Definition was included
# enhanced_cluster_comparison_with_pvals() # Definition was included
# plot_persistence() # Definition was included
# plot_landscape() # Definition was included
# generate_heatmaps_for_all_statistics() # Definition was included
# ... and any other helper functions they rely on ...


# ---------------------------
# Data Iteration Setup (Define where to find PRE-PROCESSED data)
# ---------------------------
# This list should now point to the outputs of your preprocessing script.
# The Seurat objects should already have necessary metadata (Tissue, Approach, SRA_col, Random_Group etc.)
# Matrices and lists should already be aligned to the Seurat object identifiers.

data_iterations <- list(
  list(
    name = "Raw", # Keep the identifier for the analysis type
# Path to the Seurat object *after* preprocessing/metadata merging
    seurat_obj_path = file.path(results_folder, "preprocessed", "raw_seurat_object_preprocessed.rds"),
    assay = "RNA", # Assay used for analysis within this object
# Paths to *aligned* matrices/lists from preprocessing
    bdm_matrix_aligned = file.path(results_folder, "preprocessed", "raw_bdm_matrix_aligned.rds"),
    sdm_matrix_aligned = file.path(results_folder, "preprocessed", "raw_sdm_matrix_aligned.rds"),
    pd_list_aligned = file.path(results_folder, "preprocessed", "raw_pd_list_aligned.rds"),
    landscape_matrix_aligned = file.path(results_folder, "preprocessed", "raw_landscape_matrix_aligned.rds"),
    landscape_list_aligned = file.path(results_folder, "preprocessed", "raw_landscape_list_aligned.rds"),
    variable_features = NULL # Path if needed by perform_standard_seurat_clustering
  ),
   list(
    name = "SCT_Individual",
    seurat_obj_path = file.path(results_folder, "preprocessed", "sct_individual_seurat_object_preprocessed.rds"),
    assay = "SCT_Ind",
    bdm_matrix_aligned = file.path(results_folder, "preprocessed", "sctInd_bdm_matrix_aligned.rds"),
    sdm_matrix_aligned = file.path(results_folder, "preprocessed", "sctInd_sdm_matrix_aligned.rds"),
    pd_list_aligned = file.path(results_folder, "preprocessed", "sctInd_pd_list_aligned.rds"),
    landscape_matrix_aligned = file.path(results_folder, "preprocessed", "sctInd_landscape_matrix_aligned.rds"),
    landscape_list_aligned = file.path(results_folder, "preprocessed", "sctInd_landscape_list_aligned.rds"),
     variable_features = NULL
  ),
#SCT Whole

  list(
    name = "Integrated",
     seurat_obj_path = file.path(results_folder, "preprocessed", "integrated_seurat_object_preprocessed.rds"),
    assay = "integrated",
    bdm_matrix_aligned = file.path(results_folder, "preprocessed", "integrated_bdm_matrix_aligned.rds"),
    sdm_matrix_aligned = file.path(results_folder, "preprocessed", "integrated_sdm_matrix_aligned.rds"),
    pd_list_aligned = file.path(results_folder, "preprocessed", "integrated_pd_list_aligned.rds"),
    landscape_matrix_aligned = file.path(results_folder, "preprocessed", "integrated_landscape_matrix_aligned.rds"),
    landscape_list_aligned = file.path(results_folder, "preprocessed", "integrated_landscape_list_aligned.rds"),
    variable_features = file.path(results_folder, "preprocessed", "integration_features.rds") # Example path if needed
  )
# Add other iterations similarly...
)

# ---------------------------
# Main Analysis Loop (Focus on using pre-processed data)
# ---------------------------
log_message("Starting Main Analysis Loop...")

for (data_iter in data_iterations) {
  dataset_name <- data_iter$name
  assay <- data_iter$assay
  log_message(paste("Starting Analysis for:", dataset_name))

  # --- Load Pre-processed Seurat Object ---
  if (!file.exists(data_iter$seurat_obj_path)) {
    log_message(paste("ERROR: Preprocessed Seurat object file not found:", data_iter$seurat_obj_path, "- Skipping iteration."))
    next
  }
  seurat_obj <- readRDS(data_iter$seurat_obj_path)
  log_message(paste("Loaded Preprocessed Seurat Object:", data_iter$seurat_obj_path))

  # --- Load Pre-aligned Matrices & Lists ---
  log_message("Loading pre-aligned matrices and lists...")
  bdm_matrix <- tryCatch(as.matrix(readRDS(data_iter$bdm_matrix_aligned)), error = function(e) { log_message(paste("Error loading BDM:", e$message)); NULL })
  sdm_matrix <- tryCatch(as.matrix(readRDS(data_iter$sdm_matrix_aligned)), error = function(e) { log_message(paste("Error loading SDM:", e$message)); NULL })
  pd_list <- tryCatch(readRDS(data_iter$pd_list_aligned), error = function(e) { log_message(paste("Error loading PD list:", e$message)); NULL })
  landscape_matrix <- tryCatch(as.matrix(readRDS(data_iter$landscape_matrix_aligned)), error = function(e) { log_message(paste("Error loading Landscape matrix:", e$message)); NULL })
  landscape_list <- tryCatch(readRDS(data_iter$landscape_list_aligned), error = function(e) { log_message(paste("Error loading Landscape list:", e$message)); NULL })

  # Check if essential inputs loaded correctly
  if (is.null(pd_list) || is.null(landscape_list)) {
    log_message(paste("ERROR: PD list or Landscape list failed to load for", dataset_name, "- Skipping iteration."))
    next
  }
  if (is.null(bdm_matrix) && (run_hierarchical_ph_clustering || run_spectral_clustering || run_sample_level_heatmap)) {
    log_message(paste("WARNING: BDM matrix failed to load, but is needed for requested analyses for", dataset_name))
    # Decide whether to skip or proceed without BDM-based analyses
  }
  # Add similar checks for sdm_matrix and landscape_matrix if they are essential inputs for enabled steps

  variable_features_path <- if ("variable_features" %in% names(data_iter)) data_iter$variable_features else NULL

  # --- Determine SRA column used during preprocessing ---
  # This should ideally be stored/passed from preprocessing, but we can infer it
  SRA_col <- if ("SRA" %in% colnames(seurat_obj@meta.data)) {
    "SRA"
  } else if ("SRA_Number" %in% colnames(seurat_obj@meta.data)) {
    "SRA_Number"
  } else {
    "orig.ident" # Fallback
  }
  log_message(paste("Using column '", SRA_col, "' as SRA identifier for analysis."))

  # --- Plot Persistence Diagrams (Optional) ---
  if (save_pd_plots) {
    # (Keep the tryCatch and mclapply block from your original script here)
    # Make sure plot_persistence and plot_landscape are defined/sourced
    # Use the plot_folders list for output path: plot_folders$persistence_diagrams
    # Example call inside mclapply:
    # plot_persistence(pd, file.path(plot_folders$persistence_diagrams, paste0("PD_plot_", name)), plot_title)
    # plot_landscape(landscape, file.path(plot_folders$persistence_diagrams, paste0("Landscape_plot_", name)), plot_title, grid = seq(0, 1, length.out = 100))
    log_message("Plotting Persistence Diagrams and Landscapes...")
    tryCatch({
      plot_dir <- plot_folders$persistence_diagrams
      # Ensure plot_persistence and plot_landscape are defined
      mclapply(seq_along(pd_list), function(i) {
        pd <- pd_list[[i]]
        landscape <- landscape_list[[i]]
        name <- names(pd_list)[i] %||% as.character(i)

        if (!is.matrix(pd) || ncol(pd) != 3) {
          log_message(paste("Skipping plot for", name, "as PD is not a valid matrix."))
          return()
        }

        # Use file.path to construct paths robustly
        pd_plot_file_base <- file.path(plot_dir, paste0("PD_", name))
        pl_plot_file_base <- file.path(plot_dir, paste0("Landscape_", name))

        plot_title_pd <- paste("Persistence Diagram -", name, "-", dataset_name)
        plot_title_pl <- paste("Persistence Landscape -", name, "-", dataset_name)

        # Assuming plot_persistence and plot_landscape save the files
        plot_persistence(pd, pd_plot_file_base, plot_title_pd)
        plot_landscape(landscape, pl_plot_file_base, plot_title_pl, grid = seq(0, 1, length.out = 100))

        log_message(paste("Saved PD/Landscape plots for:", name, "to", plot_dir))
      }, mc.cores = num_cores)
    }, error = function(e) {
      log_message(paste("Error plotting PDs/Landscapes for", dataset_name, ":", e$message))
    })
  }

  # --- Hierarchical PH Clustering ---
  if (run_hierarchical_ph_clustering && !is.null(bdm_matrix) && !is.null(sdm_matrix) && !is.null(landscape_matrix)) {
    log_message("Running Hierarchical PH Clustering...")
    # (Keep the relevant block from your original script here)
    # This block uses perform_hierarchical_clustering_ph and assign_ph_clusters
    # It needs k_tissue, k_sra, k_approach - calculate them from the loaded seurat_obj
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) k_tissue <- length(unique(seurat_obj$Tissue)) else k_tissue <- NULL
    if (SRA_col %in% colnames(seurat_obj@meta.data)) k_sra <- length(unique(seurat_obj@meta.data[[SRA_col]])) else k_sra <- NULL
    if ("Approach" %in% colnames(seurat_obj@meta.data)) k_approach <- length(unique(seurat_obj$Approach)) else k_approach <- NULL
    random_group_cols <- grep("^Random_Group", colnames(seurat_obj@meta.data), value = TRUE) # For logging

    # ... (rest of the hierarchical clustering block for BDM, SDM, Landscape) ...
    # Remember to handle cases where k_* is NULL if the column doesn't exist
    # Example for BDM Tissue:
    if (!is.null(k_tissue)) {
      cat(paste(Sys.time(), "- Running hierarchical clustering on BDM for tissues with k =", k_tissue, "\n"))
      ph_tissue_result_bdm <- try(perform_hierarchical_clustering_ph(bdm_matrix, k = k_tissue), silent = TRUE)
      if (!inherits(ph_tissue_result_bdm, "try-error")) {
        ph_clusters_tissue_bdm <- ph_tissue_result_bdm$clusters
        hc_tissue_bdm <- ph_tissue_result_bdm$tree # Store the tree if needed later (e.g., for heatmaps)
        seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_tissue_bdm, paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_tissue"))
        cat(paste(Sys.time(), "- Hierarchical clustering on BDM for tissues completed.\n"))
      } else { cat(paste(Sys.time(), "- Hierarchical clustering on BDM for tissues failed:", ph_tissue_result_bdm, ". Moving on.\n")) }
    } else { cat(paste(Sys.time(), "- Tissue column not found or has only one level. Skipping BDM tissue clustering.\n")) }
    # ... Repeat for SRA, Approach for BDM, SDM, Landscape ...

  } else {
    log_message("Skipping Hierarchical PH Clustering (flag disabled or input matrix missing).")
  }

  # --- Standard Seurat Clustering ---
  if (run_standard_seurat_clustering) {
    log_message("Running Standard Seurat Clustering...")
    # (Keep the relevant block from your original script here)
    # This uses perform_standard_seurat_clustering
    seurat_obj <- perform_standard_seurat_clustering(seurat_obj, assay = assay, variable_features_path = variable_features_path)
    # The function should handle FindNeighbors/FindClusters and update Idents.
    # We'll store the result in a uniquely named column.
    seurat_obj@meta.data[[paste0("seurat_cluster_", tolower(dataset_name))]] <- Idents(seurat_obj)
    # Optional: Remove the default 'seurat_clusters' column if it exists from FindClusters
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) seurat_obj$seurat_clusters <- NULL
  }

  # --- K-means Clustering ---
  if (run_kmeans_clustering) {
    log_message("Running K-means Clustering...")
    # (Keep the relevant block from your original script here)
    # This uses perform_kmeans_clustering
    # Needs k_tissue, k_sra, k_approach calculated from the loaded seurat_obj
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) k_tissue <- length(unique(seurat_obj$Tissue)) else k_tissue <- NULL
    if (SRA_col %in% colnames(seurat_obj@meta.data)) k_sra <- length(unique(seurat_obj@meta.data[[SRA_col]])) else k_sra <- NULL
    if ("Approach" %in% colnames(seurat_obj@meta.data)) k_approach <- length(unique(seurat_obj$Approach)) else k_approach <- NULL

    # ... (rest of the k-means block for Tissue, SRA, Approach) ...
    # Example for Tissue:
    if (!is.null(k_tissue)) {
      cat(paste(Sys.time(), "- Running K-means clustering with k =", k_tissue, "(number of tissues)\n"))
      seurat_obj <- perform_kmeans_clustering(seurat_obj, assay = assay, dims = 1:50, k = k_tissue)
      # Store result in unique column, then remove the temporary column
      seurat_obj@meta.data[[paste0("kmeans_cluster_", tolower(dataset_name), "_tissue")]] <- seurat_obj$kmeans_cluster
      seurat_obj$kmeans_cluster <- NULL
      cat(paste(Sys.time(), "- K-means tissue clustering completed.\n"))
    } else { cat(paste(Sys.time(), "- Tissue column not found or has only one level. Skipping K-means tissue clustering.\n")) }
    # ... Repeat for SRA, Approach ...

  }

  # --- Spectral Clustering ---
  if (run_spectral_clustering && !is.null(bdm_matrix) && !is.null(sdm_matrix) && !is.null(landscape_matrix)) {
    log_message("Running Spectral Clustering...")
    # (Keep the relevant block from your original script here)
    # This uses perform_spectral_clustering and assign_ph_clusters
    # Needs k_tissue, k_sra, k_approach calculated from the loaded seurat_obj
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) k_tissue <- length(unique(seurat_obj$Tissue)) else k_tissue <- NULL
    if (SRA_col %in% colnames(seurat_obj@meta.data)) k_sra <- length(unique(seurat_obj@meta.data[[SRA_col]])) else k_sra <- NULL
    if ("Approach" %in% colnames(seurat_obj@meta.data)) k_approach <- length(unique(seurat_obj$Approach)) else k_approach <- NULL

    # ... (rest of the spectral clustering block for BDM, SDM, Landscape) ...
    # Example for BDM Tissue:
    if (!is.null(k_tissue)) {
      cat(paste(Sys.time(), "- Running Spectral clustering on BDM for tissues with k =", k_tissue, "\n"))
      spectral_tissue_clusters_bdm <- try(perform_spectral_clustering(bdm_matrix, k = k_tissue), silent = TRUE)
      if (!inherits(spectral_tissue_clusters_bdm, "try-error")) {
        seurat_obj <- assign_ph_clusters(seurat_obj, spectral_tissue_clusters_bdm, paste0("spectral_cluster_bdm_", tolower(dataset_name), "_tissue"))
        cat(paste(Sys.time(), "- Spectral BDM tissue clustering completed.\n"))
      } else { cat(paste(Sys.time(), "- Spectral clustering on BDM for tissues failed:", spectral_tissue_clusters_bdm, ". Moving on.\n")) }
    } else { cat(paste(Sys.time(), "- Tissue column not found or has only one level. Skipping Spectral BDM tissue clustering.\n")) }
    # ... Repeat for SRA, Approach for BDM, SDM, Landscape ...

  } else {
    log_message("Skipping Spectral Clustering (flag disabled or input matrix missing).")
  }

  # --- Save Processed Seurat Object ---
  # Save the object *after* all clustering methods have potentially added metadata
  save_path_obj <- file.path(config_paths$results_folder, "seurat_objects", paste0(tolower(dataset_name), "_seurat_object_final.rds"))
  saveRDS(seurat_obj, save_path_obj)
  log_message(paste("Final Seurat object saved for:", dataset_name, "at", save_path_obj))

  # --- UMAP Visualizations ---
  if (run_visualizations) {
    log_message("Running UMAP Visualizations...")
    # (Keep the relevant block from your original script here)
    # Ensure RunUMAP is called if needed (it's inside the block)
    # Use the plot_folders list for output path: plot_folders$umap_plots (or just plots_folder)
    # Example path: file.path(plots_folder, paste0("UMAP_Plots_", dataset_name, "_All_Clusters.pdf")) 
    umap_plot_path <- file.path(plot_folders$umap_plots %||% plots_folder, paste0("UMAP_Plots_", dataset_name, "_All_Clusters.pdf"))
    ensure_directory(dirname(umap_plot_path))
    pdf(umap_plot_path)
    # ... (The loop for DimPlotting over cluster_types) ...
    dev.off()
    log_message(paste("UMAP plots saved to", umap_plot_path))
  }

  # --- Sample Level Heatmaps ---
  if (run_sample_level_heatmap) {
    log_message("Running Sample-Level Heatmaps...")
    heatmap_plot_folder <- plot_folders$heatmaps %||% plots_folder # Use specific or general plots folder
    ensure_directory(heatmap_plot_folder)

    # Reload the unfiltered metadata if needed by generate_heatmaps
    # Or better: pass the PRE-FILTERED metadata from the preprocessing step
    # Assuming metadata_filtered exists from a previous step (or reload it)

    # Call generate_heatmaps for BDM, SDM, Landscape if they exist
    # Need the hc_ objects from the hierarchical clustering step - make sure they were stored
    if (!is.null(bdm_matrix) && exists("hc_tissue_bdm") && exists("hc_sra_bdm") && exists("hc_approach_bdm")) {
      generate_heatmaps(
         dataset_name = paste0(dataset_name, "_BDM"),
         metadata = metadata_filtered, # Use the filtered metadata from preprocessing
         seurat_obj = seurat_obj, # Pass the current object
         bdm_matrix = bdm_matrix, # Pass the aligned matrix
         plots_folder = heatmap_plot_folder, # Specify output dir
         hc_tissue = hc_tissue_bdm, hc_sra = hc_sra_bdm, hc_approach = hc_approach_bdm, # Pass the hclust results
         k_tissue = k_tissue, k_sra = k_sra, k_approach = k_approach # Pass the k values
       )
    }
    # ... Repeat for SDM and Landscape if matrices and hclust objects exist ...

  }

  # --- Betti Curve Statistical Comparisons ---
  if (run_comparative_statistical_analysis_ph) {
    log_message("Running Betti Curve Statistical Comparisons...")
    # (Keep the relevant block from your original script here)
    # This uses compute_and_compare_betti_curves and generate_heatmaps_for_all_statistics
    # Ensure the necessary inputs (pd_list, landscape_list, seurat_obj) are available
    # The function handles saving plots/results to subdirectories based on results_folder
    # Example call for tissue comparison:
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
      results_list$tissue_comparison <- tryCatch(
          compute_and_compare_betti_curves(
            pd_list = pd_list, # Use aligned list
            landscape_list = landscape_list, # Use aligned list
            seurat_objects = list(seurat_obj), # Pass current object in a list
            group_by_col = "Tissue",
            dataset_name = dataset_name,
            results_folder = results_folder # Main results folder
      # Pass other params like grid_points, num_permutations, etc.
          ), error = function(e) { log_message(paste("Error in Betti tissue comparison:", e$message)); NULL })
    }
    # ... Repeat for SRA, Approach, Random_Group, and cluster columns ...

    # After running all comparisons for this iteration, generate the heatmaps
    if (exists("results_list") && length(results_list) > 0) {
      # Extract null summaries if computed (e.g., from random group comparison)
      overall_null_summary <- list(
             euler_effect_size = results_list$random_group_comparison_Random_Group_bootstrap_1$bootstrap_nulls$null_effect, # Example path, adjust based on actual output
             euler_ks = results_list$random_group_comparison_Random_Group_bootstrap_1$bootstrap_nulls$null_ks, # Example path
             landscape_effect_size = results_list$random_group_comparison_Random_Group_bootstrap_1$landscape_nulls$null_effect, # Example path
             landscape_ks = results_list$random_group_comparison_Random_Group_bootstrap_1$landscape_nulls$null_ks # Example path
         )

      generate_heatmaps_for_all_statistics(
            results_list = results_list,
            dataset_name = dataset_name,
            plots_folder = plots_folder, # Pass main plots folder
            null_summary = overall_null_summary # Pass computed nulls
          )
    } else {
      log_message("No Betti comparison results to generate heatmaps for.")
    }

  }

  # --- Enhanced Cluster Comparison Metrics ---
  if (run_cluster_comparison_metrics) {
    log_message("Running Enhanced Cluster Comparison Metrics...")
    # (Keep the relevant block from your original script here)
    # This uses enhanced_cluster_comparison_with_pvals
    # The function handles saving plots/results based on the plots_folder argument
    results_metrics <- enhanced_cluster_comparison_with_pvals(
      seurat_obj = seurat_obj,
      dataset_name = dataset_name,
      plots_folder = plot_folders$within_iteration_comparison, # Use specific subfolder
      run_comparative_metrics = TRUE,
      verbose = TRUE,
      SRA_col = SRA_col # Pass the determined SRA column name
    )
    # Optionally save the returned 'results_metrics' object
    saveRDS(results_metrics, file.path(config_paths$metric_folders$ari_nmi, paste0(dataset_name, "_cluster_metrics_results.rds")))
  }

  log_message(paste("Completed Analysis for:", dataset_name))
  rm(seurat_obj, bdm_matrix, sdm_matrix, pd_list, landscape_matrix, landscape_list, metadata_filtered) # Clean up memory
  gc() # Garbage collect
}

log_message("Main Analysis Loop Finished.")

# ---------------------------
# Cross-Iteration Comparisons (Optional - Keep if needed)
# ---------------------------
# The cross-iteration comparison functions (like cross_iteration_comparison_with_betti, 
# compare_all_iterations_with_mapped_methods_random_baseline_including_purity_pvals) 
# could remain here, but they would need to load the *final* saved Seurat objects 
# from the loop above. 
# Example:
# Load all final Seurat objects into a list first
# final_seurat_objects <- lapply(data_iterations, function(iter) {
#   readRDS(file.path(results_folder, "seurat_objects", paste0(tolower(iter$name), "_seurat_object_final.rds")))
# })
# names(final_seurat_objects) <- sapply(data_iterations, `[[`, "name")
# Now run comparisons using final_seurat_objects...

log_message("Script finished.")