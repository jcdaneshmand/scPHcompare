#VDT VERSION

# ---------------------------
# Libraries
# ---------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
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
library(dendextend) # For coloring dendrogram branches
library(digest)
library(transport) # For EMD
source("PH_Functions.R")
source("betti_utils.R")

#' Run post-processing modules
#'
#' Wrapper that selectively runs cluster comparison, Betti curve analysis,
#' and cross-iteration analysis based on provided flags.
run_modular_analysis <- function(ph_results,
                                 results_dir = "results",
                                 run_cluster = FALSE,
                                 run_betti = FALSE,
                                 run_cross_iteration = FALSE,
                                 ...) {
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  source("cluster_comparison.R")
  source("cross_iteration_functions.R")

  results <- list()
  data_iterations <- ph_results$data_iterations

  if (run_cluster && exists("run_cluster_comparison") && !is.null(data_iterations)) {
    results$cluster <- try(run_cluster_comparison(data_iterations,
                                                 results_folder = results_dir,
                                                 ...),
                           silent = TRUE)
  }

  if (run_betti && !is.null(data_iterations)) {
    betti_results <- list()
    for (iter in data_iterations) {
      pd_list <- readRDS(iter$pd_list)
      landscape_list <- if (!is.null(iter$landscape_list)) readRDS(iter$landscape_list) else NULL
      betti_results[[iter$name]] <- try(
        compute_and_compare_betti_curves(
          pd_list = pd_list,
          landscape_list = landscape_list,
          seurat_objects = list(iter$seurat_obj),
          group_by_col = "Tissue",
          dataset_name = iter$name,
          results_folder = results_dir,
          ...
        ),
        silent = TRUE
      )
    }
    results$betti <- betti_results
  }

  if (run_cross_iteration && exists("run_cross_iteration") && !is.null(data_iterations)) {
    results$cross_iteration <- try(run_cross_iteration(data_iterations,
                                                      results_folder = results_dir,
                                                      ...),
                                   silent = TRUE)
  }

  invisible(results)
}

# ---------------------------
# Set Options
# ---------------------------
# results_folder <- "results_march25_bonemarrow"
results_folder <- "results_march25"

plots_folder <- file.path(results_folder, "plots")
betti_plots_folder <- file.path(plots_folder, "betti_curves")

if (!dir.exists(results_folder)) dir.create(results_folder)
if (!dir.exists(plots_folder)) dir.create(plots_folder)
if (!dir.exists(betti_plots_folder)) dir.create(betti_plots_folder)

# Setting up directories for different types of outputs
ari_nmi_folder <- file.path(results_folder, "metrics", "ARI_NMI")
purity_folder <- file.path(results_folder, "metrics", "Purity")
silhouette_folder <- file.path(results_folder, "metrics", "Silhouette")
chi_square_folder <- file.path(results_folder, "metrics", "Chi_Square")
betti_stats_folder <- file.path(results_folder, "metrics", "Betti_Statistics")
plot_folders <- list(
  cluster_size = file.path(plots_folder, "Cluster_Size_Distribution"),
  betti_curves = file.path(plots_folder, "Betti_Curves"),
  persistence_diagrams = file.path(plots_folder, "Persistence_Diagrams")
)

# Create directories if they don't exist
dir.create(ari_nmi_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(purity_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(silhouette_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(chi_square_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(betti_stats_folder, recursive = TRUE, showWarnings = FALSE)
lapply(plot_folders, function(folder) dir.create(folder, recursive = TRUE, showWarnings = FALSE))

log_message <- function(message, task = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (!is.null(task)) {
    cat(sprintf("[%s] Task: %s - %s\n", timestamp, task, message))
  } else {
    cat(sprintf("[%s] %s\n", timestamp, message))
  }
}

# ---------------------------
# Run Options
# ---------------------------
run_standard_seurat_clustering <- TRUE
run_kmeans_clustering <- TRUE
run_hierarchical_ph_clustering <- TRUE
run_comparative_metrics <- TRUE
run_direct_permutation_tests <- TRUE
run_visualizations <- TRUE
run_sample_level_heatmap <- TRUE
run_comparative_statistical_analysis_ph <- TRUE
run_spectral_clustering <- TRUE
save_pd_plots <- TRUE
debug_tinydata <- FALSE

num_cores <- 12
# ---------------------------
# Labels to Test
# ---------------------------
labels_to_test <- c("Tissue", "SRA", "Approach", "orig.ident")

# ---------------------------
# Data Iteration Setup
# ---------------------------

# integrated_seurat <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/final_integrated_seurat.rds")
# integrated_seurat_bonemarrow <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/final_integrated_seurat_bonemarrow.rds")
# merged_seurat_unintegrated <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/merged_seurat_unintegrated.Rds")
# merged_seurat_unintegrated_bonemarrow <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/merged_seurat_unintegrated_bonemarrow.Rds")


if (debug_tinydata == TRUE) {
  # tiny_cells_unintegrated <- sample(Cells(merged_seurat_unintegrated), size = 1000)
  # merged_seurat_unintegrated_tiny <- subset(merged_seurat_unintegrated, cells = tiny_cells_unintegrated)
  # merged_seurat_unintegrated<- merged_seurat_unintegrated_tiny

  # tiny_cells_integrated <- sample(Cells(integrated_seurat), size = 1000)
  # integrated_seurat_tiny <- subset(integrated_seurat, cells = tiny_cells_integrated)
  # integrated_seurat <- integrated_seurat_tiny

  # tiny_cells_unintegrated <- sample(Cells(merged_seurat_unintegrated_bonemarrow), size = 1000)
  # merged_seurat_unintegrated_tiny <- subset(merged_seurat_unintegrated_bonemarrow, cells = tiny_cells_unintegrated)
  # merged_seurat_unintegrated_bonemarrow <- merged_seurat_unintegrated_tiny

  # tiny_cells_integrated <- sample(Cells(integrated_seurat_bonemarrow), size = 1000)
  # integrated_seurat_tiny <- subset(integrated_seurat_bonemarrow, cells = tiny_cells_integrated)
  # integrated_seurat_bonemarrow <- integrated_seurat_tiny
}

# Step 1: Create a lookup table for Approach based on orig.ident
# approach_lookup <- unique(merged_seurat_unintegrated@meta.data[, c("orig.ident", "Approach")])
# tissue_lookup <- unique(merged_seurat_unintegrated@meta.data[, c("orig.ident", "Tissue")])

# Step 2: Merge the lookup table with the meta.data of integrated_seurat
# integrated_seurat@meta.data <- merge(integrated_seurat@meta.data, approach_lookup, by = "orig.ident", all.x = TRUE)
# integrated_seurat@meta.data <- merge(integrated_seurat@meta.data, tissue_lookup, by = "orig.ident", all.x = TRUE)

# Step 3: Verify the addition of the Approach column
# integrated_seurat@meta.data$Tissue <- integrated_seurat@meta.data$Tissue.y
# str(integrated_seurat@meta.data)
updateApproach <- function(seurat_obj, metadata = NULL, new_value = NULL) {
  # Check if the 'Approach' column exists in the Seurat object's meta.data
  if (!"Approach" %in% colnames(seurat_obj@meta.data)) {
    stop("The Seurat object does not contain a column named 'Approach'.")
  }

  # Save the existing 'Approach' column as 'Inferred_Approach'
  seurat_obj@meta.data$Inferred_Approach <- seurat_obj@meta.data$Approach

  if (!is.null(metadata)) {
    # Construct a key in the metadata by concatenating SRA and SRS numbers.
    # Adjust the separator if your orig.ident uses a different format.
    metadata$key <- paste(metadata$`SRA Number`, metadata$`SRS Number`, sep = "_")

    # Get all unique orig.ident values in the Seurat object's meta.data
    unique_ids <- unique(seurat_obj@meta.data$orig.ident)

    # Loop over each unique identifier and update the Approach accordingly
    for (id in unique_ids) {
      if (id %in% metadata$key) {
        # Use the corresponding Approach value from the metadata (first one if there are multiple)
        approach_val <- metadata$Approach[metadata$key == id][1]
        seurat_obj@meta.data$Approach[seurat_obj@meta.data$orig.ident == id] <- approach_val
      } else {
        warning(paste("No metadata found for orig.ident:", id))
        # Optionally, use new_value if provided
        if (!is.null(new_value)) {
          seurat_obj@meta.data$Approach[seurat_obj@meta.data$orig.ident == id] <- new_value
        }
      }
    }
  } else if (!is.null(new_value)) {
    # If no metadata is provided, simply set the new_value for all cells
    seurat_obj@meta.data$Approach <- new_value
  } else {
    stop("Either metadata or new_value must be provided")
  }

  return(seurat_obj)
}

metadata <- read_csv("./data/VastlyDifferentTissues/metadata.csv")
# metadata <- read_csv("./data/GSE120221/metadata.csv")

# merged_seurat_unintegrated <- updateApproach(merged_seurat_unintegrated, metadata = metadata)
# integrated_seurat <- updateApproach(integrated_seurat, metadata = metadata)

# merged_seurat_unintegrated_bonemarrow <- updateApproach(merged_seurat_unintegrated_bonemarrow, new_value = 'scRNA-seq')
# integrated_seurat_bonemarrow <- updateApproach(integrated_seurat_bonemarrow, new_value = 'scRNA-seq')

# here add function to get correct approach into metadata based on assoc with orig.ident

assignRandomGroup <- function(seurat_obj, k, new_col_name = "Random_Group", seed = 123, num_bootstraps = 1) {
  # Check if the 'orig.ident' column exists in meta.data
  if (!"orig.ident" %in% colnames(seurat_obj@meta.data)) {
    stop("The Seurat object does not contain a column named 'orig.ident'.")
  }

  # Extract unique orig.ident values
  unique_ids <- unique(seurat_obj@meta.data$orig.ident)

  # Loop through the desired number of bootstrap iterations
  for (i in seq_len(num_bootstraps)) {
    # Set a reproducible seed for this bootstrap iteration (offset by iteration)
    set.seed(seed + i)

    # Assign each unique orig.ident to a random group from 1 to k in a balanced way
    random_groups <- sample(rep(1:k, length.out = length(unique_ids)))

    # Create a mapping from orig.ident to its assigned random group
    group_mapping <- setNames(random_groups, unique_ids)

    # Define the column name; if multiple bootstraps, add a suffix for each iteration
    current_col_name <- if (num_bootstraps == 1) new_col_name else paste0(new_col_name, "_bootstrap_", i)

    # Add the new group assignment column to the metadata
    seurat_obj@meta.data[[current_col_name]] <- group_mapping[seurat_obj@meta.data$orig.ident]
  }

  return(seurat_obj)
}

# Suppose 'seurat_obj' is your Seurat object and you want to split into 3 groups
# seurat_obj <- assignRandomGroup(seurat_obj, k = 3, new_col_name = "Random_Group")

# merged_seurat_unintegrated <- assignRandomGroup(merged_seurat_unintegrated, k = 5, new_col_name = "Random_Group", num_bootstraps = 10)
# integrated_seurat <- assignRandomGroup(integrated_seurat, k = 5, new_col_name = "Random_Group", num_bootstraps = 10)

# merged_seurat_unintegrated_bonemarrow <- assignRandomGroup(merged_seurat_unintegrated_bonemarrow, k = 5, new_col_name = "Random_Group", num_bootstraps = 10)
# integrated_seurat_bonemarrow <- assignRandomGroup(integrated_seurat_bonemarrow, k = 5, new_col_name = "Random_Group", num_bootstraps = 10)


data_iterations <- list(
  list(
    name = "Raw",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "raw_seurat_object.rds")),
# seurat_obj = merged_seurat_unintegrated,
    assay = "RNA",
    bdm_matrix = "BDM_unintegrated_Raw.rds",
    sdm_matrix = "SDM_unintegrated_Raw.rds",
    pd_list = "PD_list_after_retries_unintegrated_Raw.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_Raw.Rds",
    landscape_list = "landscape_list_Raw.Rds",
    expr_list = readRDS("expr_list_raw.Rds")
  )
,
# Uncomment SCT_Individual iteration if needed
  list(
    name = "SCT_Individual",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_individual_seurat_object.rds")),
# seurat_obj = merged_seurat_unintegrated,
    assay = "SCT_Ind",
    bdm_matrix = "BDM_unintegrated_sctInd.rds",
    sdm_matrix = "SDM_unintegrated_sctInd.rds",
    pd_list = "PD_list_after_retries_unintegrated_sctInd.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctInd.Rds",
    landscape_list = "landscape_list_sctInd.Rds",
    expr_list = readRDS("expr_list_sctInd.Rds")
  ),
# list(
#   name = "SCT_Whole",
#   seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_whole_seurat_object.rds")),
#   # seurat_obj = merged_seurat_unintegrated,
#   assay = "SCT",
#   bdm_matrix = "BDM_unintegrated_scTransformed_whole.rds",
#   sdm_matrix = "SDM_unintegrated_scTransformed_whole.rds",
#   pd_list = "PD_list_after_retries_unintegrated_sctWhole.rds",
#   landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctWhole.Rds",
#   landscape_list = "landscape_list_sctWhole.Rds",
#   expr_list = readRDS("expr_list_sctwhole.Rds")
# ),
  list(
    name = "Integrated",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "integrated_seurat_object.rds")),
# seurat_obj = integrated_seurat,
    assay = "integrated",
    bdm_matrix = "BDM_integrated.rds",
    sdm_matrix = "SDM_integrated.rds",
    pd_list = "PD_list_dim1_th-1_integrated.rds",
#variable_features = "./intermediate_seuratIntegration_results/integration_features.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_integrated.Rds",
    landscape_list = "landscape_list_integrated.Rds",
    expr_list = readRDS("expr_list_integrated.Rds")
  )
)

# data_iterations_bonemarrow <- list(
#   list(
#     name = "Raw_bonemarrow",
#     seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "raw_bonemarrow_seurat_object.rds")),
#     # seurat_obj = merged_seurat_unintegrated_bonemarrow,
#     assay = "RNA",
#     bdm_matrix = "BDM_unintegrated_Raw_bonemarrow.rds",
#     sdm_matrix = "SDM_unintegrated_Raw_bonemarrow.rds",
#     pd_list = "PD_list_dim1_th-1_unintegrated_Raw_bonemarrow.Rds",
#     landscape_l2_distance_matrix = "landscape_l2_distance_matrix_Raw_bonemarrow.Rds",
#     landscape_list = "landscape_list_Raw_bonemarrow.Rds",
#     expr_list = readRDS("expr_list_raw_bonemarrow.Rds")
#   ),
#   # Uncomment SCT_Individual iteration if needed
#   list(
#     name = "SCT_Individual_bonemarrow",
#     seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_individual_bonemarrow_seurat_object.rds")),
#     # seurat_obj = merged_seurat_unintegrated_bonemarrow,
#     assay = "SCT_Ind",
#     bdm_matrix = "BDM_unintegrated_sctInd_bonemarrow.rds",
#     sdm_matrix = "SDM_unintegrated_sctInd_bonemarrow.rds",
#     pd_list = "PD_list_dim1_th-1_unintegrated_sctInd_bonemarrow.rds",
#     landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctInd_bonemarrow.Rds",
#     landscape_list = "landscape_list_sctInd_bonemarrow.Rds",
#     expr_list = readRDS("expr_list_sctInd_bonemarrow.Rds")
#   ),
#   list(
#     name = "SCT_Whole_bonemarrow",
#     seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_whole_bonemarrow_seurat_object.rds")),
#     # seurat_obj = merged_seurat_unintegrated_bonemarrow,
#     assay = "SCT",
#     bdm_matrix = "BDM_unintegrated_scTransformed_whole_bonemarrow.rds",
#     sdm_matrix = "SDM_unintegrated_scTransformed_whole_bonemarrow.rds",
#     pd_list = "PD_list_dim1_th-1_unintegrated_sctWhole_bonemarrow.Rds",
#     landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctWhole_bonemarrow.Rds",
#     landscape_list = "landscape_list_sctWhole_bonemarrow.Rds",
#     expr_list = readRDS("expr_list_sctwhole_bonemarrow.Rds")
#   ),
#   list(
#     name = "Integrated_bonemarrow",
#     seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "integrated_bonemarrow_seurat_object.rds")),
#     # seurat_obj = merged_seurat_unintegrated,
#     assay = "integrated",
#     bdm_matrix = "BDM_integrated_bonemarrow.rds",
#     sdm_matrix = "SDM_integrated_bonemarrow.rds",
#     pd_list = "PD_list_dim1_th-1_integrated_bonemarrow.rds",
#     #variable_features = "./intermediate_seuratIntegration_results/integration_features_bonemarrow.rds",
#     landscape_l2_distance_matrix = "landscape_l2_distance_matrix_integrated_bonemarrow.Rds",
#     landscape_list = "landscape_list_integrated_bonemarrow.Rds",
#     expr_list = readRDS("expr_list_integrated_bonemarrow.Rds")
#   )
# )


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
  library("TDA")
  # Wrap the computation in tryCatch so that errors are caught
  res <- tryCatch({
    landscape0 <- landscape(Diag = pd, dimension = 0, tseq = grid)
    landscape1 <- landscape(Diag = pd, dimension = 1, tseq = grid)
    list(dim0 = landscape0, dim1 = landscape1)
  },
    error = function(e) {
      log_message(paste("Error computing persistence landscape for a PD:", e$message))
      return(NULL)
    }
  )
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
validate_matrix <- function(mat, expected_dim = NULL, expected_names = NULL) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    cat("Validation failed: Object is not a matrix or data.frame.\n")
    return(FALSE)
  }

  dims <- dim(mat)
  cat("Matrix dimensions:", dims, "\n")
  if (!is.null(expected_dim) && !all(dims == expected_dim)) {
    cat("Dimension mismatch. Expected:", expected_dim, "\n")
    return(FALSE)
  }

  if (!is.null(expected_names)) {
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      cat("Row or column names are missing.\n")
      return(FALSE)
    }
    if (!all(rownames(mat) == expected_names)) {
      cat("Row names do not match expected names.\n")
      return(FALSE)
    }
    if (!all(colnames(mat) == expected_names)) {
      cat("Column names do not match expected names.\n")
      return(FALSE)
    }
  }

  if (!isTRUE(all.equal(mat, t(mat)))) {
    cat("Warning: Matrix is not symmetric.\n")
  } else {
    cat("Matrix is symmetric.\n")
  }

  return(TRUE)
}

# Validate that the landscape list is structured properly.
validate_landscape_list <- function(landscape_list, expected_length = NULL, grid_length = 100) {
  if (!is.list(landscape_list)) {
    cat("Validation failed: Landscape object is not a list.\n")
    return(FALSE)
  }

  n <- length(landscape_list)
  cat("Length of landscape list:", n, "\n")
  if (!is.null(expected_length) && n != expected_length) {
    cat("Landscape list length does not match expected value:", expected_length, "\n")
    return(FALSE)
  }

  valid <- TRUE
  for (i in seq_along(landscape_list)) {
    elem <- landscape_list[[i]]
    if (!is.list(elem) || is.null(elem$dim0) || is.null(elem$dim1)) {
      cat("Element", i, "is not properly structured (missing dim0 or dim1).\n")
      valid <- FALSE
    } else {
      len0 <- length(elem$dim0)
      len1 <- length(elem$dim1)
      cat("Element", i, "has dim0 length", len0, "and dim1 length", len1, "\n")
      if (len0 != grid_length || len1 != grid_length) {
        cat("Element", i, "does not have the expected grid length of", grid_length, "\n")
        valid <- FALSE
      }
    }
  }
  return(valid)
}

# ---------------------------
# Run the Iterations in Parallel
# ---------------------------
num_workers <- min(length(data_iterations), detectCores() - 1)
mclapply(data_iterations, function(iteration) {
  process_iteration_calculate_matrices(iteration, num_cores = 6, log_message = log_message)
}, mc.cores = num_workers)

# # Suppose your expected number of identifiers is 25 and your names vector is 'ordered_orig_idents'
# bdm <- readRDS("BDM_unintegrated_Raw_bonemarrow.rds")
# is_valid_bdm <- validate_matrix(bdm, expected_dim = c(25, 25), expected_names = ordered_orig_idents)
# if (is_valid_bdm) {
#   cat("BDM validated successfully.\n")
# }

# landscape_list <- readRDS("landscape_list_Raw_bonemarrow.Rds")
# is_valid_landscape <- validate_landscape_list(landscape_list, expected_length = 25, grid_length = 100)
# if (is_valid_landscape) {
#   cat("Landscape list validated successfully.\n")
# }


# landscape_matrix <- readRDS("landscape_l2_distance_matrix_Raw_bonemarrow.Rds")
# is_valid_landscape_mat <- validate_matrix(landscape_matrix, expected_dim = c(25, 25), expected_names = ordered_orig_idents)
# if (is_valid_landscape_mat) {
#   cat("Combined landscape distance matrix validated successfully.\n")
# }

# log_message("Distance matrix computation for all iterations completed.")

# ---------------------------
# Functions for Modular Analysis
# ---------------------------


# Function: Perform standard Seurat clustering, handling various assays and precomputed variable features
perform_standard_seurat_clustering <- function(seurat_obj,
                                               assay = "integrated",
                                               resolution = 0.5,
                                               variable_features_path = NULL) {
  # Set active assay
  DefaultAssay(seurat_obj) <- assay

  # Logging helper function
  log_message <- function(message) {
    cat(Sys.time(), "-", message, "\n")
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
perform_kmeans_clustering <- function(seurat_obj, assay, dims = 1:50, k = 5) {
  log_message <- function(message) {
    cat(Sys.time(), "-", message, "\n")
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
perform_hierarchical_clustering_ph <- function(bdm_matrix, k) {
  log_message <- function(message) {
    cat(Sys.time(), "-", message, "\n")
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
assign_ph_clusters <- function(seurat_obj, clusters_ph, new_cluster_col) {
  log_message <- function(message) {
    cat(Sys.time(), "-", message, "\n")
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
  if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
    seurat_sample_ids <- seurat_obj@meta.data$orig.ident
  } else if (any(grepl("__", rownames(seurat_obj@meta.data)))) {
    seurat_sample_ids <- sapply(strsplit(rownames(seurat_obj@meta.data), "__"), `[`, 1)
  } else {
    stop("No suitable sample identifier found in seurat_obj@meta.data.")
  }

  log_message("Sample IDs extracted from Seurat object:")
  print(head(seurat_sample_ids))

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



# Function to plot persistence diagrams and barcodes
plot_persistence <- function(pd, output_file_base, plot_title) {
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

    # Create the persistence diagram plot
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

    # Create the persistence barcode plot
    barcode_plot <- ggplot(pd_df, aes(y = as.factor(Dimension), x = Birth, xend = Death, color = as.factor(Dimension))) +
        geom_segment(size = 2) +
        scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
        labs(
          title = paste("Persistence Barcode -", plot_title),
          x = "Time",
          y = "Dimension",
    dim0_df$t <- grid
    df0 <- melt(dim0_df, id.vars = "t", variable.name = "Level", value.name = "Value")
    df0$Dimension <- "0"

    # Convert the landscape for dimension 1
    dim1_df <- as.data.frame(landscape$dim1)
    dim1_df$t <- grid
    df1 <- melt(dim1_df, id.vars = "t", variable.name = "Level", value.name = "Value")
    df1$Dimension <- "1"

    # Combine data from both dimensions
    df <- rbind(df0, df1)
    df$Level <- as.factor(df$Level)

    # Create the plot: one facet per dimension
    p <- ggplot(df, aes(x = t, y = Value, color = Level)) +
      geom_line() +
      facet_wrap(~Dimension, scales = "free_y") +
      labs(title = paste("Persistence Landscape -", plot_title),
           x = "t", y = "Landscape Value") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "top")

    # Save the plot in multiple formats
    ggsave(paste0(output_file_base, ".png"), p, width = 8, height = 6, dpi = 300)
    ggsave(paste0(output_file_base, ".svg"), p, width = 8, height = 6)
    ggsave(paste0(output_file_base, ".pdf"), p, width = 8, height = 6)
    log_message(paste("Saved landscape plot to", output_file_base))

  }, error = function(e) {
    log_message(paste("Error in plotting landscape data:", e$message))
  })
}


# ---------------------------
# Main Analysis Loop for Data Iterations with BDM and PD Loading
# ---------------------------
for (data_iter in data_iterations) {
  dataset_name <- data_iter$name
  seurat_obj <- data_iter$seurat_obj
  assay <- data_iter$assay

  variable_features_path <- if ("variable_features" %in% names(data_iter)) data_iter$variable_features else NULL

  source("PH_Functions.R")
  metadata <- read_csv("./data/VastlyDifferentTissues/metadata.csv")
  # metadata <- read_csv("./data/GSE120221/metadata.csv")

  colnames(metadata) <- gsub(" ", "_", colnames(metadata))

  # ---------------------------
  # Load BDMs and PD Lists
  # ---------------------------
  bdm_matrix <- as.matrix(readRDS(data_iter$bdm_matrix))
  sdm_matrix <- as.matrix(readRDS(data_iter$sdm_matrix))
  pd_list <- readRDS(data_iter$pd_list)

  # Load the Landscape L2 Distance Matrix (computed earlier in the pipeline)
  landscape_matrix <- as.matrix(readRDS(data_iter$landscape_l2_distance_matrix))
  landscape_list <- readRDS(data_iter$landscape_list)

  # Handle NA values in BDM
  bdm_matrix[is.na(bdm_matrix)] <- 0
  sdm_matrix[is.na(sdm_matrix)] <- 0
  landscape_matrix[is.na(landscape_matrix)] <- 0

    if ("SRA_Number" %in% colnames(metadata)) {
      SRA_col = "SRA_Number"
    } else {
      SRA_col = "orig.ident"
    }
    
  # Load and prepare the metadata
  metadata <- metadata %>%
    mutate(
      orig.ident = if ("SRA_Number" %in% colnames(metadata) & "SRS_Number" %in% colnames(metadata)) {
        paste(SRA_Number, SRS_Number, sep = "_")
      } else if ("File_Path" %in% colnames(metadata)) {
        basename(File_Path)
      } else {
    stop("Neither SRA_Number & SRS_Number nor File_Path columns exist in metadata.")
  }
    )

  # Load and prepare the metadata
  metadata <- metadata %>%
    mutate(
      Tissue = if ("Tissue" %in% colnames(metadata)) Tissue else if ("source_name" %in% colnames(metadata)) source_name else NA_character_
    )

  # Get the unique orig.idents from seurat_obj@meta.data
  orig_idents_in_seurat <- unique(seurat_obj@meta.data$orig.ident)

  # Filter the metadata to keep only the orig.idents that appear in the seurat_obj
  metadata_filtered <- metadata %>%
    filter(orig.ident %in% orig_idents_in_seurat)

  # Extract the ordered list of orig.idents after filtering
  ordered_orig_idents <- as.character(metadata_filtered$orig.ident)

  # Set the dimnames of bdm_matrix to match the filtered orig.ident list
  dimnames(bdm_matrix) <- list(ordered_orig_idents, ordered_orig_idents)
  dimnames(sdm_matrix) <- list(ordered_orig_idents, ordered_orig_idents)
  dimnames(landscape_matrix) <- list(ordered_orig_idents, ordered_orig_idents)

  names(pd_list) <- ordered_orig_idents
  names(landscape_list) <- ordered_orig_idents

  # Verify the updated rownames and colnames
  print(rownames(bdm_matrix))
  print(colnames(bdm_matrix))
  print(rownames(sdm_matrix))
  print(colnames(sdm_matrix))
  print(rownames(landscape_matrix))
  print(colnames(landscape_matrix))

  if (save_pd_plots) {
    tryCatch({
      # Create folder for plots if it doesn't exist
      plot_dir <- plot_folders$persistence_diagrams
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE)
        log_message(paste("Created directory for plots:", plot_dir))
      }

      # Plot and save each persistence diagram
      mclapply(seq_along(pd_list), function(i) {
        pd <- pd_list[[i]]
        landscape <- landscape_list[[i]]
        name <- names(pd_list)[i] %||% as.character(i) # Use name or default to index

        if (!is.matrix(pd) || ncol(pd) != 3) {
          log_message(paste("Skipping plot for", name, "as PD is not a valid matrix."))
          return()
        }

        pd_plot_file <- file.path(plot_dir, paste0("PD_plot_", name, ".png"))
        plot_title <- paste("Persistence Diagram -", name)
        plot_persistence(pd, pd_plot_file, plot_title)
        log_message(paste("Saved plot for PD:", data_iter$name, "- Name:", name, "to", pd_plot_file))

        # Define output file base (you can adjust naming as needed)
        pl_plot_file <- file.path(plot_dir, paste0("Landscape_plot_", name))
        plot_title <- paste("Persistence Landscape -", name)
        plot_landscape(landscape, pl_plot_file, plot_title, grid = seq(0, 1, length.out = 100))
        log_message(paste("Saved landscape plot for:", name, "to", pl_plot_file))
      }, mc.cores = detectCores() - 1) # Use available cores minus one
    },
      error = function(e) {
        log_message(paste("Error in processing iteration", data_iter$name, ":", e$message))
      }
    )
  }

  # ---------------------------
  # Hierarchical PH Clustering
  # ---------------------------
  if (run_hierarchical_ph_clustering) {
    # Only define grouping variables if the corresponding column exists
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
      k_tissue <- length(unique(seurat_obj$Tissue))
    }
    if (SRA_col %in% colnames(seurat_obj@meta.data)) {
      k_sra <- length(unique(seurat_obj@meta.data[[SRA_col]]))
    }
    if ("Approach" %in% colnames(seurat_obj@meta.data)) {
      k_approach <- length(unique(seurat_obj$Approach))
    }

    # For random groups, look for any column starting with "Random_Group"
    random_group_cols <- grep("^Random_Group", colnames(seurat_obj@meta.data), value = TRUE)

    # ---------------------------
    # Hierarchical Clustering on BDM
    # ---------------------------
    if (!is.null(bdm_matrix)) {
      # Tissue clustering
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on BDM for tissues with k =", k_tissue, "\n"))
        ph_tissue_result_bdm <- try(perform_hierarchical_clustering_ph(bdm_matrix, k = k_tissue), silent = TRUE)
        if (!inherits(ph_tissue_result_bdm, "try-error")) {
          ph_clusters_tissue_bdm <- ph_tissue_result_bdm$clusters
          hc_tissue_bdm <- ph_tissue_result_bdm$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_tissue_bdm,
            paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_tissue"))
          cat(paste(Sys.time(), "- Hierarchical clustering on BDM for tissues completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on BDM for tissues failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Tissue column not found in seurat_obj. Skipping BDM tissue clustering.\n"))
      }

      # SRA clustering
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on BDM for SRAs with k =", k_sra, "\n"))
        ph_sra_result_bdm <- try(perform_hierarchical_clustering_ph(bdm_matrix, k = k_sra), silent = TRUE)
        if (!inherits(ph_sra_result_bdm, "try-error")) {
          ph_clusters_sra_bdm <- ph_sra_result_bdm$clusters
          hc_sra_bdm <- ph_sra_result_bdm$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_sra_bdm,
            paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_sra"))
          cat(paste(Sys.time(), "- Hierarchical clustering on BDM for SRAs completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on BDM for SRAs failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- orig.ident column not found in seurat_obj. Skipping BDM SRA clustering.\n"))
      }

      # Approach clustering
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on BDM for Approach with k =", k_approach, "\n"))
        ph_approach_result_bdm <- try(perform_hierarchical_clustering_ph(bdm_matrix, k = k_approach), silent = TRUE)
        if (!inherits(ph_approach_result_bdm, "try-error")) {
          ph_clusters_approach_bdm <- ph_approach_result_bdm$clusters
          hc_approach_bdm <- ph_approach_result_bdm$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_approach_bdm,
            paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_approach"))
          cat(paste(Sys.time(), "- Hierarchical clustering on BDM for Approach completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on BDM for Approach failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Approach column not found in seurat_obj. Skipping BDM Approach clustering.\n"))
      }

      # # Random_Group clustering (looping over all bootstrapped columns)
      # if (length(random_group_cols) > 0) {
      #   for (random_col in random_group_cols) {
      #     k_random_group <- length(unique(seurat_obj@meta.data[[random_col]]))
      #     cat(paste(Sys.time(), "- Running hierarchical clustering on BDM for", random_col, "with k =", k_random_group, "\n"))
      #     ph_random_result_bdm <- try(perform_hierarchical_clustering_ph(bdm_matrix, k = k_random_group), silent = TRUE)
      #     if (!inherits(ph_random_result_bdm, "try-error")) {
      #       ph_clusters_random_bdm <- ph_random_result_bdm$clusters
      #       hc_random_bdm <- ph_random_result_bdm$tree
      #       # Append the bootstrap column name to the output column name
      #       new_cluster_name <- paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_", random_col)
      #       seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_random_bdm, new_cluster_name)
      #       cat(paste(Sys.time(), "- Hierarchical clustering on BDM for", random_col, "completed.\n"))
      #     } else {
      #       cat(paste(Sys.time(), "- Hierarchical clustering on BDM for", random_col, "failed. Moving on.\n"))
      #     }
      #   }
      # } else {
      #   cat(paste(Sys.time(), "- No Random_Group column found in seurat_obj. Skipping BDM Random_Group clustering.\n"))
      # }
    }

    # ---------------------------
    # Hierarchical Clustering on SDM
    # ---------------------------
    if (!is.null(sdm_matrix)) {
      # Tissue clustering
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on SDM for tissues with k =", k_tissue, "\n"))
        ph_tissue_result_sdm <- try(perform_hierarchical_clustering_ph(sdm_matrix, k = k_tissue), silent = TRUE)
        if (!inherits(ph_tissue_result_sdm, "try-error")) {
          ph_clusters_tissue_sdm <- ph_tissue_result_sdm$clusters
          hc_tissue_sdm <- ph_tissue_result_sdm$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_tissue_sdm,
            paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_tissue"))
          cat(paste(Sys.time(), "- Hierarchical clustering on SDM for tissues completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on SDM for tissues failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Tissue column not found in seurat_obj. Skipping SDM tissue clustering.\n"))
      }

      # SRA clustering
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on SDM for SRAs with k =", k_sra, "\n"))
        ph_sra_result_sdm <- try(perform_hierarchical_clustering_ph(sdm_matrix, k = k_sra), silent = TRUE)
        if (!inherits(ph_sra_result_sdm, "try-error")) {
          ph_clusters_sra_sdm <- ph_sra_result_sdm$clusters
          hc_sra_sdm <- ph_sra_result_sdm$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_sra_sdm,
            paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_sra"))
          cat(paste(Sys.time(), "- Hierarchical clustering on SDM for SRAs completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on SDM for SRAs failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- orig.ident column not found in seurat_obj. Skipping SDM SRA clustering.\n"))
      }

      # Approach clustering
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on SDM for Approach with k =", k_approach, "\n"))
        ph_approach_result_sdm <- try(perform_hierarchical_clustering_ph(sdm_matrix, k = k_approach), silent = TRUE)
        if (!inherits(ph_approach_result_sdm, "try-error")) {
          ph_clusters_approach_sdm <- ph_approach_result_sdm$clusters
          hc_approach_sdm <- ph_approach_result_sdm$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_approach_sdm,
            paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_approach"))
          cat(paste(Sys.time(), "- Hierarchical clustering on SDM for Approach completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on SDM for Approach failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Approach column not found in seurat_obj. Skipping SDM Approach clustering.\n"))
      }

      # # Random_Group clustering (looping over all bootstrapped columns)
      # if (length(random_group_cols) > 0) {
      #   for (random_col in random_group_cols) {
      #     k_random_group <- length(unique(seurat_obj@meta.data[[random_col]]))
      #     cat(paste(Sys.time(), "- Running hierarchical clustering on SDM for", random_col, "with k =", k_random_group, "\n"))
      #     ph_random_result_sdm <- try(perform_hierarchical_clustering_ph(sdm_matrix, k = k_random_group), silent = TRUE)
      #     if (!inherits(ph_random_result_sdm, "try-error")) {
      #       ph_clusters_random_sdm <- ph_random_result_sdm$clusters
      #       hc_random_sdm <- ph_random_result_sdm$tree
      #       new_cluster_name <- paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_", random_col)
      #       seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_random_sdm, new_cluster_name)
      #       cat(paste(Sys.time(), "- Hierarchical clustering on SDM for", random_col, "completed.\n"))
      #     } else {
      #       cat(paste(Sys.time(), "- Hierarchical clustering on SDM for", random_col, "failed. Moving on.\n"))
      #     }
      #   }
      # } else {
      #   cat(paste(Sys.time(), "- No Random_Group column found in seurat_obj. Skipping SDM Random_Group clustering.\n"))
      # }
    }

    # ---------------------------
    # Hierarchical Clustering on Landscape L2 Distance Matrix
    # ---------------------------
    if (!is.null(landscape_matrix)) {
      # Tissue clustering on Landscape distances
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on Landscape distances for tissues with k =", k_tissue, "\n"))
        ph_tissue_result_landscape <- try(perform_hierarchical_clustering_ph(landscape_matrix, k = k_tissue), silent = TRUE)
        if (!inherits(ph_tissue_result_landscape, "try-error")) {
          ph_clusters_tissue_landscape <- ph_tissue_result_landscape$clusters
          hc_tissue_landscape <- ph_tissue_result_landscape$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_tissue_landscape,
                           paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_tissue"))
          cat(paste(Sys.time(), "- Hierarchical clustering on Landscape distances for tissues completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on Landscape distances for tissues failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Tissue column not found in seurat_obj. Skipping Landscape tissue clustering.\n"))
      }

      # SRA clustering on Landscape distances
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on Landscape distances for SRAs with k =", k_sra, "\n"))
        ph_sra_result_landscape <- try(perform_hierarchical_clustering_ph(landscape_matrix, k = k_sra), silent = TRUE)
        if (!inherits(ph_sra_result_landscape, "try-error")) {
          ph_clusters_sra_landscape <- ph_sra_result_landscape$clusters
          hc_sra_landscape <- ph_sra_result_landscape$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_sra_landscape,
                           paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_sra"))
          cat(paste(Sys.time(), "- Hierarchical clustering on Landscape distances for SRAs completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on Landscape distances for SRAs failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- orig.ident column not found in seurat_obj. Skipping Landscape SRA clustering.\n"))
      }

      # Approach clustering on Landscape distances
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running hierarchical clustering on Landscape distances for Approach with k =", k_approach, "\n"))
        ph_approach_result_landscape <- try(perform_hierarchical_clustering_ph(landscape_matrix, k = k_approach), silent = TRUE)
        if (!inherits(ph_approach_result_landscape, "try-error")) {
          ph_clusters_approach_landscape <- ph_approach_result_landscape$clusters
          hc_approach_landscape <- ph_approach_result_landscape$tree
          seurat_obj <- assign_ph_clusters(seurat_obj, ph_clusters_approach_landscape,
                           paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_approach"))
          cat(paste(Sys.time(), "- Hierarchical clustering on Landscape distances for Approach completed.\n"))
        } else {
          cat(paste(Sys.time(), "- Hierarchical clustering on Landscape distances for Approach failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Approach column not found in seurat_obj. Skipping Landscape Approach clustering.\n"))
      }

      # Optional: Add Random_Group clustering for Landscape distances here if needed.
    }

    # ---------------------------
    # Log Clustering Outputs
    # ---------------------------
    cat("Clustering completed:\n")
    if (!is.null(bdm_matrix)) {
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat("BDM - Tissue k:", k_tissue, "\n")
      }
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat("BDM - SRA k:", k_sra, "\n")
      }
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat("BDM - Approach k:", k_approach, "\n")
      }
      if (length(random_group_cols) > 0) {
        for (random_col in random_group_cols) {
          k_random_group <- length(unique(seurat_obj@meta.data[[random_col]]))
          cat("BDM -", random_col, "k:", k_random_group, "\n")
        }
      }
    }
    if (!is.null(sdm_matrix)) {
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat("SDM - Tissue k:", k_tissue, "\n")
      }
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat("SDM - SRA k:", k_sra, "\n")
      }
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat("SDM - Approach k:", k_approach, "\n")
      }
      if (length(random_group_cols) > 0) {
        for (random_col in random_group_cols) {
          k_random_group <- length(unique(seurat_obj@meta.data[[random_col]]))
          cat("SDM -", random_col, "k:", k_random_group, "\n")
        }
      }
    }
    if (!is.null(landscape_matrix)) {
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat("Landscape - Tissue k:", k_tissue, "\n")
      }
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat("Landscape - SRA k:", k_sra, "\n")
      }
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat("Landscape - Approach k:", k_approach, "\n")
      }
      if (length(random_group_cols) > 0) {
        for (random_col in random_group_cols) {
          k_random_group <- length(unique(seurat_obj@meta.data[[random_col]]))
          cat("Landscape -", random_col, "k:", k_random_group, "\n")
        }
      }
    }
  }

  # ---------------------------
  # Standard Seurat Clustering
  # ---------------------------
  if (run_standard_seurat_clustering) {
    seurat_obj <- perform_standard_seurat_clustering(seurat_obj, assay = assay, variable_features_path = variable_features_path)
    seurat_obj@meta.data[[paste0("seurat_cluster_", tolower(dataset_name))]] <- Idents(seurat_obj)
    Idents(seurat_obj) <- seurat_obj@meta.data[[paste0("seurat_cluster_", tolower(dataset_name))]]
    seurat_obj@meta.data$seurat_clusters <- NULL
  }

  # ---------------------------
  # K-means Clustering by Tissues, SRAs, Approach, and Random_Group
  # ---------------------------
  if (run_kmeans_clustering) {
    # K-means clustering for Tissue
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
      k_tissue <- length(unique(seurat_obj$Tissue))
      cat(paste(Sys.time(), "- Running K-means clustering with k =", k_tissue, "(number of tissues)\n"))
      seurat_obj <- perform_kmeans_clustering(seurat_obj, assay = assay, dims = 1:50, k = k_tissue)
      seurat_obj@meta.data[[paste0("kmeans_cluster_", tolower(dataset_name), "_tissue")]] <- seurat_obj$kmeans_cluster
      cat(paste(Sys.time(), "- K-means clustering with k =", k_tissue, "completed and stored in",
                paste0("kmeans_cluster_", tolower(dataset_name), "_tissue"), "\n"))
    } else {
      cat(paste(Sys.time(), "- Tissue column not found in seurat_obj. Skipping K-means tissue clustering.\n"))
    }

    # K-means clustering for SRA
    if (SRA_col %in% colnames(seurat_obj@meta.data)) {
      k_sra <- length(unique(seurat_obj@meta.data[[SRA_col]]))
      cat(paste(Sys.time(), "- Running K-means clustering with k =", k_sra, "(number of SRAs)\n"))
      seurat_obj <- perform_kmeans_clustering(seurat_obj, assay = assay, dims = 1:50, k = k_sra)
      seurat_obj@meta.data[[paste0("kmeans_cluster_", tolower(dataset_name), "_sra")]] <- seurat_obj$kmeans_cluster
      cat(paste(Sys.time(), "- K-means clustering with k =", k_sra, "completed and stored in",
                paste0("kmeans_cluster_", tolower(dataset_name), "_sra"), "\n"))
    } else {
      cat(paste(Sys.time(), "- orig.ident column not found in seurat_obj. Skipping K-means SRA clustering.\n"))
    }

    # K-means clustering for Approach
    if ("Approach" %in% colnames(seurat_obj@meta.data)) {
      k_approach <- length(unique(seurat_obj$Approach))
      cat(paste(Sys.time(), "- Running K-means clustering with k =", k_approach, "(number of Approaches)\n"))
      seurat_obj <- perform_kmeans_clustering(seurat_obj, assay = assay, dims = 1:50, k = k_approach)
      seurat_obj@meta.data[[paste0("kmeans_cluster_", tolower(dataset_name), "_approach")]] <- seurat_obj$kmeans_cluster
      cat(paste(Sys.time(), "- K-means clustering with k =", k_approach, "completed and stored in",
                paste0("kmeans_cluster_", tolower(dataset_name), "_approach"), "\n"))
    } else {
      cat(paste(Sys.time(), "- Approach column not found in seurat_obj. Skipping K-means Approach clustering.\n"))
    }

    # K-means clustering for Random_Group
    # # Look for all columns starting with "Random_Group" (this includes bootstrapped columns)
    # random_group_cols <- grep("^Random_Group", colnames(seurat_obj@meta.data), value = TRUE)
    # if (length(random_group_cols) > 0) {
    #   for (rg in random_group_cols) {
    #     k_random_group <- length(unique(seurat_obj@meta.data[[rg]]))
    #     cat(paste(Sys.time(), "- Running K-means clustering with k =", k_random_group, "for", rg, "\n"))
    #     seurat_obj <- perform_kmeans_clustering(seurat_obj, assay = assay, dims = 1:50, k = k_random_group)
    #     new_col_name <- paste0("kmeans_cluster_", tolower(dataset_name), "_", rg)
    #     seurat_obj@meta.data[[new_col_name]] <- seurat_obj$kmeans_cluster
    #     cat(paste(Sys.time(), "- K-means clustering for", rg, "with k =", k_random_group, "completed and stored in", new_col_name, "\n"))
    #   }
    # } else {
    #   cat(paste(Sys.time(), "- Random_Group column not found in seurat_obj. Skipping K-means Random_Group clustering.\n"))
    # }

    seurat_obj$kmeans_cluster <- NULL
  }

  # ---------------------------
  # Spectral Clustering by Tissues, SRAs, Approach, and Random_Group (BDM and SDM)
  # ---------------------------
  if (run_spectral_clustering) {
    # ---------------------------
    # Spectral Clustering Function
    # ---------------------------
    perform_spectral_clustering <- function(distance_matrix, k) {
      # Convert distance matrix to similarity matrix
      epsilon <- 1e-8
      similarity_matrix <- exp(-distance_matrix ^ 2 / (2 * (mean(distance_matrix) + epsilon) ^ 2))
      # Perform spectral clustering
      spectral_result <- specc(as.kernelMatrix(similarity_matrix), centers = k)
      # Return cluster assignments as a named vector using sample names
      cluster_assignments <- spectral_result@.Data
      names(cluster_assignments) <- rownames(distance_matrix)
      return(cluster_assignments)
    }

    library(kernlab) # For spectral clustering

    # Determine the number of clusters for each case, if the column exists
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
      k_tissue <- length(unique(seurat_obj$Tissue))
    }
    if (SRA_col %in% colnames(seurat_obj@meta.data)) {
      k_sra <- length(unique(seurat_obj@meta.data[[SRA_col]]))
    }
    if ("Approach" %in% colnames(seurat_obj@meta.data)) {
      k_approach <- length(unique(seurat_obj$Approach))
    }

    # For random groups, get all columns starting with "Random_Group"
    random_group_cols <- grep("^Random_Group", colnames(seurat_obj@meta.data), value = TRUE)

    # ---------------------------
    # Spectral Clustering on BDM
    # ---------------------------
    if (!is.null(bdm_matrix)) {
      # Spectral Clustering for Tissue (BDM)
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on BDM for tissues with k =", k_tissue, "\n"))
        spectral_tissue_clusters_bdm <- try(perform_spectral_clustering(bdm_matrix, k = k_tissue), silent = TRUE)
        if (!inherits(spectral_tissue_clusters_bdm, "try-error")) {
          seurat_obj <- assign_ph_clusters(seurat_obj, spectral_tissue_clusters_bdm,
            paste0("spectral_cluster_bdm_", tolower(dataset_name), "_tissue"))
          cat(paste(Sys.time(), "- Spectral clustering on BDM for tissues completed and stored in",
              paste0("spectral_cluster_bdm_", tolower(dataset_name), "_tissue"), "\n"))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on BDM for tissues failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Tissue column not found in seurat_obj. Skipping Spectral clustering for tissues (BDM).\n"))
      }

      # Spectral Clustering for SRA (BDM)
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on BDM for SRAs with k =", k_sra, "\n"))
        spectral_sra_clusters_bdm <- try(perform_spectral_clustering(bdm_matrix, k = k_sra), silent = TRUE)
        if (!inherits(spectral_sra_clusters_bdm, "try-error")) {
          seurat_obj <- assign_ph_clusters(seurat_obj, spectral_sra_clusters_bdm,
            paste0("spectral_cluster_bdm_", tolower(dataset_name), "_sra"))
          cat(paste(Sys.time(), "- Spectral clustering on BDM for SRAs completed and stored in",
              paste0("spectral_cluster_bdm_", tolower(dataset_name), "_sra"), "\n"))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on BDM for SRAs failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- orig.ident column not found in seurat_obj. Skipping Spectral clustering for SRAs (BDM).\n"))
      }

      # Spectral Clustering for Approach (BDM)
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on BDM for Approach with k =", k_approach, "\n"))
        spectral_approach_clusters_bdm <- try(perform_spectral_clustering(bdm_matrix, k = k_approach), silent = TRUE)
        if (!inherits(spectral_approach_clusters_bdm, "try-error")) {
          seurat_obj <- assign_ph_clusters(seurat_obj, spectral_approach_clusters_bdm,
            paste0("spectral_cluster_bdm_", tolower(dataset_name), "_approach"))
          cat(paste(Sys.time(), "- Spectral clustering on BDM for Approach completed and stored in",
              paste0("spectral_cluster_bdm_", tolower(dataset_name), "_approach"), "\n"))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on BDM for Approach failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Approach column not found in seurat_obj. Skipping Spectral clustering for Approach (BDM).\n"))
      }

      # # Spectral Clustering for Random_Group (BDM) across all bootstrapped columns
      # if (length(random_group_cols) > 0) {
      #   for (random_col in random_group_cols) {
      #     k_random_group <- length(unique(seurat_obj@meta.data[[random_col]]))
      #     cat(paste(Sys.time(), "- Running Spectral clustering on BDM for", random_col, "with k =", k_random_group, "\n"))
      #     spectral_random_clusters_bdm <- try(perform_spectral_clustering(bdm_matrix, k = k_random_group), silent = TRUE)
      #     if (!inherits(spectral_random_clusters_bdm, "try-error")) {
      #       new_cluster_name <- paste0("spectral_cluster_bdm_", tolower(dataset_name), "_", random_col)
      #       seurat_obj <- assign_ph_clusters(seurat_obj, spectral_random_clusters_bdm, new_cluster_name)
      #       cat(paste(Sys.time(), "- Spectral clustering on BDM for", random_col, "completed and stored in", new_cluster_name, "\n"))
      #     } else {
      #       cat(paste(Sys.time(), "- Spectral clustering on BDM for", random_col, "failed. Moving on.\n"))
      #     }
      #   }
      # } else {
      #   cat(paste(Sys.time(), "- No Random_Group column found in seurat_obj. Skipping Spectral clustering for Random_Group (BDM).\n"))
      # }
    }

    # ---------------------------
    # Spectral Clustering on SDM
    # ---------------------------
    if (!is.null(sdm_matrix)) {
      # Spectral Clustering for Tissue (SDM)
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on SDM for tissues with k =", k_tissue, "\n"))
        spectral_tissue_clusters_sdm <- try(perform_spectral_clustering(sdm_matrix, k = k_tissue), silent = TRUE)
        if (!inherits(spectral_tissue_clusters_sdm, "try-error")) {
          seurat_obj <- assign_ph_clusters(
            seurat_obj, spectral_tissue_clusters_sdm,
            paste0("spectral_cluster_sdm_", tolower(dataset_name), "_tissue")
          )
          cat(paste(
            Sys.time(), "- Spectral clustering on SDM for tissues completed and stored in",
            paste0("spectral_cluster_sdm_", tolower(dataset_name), "_tissue"), "\n"
          ))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on SDM for tissues failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Tissue column not found in seurat_obj. Skipping Spectral clustering for tissues (SDM).\n"))
      }

      # Spectral Clustering for SRA (SDM)
      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on SDM for SRAs with k =", k_sra, "\n"))
        spectral_sra_clusters_sdm <- try(perform_spectral_clustering(sdm_matrix, k = k_sra), silent = TRUE)
        if (!inherits(spectral_sra_clusters_sdm, "try-error")) {
          seurat_obj <- assign_ph_clusters(
            seurat_obj, spectral_sra_clusters_sdm,
            paste0("spectral_cluster_sdm_", tolower(dataset_name), "_sra")
          )
          cat(paste(
            Sys.time(), "- Spectral clustering on SDM for SRAs completed and stored in",
            paste0("spectral_cluster_sdm_", tolower(dataset_name), "_sra"), "\n"
          ))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on SDM for SRAs failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- orig.ident column not found in seurat_obj. Skipping Spectral clustering for SRAs (SDM).\n"))
      }

      # Spectral Clustering for Approach (SDM)
      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on SDM for Approach with k =", k_approach, "\n"))
        spectral_approach_clusters_sdm <- try(perform_spectral_clustering(sdm_matrix, k = k_approach), silent = TRUE)
        if (!inherits(spectral_approach_clusters_sdm, "try-error")) {
          seurat_obj <- assign_ph_clusters(
            seurat_obj, spectral_approach_clusters_sdm,
            paste0("spectral_cluster_sdm_", tolower(dataset_name), "_approach")
          )
          cat(paste(
            Sys.time(), "- Spectral clustering on SDM for Approach completed and stored in",
            paste0("spectral_cluster_sdm_", tolower(dataset_name), "_approach"), "\n"
          ))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on SDM for Approach failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Approach column not found in seurat_obj. Skipping Spectral clustering for Approach (SDM).\n"))
      }

      # # Spectral Clustering for Random_Group (SDM) across all bootstrapped columns
      # if (length(random_group_cols) > 0) {
      #   for (random_col in random_group_cols) {
      #     k_random_group <- length(unique(seurat_obj@meta.data[[random_col]]))
      #     cat(paste(Sys.time(), "- Running Spectral clustering on SDM for", random_col, "with k =", k_random_group, "\n"))
      #     spectral_random_clusters_sdm <- try(perform_spectral_clustering(sdm_matrix, k = k_random_group), silent = TRUE)
      #     if (!inherits(spectral_random_clusters_sdm, "try-error")) {
      #       new_cluster_name <- paste0("spectral_cluster_sdm_", tolower(dataset_name), "_", random_col)
      #       seurat_obj <- assign_ph_clusters(seurat_obj, spectral_random_clusters_sdm, new_cluster_name)
      #       cat(paste(Sys.time(), "- Spectral clustering on SDM for", random_col, "completed and stored in", new_cluster_name, "\n"))
      #     } else {
      #       cat(paste(Sys.time(), "- Spectral clustering on SDM for", random_col, "failed. Moving on.\n"))
      #     }
      #   }
      # } else {
      #   cat(paste(Sys.time(), "- No Random_Group column found in seurat_obj. Skipping Spectral clustering for Random_Group (SDM).\n"))
      # }
    }

    # ---------------------------
    # Spectral Clustering on Landscape Distance Matrix
    # ---------------------------
    # Make sure the landscape matrix is loaded and its dimnames are set appropriately
    if (!is.null(landscape_matrix)) {
      if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on Landscape distances for tissues with k =", k_tissue, "\n"))
        spectral_tissue_clusters_landscape <- try(perform_spectral_clustering(landscape_matrix, k = k_tissue), silent = TRUE)
        if (!inherits(spectral_tissue_clusters_landscape, "try-error")) {
          seurat_obj <- assign_ph_clusters(seurat_obj, spectral_tissue_clusters_landscape,
            paste0("spectral_cluster_landscape_", tolower(dataset_name), "_tissue"))
          cat(paste(Sys.time(), "- Spectral clustering on Landscape distances for tissues completed and stored in",
              paste0("spectral_cluster_landscape_", tolower(dataset_name), "_tissue"), "\n"))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on Landscape distances for tissues failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Tissue column not found in seurat_obj. Skipping Spectral clustering for Landscape (tissues).\n"))
      }

      if (SRA_col %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on Landscape distances for SRAs with k =", k_sra, "\n"))
        spectral_sra_clusters_landscape <- try(perform_spectral_clustering(landscape_matrix, k = k_sra), silent = TRUE)
        if (!inherits(spectral_sra_clusters_landscape, "try-error")) {
          seurat_obj <- assign_ph_clusters(seurat_obj, spectral_sra_clusters_landscape,
            paste0("spectral_cluster_landscape_", tolower(dataset_name), "_sra"))
          cat(paste(Sys.time(), "- Spectral clustering on Landscape distances for SRAs completed and stored in",
              paste0("spectral_cluster_landscape_", tolower(dataset_name), "_sra"), "\n"))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on Landscape distances for SRAs failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- orig.ident column not found in seurat_obj. Skipping Spectral clustering for Landscape (SRAs).\n"))
      }

      if ("Approach" %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Running Spectral clustering on Landscape distances for Approach with k =", k_approach, "\n"))
        spectral_approach_clusters_landscape <- try(perform_spectral_clustering(landscape_matrix, k = k_approach), silent = TRUE)
        if (!inherits(spectral_approach_clusters_landscape, "try-error")) {
          seurat_obj <- assign_ph_clusters(seurat_obj, spectral_approach_clusters_landscape,
            paste0("spectral_cluster_landscape_", tolower(dataset_name), "_approach"))
          cat(paste(Sys.time(), "- Spectral clustering on Landscape distances for Approach completed and stored in",
              paste0("spectral_cluster_landscape_", tolower(dataset_name), "_approach"), "\n"))
        } else {
          cat(paste(Sys.time(), "- Spectral clustering on Landscape distances for Approach failed. Moving on.\n"))
        }
      } else {
        cat(paste(Sys.time(), "- Approach column not found in seurat_obj. Skipping Spectral clustering for Landscape (Approach).\n"))
      }

      # Optionally, add Random_Group clustering for Landscape if desired.
    }
  }

  # Save the Seurat Object
  save_path <- paste0(results_folder, "/seurat_objects/", tolower(dataset_name), "_seurat_object.rds")
  if (!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)
  saveRDS(seurat_obj, save_path)
  cat(paste(Sys.time(), "- Seurat object saved for dataset:", dataset_name, "at", save_path, "\n"))

  # ---------------------------
  # UMAP Calculations and Visualizations with Annotations
  # ---------------------------
  if (run_visualizations) {
    # Check if UMAP has already been run
    if (!"umap" %in% names(seurat_obj@reductions)) {
      cat(paste(Sys.time(), "- UMAP embedding not found. Calculating UMAP.\n"))
      seurat_obj <- RunUMAP(seurat_obj,
        dims = 1:50, reduction = "pca", assay = assay,
        verbose = TRUE, umap.method = "uwot", n.neighbors = 200L,
        min.dist = 0.001, seed.use = 1L
      )
    } else {
      cat(paste(Sys.time(), "- UMAP embedding already present. Skipping UMAP calculation.\n"))
    }

    # Identify any Random_Group cluster assignment columns (including bootstrapped ones)
    random_group_cluster_cols <- grep("random_group", colnames(seurat_obj@meta.data),
      ignore.case = TRUE, value = TRUE
    )

    # Define cluster types for visualization and analysis
    cluster_types <- c(
      paste0("seurat_cluster_", tolower(dataset_name)),
      paste0("kmeans_cluster_", tolower(dataset_name), "_tissue"),
      paste0("kmeans_cluster_", tolower(dataset_name), "_sra"),
      paste0("kmeans_cluster_", tolower(dataset_name), "_approach"),
      paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_tissue"),
      paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_sra"),
      paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_approach"),
      paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_tissue"),
      paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_sra"),
      paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_approach"),
      paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_tissue"),
      paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_sra"),
      paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_approach"),
      paste0("spectral_cluster_bdm_", tolower(dataset_name), "_tissue"),
      paste0("spectral_cluster_bdm_", tolower(dataset_name), "_sra"),
      paste0("spectral_cluster_bdm_", tolower(dataset_name), "_approach"),
      paste0("spectral_cluster_sdm_", tolower(dataset_name), "_tissue"),
      paste0("spectral_cluster_sdm_", tolower(dataset_name), "_sra"),
      paste0("spectral_cluster_sdm_", tolower(dataset_name), "_approach"),
      paste0("spectral_cluster_landscape_", tolower(dataset_name), "_tissue"),
      paste0("spectral_cluster_landscape_", tolower(dataset_name), "_sra"),
      paste0("spectral_cluster_landscape_", tolower(dataset_name), "_approach"),
      random_group_cluster_cols, # include all Random_Group (and bootstrapped) columns
      "Tissue",
      "SRA",
      "Approach",
      "orig.ident"
    )

    # --- Filter out duplicate bootstrapped random group clusters ---
    # For any cluster name that contains "_Random_Group_bootstrap_" (case-insensitive),
    # keep only the first occurrence for each base method.
    filtered_cluster_types <- c()
    boot_base_seen <- c()
    for (clust in cluster_types) {
      if (grepl("random_group_bootstrap", clust, ignore.case = TRUE)) {
        # Split the string at the pattern and take the base portion.
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

    # Open multi-page PDF for UMAP plots of each cluster type
    pdf(file.path(plots_folder, paste0("UMAP_Plots_", dataset_name, "_All_Clusters.pdf")))

    for (cluster_col in cluster_types) {
      if (cluster_col %in% colnames(seurat_obj@meta.data)) {
        cat(paste(Sys.time(), "- Plotting UMAP for", cluster_col, "\n"))
        umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = cluster_col) +
          ggtitle(paste("UMAP -", dataset_name, "(", cluster_col, ")")) +
          xlab("UMAP 1") + ylab("UMAP 2") +
          theme_minimal() +
          theme(
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 10, face = "bold")
          )
        print(umap_plot)
      } else {
        cat(paste(Sys.time(), "- Cluster type", cluster_col, "not found in Seurat object metadata. Skipping.\n"))
      }
    }

    # Close PDF after saving all pages
    dev.off()
    cat(paste(Sys.time(), "- UMAP plots saved to", file.path(plots_folder, paste0("UMAP_Plots_", dataset_name, "_All_Clusters.pdf")), "\n"))
  }

  if (run_sample_level_heatmap) {
    library(ComplexHeatmap)
    library(viridis)
    library(dendextend)

    generate_heatmaps <- function(dataset_name, metadata, seurat_obj, bdm_matrix, plots_folder,
                                  hc_tissue, hc_sra, hc_approach,
                                  k_tissue, k_sra, k_approach) {
      # Validate BDM matrix
      if (!is.matrix(bdm_matrix) || nrow(bdm_matrix) != ncol(bdm_matrix)) {
        stop("BDM must be a square matrix with equal rows and columns.")
      }
      if (is.null(rownames(bdm_matrix)) || is.null(colnames(bdm_matrix))) {
        stop("BDM matrix must have row and column names.")
      }

      # Check if metadata has multiple SRA values
      if ("SRA_Number" %in% colnames(metadata) && length(unique(metadata$SRA_Number)) > 1) {
        # Use the multi-SRA approach
        metadata2 <- metadata %>% dplyr::mutate(Sample = orig.ident)
        seurat_obj@meta.data <- seurat_obj@meta.data %>%
          dplyr::mutate(Sample = orig.ident) %>%
          dplyr::left_join(
            metadata2 %>% dplyr::select(Sample, SRA_Number, Tissue, Approach),
            by = "Sample",
            suffix = c(".seurat", ".meta")
          ) %>%
          dplyr::mutate(
              Tissue = coalesce(Tissue.meta, Tissue.seurat),
              Approach = Approach.seurat # explicitly choose the seurat's Approach
            ) %>%
          dplyr::select(-Tissue.meta, - Tissue.seurat, - Approach.meta, - Approach.seurat)

        # Aggregate sample-level metadata including all cluster columns
        sample_level_metadata <- seurat_obj@meta.data %>%
          dplyr::group_by(Sample) %>%
          dplyr::summarize(
            SRA_Number = dplyr::first(SRA_Number),
            Tissue = dplyr::first(Tissue),
            Approach = dplyr::first(Approach),
            across(contains("cluster"), ~ dplyr::first(.)),
            .groups = "drop"
          )
        annot_fields <- c("SRA_Number", "Tissue", "Approach")
      } else {
        # Use the simple approach
        metadata2 <- metadata %>% dplyr::mutate(Sample = orig.ident)
        seurat_obj@meta.data <- seurat_obj@meta.data %>%
          dplyr::mutate(Sample = orig.ident) %>%
          dplyr::left_join(
            metadata2 %>% dplyr::select(Sample, Tissue),
            by = "Sample",
            suffix = c(".seurat", ".meta")
          ) %>%
          dplyr::mutate(Tissue = coalesce(Tissue.meta, Tissue.seurat)) %>%
          dplyr::select(-Tissue.meta, - Tissue.seurat)

        sample_level_metadata <- seurat_obj@meta.data %>%
          dplyr::group_by(Sample) %>%
          dplyr::summarize(
            Tissue = dplyr::first(Tissue),
            SRA_Number = dplyr::first(SRA),
            Approach = dplyr::first(Approach),
            Sample = dplyr::first(Sample),
            across(contains("cluster"), ~ dplyr::first(.)),
            .groups = "drop"
          )
        annot_fields <- c("Tissue", "SRA_Number", "Approach", "Sample")
      }

      # Identify cluster columns from the aggregated sample metadata:
      # any column name that contains "cluster" (case-insensitive),
      # excluding any that contain "random_group".
      cluster_cols <- colnames(sample_level_metadata)[
        grepl("cluster", colnames(sample_level_metadata), ignore.case = TRUE) &
        !grepl("random_group", colnames(sample_level_metadata), ignore.case = TRUE)
      ]
      if (length(cluster_cols) == 0) {
        stop("No cluster columns (excluding Random_Group) found in the sample-level metadata.")
      }

      # Partition cluster columns into hierarchical vs. non-hierarchical
      hierarchical_cols <- cluster_cols[grepl("hierarchical", cluster_cols, ignore.case = TRUE)]
      non_hierarchical_cols <- setdiff(cluster_cols, hierarchical_cols)

      # Filter and order sample-level metadata to match the BDM matrix samples
      sample_annotations <- sample_level_metadata %>%
        dplyr::filter(Sample %in% rownames(bdm_matrix)) %>%
        dplyr::arrange(match(Sample, rownames(bdm_matrix)))

      if (!all(sample_annotations$Sample %in% rownames(bdm_matrix))) {
        stop("Sample identifiers in the metadata do not match the BDM matrix.")
      }
      bdm_matrix <- bdm_matrix[sample_annotations$Sample, sample_annotations$Sample, drop = FALSE]

      # Generate annotation colors for each annotation column
      annotation_colors <- lapply(sample_annotations, function(col) {
        if (is.factor(col) || is.character(col)) {
          unique_vals <- unique(col)
          setNames(viridis::viridis(length(unique_vals)), unique_vals)
        } else {
          NULL
        }
      })
      annotation_colors <- annotation_colors[!sapply(annotation_colors, is.null)]

      # Build the heatmap annotation; include extra fields if available
      ha <- HeatmapAnnotation(
        df = sample_annotations %>% dplyr::select(all_of(c(annot_fields, cluster_cols))),
        col = annotation_colors,
        annotation_name_side = "left"
      )

      # Overwrite any existing PDF file
      pdf_path <- file.path(plots_folder, paste0("BDM_Heatmaps_", dataset_name, ".pdf"))
      if (file.exists(pdf_path)) {
        file.remove(pdf_path)
      }
      pdf(pdf_path, height = 20, width = 25)

      # --- Hierarchical clusters: produce dendrogram-based and natural ordering heatmaps ---
      for (cluster_col in hierarchical_cols) {
        cluster_annotations <- sample_annotations %>% dplyr::pull(!!rlang::sym(cluster_col))
        cluster_colors <- annotation_colors[[cluster_col]]
        # Select the appropriate hierarchical clustering object:
        if (grepl("tissue", cluster_col, ignore.case = TRUE)) {
          hc <- hc_tissue
          k_val <- k_tissue
        } else if (grepl("approach", cluster_col, ignore.case = TRUE)) {
          hc <- hc_approach
          k_val <- k_approach
        } else if (grepl("sra", cluster_col, ignore.case = TRUE)) {
          hc <- hc_sra
          k_val <- k_sra
        } else {
          hc <- hc_sra
          k_val <- NA
        }

        sample_order <- order.dendrogram(as.dendrogram(hc))
        ordered_annotations <- cluster_annotations[sample_order]
        dend <- as.dendrogram(hc)
        colored_dendrogram <- dendextend::set(dend, labels_col = cluster_colors[ordered_annotations])
        # Color branches using the unique colors in the ordered annotations
        unique_colors <- unique(cluster_colors[ordered_annotations])
        colored_dendrogram <- dendextend::color_branches(colored_dendrogram, k = length(unique_colors),
                                                        col = unique_colors)

        # Plot dendrogram
        plot(colored_dendrogram,
            main = paste("Dendrogram -", cluster_col, dataset_name),
            xlab = "", ylab = "Height", sub = ""
        )
        text(
          x = 1, y = 0, labels = paste("Number of Clusters (k) =", k_val),
          pos = 4, col = "blue", font = 2
        )

        # Reorder matrix based on the dendrogram ordering
        reordered_bdm_matrix <- bdm_matrix[sample_order, sample_order]

        # Draw hierarchical clustering heatmaps (raw and z-score)
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

        # Also produce natural ordering heatmaps for hierarchical clusters
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

      # --- Non-hierarchical clusters: consolidate into one set of natural ordering heatmaps ---
      if (length(non_hierarchical_cols) > 0) {
        # We simply choose the first non-hierarchical column as representative
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
      cat("All heatmaps and dendrograms saved to:", pdf_path, "\n")
    }

    # Call the consolidated function for BDM
    if (!is.null(bdm_matrix)) {
      if (exists("metadata_filtered") && !is.null(metadata_filtered)) {
        generate_heatmaps(
          dataset_name = paste0(dataset_name, "_BDM"),
          metadata = metadata_filtered,
          seurat_obj = seurat_obj,
          bdm_matrix = bdm_matrix,
          plots_folder = plots_folder,
          hc_tissue = hc_tissue_bdm,
          hc_sra = hc_sra_bdm,
          hc_approach = hc_approach_bdm,
          k_tissue = k_tissue,
          k_sra = k_sra,
          k_approach = k_approach
        )
      } else {
        generate_heatmaps(
          dataset_name = paste0(dataset_name, "_BDM"),
          metadata = metadata,
          seurat_obj = seurat_obj,
          bdm_matrix = bdm_matrix,
          plots_folder = plots_folder,
          hc_tissue = hc_tissue_bdm,
          hc_sra = hc_sra_bdm,
          hc_approach = hc_approach_bdm,
          k_tissue = k_tissue,
          k_sra = k_sra,
          k_approach = k_approach
        )
      }
    }

    # Similarly, call the same function for SDM if needed.
    if (!is.null(sdm_matrix)) {
      if (exists("metadata_filtered") && !is.null(metadata_filtered)) {
        generate_heatmaps(
          dataset_name = paste0(dataset_name, "_SDM"),
          metadata = metadata_filtered,
          seurat_obj = seurat_obj,
          bdm_matrix = sdm_matrix,
          plots_folder = plots_folder,
          hc_tissue = hc_tissue_sdm,
          hc_sra = hc_sra_sdm,
          hc_approach = hc_approach_sdm,
          k_tissue = k_tissue,
          k_sra = k_sra,
          k_approach = k_approach
        )
      } else {
        generate_heatmaps(
          dataset_name = paste0(dataset_name, "_SDM"),
          metadata = metadata,
          seurat_obj = seurat_obj,
          bdm_matrix = sdm_matrix,
          plots_folder = plots_folder,
          hc_tissue = hc_tissue_sdm,
          hc_sra = hc_sra_sdm,
          hc_approach = hc_approach_sdm,
          k_tissue = k_tissue,
          k_sra = k_sra,
          k_approach = k_approach
        )
      }
    }


    # ---------------------------
    # Heatmap Generation for Landscape Matrix
    # ---------------------------
    if (!is.null(landscape_matrix)) {
      # Before passing to generate_heatmaps(), check the range.
      matrix_range <- range(landscape_matrix, na.rm = TRUE)
      if (abs(matrix_range[2] - matrix_range[1]) < .Machine$double.eps) {
        # Add a tiny offset to create a non-zero range.
        landscape_matrix <- landscape_matrix + runif(length(landscape_matrix), 0, 1e-6)
      }

      # Landscape-based clusters (both hierarchical and spectral) should be stored in seurat_obj@meta.data.
      # We use a separate heatmap generation call for the landscape matrix.
      if (exists("metadata_filtered") && !is.null(metadata_filtered)) {
        generate_heatmaps(
          dataset_name = paste0(dataset_name, "_Landscape"),
          metadata = metadata_filtered,
          seurat_obj = seurat_obj,
          bdm_matrix = landscape_matrix, # Here, landscape_matrix acts like the distance matrix.
          plots_folder = plots_folder,
          hc_tissue = hc_tissue_landscape,
          hc_sra = hc_sra_landscape,
          hc_approach = hc_approach_landscape,
          k_tissue = k_tissue,
          k_sra = k_sra,
          k_approach = k_approach
        )
      } else {
        generate_heatmaps(
          dataset_name = paste0(dataset_name, "_Landscape"),
          metadata = metadata,
          seurat_obj = seurat_obj,
          bdm_matrix = landscape_matrix,
          plots_folder = plots_folder,
          hc_tissue = hc_tissue_landscape,
          hc_sra = hc_sra_landscape,
          hc_approach = hc_approach_landscape,
          k_tissue = k_tissue,
          k_sra = k_sra,
          k_approach = k_approach
        )
      }
    }

  }

  # ---------------------------
  # Improved Betti Curve Statistical Comparisons: Between Tissues and Clusters
  # ---------------------------
  if (run_comparative_statistical_analysis_ph) {

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

                data.frame(Group1 = group_names[1], Group2 = group_names[2], Value = value, Label = label)
              }))
            }, error = function(e) {
              log_message(paste("Error collecting data for", stat, "in", dimension, "for", comparison_type, ":", conditionMessage(e)))
              NULL
            })

            if (is.null(pairwise_data) || nrow(pairwise_data) == 0) {
              log_message(paste("No valid data for statistic", stat, "in", dimension, "for", comparison_type, "- Skipping heatmap."))
              next
            }

            # Accumulate stats
            pairwise_stat <- pairwise_data[, c("Group1", "Group2", "Value")]
            colnames(pairwise_stat)[3] <- stat
            pairwise_stats_df[[stat]] <- pairwise_stat

            # Build heatmap matrix
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
              ggplot(heatmap_data, aes(x = Group1, y = Group2, fill = Value)) +
            geom_tile(color = "white") +
            fill_scale +
            labs(
              title = paste("Pairwise Heatmap for", comparison_type, "-", dimension, "-", stat),
              subtitle = subtitle_text,
              x = "Group 1", y = "Group 2", fill = fill_label
            ) +
            theme_minimal(base_size = 10) +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 7),
              plot.subtitle = element_text(size = 5),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8)
            ) + {
                if (tolower(stat) %in% c("ks_result", "ks_stat", "perm_p")) {
                  geom_text(aes(label = Label), color = "black", size = 3)
                }
              }
            }, error = function(e) {
              log_message(paste("Error generating plot for", stat, "in", dimension, "for", comparison_type, ":", conditionMessage(e)))
              NULL
            })

            if (!is.null(heatmap_plot)) {
              heatmap_list[[paste(dimension, stat, sep = "_")]] <- heatmap_plot
            }
          }

          # Save wide-format CSV for this dimension
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

        # Save combined heatmaps for this comparison type
        if (length(heatmap_list) > 0) {
          combined <- tryCatch({
            gridExtra::marrangeGrob(grobs = heatmap_list, nrow = 2, ncol = 2)
          }, error = function(e) {
            log_message(paste("Error combining heatmaps for", comparison_type, ":", conditionMessage(e)))
            NULL
          })
          if (!is.null(combined)) {
            combined_file <- file.path(plot_output_dir, paste0(dataset_name, "_", comparison_type, "_combined_heatmaps.pdf"))
            ggsave(combined_file, combined, width = 12, height = 9)
            log_message(paste("Saved combined heatmaps for", comparison_type, "to", combined_file))
          }
        }

        log_message(paste("Completed heatmap generation for", comparison_type))
      }

      # Save final master combined CSV
      if (length(all_combined_csvs) > 0) {
        all_stats_df <- do.call(rbind, all_combined_csvs)
        master_file <- file.path(plot_output_dir, paste0(dataset_name, "_all_pairwise_statistics_combined.csv"))
        write.csv(all_stats_df, master_file, row.names = FALSE)
        log_message(paste("Saved final combined CSV for all stats to", master_file))
      }

      log_message("Completed heatmap generation for all statistics.")
    }

    # dataset_name <- data_iterations_bonemarrow[[1]]$name
    # pd_list <- readRDS(data_iterations_bonemarrow[[1]]$pd_list)
    # seurat_obj <- data_iterations_bonemarrow[[1]]$seurat_obj

    results_list <- list()

    # Existing cluster comparisons
    cluster_columns <- c(
      paste0("seurat_cluster_", tolower(dataset_name)),
      paste0("kmeans_cluster_", tolower(dataset_name), "_tissue"),
      paste0("kmeans_cluster_", tolower(dataset_name), "_sra"),
      paste0("kmeans_cluster_", tolower(dataset_name), "_approach"),
      paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_tissue"),
      paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_sra"),
      paste0("hierarchical_cluster_bdm_ph_", tolower(dataset_name), "_approach"),
      paste0("spectral_cluster_bdm_", tolower(dataset_name), "_tissue"),
      paste0("spectral_cluster_bdm_", tolower(dataset_name), "_sra"),
      paste0("spectral_cluster_bdm_", tolower(dataset_name), "_approach"),
      paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_tissue"),
      paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_sra"),
      paste0("hierarchical_cluster_sdm_ph_", tolower(dataset_name), "_approach"),
      paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_tissue"),
      paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_sra"),
      paste0("hierarchical_cluster_landscape_ph_", tolower(dataset_name), "_approach"),
      paste0("spectral_cluster_sdm_", tolower(dataset_name), "_tissue"),
      paste0("spectral_cluster_sdm_", tolower(dataset_name), "_sra"),
      paste0("spectral_cluster_sdm_", tolower(dataset_name), "_approach"),
      paste0("spectral_cluster_landscape_", tolower(dataset_name), "_tissue"),
      paste0("spectral_cluster_landscape_", tolower(dataset_name), "_sra"),
      paste0("spectral_cluster_landscape_", tolower(dataset_name), "_approach")
    )

    # Main comparisons
    if ("Tissue" %in% colnames(seurat_obj@meta.data)) {
      results_list$tissue_comparison <- tryCatch(
        compute_and_compare_betti_curves(
          pd_list = pd_list,
          landscape_list = landscape_list,
          seurat_obj = list(seurat_obj),
          group_by_col = "Tissue",
          grid_points = 500,
          num_permutations = 10000,
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
          num_permutations = 10000,
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
          num_permutations = 10000,
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
          num_permutations = 10000,
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
            num_permutations = 10000,
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

    cluster_cols <- grep("hierarchical", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
    if (length(cluster_cols) > 0) {
      for (cc in cluster_cols) {
        results_list[[paste0("cluster_comparison_", cc)]] <- tryCatch(
          compute_and_compare_betti_curves(
            pd_list = pd_list,
            landscape_list = landscape_list,
            seurat_obj = list(seurat_obj),
            group_by_col = cc,
            grid_points = 500,
            num_permutations = 10000,
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
      cat("Warning: No 'cluster' column found in seurat_obj metadata. Skipping cluster-based comparison.\n")
    }

    cluster_cols <- grep("kmeans", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
    if (length(cluster_cols) > 0) {
      for (cc in cluster_cols) {
        results_list[[paste0("cluster_comparison_", cc)]] <- tryCatch(
          compute_and_compare_betti_curves(
            pd_list = pd_list,
            landscape_list = landscape_list,
            seurat_obj = list(seurat_obj),
            group_by_col = cc,
            grid_points = 500,
            num_permutations = 10000,
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
      cat("Warning: No 'cluster' column found in seurat_obj metadata. Skipping cluster-based comparison.\n")
    }

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

    # Remove NA values for each metric
    euler_null_stats_effect <- euler_null_stats_effect[!is.na(euler_null_stats_effect)]
    euler_null_stats_ks <- euler_null_stats_ks[!is.na(euler_null_stats_ks)]
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
      generate_heatmaps_for_all_statistics(
        results_list = results_list,
        dataset_name = dataset_name,
        plots_folder = plots_folder,
        null_summary = overall_null_summary
      ),
      error = function(e) {
        cat("Error in generating heatmaps:", conditionMessage(e), "\n")
      }
    )

  }

  # ---------------------------
  # Enhanced Cluster Comparison with P-values
  # - Degenerate => forcibly set metric to 0 for bar charts, skip significance (p-value = NA).
  # - Uses "∅" in the plots to mark degenerate (not testable).
  # - VI can exceed 1 on the y-axis (auto-scale).
  # - Metrics: ARI, NMI, Jaccard, VI, Purity.
  # - For seurat_cluster methods, new metadata columns are created with suffixes _sra, _tissue, and _approach.
  #   These new columns are used as separate methods throughout the analysis.
  # ---------------------------
enhanced_cluster_comparison_with_pvals <- function(
  seurat_obj,
  dataset_name = "MyData",
  plots_folder = "plots",
  run_comparative_metrics = TRUE,
  include_silhouette = FALSE,
  verbose = TRUE,
  num_cores = 16,
  SRA_col = "SRA_Number"
  ) {
  # -- Initialization & Logging helper --
  if (!run_comparative_metrics) {
    message("Comparative metrics not requested. Exiting function.")
    return(NULL)
  }
  log_message <- function(msg) {
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
    
  if (!dir.exists(plots_folder)) {
    dir.create(plots_folder, recursive = TRUE)
    log_message(sprintf("Created plots folder: %s", plots_folder))
  }
  log_message("Starting enhanced cluster comparison analysis with p-values.")

  # -- Load libraries --
  log_message("Loading libraries …")
  library(mclust); library(aricode); library(clusterSim)
  library(ggplot2); library(entropy); library(reshape2)
  library(tidyr); library(dplyr); library(grid); library(parallel)
  if (include_silhouette) {
    library(cluster)
    log_message("Loaded silhouette library.")
  }
  log_message("All libraries loaded.")

  # -- Identify random-group columns --
  random_group_cluster_cols <- grep(
    "random_group",
    colnames(seurat_obj@meta.data),
    ignore.case = TRUE,
    value = TRUE
  )
  log_message(sprintf(
    "Found %d random-group columns for baseline distribution.",
    length(random_group_cluster_cols)
  ))

  # -- Silhouette setup --
  if (include_silhouette) {
    log_message("Setting up silhouette distance matrix …")
    if (!"pca" %in% names(seurat_obj@reductions)) {
      stop("PCA reduction not found—needed for silhouette.")
    }
    emb <- Embeddings(seurat_obj, "pca")[, 1:10, drop = FALSE]
    dist_mat <- dist(emb)
    metric_silhouette <- function(clusters, unused) {
      labs <- as.integer(as.factor(clusters))
      if (length(unique(labs)) < 2) return(NA_real_)
      mean(silhouette(labs, dist_mat)[, "sil_width"])
    }
    log_message("Silhouette distance matrix ready.")
  }

  # -- Define clustering method names --
  log_message("Defining original clustering methods …")
  dataset_lower <- tolower(dataset_name)
  original_methods <- c(
    paste0("seurat_cluster_", dataset_lower),
    paste0("kmeans_cluster_", dataset_lower, c("_tissue", "_sra", "_approach")),
    paste0("hierarchical_cluster_", c("bdm_ph","sdm_ph","landscape_ph"), "_", dataset_lower, c("_tissue","_sra","_approach")),
    paste0("spectral_cluster_", c("bdm","sdm","landscape"), "_", dataset_lower, c("_tissue","_sra","_approach"))
  ) %>% unlist()
  log_message(sprintf("Original methods: %s", paste(original_methods, collapse = ", ")))

  # -- Expand seurat_cluster into three suffixes --
  log_message("Expanding seurat_cluster methods to include _sra/_tissue/_approach …")
  final_methods <- c()
  for (meth in original_methods) {
    if (!meth %in% colnames(seurat_obj@meta.data)) next
    if (grepl("^seurat_cluster_", meth, ignore.case = TRUE)) {
      for (suf in c("_sra", "_tissue", "_approach")) {
        newcol <- paste0(meth, suf)
        seurat_obj@meta.data[[newcol]] <- seurat_obj@meta.data[[meth]]
        final_methods <- c(final_methods, newcol)
      }
    } else {
      final_methods <- c(final_methods, meth)
    }
  }
  log_message(sprintf("Final methods to evaluate: %s", paste(final_methods, collapse = ", ")))

  # -- Reference column helper --
  get_reference_column <- function(m) {
    if      (grepl("_sra$", m))    SRA_col
    else if (grepl("_tissue$", m)) "Tissue"
    else if (grepl("_approach$", m)) "Approach"
    else                           "Tissue"
  }

  # -- Metric functions & normalization --
  DEGENERATE_SYMBOL <- "X"
  safe_metric <- function(x, y, fun) {
    if (length(unique(x)) == 1 || length(unique(y)) == 1) {
      list(val = 0,
    else {
      list(val = fun(x, y), label = "", is_deg = FALSE)
    }
  }
  metric_ari    <- function(x, y) adjustedRandIndex(x, y)
  metric_nmi    <- function(x, y) NMI(x, y)
  metric_jac    <- function(x, y) { ctbl <- table(x, y); a <- sum(choose(rowSums(ctbl),2)); b <- sum(choose(colSums(ctbl),2)); c <- sum(choose(ctbl,2)); if(a+b-c==0)0 else c/(a+b-c) }
  metric_vi     <- function(x, y) { xn <- as.numeric(as.factor(x)); yn <- as.numeric(as.factor(y)); joint <- table(xn,yn)/length(xn); H1 <- entropy.empirical(rowSums(joint)); H2 <- entropy.empirical(colSums(joint)); HJ <- entropy.empirical(as.vector(joint)); 2*HJ-(H1+H2) }
  metric_purity <- function(x, y) { if(length(x)==0) return(NA); s <- sapply(unique(x), function(cl) max(table(y[x==cl]))); sum(s)/length(x) }
  safe_normalize <- function(obs, rm, bigger=TRUE, eps=1e-8) {
    if (!is.finite(rm)) return(NA)
    if (bigger) { d <- 1 - rm; if(abs(d)<eps)d<-eps; (obs - rm)/d }
    else        { if(abs(rm)<eps)rm<-eps; (rm - obs)/rm }
  }
  significance_code <- function(p) {
    if      (is.na(p)) "" 
    else if (p < 0.001) "***" 
    else if (p < 0.01)  "**" 
    else if (p < 0.05)  "*" 
    else if (p < 0.1)   "." 
    se                 ""
  }

  # -- Compute real vs random metrics --
  log_message("Beginning metric computation loop …")
  raw_comparison <- data.frame()
  for (method in final_methods) {
    log_messa
     (!method %in% colnames(seurat_obj@meta.data)) {
      log_message(" → column not found, skipping.")
      next
    }
    ref_col <- get_reference_column(method)
    if (!ref_col %in% colnames(seurat_obj@meta.data)) {
      log_message(sprintf(" → reference '%s' missing, skipping.", ref_col))
      next
    }
    real_cluster <- as.character(seurat_obj@meta.data[[method]])
    truth        <- as.character(seurat_obj@meta.data[[ref_col]])
    if (length(real_cluster) != length(truth)) {
      log_message(" → length mismatch, skipping.")
      next
    }

    # real metrics
    ari_obj <- safe_metric(real_cluster, truth, metric_ari)
    nmi_obj <- safe_metric(real_cluster, truth, metric_nmi)
    jac_obj <- safe_metric(real_cluster, truth, metric_jac)
    vi_obj  <- safe_metric(real_cluster, truth, metric_vi)
    pur_obj <- safe_metric(real_cluster, truth, metric_purity)

    # silhouette real
    if (include_silhouette) {
      sil_obj <- if (length(unique(real_cluster)) < 2) {
        list(val=0,label=DEGENERATE_SYMBOL,is_deg=TRUE)
      } else {
        list(val=metric_silhouette(real_cluster,NULL), label="", is_deg=FALSE)
      }
    }

    # random baseline
    rand_list <- mclapply(random_group_cluster_cols, function(rg) {
      if (!rg%in%colnames(seurat_obj@meta.data)) return(NULL)
      rc <- as.character(seurat_obj@meta.data[[rg]])
      if (length(rc)!=length(real_cluster)) return(NULL)
      out <- list(
        ari=metric_ari(real_cluster,rc),
        nmi={v<-metric_nmi(real_cluster,rc);if(is.na(v))0 else v},
        jac={v<-metric_jac(real_cluster,rc);if(is.na(v))0 else v},
        vi={v<-metric_vi(real_cluster,rc);if(is.na(v))0 else v},
        pur={v<-metric_purity(real_cluster,rc);if(is.na(v))0 else v}
      )
      if (include_silhouette) {
        out$sil <- {v<-metric_silhouette(real_cluster,rc); if(is.na(v))0 else v}
      }
      out
    }, mc.cores = num_cores)
    rand_list <- Filter(Negate(is.null), rand_list)

    # helper to summarise
    stat <- function(name, robj, bigger=TRUE) {
      if (length(rand_list)==0) return(list(mean=NA,pval=NA))
      vals <- sapply(rand_list, `[[`, name)
      mean_val <- mean(vals, na.rm=TRUE)
      p_val <- if (robj$is_deg) NA else if (bigger) mean(vals>=robj$val, na.rm=TRUE) else mean(vals<=robj$val, na.rm=TRUE)
      list(mean=mean_val, pval=p_val)
    }

    ari_stat <- stat("ari", ari_obj, TRUE)
    nmi_stat <- stat("nmi", nmi_obj, TRUE)
    jac_stat <- stat("jac", jac_obj, TRUE)
    vi_stat  <- stat("vi",  vi_obj,  FALSE)
    pur_stat <- stat("pur", pur_obj, TRUE)
    if (include_silhouette) sil_stat <- stat("sil", sil_obj, TRUE)

    # assemble row
    rowdf <- data.frame(
      Method           = method,
      ReferenceCol     = ref_col,
      ARI_Real         = ari_obj$val,
      ARI_RandMean     = ari_stat$mean,
      ARI_pval         = ari_stat$pval,
      ARI_flag         = ari_obj$label,

      NMI_Real         = nmi_obj$val,
      NMI_RandMean     = nmi_stat$mean,
      NMI_pval         = nmi_stat$pval,
      NMI_flag         = nmi_obj$label,

      Jaccard_Real     = jac_obj$val,
      Jaccard_RandMean = jac_stat$mean,
      Jaccard_pval     = jac_stat$pval,
      Jaccard_flag     = jac_obj$label,

      VI_Real          = vi_obj$val,
      VI_RandMean      = vi_stat$mean,
      VI_pval          = vi_stat$pval,
      VI_flag          = vi_obj$label,

      Purity_Real      = pur_obj$val,
      Purity_RandMean  = pur_stat$mean,
      Purity_pval      = pur_stat$pval,
      Purity_flag      = pur_obj$label,
      stringsAsFactors = FALSE
    )
    if (include_silhouette) {
      rowdf$Silhouette_Real     <- sil_obj$val
      rowdf$Silhouette_RandMean <- sil_stat$mean
      rowdf$Silhouette_pval     <- sil_stat$pval
      rowdf$Silhouette_flag     <- sil_obj$label
    }
    raw_comparison <- rbind(raw_comparison, rowdf)
  }
  log_message("Finished computing raw metrics for all methods.")

  # -- Adjust p-values --
  if (nrow(raw_comparison)>0) {
    pv_cols <- c("ARI_pval","NMI_pval","Jaccard_pval","VI_pval","Purity_pval")
    if (include_silhouette) pv_cols <- c(pv_cols,"Silhouette_pval")
    for (col in pv_cols) {
      raw_comparison[[paste0(col,"_adj_fdr")]]      <- p.adjust(raw_comparison[[col]], method="fdr")
      raw_comparison[[paste0(col,"_adj_bh")]]       <- p.adjust(raw_comparison[[col]], method="BH")
      raw_comparison[[paste0(col,"_adj_bonferroni")]] <- p.adjust(raw_comparison[[col]], method="bonferroni")
    }
    log_message("Adjusted p-values for multiple testing.")
  }

  # -- Write raw comparison CSV --
  raw_csv <- file.path(plots_folder, paste0(dataset_name,"_raw_comparison_metrics_with_pvals.csv"))
  write.csv(raw_comparison, raw_csv, row.names=FALSE)
  log_message(sprintf("Wrote raw metrics CSV: %s", raw_csv))

  # -- Compute and write normalized metrics CSV --
  norm_list <- lapply(seq_len(nrow(raw_comparison)), function(i) {
    rowi <- raw_comparison[i,]
    df <- data.frame(
      Method           = rowi$Method,
      ReferenceCol     = rowi$ReferenceCol,
      ARI_Observed     = rowi$ARI_Real,
      ARI_RandMean     = rowi$ARI_RandMean,
      ARI_Norm         = safe_normalize(rowi$ARI_Real, rowi$ARI_RandMean, TRUE),
      ARI_pval         = rowi$ARI_pval,
      ARI_pval_adj_fdr = rowi$ARI_pval_adj_fdr,
      ARI_flag         = rowi$ARI_flag,
      NMI_Observed     = rowi$NMI_Real,
      NMI_RandMean     = rowi$NMI_RandMean,
      NMI_Norm         = safe_normalize(rowi$NMI_Real, rowi$NMI_RandMean, TRUE),
      NMI_pval         = rowi$NMI_pval,
      NMI_pval_adj_fdr = rowi$NMI_pval_adj_fdr,
      NMI_flag         = rowi$NMI_flag,
      Jaccard_Observed = rowi$Jaccard_Real,
      Jaccard_RandMean = rowi$Jaccard_RandMean,
      Jaccard_Norm     = safe_normalize(rowi$Jaccard_Real, rowi$Jaccard_RandMean, TRUE),
      Jaccard_pval     = rowi$Jaccard_pval,
      Jaccard_pval_adj_fdr = rowi$Jaccard_pval_adj_fdr,
      Jaccard_flag     = rowi$Jaccard_flag,
      VI_Observed      = rowi$VI_Real,
      VI_RandMean      = rowi$VI_RandMean,
      VI_Norm          = safe_normalize(rowi$VI_Real, rowi$VI_RandMean, FALSE),
      VI_pval          = rowi$VI_pval,
      VI_pval_adj_fdr  = rowi$VI_pval_adj_fdr,
      VI_flag          = rowi$VI_flag,
      Purity_Observed  = rowi$Purity_Real,
      Purity_RandMean  = rowi$Purity_RandMean,
      Purity_Norm      = safe_normalize(rowi$Purity_Real, rowi$Purity_RandMean, TRUE),
      Purity_pval      = rowi$Purity_pval,
      Purity_pval_adj_fdr = rowi$Purity_pval_adj_fdr,
      Purity_flag      = rowi$Purity_flag,
      stringsAsFactors = FALSE
    )
    if (include_silhouette) {
      df$Silhouette_Observed    <- rowi$Silhouette_Real
      df$Silhouette_RandMean    <- rowi$Silhouette_RandMean
      df$Silhouette_Norm        <- safe_normalize(rowi$Silhouette_Real,rowi$Silhouette_RandMean, TRUE)
      df$Silhouette_pval        <- rowi$Silhouette_pval
      df$Silhouette_pval_adj_fdr<- rowi$Silhouette_pval_adj_fdr
      df$Silhouette_flag        <- rowi$Silhouette_flag
    }
    df
  })
  norm_df <- do.call(rbind, norm_list)
  norm_csv <- file.path(plots_folder, paste0(dataset_name,"_normalized_metrics_with_pvals.csv"))
  write.csv(norm_df, norm_csv, row.names=FALSE)
  log_message(sprintf("Wrote normalized metrics CSV: %s", norm_csv))

  # -- Generate PDF + individual PNG/SVG plots --
  pdf_file <- file.path(plots_folder, paste0("Combined_Cluster_Metrics_",dataset_name,"_with_pvals.pdf"))
  pdf(pdf_file, width=11, height=8.5)
  log_message(sprintf("Opened PDF device: %s", pdf_file))

  grid.newpage()
    text_explanation <- paste0(
     "Enhanced Cluster Comparison with P-values\n\n",
     "Metrics (ARI, NMI, Jaccard, VI, Purity, Silhouette):\n",
     " - ARI, NMI, Jaccard, Purity: bigger is better.\n",
     " - VI: smaller is better (can exceed 1).\n\n",
     "Random Baseline:\n",
     " - 'RandMean' is the average from 'random_group' replicates.\n",
     " - p-value: fraction of random replicates that exceed (or for VI, fall below) the observed.\n\n",
     "Degenerate Cases:\n",
     " - If method or reference has only one cluster, the metric is set to 0 and flagged as degenerate ('",
      DEGENERATE_SYMBOL, "').\n",
     " - p-value is set to NA (not testable).\n\n",
     "For seurat_cluster methods, separate columns were created for SRA, Tissue, and Approach.\n",
     "Normalized Score:\n",
     " - Bigger-is-better: (Obs - RandMean)/(1 - RandMean).\n",
     " - For VI: (RandMean - Obs)/RandMean.\n\n",
     "Significance Stars:\n",
     " *** p<0.001, ** p<0.01, * p<0.05, . p<0.1, '' if p>=0.1 or degenerate.\n\n",
     "Unified Method names indicate the reference used (displayed in parentheses).\n\n",
     "This PDF includes:\n",
     "1) Heatmaps of real metrics (with degenerate => 0, p-value=NA => '", DEGENERATE_SYMBOL, "').\n",
     "2) Normalized bar & bubble charts, and Observed vs. Random bar charts.\n",
     "   (Note: Negative normalized values appear if real < random baseline; VI can exceed 1.)"
     )
    grid.text(text_explanation, x = 0.05, y = 0.95, just = c("left", "top"),
            gp = gpar(fontfamily = "sans"))
  # Heatmap plotting function
  plot_metric_heatmap <- function(df, rc, pc, rm, fg, mname) {
    if (!all(c("Method",rc,pc,rm,fg) %in% colnames(df))) return(NULL)
    vals <- df[[rc]]; if (all(is.na(vals))) return(NULL)
    mat <- matrix(vals, nrow=nrow(df), ncol=1, dimnames=list(df$Method,mname))
    molten <- reshape2::melt(mat,varnames=c("Method","Metric"),value.name="RealValue")
    extra <- df[,c("Method",pc,rm,fg)]
    colnames(extra) <- c("Method","pval","RandMean","flag")
    molten <- merge(molten, extra, by="Method")
    molten$star <- significance_code(molten$pval)
    molten$label<-paste0(
      formatC(molten$RealValue,digits=3,format="f"), " ", molten$star,
      ifelse(molten$flag=="X"," X",""),
      "\nR=",formatC(molten$RandMean,digits=3,format="f"),
      "\nP=",formatC(molten$pval,digits=3,format="f")
    )
    p <- ggplot(molten,aes(x=Metric,y=Method,fill=RealValue))+
      geom_tile(color="grey70")+
      geom_text(aes(label=label),size=2)+
      scale_fill_gradient(low="white",high="green",na.value="grey50")+
      theme_minimal(base_size=8)+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
            axis.text.y=element_text(size=6))+
      labs(title=paste(mname,"Heatmap -",dataset_name),x=NULL,y=NULL)

    log_message(sprintf("Plotting heatmap for %s", mname))
    print(p)
    ggsave(file.path(plots_folder,paste0(dataset_name,"_heatmap_",mname,".png")),
           p, width=6, height=4, dpi=300)
    ggsave(file.path(plots_folder,paste0(dataset_name,"_heatmap_",mname,".svg")),
           p, width=6, height=4, device="svg")
    log_message(sprintf("Saved heatmap_%s.png/svg", mname))
    grid.newpage()
  }

  metrics <- c("ARI","NMI","Jaccard","VI","Purity")
  if (include_silhouette) metrics <- c(metrics,"Silhouette")
  for (m in metrics) {
    plot_metric_heatmap(
      raw_comparison,
      paste0(m,"_Real"),
      paste0(m,"_pval_adj_fdr"),
      paste0(m,"_RandMean"),
      paste0(m,"_flag"),
      m
    )
  }

  # Gather normalized long for bar & bubble
  gather_norm_metrics <- function(row, m) {
    df <- data.frame(
      Method     = row$Method,
      ReferenceCol = row$ReferenceCol,
      Metric     = m,
      Observed   = row[[paste0(m,"_Observed")]],
      RandMean   = row[[paste0(m,"_RandMean")]],
      Normalized = row[[paste0(m,"_Norm")]],
      p_adj      = row[[paste0(m,"_pval_adj_fdr")]],
      flag       = row[[paste0(m,"_flag")]],
      stringsAsFactors = FALSE
    )
    df$annot   <- paste0(ifelse(df$flag=="X","X",""), significance_code(df$p_adj))
    df$annot_y <- df$Normalized + ifelse(df$flag=="X", sign(df$Normalized)*0.05, 0)
    df
  }
  norm_long <- do.call(rbind,
    lapply(metrics, function(m)
      do.call(rbind,
        lapply(seq_len(nrow(norm_df)), function(i)
          gather_norm_metrics(norm_df[i,], m)))))

  norm_long <- norm_long %>%
    mutate(MethodRef = factor(paste0(Method,"\n(Ref=",ReferenceCol,")"),
                              levels=unique(paste0(Method,"\n(Ref=",ReferenceCol,")"))))

  # -- 1) Normalized bar chart --
  log_message("Plotting normalized bar chart …")
  p_bar <- ggplot(norm_long, aes(x=MethodRef,y=Normalized,fill=Metric))+
    geom_vline(xintercept=seq(1.5,length(unique(norm_long$MethodRef))-0.5,by=1),
               color="grey40")+
    geom_bar(stat="identity",position=position_dodge(width=0.9))+
    geom_text(aes(label=annot,y=annot_y),position=position_dodge(width=0.9),size=3)+
    geom_hline(yintercept=0,linetype="dotted")+
    theme_minimal(base_size=8)+
    theme(axis.text.x=element_text(angle=75,hjust=1,size=6),
          axis.text.y=element_text(size=6))+
    labs(title=paste("Normalized Metrics by Method -",dataset_name),
         x=NULL,y="Normalized Value",fill="Metric")
  print(p_bar)
  ggsave(file.path(plots_folder,paste0(dataset_name,"_normalized_bar.png")),
         p_bar,width=8,height=4,dpi=300)
  ggsave(file.path(plots_folder,paste0(dataset_name,"_normalized_bar.svg")),
         p_bar,width=8,height=4,device="svg")
  log_message("Saved normalized_bar.png/svg")
  grid.newpage()

  # -- 2) Bubble chart --
  log_message("Plotting bubble chart …")
  p_bubble <- ggplot(norm_long, aes(x=Metric,y=MethodRef,size=Normalized,color=p_adj))+
    geom_point()+ geom_text(aes(label=annot),vjust=-2,size=2)+
    scale_color_gradient(low="red",high="blue")+
    theme_minimal(base_size=8)+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
          axis.text.y=element_text(size=6))+
    labs(title=paste("Bubble Plot of Normalized Metrics -",dataset_name),
         x=NULL,y=NULL,size="Norm Val",color="Adj P")
  print(p_bubble)
  ggsave(file.path(plots_folder,paste0(dataset_name,"_bubble.png")),
         p_bubble,width=8,height=6,dpi=300)
  ggsave(file.path(plots_folder,paste0(dataset_name,"_bubble.svg")),
         p_bubble,width=8,height=6,device="svg")
  log_message("Saved bubble.png/svg")
  grid.newpage()

  # -- 3) Observed vs Random bar chart --
  log_message("Plotting Observed vs Random chart …")
  obs_rand <- norm_long %>%
    select(MethodRef,Metric,Observed,RandMean,annot) %>%
    pivot_longer(c(Observed,RandMean),names_to="Type",values_to="Value")
  p_obs <- ggplot(obs_rand, aes(x=MethodRef,y=Value,fill=Type))+
    geom_bar(stat="identity",position="dodge")+
    geom_text(data=filter(obs_rand,Type=="Observed"),
              aes(label=annot),position=position_dodge(width=0),size=2)+
    facet_wrap(~Metric,scales="free_y")+
    theme_minimal(base_size=8)+
    theme(axis.text.x=element_text(angle=90,hjust=1,size=6),
          axis.text.y=element_text(size=6))+
    labs(title=paste("Observed vs Random Baseline -",dataset_name),
         x=NULL,y=NULL,fill=NULL)
  print(p_obs)
  ggsave(file.path(plots_folder,paste0(dataset_name,"_obs_vs_random.png")),
         p_obs,width=10,height=6,dpi=300)
  ggsave(file.path(plots_folder,paste0(dataset_name,"_obs_vs_random.svg")),
         p_obs,width=10,height=6,device="svg")
  log_message("Saved obs_vs_random.png/svg")

  # -- Close PDF and finish --
  dev.off()
  log_message("Closed PDF device and completed all plot saving.")
  list(raw_comparison = raw_comparison, norm_metri
    

  results <- enhanced_cluster_comparison_with_pvals(
  seurat_obj = seurat_obj,
  dataset_name = dataset_name,
  plots_folder = paste0(results_folder, "/plots/withinIterationClusterComparison"),
  run_comparative_metrics = TRUE,
  verbose = TRUE,
  S
   )

  # # to run outside of loop if needed (ie, to adjust plots, debugging)
  # results <- enhanced_cluster_comparison_with_pvals(
  #   seurat_obj         = data_iterations_bonemarrow[[1]]$seurat_obj,
  #   dataset_name       = data_iterations_bonemarrow[[1]]$name,
  #   plots_folder       = paste0(results_folder, "/plots/withinIterationClusterComparison"),
  #   run_comparative_metrics = TRUE,
  #   verbose            = TRUE,
  #   SRA_col = SRA_col
  # )

  # results <- enhanced_cluster_comparison_with_pvals(
  #   seurat_obj         = data_iterations_bonemarrow[[2]]$seurat_obj,
  #   dataset_name       = data_iterations_bonemarrow[[2]]$name,
  #   plots_folder       = paste0(results_folder, "/plots/withinIterationClusterComparison"),
  #   run_comparative_metrics = TRUE,
  #   verbose            = TRUE,
  #   SRA_col = SRA_col
  # )

  # results <- enhanced_cluster_comparison_with_pvals(
  #   seurat_obj         = data_iterations_bonemarrow[[3]]$seurat_obj,
  #   dataset_name       = data_iterations_bonemarrow[[3]]$name,
  #   plots_folder       = paste0(results_folder, "/plots/withinIterationClusterComparison"),
  #   run_comparative_metrics = TRUE,
  #   verbose            = TRUE,
  #   SRA_col = SRA_col
  # )

  # results <- enhanced_cluster_comparison_with_pvals(
  #   seurat_obj         = data_iterations_bonemarrow[[4]]$seurat_obj,
  #   dataset_name       = data_iterations_bonemarrow[[4]]$name,
  #   plots_folder       = paste0(results_folder, "/plots/withinIterationClusterComparison"),
  #   run_comparative_metrics = TRUE,
  #   verbose            = TRUE,
  #   SRA_col = SRA_col
  # )

}

# ---------------------------
# Cross-Iteration, Multi-Cluster-Type Betti Curve and Cluster Metric Comparisons with Visualizations
# ---------------------------

# integrated_seurat <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/final_integrated_seurat.rds")
# integrated_seurat_bonemarrow <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/final_integrated_seurat_bonemarrow.rds")
# merged_seurat_unintegrated <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/merged_seurat_unintegrated.Rds")
# merged_seurat_unintegrated_bonemarrow <- readRDS("/mnt/e/Repositories/Jonah/PH_ClusteringApp/merged_seurat_unintegrated_bonemarrow.Rds")

# raw_bonemarrow_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "raw_bonemarrow_seurat_object.rds"))
# sct_individual_bonemarrow_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "sct_individual_bonemarrow_seurat_object.rds"))
# sct_whole_bonemarrow_seurat_objectsct_whole_bonemarrow_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "sct_whole_bonemarrow_seurat_object.rds"))
# integrated_bonemarrow_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "integrated_bonemarrow_seurat_object.rds"))

# raw_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "raw_seurat_object.rds"))
# sct_individual_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "sct_individual_seurat_object.rds"))
# sct_whole_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "sct_whole_seurat_object.rds"))
# integrated_seurat_object <- readRDS(paste0(results_folder, "/seurat_objects/", "integrated_seurat_object.rds"))

data_iterations <- list(
  list(
    name = "Raw",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "raw_seurat_object.rds")),
# seurat_obj = merged_seurat_unintegrated,
    assay = "RNA",
    bdm_matrix = "BDM_unintegrated_Raw.rds",
    sdm_matrix = "SDM_unintegrated_Raw.rds",
    pd_list = "PD_list_after_retries_unintegrated_Raw.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_Raw.Rds",
    landscape_list = "landscape_list_Raw.Rds"
    # expr_list = readRDS("expr_list_raw.Rds")
  ),
# Uncomment SCT_Individual iteration if needed
  list(
    name = "SCT_Individual",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_individual_seurat_object.rds")),
# seurat_obj = merged_seurat_unintegrated,
    assay = "SCT_Ind",
     bdm_matrix = "BDM_unintegrated_sctInd.rds",
    sdm_matrix = "SDM_unintegrated_sctInd.rds",
    pd_list = "PD_list_after_retries_unintegrated_sctInd.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctInd.Rds",
    landscape_list = "landscape_list_sctInd.Rds"
    # expr_list = readRDS("expr_list_sctInd.Rds")
  ),
# list(
#   name = "SCT_Whole",
#   seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_whole_seurat_object.rds")),
#   # seurat_obj = merged_seurat_unintegrated,
#   assay = "SCT",
#   bdm_matrix = "BDM_unintegrated_scTransformed_whole.rds",
#    sdm_matrix = "SDM_unintegrated_scTransformed_whole.rds",
#   pd_list = "PD_list_after_retries_unintegrated_sctWhole.rds",
#   landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctWhole.Rds",
#   landscape_list = "landscape_list_sctWhole.Rds",
#   expr_list = readRDS("expr_list_sctwhole.Rds")
# ),
  list(
    name = "Integrated",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "integrated_seurat_object.rds")),
# seurat_obj = integrated_seurat,
    assay = "integrated",
    bdm_matrix = "BDM_integrated.rds",
    sdm_matrix = "SDM_integrated.rds",
    pd_list = "PD_list_dim1_th-1_integrated.rds",
#variable_features = "./intermediate_seuratIntegration_results/integration_features.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_integrated.Rds",
    landscape_list = "landscape_list_integrated.Rds"
    # expr_list = readRDS("expr_list_integrated.Rds")
  )
)


data_iterations_bonemarrow <- list(
  list(
    name = "Raw_bonemarrow",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "raw_bonemarrow_seurat_object.rds")),
    assay = "RNA",
    bdm_matrix = "BDM_unintegrated_Raw_bonemarrow.rds",
    sdm_matrix = "SDM_unintegrated_Raw_bonemarrow.rds",
    pd_list = "PD_list_dim1_th-1_unintegrated_Raw_bonemarrow.Rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_Raw_bonemarrow.Rds",
    landscape_list = "landscape_list_Raw_bonemarrow.Rds"
  ),
  # Uncomment SCT_Individual iteration if needed
  list(
    name = "SCT_Individual_bonemarrow",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_individual_bonemarrow_seurat_object.rds")),
    assay = "SCT_Ind",
    bdm_matrix = "BDM_unintegrated_sctInd_bonemarrow.rds",
    
  pd_list = "PD_list_dim1_th-1_unintegrated_sctInd_bonemarrow.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctInd_bonemarrow.Rds",
    landscape_list = "landscape_list_sctInd_bonemarrow.Rds"
  ),
  list(
    name = "SCT_Whole_bonemarrow",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "sct_whole_bonemarrow_seurat_object.rds")),
    assay = "SCT",
    bdm_matrix = "BDM_unintegrated_scTransformed_whole_bonemarrow.rds",
    sdm_matrix = "SDM_unintegrated_scTransformed_whole_bonemarrow.rds",
    pd_list = "PD_list_dim1_th-1_unintegrated_sctWhole_bonemarrow.Rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_sctWhole_bonemarrow.Rds",
    landscape_list = "landscape_list_sctWhole_bonemarrow.Rds"
  ),
  list(
    name = "Integrated_bonemarrow",
    seurat_obj = readRDS(paste0(results_folder, "/seurat_objects/", "integrated_bonemarrow_seurat_object.rds")),
    assay = "integrated",
    bdm_matrix = "BDM_integrated_bonemarrow.rds",
    sdm_matrix = "SDM_integrated_bonemarrow.rds",
    pd_list = "PD_list_dim1_th-1_integrated_bonemarrow.rds",
    #variable_features = "./intermediate_seuratIntegration_results/integration_features_bonemarrow.rds",
    landscape_l2_distance_matrix = "landscape_l2_distance_matrix_integrated_bonemarrow.Rds",
    landscape_list = "landscape_list_integrated_bonemarrow.Rds"
  )
)

# ---------------------------
oss-Iteration Betti Curve and Cluster Metric Comparisons
# ---------------------------
 
  # Simple logging function.
log_message <- function(msg) {
   cat(sprintf("[%s] %s\n", Sys.time(), msg))
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
  valid_pd <- pd[pd[, 1] == dimension & pd[, 2] < pd[, 3], ]
  if (nrow(valid_pd) == 0) return(rep(0, grid_points))
  tau_grid <- seq(0, tau_max, length.out = grid_points)
  gaussian <- function(x, mu, sigma) {
    exp(-((x - mu)^2) / (2 * sigma^2)) / (sqrt(2 * pi) * sigma)
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
      any(sapply(betti_curves, function(x) !is.list(x) || !"mean" %in% names(x)))) {
    stop("Invalid Betti curves structure. Ensure each dimension has a 'mean' component.")
  }
  tryCatch({
    consistency_results <- lapply(seq_along(tau_vals), function(tau_idx) {
      tau <- tau_vals[tau_idx]
      euler_at_tau <- euler_curve$mean[tau_idx]
      betti_sum_at_tau <- sum(sapply(seq_along(dimensions), function(dim_idx) {
        if (!is.null(betti_curves[[dim_idx]]$mean))
          betti_curves[[dim_idx]]$mean[tau_idx] * (-1)^dimensions[dim_idx]
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
  # Load required package
  library(parallel)

  # Helper function to compute the integrated difference
  integrated_diff <- function(curve_diff, tau_vals) {
    delta_tau <- diff(tau_vals)
    integrated_sq <- sum(delta_tau * (curve_diff[-length(curve_diff)] ^ 2 +
                                        curve_diff[-1] ^ 2) / 2, na.rm = TRUE)
    sqrt(integrated_sq)
  }

  # Optional logging
  log_message <- function(msg) {
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
    }

  # Combine metadata and identify bootstrap columns
  combined_meta <- do.call(rbind, metadata_list)
  bs_cols <- grep("^random_group_bootstrap_", colnames(combined_meta),
                  value = TRUE, ignore.case = TRUE)

  # Parallelize over bootstrap columns using mclapply
  results <- mclapply(bs_cols, function(col) {
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

  # Filter out NULL results
  results <- results[!sapply(results, is.null)]
  null_effects <- sapply(results, function(x) x$eff_size)
  null_ks <- sapply(results, function(x) x$ks_stat)

  # Return summary statistics
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

# --- Hashing for Caching ---
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

# --- Bootstrapping for PD Curves ---
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

# --- Euler Curve Computation (bootstrapping per dimension then aggregating) ---
compute_euler_curve <- function(pd_subset, dimensions, group_name, dataset_name, grid_points, tau_max, base_sigma, n_bootstrap, results_folder) {
  log_message("Computing Euler curve with bootstrapping.")
  individual_euler <- lapply(seq_along(pd_subset), function(i) {
    pd <- pd_subset[[i]]
    if (is.matrix(pd) && nrow(pd) > 0) {
      euler <- rep(0, grid_points)
      for (d in dimensions) {
        curv <- get_or_compute_curve(pd, d, dataset_name, base_sigma, grid_points, tau_max, results_folder)
        euler <- euler + ((-1)^d) * curv
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
    euler_mean <- euler_mean + ((-1)^dimensions[i]) * bootstrapped[[i]]$mean
    euler_lower <- euler_lower + ((-1)^dimensions[i]) * bootstrapped[[i]]$lower
    euler_upper <- euler_upper + ((-1)^dimensions[i]) * bootstrapped[[i]]$upper
  }
  list(mean = euler_mean, lower = euler_lower, upper = euler_upper, individual_euler_curves = individual_euler)
}

# --- Integrated Difference Helper (Trapezoidal Rule) ---
integrated_diff <- function(curve_diff, tau_vals) {
  if (length(curve_diff) < 2) return(0)
  delta_tau <- diff(tau_vals)
  sqrt(sum(delta_tau * (curve_diff[-length(curve_diff)]^2 + curve_diff[-1]^2) / 2, na.rm = TRUE))
}

# --- Summary Statistic for Individual Landscape Curve (AUC) ---
summary_landscape <- function(landscape_obj, grid = seq(0, 1, length.out = grid_points)) {
  curve <- compute_landscape_curve(landscape_obj, grid = grid)
  delta <- diff(grid)
  auc <- sum(((curve[-length(curve)] + curve[-1]) / 2) * delta)
  auc
}

bootstrap_null_stats_landscape <- function(landscape_list, n_bootstrap = 50, grid_points = 500) {
  # We assume summary_landscape() is defined. It takes a landscape object and returns a summary statistic
  # (for example, the AUC of the landscape curve computed over a grid of length grid_points).
  # If summary_landscape() is not defined, here is a simple version:
  if (!exists("summary_landscape")) {
    summary_landscape <- function(landscape_obj, grid) {
      # For demonstration, we return the area under the aggregated curve (computed using compute_landscape_curve)
      agg <- compute_landscape_curve(landscape_obj, grid = grid)
      delta <- diff(grid)
      sum(((agg[-length(agg)] + agg[-1]) / 2) * delta)
    }
  }
  
  # Compute a summary statistic (e.g., AUC) for each landscape object.
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
    
    # Compute KS statistic between two groups’ empirical distributions.
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


# --- Null Stats for Individual Landscape Curves via Bootstrapping ---
compute_null_stats_individual_landscapes <- function(landscape_list, n_bootstrap = 50, grid_points = 500) {
  N <- length(landscape_list)
  if (N < 2) return(NULL)
  
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

# Compute the Landscape Curve for a Single Persistence Diagram
compute_landscape_curve <- function(landscape_obj, grid = seq(0, 1, length.out = grid_points)) {
  # If no landscape object is provided, return a zero curve.
  if (is.null(landscape_obj)) return(rep(0, length(grid)))
  
  # Extract dimension 0 and dimension 1 components.
  # If the components are matrices (each row representing a level), compute the row means.
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
  
  # Compute the aggregated landscape curve by taking the Euclidean norm.
  # This combines the information from both dimensions.
  aggregated_curve <- sqrt(curve0^2 + curve1^2)
  
  # If the computed curve length does not match the desired grid length,
  # interpolate the curve to the grid.
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


# --- Aggregated Landscape Curve for a Group ---
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

# --- Null Stats for Aggregated Landscape Curves ---
compute_null_stats_aggregated_landscapes <- function(pd_list, landscape_list, tau_vals, n_bootstrap = 50, num_cores = 8, grid_points = 500) {
  integrated_diff_local <- function(curve_diff, grid) {
    delta <- diff(grid)
    sqrt(sum(delta * (curve_diff[-length(curve_diff)]^2 + curve_diff[-1]^2) / 2, na.rm = TRUE))
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

# --- Bootstrap Aggregated Landscape Curve (for deriving null stat annotations) ---
bootstrap_aggregated_landscape_curve <- function(individual_curves, grid_points, n_bootstrap = 100) {
  boot_samples <- replicate(n_bootstrap, {
    sampled <- sample(individual_curves, replace = TRUE)
    rowMeans(do.call(cbind, sampled), na.rm = TRUE)
  })
  mean_curve <- rowMeans(boot_samples, na.rm = TRUE)
  lower_curve <- apply(boot_samples, 1, quantile, probs = 0.025, na.rm = TRUE)
  upper_curve <- apply(boot_samples, 1, quantile, probs = 0.975, na.rm = TRUE)
  list(
    mean = mean_curve,
    lower = lower_curve,
    upper = upper_curve,
    bootstrap_curves = boot_samples
  )
}

analyze_and_plot_curves <- function(curve1, curve2, grp1, grp2, curve_type, tau_vals,
                                      alpha, effect_threshold, num_permutations,
                                      null_effect_stats = NULL, null_ks_stats = NULL,
                                      verbose = TRUE) {
  # Local logging helper.
  log_message <- function(msg) { 
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg)) 
  }
  
  # Integrated difference helper using trapezoidal rule.
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
    sqrt(sum(delta_tau * (curve_diff[-length(curve_diff)]^2 + curve_diff[-1]^2) / 2, na.rm = TRUE))
  }
  
  # Build a data frame for plotting.
  build_plot_data <- function(curve, grp) {
    data.frame(Tau = tau_vals, Mean = curve$mean, 
               Lower = curve$lower, Upper = curve$upper, Group = grp)
  }
  
  pd1 <- build_plot_data(curve1, grp1)
  pd2 <- build_plot_data(curve2, grp2)
  plot_data <- rbind(pd1, pd2)
  
  # --- Select individual curves based on curve type ---
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
  
  # --- Compute aggregated mean curves ---
  if (!is.null(inds1) && !is.null(inds2) && is.matrix(inds1) && is.matrix(inds2) &&
      ncol(inds1) > 0 && ncol(inds2) > 0) {
    m1 <- rowMeans(inds1, na.rm = TRUE)
    m2 <- rowMeans(inds2, na.rm = TRUE)
  } else {
    log_message("Falling back to using the aggregated 'mean' fields.")
    m1 <- curve1$mean
    m2 <- curve2$mean
  }
  
  # Defensive: if too short, replace with zeros.
  if (length(m1) < 2) m1 <- rep(0, length(tau_vals))
  if (length(m2) < 2) m2 <- rep(0, length(tau_vals))
  
  # Interpolate m1 and m2 if necessary.
  if (length(m1) != length(tau_vals)) {
    log_message("Interpolating m1 to match tau_vals length.")
    m1 <- approx(x = seq(0, 1, length.out = length(m1)), y = m1,
                 xout = seq(0, 1, length.out = length(tau_vals)), rule = 2)$y
  }
  if (length(m2) != length(tau_vals)) {
    log_message("Interpolating m2 to match tau_vals length.")
    m2 <- approx(x = seq(0, 1, length.out = length(m2)), y = m2,
                 xout = seq(0, 1, length.out = length(tau_vals)), rule = 2)$y
  }
  
  # --- Compute Wasserstein Distance ---
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
  
  # --- Permutation Test and KS Testing ---
  if (!is.null(inds1) && !is.null(inds2) && is.matrix(inds1) && is.matrix(inds2) &&
      ncol(inds1) > 0 && ncol(inds2) > 0) {
    all_inds <- cbind(inds1, inds2)
    n1 <- ncol(inds1)
    n_total <- ncol(all_inds)
    perm_diffs <- replicate(num_permutations, {
      perm <- sample(n_total)
      m1p <- rowMeans(all_inds[, perm[1:n1], drop = FALSE], na.rm = TRUE)
      m2p <- rowMeans(all_inds[, perm[(n1+1):n_total], drop = FALSE], na.rm = TRUE)
      integrated_diff_local(m1p - m2p, tau_vals)
    })
    perm_p <- mean(perm_diffs >= obs_diff)
    
    ks_res <- lapply(1:nrow(inds1), function(i) {
      res <- tryCatch(ks.test(inds1[i, ], inds2[i, ]), error = function(e) NULL)
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
  
  # --- Append Null Stats for Euler and Landscape (if provided) ---
  if (grepl("Euler", curve_type, ignore.case = TRUE) || grepl("Landscape", curve_type, ignore.case = TRUE)) {
    # Try to extract the median, lower, and upper values.
    extract_null_values <- function(null_stats) {
      if (!is.null(null_stats$median)) {
        list(median = null_stats$median, lower = null_stats$lower, upper = null_stats$upper)
      } else if (!is.null(null_stats$quantiles)) {
        # Extract values from the quantiles vector.
        list(median = as.numeric(null_stats$quantiles["50%"]),
             lower  = as.numeric(null_stats$quantiles["2.5%"]),
             upper  = as.numeric(null_stats$quantiles["97.5%"]))
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
  
  x_scale <- if (all(tau_vals > 0)) scale_x_continuous(trans = "log10") else scale_x_continuous()
  plot_data$Group <- as.factor(plot_data$Group)
  max_x <- max(tau_vals)
  p <- ggplot(plot_data, aes(x = Tau, y = Mean, color = Group, group = Group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group), alpha = 0.2) +
    labs(title = paste(curve_type, "Comparison:", grp1, "vs.", grp2), 
         y = paste(curve_type, "Value")) +
    x_scale +
    annotate("text", x = max_x * 0.99, y = Inf, label = annot, hjust = 0, vjust = 1.2, size = 2) +
    coord_cartesian(clip = "off") +
    theme_minimal()
  
  list(stats = list(wasserstein = wass,
                      perm_p = perm_p,
                      effect_size = eff_size,
                      ks_stat = ks_combined_stat,
                      ks_adj_p = ks_adj_pval),
       plot = p)
}

# --- Pairwise Comparisons Function ---
perform_pairwise_comparisons <- function(curves, euler_curves, groups, tau_vals, dimensions = c(0,1),
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
  use_log_scale <- all(tau_vals > 0)
  
  build_plot_data <- function(curve, grp) {
    data.frame(Tau = tau_vals, Mean = curve$mean, Lower = curve$lower, Upper = curve$upper, Group = grp)
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
        ar <- analyze_and_plot_curves(c1, c2, grp1, grp2, "Landscape", tau_vals,
                                      alpha, effect_threshold, num_permutations,
                                      null_effect_stats, null_ks_stats)
        if (!is.null(ar$stats)) comp[["landscape"]] <- ar$stats
        if (!is.null(ar$plot)) plots <- c(plots, list(ar$plot))
      } else {
        for (d in dimensions) {
          c1 <- tryCatch(curves[[grp1]][[d+1]], error = function(e) NULL)
          c2 <- tryCatch(curves[[grp2]][[d+1]], error = function(e) NULL)
          if (is.null(c1) || is.null(c2)) next
          ar <- analyze_and_plot_curves(c1, c2, grp1, grp2, paste("Betti Dim", d), tau_vals,
                                        alpha, effect_threshold, num_permutations)
          if (!is.null(ar$stats)) comp[[paste0("dimension_", d)]] <- ar$stats
          if (!is.null(ar$plot)) plots <- c(plots, list(ar$plot))
        }
        e1 <- tryCatch(euler_curves[[grp1]], error = function(e) NULL)
        e2 <- tryCatch(euler_curves[[grp2]], error = function(e) NULL)
        if (!is.null(e1) && !is.null(e2)) {
          er <- analyze_and_plot_curves(e1, e2, grp1, grp2, "Euler", tau_vals,
                                        alpha, effect_threshold, num_permutations,
                                        null_effect_stats, null_ks_stats)
          comp[["euler"]] <- er$stats
          plots <- c(plots, list(er$plot))
        }
      }
      if (length(plots) > 0) {
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
        cairo_pdf(filename = plot_file_pdf, width = 10, height = 8 * length(plots))
        gridExtra::grid.arrange(grobs = plots, ncol = 1)
        dev.off()
        log_message(paste("Saved plot to", plot_file_pdf))
        
        plot_file_svg <- file.path(plot_output_dir, group_by_col, file_name_svg)
        ensure_directory(dirname(plot_file_svg))
        svglite::svglite(plot_file_svg, width = 10, height = 8 * length(plots))
        gridExtra::grid.arrange(grobs = plots, ncol = 1)
        dev.off()
        log_message(paste("Saved plot to", plot_file_svg))
      }
      res[[paste0(grp1, "_vs_", grp2)]] <- comp
    }
  }
  res
}

# --- Plot Aggregated Curves ---
plot_aggregated_curves <- function(group_curves, groups, dimensions, tau_vals, dataset_name, group_by_col, euler_curves = NULL, results_folder) {
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
      geom_line(size = 1, alpha = 0.5) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group), alpha = 0.2) +
      labs(title = paste(curve_type, "Curve", if(!is.null(dim)) paste("(Dim", dim, ")") else ""),
           x = "Tau", y = paste(curve_type, "Value")) +
      theme_minimal() +
      scale_x_continuous(trans = "log10")
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
      plots <- list(gen_plot(all_data, curve_label))
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
    if (is.null(dimensions)) {
      file_name_pdf <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, ".pdf")
      file_name_svg <- paste0(dataset_name, "_landscape_curves_by_", group_by_col, ".svg")
    } else {
      file_name_pdf <- paste0(dataset_name, "_betti_and_euler_curves_by_", group_by_col, ".pdf")
      file_name_svg <- paste0(dataset_name, "_betti_and_euler_curves_by_", group_by_col, ".svg")
    }
    plot_file_pdf <- file.path(out_dir, file_name_pdf)
    ensure_directory(dirname(plot_file_pdf))
    pdf(plot_file_pdf, width = 10, height = 8 * length(plots))
    gridExtra::grid.arrange(grobs = plots, ncol = 1)
    dev.off()
    log_message(paste("Saved aggregated plot to", plot_file_pdf))
    
    plot_file_svg <- file.path(out_dir, file_name_svg)
    ensure_directory(dirname(plot_file_svg))
    svglite::svglite(plot_file_svg, width = 10, height = 8 * length(plots))
    gridExtra::grid.arrange(grobs = plots, ncol = 1)
    dev.off()
    log_message(paste("Saved aggregated plot to", plot_file_svg))
  } else {
    log_message("No plots to save.")
  }
}


compute_null_stats_cross_iterations <- function(pd_list, metadata_list, tau_vals, 
                                                  dimensions = c(0,1),
                                                  dataset_name, base_sigma, grid_points, tau_max,
                                                  results_folder, verbose = TRUE, num_cores = 16,
                                                  landscape_list = NULL) {
  # Compute the integrated difference between curves using the trapezoidal rule.
  integrated_diff <- function(curve_diff, tau_vals) {
    delta_tau <- diff(tau_vals)
    integrated_sq <- sum(delta_tau * (curve_diff[-length(curve_diff)]^2 + curve_diff[-1]^2) / 2)
    sqrt(integrated_sq)
  }
  
  log_message <- function(msg) {
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
  }
  
  combined_meta <- do.call(rbind, metadata_list)
  bs_cols <- grep("^random_group_bootstrap_", colnames(combined_meta), value = TRUE, ignore.case = TRUE)
  
  # Compute Betti/Euler null statistics.
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
      test_result <- tryCatch(ks.test(individual_curves1[i, ], individual_curves2[i, ]), error = function(e) NULL)
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
  
  # Compute landscape null statistics if landscape_list is provided.
  landscape_nulls <- NULL
  if (!is.null(landscape_list)) {
    landscape_nulls <- bootstrap_null_stats_landscape(landscape_list, n_bootstrap = 50, grid_points = grid_points)
  }
  
  list(betti_euler = betti_euler_nulls, landscape = landscape_nulls)
}



#=== Function: compute_and_compare_betti_iterations() ===#
compute_and_compare_betti_iterations <- function(
  pd_list, metadata_list, group_by_col = NULL,
  dimensions = c(0, 1), grid_points = 500, num_permutations = 10000,
  bootstrap_samples = 100, dataset_name = NULL, comparison_type = "group",
  base_sigma = 1, num_cores = 20, results_folder = "results", verbose = TRUE,
  landscape_list = NULL
) {
  library(ggplot2)
  library(parallel)
  library(cowplot)
  library(dplyr)
  library(digest)
  library(transport)
  
  ensure_directory <- function(path) {
    if (length(path) == 0 || is.na(path) || path == "")
      stop("Invalid directory path provided: ", path)
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      log_message(paste("Created directory:", path))
    }
  }
  
  log_message <- function(msg) {
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
  }
  
  # Caching functions for PD-based curves.
  generate_cache_filename <- function(dataset_name, hash, results_folder) {
    path <- file.path(results_folder, "plots", "betti_plots", "betti_cache", dataset_name, paste0("cache_", hash, ".rds"))
    ensure_directory(dirname(path))
    path
  }
  
  compute_betti_curve <- function(pd, dimension, grid_points, tau_max, base_sigma) {
    pd <- as.data.frame(pd)
    valid_pd <- pd[pd[, 1] == dimension & pd[, 2] < pd[, 3], ]
    if (nrow(valid_pd) == 0) return(rep(0, grid_points))
    tau_grid <- seq(0, tau_max, length.out = grid_points)
    gaussian <- function(x, mu, sigma) {
      exp(-((x - mu)^2) / (2 * sigma^2)) / (sqrt(2 * pi) * sigma)
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
  
  # Landscape helper functions (e.g., compute_landscape_curve, bootstrap_aggregated_landscape_curve)
  # are assumed to be defined externally.
  
  tau_max <- max(unlist(lapply(pd_list, function(pd) max(pd[, 3], na.rm = TRUE))), na.rm = TRUE)
  tau_vals <- seq(0.0001, tau_max, length.out = grid_points)
  
  # Align metadata using "orig.ident.Iter" keys (which are assumed to match pd_list names).
  group_metadata_mapping <- do.call(rbind, lapply(seq_along(metadata_list), function(i) {
    metadata <- metadata_list[[i]]
    pd_names <- names(pd_list)[names(pd_list) %in% metadata$orig.ident.Iter]
    metadata[metadata$orig.ident.Iter %in% pd_names, ]
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
  
  # --- LANDSCAPE ANALYSIS ---  
  aggregated_landscape_curves_by_group <- NULL
  if (!is.null(landscape_list)) {
    aggregated_landscape_curves_by_group <- setNames(lapply(as.character(groups), function(group) {
      group_ids <- group_metadata_mapping$orig.ident.Iter[group_metadata_mapping[[group_by_col]] == group]
      grid <- seq(0, 1, length.out = grid_points)
      if (length(group_ids) == 0) {
        return(list(
          mean = rep(0, grid_points), lower = rep(0, grid_points),
          upper = rep(0, grid_points), individual_landscape_curves = list()
        ))
      }
      # Because landscape_list corresponds 1:1 with pd_list, we assume that the keys are identical.
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
  
  #--- Pairwise Comparisons for Betti/Euler ---  
  pairwise_results <- perform_pairwise_comparisons(
      curves = betti_curves_by_group, groups = groups, tau_vals = tau_vals, dimensions = dimensions,
      euler_curves = euler_curves_by_group, dataset_name = dataset_name, group_by_col = group_by_col,
      plot_output_dir = file.path(results_folder, "plots", "betti_plots", "cross_iteration_comparisons", dataset_name, "pairwise_comparisons"),
      null_effect_stats = bootstrap_nulls$betti_euler$null_effect, 
      null_ks_stats = bootstrap_nulls$betti_euler$null_ks,
      alpha = 0.05, effect_threshold = 0.5, num_permutations = 1000
  )

  
  #--- Pairwise Comparisons for Landscapes ---  
  landscape_pairwise_results <- NULL
  if (!is.null(aggregated_landscape_curves_by_group)) {
    landscape_nulls <- bootstrap_null_stats_landscape(landscape_list, n_bootstrap = 50, grid_points = grid_points)
    landscape_pairwise_results <- perform_pairwise_comparisons(
      curves = aggregated_landscape_curves_by_group, groups = groups, tau_vals = tau_vals,
      dimensions = NULL,  # Aggregated landscape curves are one-dimensional.
      dataset_name = dataset_name, group_by_col = group_by_col,
      plot_output_dir = file.path(results_folder, "plots", "landscape_curves", dataset_name, "pairwise_comparisons"),
      null_effect_stats = landscape_nulls$null_effect, null_ks_stats = landscape_nulls$null_ks,
      alpha = 0.05, effect_threshold = 0.5, num_permutations = 1000
    )
  }
  
  log_message("Generating aggregated plots for Betti/Euler curves.")
  plot_aggregated_curves(
    group_curves = betti_curves_by_group, groups = groups, dimensions = dimensions,
    tau_vals = tau_vals, dataset_name = dataset_name, group_by_col = group_by_col,
    euler_curves = euler_curves_by_group, results_folder = results_folder
  )
  
  if (!is.null(aggregated_landscape_curves_by_group)) {
    log_message("Generating aggregated landscape plots.")
    plot_aggregated_curves(
      group_curves = aggregated_landscape_curves_by_group, groups = groups, dimensions = NULL,
      tau_vals = seq(0, 1, length.out = grid_points), dataset_name = paste0(dataset_name, "_landscape"),
      group_by_col = group_by_col, euler_curves = NULL, results_folder = results_folder
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


#=== Wrapper: cross_iteration_comparison_with_betti() ===#

# The revised wrapper:
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
  library(dplyr)
  library(readr)
  library(parallel)
  
  group_values <- unique(unlist(lapply(data_iterations, function(iter) 
    unique(iter$seurat_obj@meta.data[[group_by_col]]))))
  
  log_message <- function(msg) {
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
  }
  
  reference_cell_names <- colnames(data_iterations[[1]]$seurat_obj)
  log_message("Matching cell names across iterations")
  for (i in seq_along(data_iterations)) {
    seurat_obj <- data_iterations[[i]]$seurat_obj
    if (!all(colnames(seurat_obj) == reference_cell_names)) {
      if (ncol(seurat_obj) == length(reference_cell_names)) {
        colnames(seurat_obj) <- reference_cell_names
        data_iterations[[i]]$seurat_obj <- seurat_obj
        log_message(paste0("Iteration ", data_iterations[[i]]$name, " renamed to match reference cell names."))
      } else {
        log_message(paste0("Iteration ", data_iterations[[i]]$name, " has a different number of cells; skipping renaming."))
      }
    }
  }
  
  process_group <- function(group_value) {
    log_message(paste("Processing", group_by_col, ":", group_value))
    
    combined_pd_list <- list()
    metadata_list <- list()
    combined_landscape_list <- list()  # To hold landscape data from each iteration
    
    iter_results <- mclapply(data_iterations, function(iter) {
      # Load pd_list from file.
      pd_list <- readRDS(iter$pd_list)
      if (!all(startsWith(names(pd_list), paste0(iter$name, "_")))) {
        valid_names <- paste0(iter$name, "_", names(pd_list))
        names(pd_list) <- valid_names
      } else {
        valid_names <- names(pd_list)
      }
      
      # Load landscape_list from file.
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
      
      group_metadata <- metadata_local[metadata_local[[group_by_col]] == group_value, ]
      group_pd_list <- pd_list[valid_names %in% group_metadata$orig.ident.Iter]
      
      if (length(group_pd_list) == 0) {
        log_message(paste("No persistence diagrams for", group_by_col, ":", group_value, "in iteration:", iter$name))
        return(NULL)
      }
      
      pd_names_filtered <- valid_names[valid_names %in% group_metadata$orig.ident.Iter]
      names(group_pd_list) <- pd_names_filtered
      
      aligned_metadata <- lapply(pd_names_filtered, function(nm) {
        rows <- group_metadata[group_metadata$orig.ident.Iter == nm, ]
        if (nrow(rows) > 0) return(rows[1, , drop = FALSE]) else return(NULL)
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
  
  results <- mclapply(group_values, process_group, mc.cores = num_cores)
  log_message("All group-specific comparisons completed.")
}


compare_all_iterations_with_mapped_methods_random_baseline_including_purity_pvals <- function(
  data_iterations,         # list of iteration objects, each with $name, $seurat_obj, and $pd_list (RDS path)
  method_mapping_list,     # a named list; each element is a mapping for one unified method (names are iteration names, values are column names)
  dims_use = 1:10,         # PCA dims for centroid matching
  reduction_name = "pca",
  n_random = 100,          # (for additional uses)
  results_folder = file.path(results_folder, "cross_iteration_comparisons"),
  seed = 42,
  dataset_name = "CrossIter",
  verbose = TRUE,
  n_cores = 8             # number of cores for parallel processing
) {
  ###############################
  # Helper: Logging
  ###############################
  log_message <- function(msg) {
    if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
  }
  log_message("Function started: compare_all_iterations_with_mapped_methods_random_baseline_including_purity_pvals")

  ###############################
  # Ensure results folder exists
  ###############################
  if (!dir.exists(results_folder)) {
    dir.create(results_folder, recursive = TRUE)
    log_message(sprintf("Created results folder: %s", results_folder))
  } else {
    log_message(sprintf("Results folder already exists: %s", results_folder))
  }

  ###############################
  # Load required Libraries
  ###############################
  library(mclust)     # ARI
  library(aricode)    # NMI
  library(entropy)    # Variation of information
  library(dplyr)
  library(parallel)   # mclapply
  library(reshape2)
  library(ggplot2)
  library(grid)
  library(tidyr)

  ###############################
  # Helper: Global Linked Cluster Mapping (method-specific)
  ###############################
  global_linked_cluster_mapping <- function(iteration_map, reference_iter_name, dims_use, reduction_name = "pca", mapping_sub, verbose = TRUE) {
    log_message <- function(msg) {
      if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
    }
    log_message(sprintf("Starting global linked cluster mapping for unified method using '%s' as reference.", reference_iter_name))

    reference_iter <- iteration_map[[reference_iter_name]]
    so_ref <- reference_iter$seurat_obj
    ref_col <- mapping_sub[[reference_iter_name]]
    if (is.null(ref_col) || !(ref_col %in% colnames(so_ref@meta.data))) {
      stop(sprintf("Reference iteration '%s' does not have the specified mapping column '%s'.", reference_iter_name, ref_col))
    }
    ref_clusters <- unique(so_ref@meta.data[[ref_col]])
    centroids_ref <- lapply(ref_clusters, function(cl) {
      cells <- rownames(so_ref@meta.data)[so_ref@meta.data[[ref_col]] == cl]
      if (length(cells) == 0) return(NULL)
      emb <- Embeddings(so_ref, reduction = reduction_name)[cells, dims_use, drop = FALSE]
      colMeans(emb, na.rm = TRUE)
    })
    names(centroids_ref) <- as.character(ref_clusters)
    log_message(sprintf("Reference iteration '%s' has %d clusters based on column '%s'.", reference_iter_name, length(centroids_ref), ref_col))
    so_ref@meta.data$LinkedCluster <- as.character(so_ref@meta.data[[ref_col]])
    iteration_map[[reference_iter_name]]$seurat_obj <- so_ref
    rm(so_ref, ref_clusters); gc()

    for (iter_name in names(iteration_map)) {
      if (iter_name == reference_iter_name) next
      iter <- iteration_map[[iter_name]]
      so <- iter$seurat_obj
      cur_col <- mapping_sub[[iter_name]]
      if (is.null(cur_col) || !(cur_col %in% colnames(so@meta.data))) {
        log_message(sprintf("Iteration '%s' missing mapping column '%s'; skipping global mapping for this iteration.", iter_name, cur_col))
        next
      }
      clusters <- unique(so@meta.data[[cur_col]])
      new_linked <- rep(NA, nrow(so@meta.data))
      names(new_linked) <- rownames(so@meta.data)
      for (cl in clusters) {
        cells <- rownames(so@meta.data)[so@meta.data[[cur_col]] == cl]
        if (length(cells) == 0) next
        emb <- Embeddings(so, reduction = reduction_name)[cells, dims_use, drop = FALSE]
        centroid_iter <- colMeans(emb, na.rm = TRUE)
        dist_vec <- sapply(centroids_ref, function(ref_centroid) {
          sqrt(sum((ref_centroid - centroid_iter)^2))
        })
        matched_label <- names(which.min(dist_vec))
        new_linked[cells] <- matched_label
      }
      so@meta.data$LinkedCluster <- new_linked
      iter$seurat_obj <- so
      iteration_map[[iter_name]] <- iter
      rm(so, clusters, new_linked); gc()
    }
    log_message("Global linked cluster mapping completed for this unified method.")
    return(iteration_map)
  }

  ###############################
  # Helper: Process Linked Clusters
  ###############################
  process_linked_clusters <- function(data_iterations, linked_cluster_col = "LinkedCluster", verbose = TRUE) {
    log_message <- function(msg) {
      if (verbose) cat(sprintf("[%s] %s\n", Sys.time(), msg))
    }
    new_iterations <- lapply(data_iterations, function(iter) {
      so <- iter$seurat_obj
      if (!(linked_cluster_col %in% colnames(so@meta.data))) {
        log_message(sprintf("Iteration '%s' does not have '%s'. Skipping linked cluster processing.", iter$name, linked_cluster_col))
        return(NULL)
      }
      valid_cells <- rownames(so@meta.data)[!is.na(so@meta.data[[linked_cluster_col]]) & so@meta.data[[linked_cluster_col]] != ""]
      if (length(valid_cells) == 0) {
        log_message(sprintf("Iteration '%s' has no cells with linked cluster assignments. Skipping.", iter$name))
        return(NULL)
      }
      pd_list <- readRDS(iter$pd_list)
      pd_list_subset <- pd_list[names(pd_list) %in% valid_cells]
      if (length(pd_list_subset) == 0) {
        log_message(sprintf("Iteration '%s' has no persistence diagrams for linked cells. Skipping.", iter$name))
        return(NULL)
      }
      iter$pd_list <- pd_list_subset
      iter$seurat_obj@meta.data <- so@meta.data[valid_cells, , drop = FALSE]
      rm(pd_list, pd_list_subset, valid_cells); gc()
      log_message(sprintf("Iteration '%s' processed for linked clusters.", iter$name))
      return(iter)
    })
    new_iterations <- Filter(Negate(is.null), new_iterations)
    return(new_iterations)
  }

  ###############################
  # Metric Functions
  ###############################
  variation_of_information <- function(c1, c2) {
    c1 <- as.numeric(as.factor(c1))
    c2 <- as.numeric(as.factor(c2))
    joint_dist <- table(c1, c2) / length(c1)
    H1 <- entropy.empirical(rowSums(joint_dist))
    H2 <- entropy.empirical(colSums(joint_dist))
    H_joint <- entropy.empirical(as.vector(joint_dist))
    H1 + H2 - 2 * H_joint
  }

  jaccard_index <- function(c1, c2) {
    mat1 <- outer(c1, c1, FUN = "==")
    mat2 <- outer(c2, c2, FUN = "==")
    intersect_count <- sum(mat1 & mat2)
    union_count <- sum(mat1 | mat2)
    if (union_count == 0) return(0)
    intersect_count / union_count
  }

  purity_single <- function(pred, truth) {
    total <- length(pred)
    if (total < 1) return(NA)
    s <- sapply(unique(pred), function(cl) {
      max(table(truth[pred == cl]))
    })
    sum(s) / total
  }

  safe_normalize <- function(obs, rand_mean, bigger_better = TRUE, epsilon = 1e-8) {
    if (is.na(rand_mean) || !is.finite(rand_mean)) return(NA)
    if (bigger_better) {
      denom <- 1 - rand_mean
      if (abs(denom) < epsilon) denom <- epsilon
      (obs - rand_mean) / denom
    } else {
      if (abs(rand_mean) < epsilon) rand_mean <- epsilon
      (rand_mean - obs) / rand_mean
    }
  }

  DEGENERATE_SYMBOL <- "X"
  safe_metric <- function(x, y, fun) {
    if (length(unique(x)) == 1 || length(unique(y)) == 1) {
      list(val = 0, label = DEGENERATE_SYMBOL, is_deg = TRUE)
    } else {
      list(val = fun(x, y), label = "", is_deg = FALSE)
    }
  }

  metric_ari <- function(x, y) adjustedRandIndex(x, y)
  metric_nmi <- function(x, y) NMI(x, y)
  metric_jac <- function(x, y) jaccard_index(x, y)
  metric_vi <- function(x, y) variation_of_information(x, y)
  metric_purity <- function(x, y) purity_single(x, y)

  ###############################
  # Build iteration map
  ###############################
  iteration_map <- setNames(data_iterations, vapply(data_iterations, `[[`, "", "name"))
  iter_names <- names(iteration_map)
  log_message(sprintf("Found %d iterations: %s", length(iter_names), paste(iter_names, collapse = ", ")))
  if (length(iter_names) < 2) stop("Need at least two iterations in data_iterations.")

  ###############################
  # Loop over each unified method in the mapping list
  ###############################
  results_list <- list()
  result_counter <- 1

  for (unified_method in names(method_mapping_list)) {
    log_message(sprintf("Processing unified method: '%s'", unified_method))
    mapping_sub <- method_mapping_list[[unified_method]]
    mapped_iter_names <- intersect(iter_names, names(mapping_sub))
    log_message(sprintf("Unified method '%s' is mapped to iterations: %s", unified_method, paste(mapped_iter_names, collapse = ", ")))
    if (length(mapped_iter_names) < 2) {
      log_message(sprintf("Warning: Unified method '%s' found in fewer than 2 iterations. Skipping.", unified_method))
      next
    }
    # Global linked cluster mapping for this unified method.
    reference_iter_name <- mapped_iter_names[1]
    log_message(sprintf("For unified method '%s', using iteration '%s' as reference for global mapping.", unified_method, reference_iter_name))
    iteration_map_method <- global_linked_cluster_mapping(iteration_map, reference_iter_name, dims_use, reduction_name, mapping_sub, verbose)
    data_iterations_method <- unname(iteration_map_method)
    rm(iteration_map_method); gc()

    # Pairwise comparisons for this unified method.
    for (i in seq_len(length(mapped_iter_names) - 1)) {
      for (j in (i + 1):length(mapped_iter_names)) {
        iterA_name <- mapped_iter_names[i]
        iterB_name <- mapped_iter_names[j]
        colA <- mapping_sub[[iterA_name]]
        colB <- mapping_sub[[iterB_name]]
        log_message(sprintf("Unified method '%s': Comparing iteration '%s' (col='%s') vs '%s' (col='%s')", 
                            unified_method, iterA_name, colA, iterB_name, colB))
        itA <- iteration_map[[iterA_name]]
        itB <- iteration_map[[iterB_name]]
        soA <- itA$seurat_obj
        soB <- itB$seurat_obj
        if (!(colA %in% colnames(soA@meta.data))) {
          log_message(sprintf("Warning: Column '%s' not found in iteration '%s'. Skipping.", colA, iterA_name))
          next
        }
        if (!(colB %in% colnames(soB@meta.data))) {
          log_message(sprintf("Warning: Column '%s' not found in iteration '%s'. Skipping.", colB, iterB_name))
          next
        }
        # Ensure PCA reduction exists; run PCA if needed.
        if (!(reduction_name %in% names(soA@reductions))) {
          log_message(sprintf("PCA reduction '%s' not found in iteration '%s'. Running PCA...", reduction_name, iterA_name))
          tryCatch({
            soA <- RunPCA(soA, npcs = max(dims_use), verbose = FALSE)
          }, error = function(e) {
            log_message(sprintf("Error running PCA on iteration '%s': %s. Skipping pair.", iterA_name, e$message))
            next
          })
        } else {
          log_message(sprintf("Using existing PCA reduction in iteration '%s'", iterA_name))
        }
        if (!(reduction_name %in% names(soB@reductions))) {
          log_message(sprintf("PCA reduction '%s' not found in iteration '%s'. Running PCA...", reduction_name, iterB_name))
          tryCatch({
            soB <- RunPCA(soB, npcs = max(dims_use), verbose = FALSE)
          }, error = function(e) {
            log_message(sprintf("Error running PCA on iteration '%s': %s. Skipping pair.", iterB_name, e$message))
            next
          })
        } else {
          log_message(sprintf("Using existing PCA reduction in iteration '%s'", iterB_name))
        }
        common_cells <- intersect(rownames(soA@meta.data), rownames(soB@meta.data))
        log_message(sprintf("Found %d common cells between iterations '%s' and '%s'", length(common_cells), iterA_name, iterB_name))
        if (length(common_cells) < 2) {
          log_message(sprintf("Warning: Less than 2 overlapping cells for '%s' and '%s'. Skipping.", iterA_name, iterB_name))
          next
        }
        embA <- Embeddings(soA, reduction = reduction_name)[common_cells, dims_use, drop = FALSE]
        embB <- Embeddings(soB, reduction = reduction_name)[common_cells, dims_use, drop = FALSE]
        metaA <- soA@meta.data[common_cells, , drop = FALSE]
        metaB <- soB@meta.data[common_cells, , drop = FALSE]
        metaA[[colA]] <- as.character(metaA[[colA]])
        metaB[[colB]] <- as.character(metaB[[colB]])
        rm(soA, soB); gc()

        # Centroid matching between iterations.
        log_message("Performing centroid matching between iterations")
        dfA <- data.frame(Cell = common_cells, Cluster = metaA[[colA]], embA, stringsAsFactors = FALSE)
        dfB <- data.frame(Cell = common_cells, Cluster = metaB[[colB]], embB, stringsAsFactors = FALSE)
        centroidsA <- dfA %>% group_by(Cluster) %>% summarize(across(where(is.numeric), mean), .groups = "drop")
        centroidsA$Cluster <- as.character(centroidsA$Cluster)
        centroidsB <- dfB %>% group_by(Cluster) %>% summarize(across(where(is.numeric), mean), .groups = "drop")
        centroidsB$Cluster <- as.character(centroidsB$Cluster)
        log_message(sprintf("Centroids: %d clusters in iteration A, %d in iteration B", nrow(centroidsA), nrow(centroidsB)))
        matched_B_labels <- list()
        for (kk in seq_len(nrow(centroidsA))) {
          ref_cluster_label <- centroidsA$Cluster[kk]
          ref_coords <- as.numeric(centroidsA[kk, -1])
          if (nrow(centroidsB) == 0) next
          dist_vec <- apply(centroidsB[, -1, drop = FALSE], 1, function(x) sum((x - ref_coords)^2))
          min_idx <- which.min(dist_vec)
          matched_cluster_label <- centroidsB$Cluster[min_idx]
          matched_B_labels[[matched_cluster_label]] <- ref_cluster_label
        }
        old_B_clusters <- metaB[[colB]]
        new_B_clusters <- sapply(old_B_clusters, function(b_val) {
          if (b_val %in% names(matched_B_labels)) {
            matched_B_labels[[b_val]]
          } else {
            b_val
          }
        })
        metaB[[colB]] <- new_B_clusters
        rm(dfA, dfB, centroidsA, centroidsB, matched_B_labels, dist_vec); gc()

        # Compute metrics.
        log_message("Computing metrics for the current iteration pair")
        clusterA <- metaA[[colA]]
        clusterB <- metaB[[colB]]
        ari_obj <- safe_metric(clusterA, clusterB, metric_ari)
        nmi_obj <- safe_metric(clusterA, clusterB, metric_nmi)
        jac_obj <- safe_metric(clusterA, clusterB, metric_jac)
        vi_obj  <- safe_metric(clusterA, clusterB, metric_vi)
        pur_obj <- safe_metric(clusterA, clusterB, metric_purity)
        ari_real <- ari_obj$val
        nmi_real <- nmi_obj$val
        jac_real <- jac_obj$val
        vi_real  <- vi_obj$val
        pur_real <- pur_obj$val

        # Compute null baseline from random_group columns.
        rand_group_cols_A <- grep("random_group", colnames(metaA), ignore.case = TRUE, value = TRUE)
        if (length(rand_group_cols_A) > 0) {
          rand_list_A <- mclapply(rand_group_cols_A, function(rg) {
            rc <- as.character(metaA[[rg]])
            if (length(rc) != length(clusterA)) return(NULL)
            list(
              ari = metric_ari(clusterA, rc),
              nmi = metric_nmi(clusterA, rc),
              jac = metric_jac(clusterA, rc),
              vi  = metric_vi(clusterA, rc),
              pur = metric_purity(clusterA, rc)
            )
          }, mc.cores = n_cores)
          rand_list_A <- Filter(Negate(is.null), rand_list_A)
          if (length(rand_list_A) > 0) {
            all_ari_A <- unlist(lapply(rand_list_A, function(x) x$ari))
            all_nmi_A <- unlist(lapply(rand_list_A, function(x) x$nmi))
            all_jac_A <- unlist(lapply(rand_list_A, function(x) x$jac))
            all_vi_A  <- unlist(lapply(rand_list_A, function(x) x$vi))
            all_pur_A <- unlist(lapply(rand_list_A, function(x) x$pur))
          } else {
            all_ari_A <- all_nmi_A <- all_jac_A <- all_vi_A <- all_pur_A <- numeric(0)
          }
        } else {
          all_ari_A <- all_nmi_A <- all_jac_A <- all_vi_A <- all_pur_A <- numeric(0)
        }
        rand_group_cols_B <- grep("random_group", colnames(metaB), ignore.case = TRUE, value = TRUE)
        if (length(rand_group_cols_B) > 0) {
          rand_list_B <- mclapply(rand_group_cols_B, function(rg) {
            rc <- as.character(metaB[[rg]])
            if (length(rc) != length(clusterB)) return(NULL)
            list(
              ari = metric_ari(clusterB, rc),
              nmi = metric_nmi(clusterB, rc),
              jac = metric_jac(clusterB, rc),
              vi  = metric_vi(clusterB, rc),
              pur = metric_purity(clusterB, rc)
            )
          }, mc.cores = n_cores)
          rand_list_B <- Filter(Negate(is.null), rand_list_B)
          if (length(rand_list_B) > 0) {
            all_ari_B <- unlist(lapply(rand_list_B, function(x) x$ari))
            all_nmi_B <- unlist(lapply(rand_list_B, function(x) x$nmi))
            all_jac_B <- unlist(lapply(rand_list_B, function(x) x$jac))
            all_vi_B  <- unlist(lapply(rand_list_B, function(x) x$vi))
            all_pur_B <- unlist(lapply(rand_list_B, function(x) x$pur))
          } else {
            all_ari_B <- all_nmi_B <- all_jac_B <- all_vi_B <- all_pur_B <- numeric(0)
          }
        } else {
          all_ari_B <- all_nmi_B <- all_jac_B <- all_vi_B <- all_pur_B <- numeric(0)
        }
        all_ari <- c(all_ari_A, all_ari_B)
        all_nmi <- c(all_nmi_A, all_nmi_B)
        all_jac <- c(all_jac_A, all_jac_B)
        all_vi  <- c(all_vi_A, all_vi_B)
        all_pur <- c(all_pur_A, all_pur_B)
        ari_randgroup <- if (length(all_ari) > 0) mean(all_ari, na.rm = TRUE) else NA
        nmi_randgroup <- if (length(all_nmi) > 0) mean(all_nmi, na.rm = TRUE) else NA
        jac_randgroup <- if (length(all_jac) > 0) mean(all_jac, na.rm = TRUE) else NA
        vi_randgroup  <- if (length(all_vi) > 0) mean(all_vi, na.rm = TRUE) else NA
        pur_randgroup <- if (length(all_pur) > 0) mean(all_pur, na.rm = TRUE) else NA
        rm(rand_list_A, rand_list_B, all_ari_A, all_nmi_A, all_jac_A, all_vi_A, all_pur_A, 
           all_ari_B, all_nmi_B, all_jac_B, all_vi_B, all_pur_B); gc()

        ari_pval <- if (ari_obj$is_deg || length(all_ari) == 0) NA_real_ else mean(all_ari >= ari_real, na.rm = TRUE)
        nmi_pval <- if (nmi_obj$is_deg || length(all_nmi) == 0) NA_real_ else mean(all_nmi >= nmi_real, na.rm = TRUE)
        jac_pval <- if (jac_obj$is_deg || length(all_jac) == 0) NA_real_ else mean(all_jac >= jac_real, na.rm = TRUE)
        vi_pval  <- if (vi_obj$is_deg  || length(all_vi)  == 0) NA_real_ else mean(all_vi <= vi_real, na.rm = TRUE)
        pur_pval <- if (pur_obj$is_deg || length(all_pur) == 0) NA_real_ else mean(all_pur >= pur_real, na.rm = TRUE)

        iter_pair <- paste0(iterA_name, "_vs_", iterB_name)
        row_res <- data.frame(
          UnifiedMethod    = unified_method,
          IterationA       = iterA_name,
          IterationB       = iterB_name,
          IterPair         = iter_pair,
          ARI_Real         = ari_real,
          ARI_RandMean     = ari_randgroup,
          ARI_pval         = ari_pval,
          ARI_flag         = ari_obj$label,
          NMI_Real         = nmi_real,
          NMI_RandMean     = nmi_randgroup,
          NMI_pval         = nmi_pval,
          NMI_flag         = nmi_obj$label,
          Jaccard_Real     = jac_real,
          Jaccard_RandMean = jac_randgroup,
          Jaccard_pval     = jac_pval,
          Jaccard_flag     = jac_obj$label,
          VI_Real          = vi_real,
          VI_RandMean      = vi_randgroup,
          VI_pval          = vi_pval,
          VI_flag          = vi_obj$label,
          Purity_Real      = pur_real,
          Purity_RandMean  = pur_randgroup,
          Purity_pval      = pur_pval,
          Purity_flag      = pur_obj$label,
          stringsAsFactors = FALSE
        )
        results_list[[result_counter]] <- row_res
        result_counter <- result_counter + 1
        rm(metaA, metaB, embA, embB, common_cells); gc()
        log_message(sprintf("Stored results for iteration pair: %s", iter_pair))
      }
    }

    if (result_counter == 1) stop("No valid cross-iteration comparisons found.")

    ###############################
    # Adjust p-values for multiple comparisons
    ###############################
    log_message("Adjusting p-values for multiple comparisons")
    all_results <- do.call(rbind, results_list)
    for (col in c("ARI_pval", "NMI_pval", "Jaccard_pval", "VI_pval", "Purity_pval")) {
      all_results[[paste0(col, "_adj_fdr")]] <- p.adjust(all_results[[col]], method = "fdr")
      all_results[[paste0(col, "_adj_bh")]]  <- p.adjust(all_results[[col]], method = "BH")
      all_results[[paste0(col, "_adj_bonferroni")]] <- p.adjust(all_results[[col]], method = "bonferroni")
    }

    raw_csv <- file.path(results_folder, paste0(dataset_name, "_raw_comparison_metrics_with_pvals.csv"))
    write.csv(all_results, raw_csv, row.names = FALSE)
    log_message(sprintf("Wrote raw comparison metrics to CSV: %s", raw_csv))

    ###############################
    # Compute normalized improvements
    ###############################
    norm_list <- lapply(seq_len(nrow(all_results)), function(i) {
      rowi <- all_results[i, ]
      data.frame(
        UnifiedMethod    = rowi$UnifiedMethod,
        IterationA       = rowi$IterationA,
        IterationB       = rowi$IterationB,
        IterPair         = rowi$IterPair,
        ARI_Observed     = rowi$ARI_Real,
        ARI_RandMean     = rowi$ARI_RandMean,
        ARI_Norm         = safe_normalize(rowi$ARI_Real, rowi$ARI_RandMean, bigger_better = TRUE),
        ARI_pval         = rowi$ARI_pval,
        ARI_flag         = rowi$ARI_flag,

        NMI_Observed     = rowi$NMI_Real,
        NMI_RandMean     = rowi$NMI_RandMean,
        NMI_Norm         = safe_normalize(rowi$NMI_Real, rowi$NMI_RandMean, bigger_better = TRUE),
        NMI_pval         = rowi$NMI_pval,
        NMI_flag         = rowi$NMI_flag,

        Jaccard_Observed = rowi$Jaccard_Real,
        Jaccard_RandMean = rowi$Jaccard_RandMean,
        Jaccard_Norm     = safe_normalize(rowi$Jaccard_Real, rowi$Jaccard_RandMean, bigger_better = TRUE),
        Jaccard_pval     = rowi$Jaccard_pval,
        Jaccard_flag     = rowi$Jaccard_flag,

        VI_Observed      = rowi$VI_Real,
        VI_RandMean      = rowi$VI_RandMean,
        VI_Norm          = safe_normalize(rowi$VI_Real, rowi$VI_RandMean, bigger_better = FALSE),
        VI_pval          = rowi$VI_pval,
        VI_flag          = rowi$VI_flag,

        Purity_Observed  = rowi$Purity_Real,
        Purity_RandMean  = rowi$Purity_RandMean,
        Purity_Norm      = safe_normalize(rowi$Purity_Real, rowi$Purity_RandMean, bigger_better = TRUE),
        Purity_pval      = rowi$Purity_pval,
        Purity_flag      = rowi$Purity_flag,
        stringsAsFactors = FALSE
      )
    })
    norm_df <- do.call(rbind, norm_list)
    norm_csv <- file.path(results_folder, paste0(dataset_name, "_normalized_metrics_with_pvals.csv"))
    write.csv(norm_df, norm_csv, row.names = FALSE)
    log_message(sprintf("Wrote normalized metrics to CSV: %s", norm_csv))

    ###############################
    # Visualization
    ###############################
    pdf_file <- file.path(results_folder, paste0("Combined_Iteration_Metrics_", dataset_name, "_with_pvals.pdf"))
    pdf(pdf_file, width = 11, height = 8.5)

    text_explanation <- paste0(
      "Enhanced Cross-Iteration Comparison with P-values\n\n",
      "Metrics (ARI, NMI, Jaccard, VI, Purity):\n",
      " - ARI, NMI, Jaccard, Purity: bigger is better.\n",
      " - VI: smaller is better (can exceed 1).\n\n",
      "Random Baseline (Null) from random_group columns:\n",
      " - 'RandMean' is the average from random_group columns.\n",
      " - p-value: fraction of random group null values meeting the criterion.\n\n",
      "Normalized Score:\n",
      " - For ARI, NMI, Jaccard, Purity: (Obs - RandMean)/(1 - RandMean).\n",
      " - For VI: (RandMean - Obs)/RandMean.\n\n",
      "Significance Stars:\n",
      " *** p<0.001, ** p<0.01, * p<0.05, . p<0.1, '' if p>=0.1 or degenerate.\n\n",
      "Unified Method names (displayed in IterPair) indicate the mapping.\n\n",
      "This PDF includes:\n",
      "1) Heatmaps of real metrics (annotated with null baseline & adjusted p-value),\n",
      "2) Normalized bar charts, bubble charts, and Observed vs. Random bar charts."
    )
    grid.text(text_explanation, x = 0.05, y = 0.95, just = c("left", "top"),
              gp = gpar(fontfamily = "sans"))

    plot_metric_heatmap <- function(df, real_col, randmean_col, pval_col, metric_name = "Metric", digits = 3) {
      molten <- data.frame(
        IterPair = df$IterPair,
        UnifiedMethod = df$UnifiedMethod,
        RealValue = df[[real_col]],
        RandMean = df[[randmean_col]],
        Pval = df[[pval_col]],
        stringsAsFactors = FALSE
      )
      p <- ggplot(molten, aes(x = UnifiedMethod, y = IterPair, fill = RealValue)) +
        geom_tile(color = "grey70") +
        geom_text(aes(label = paste0(
          formatC(RealValue, digits = digits, format = "f"), "\nR=",
          formatC(RandMean, digits = digits, format = "f"), "\nAdjP=",
          formatC(Pval, digits = digits, format = "f")
        )), size = 2.8) +
        scale_fill_gradient(low = "white", high = "blue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste(metric_name, "Heatmap -", dataset_name),
            subtitle = "Fill = Real, annotated with null baseline & adjusted p-value",
            x = "Unified Method", y = "Iteration Pair")
      return(p)
    }

    metrics <- c("ARI", "NMI", "Jaccard", "VI", "Purity")
    for (mname in metrics) {
      real_col   <- paste0(mname, "_Real")
      rand_col   <- paste0(mname, "_RandMean")
      pval_col   <- paste0(mname, "_pval_adj_fdr")
      p_hm <- plot_metric_heatmap(all_results, real_col, rand_col, pval_col, mname)
      if (!is.null(p_hm)) {
        print(p_hm)
        svg_file <- file.path(results_folder, paste0(dataset_name, "_Heatmap_", mname, ".svg"))
        ggsave(svg_file, plot = p_hm, device = "svg", width = 8, height = 6)
      }
    }

    # Prepare normalized data for bar and bubble charts.
    # Here we create an annotation column using the degenerate flag.
    norm_long <- norm_df %>% mutate(IterPair = as.character(IterPair))
    norm_long <- norm_long %>% mutate(
      # For ARI we check the flag; if degenerate, prepend "X" to the significance code.
      annot = ifelse(ARI_flag == DEGENERATE_SYMBOL,
                     paste0("X", significance_code(ARI_pval)),
                     significance_code(ARI_pval))
    )
    norm_long <- norm_long %>% mutate(MethodRef = UnifiedMethod)
    unique_methods <- unique(norm_long$MethodRef)
    norm_long$MethodRef <- factor(norm_long$MethodRef, levels = unique_methods)
    y_breaks <- seq(-1, 1.2, 0.2)

    # Normalized Metrics Bar Chart for ARI with adjusted text position.
    p_bar <- ggplot(norm_long, aes(x = MethodRef, y = ARI_Norm, fill = "ARI")) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = annot, vjust = ifelse(ARI_Norm >= 0, -0.2, 1.2)),
                position = position_dodge(width = 0.9), size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0.2, 0.2))) +
      facet_wrap(~IterPair, scales = "free_x") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Normalized Metrics by Unified Method & IterPair -", dataset_name),
          x = "Unified Method", y = "Normalized Value", fill = "Metric")
    print(p_bar)
    ggsave(file.path(results_folder, paste0("Normalized_Metrics_Bar_", dataset_name, ".svg")),
           p_bar, device = "svg", width = 9, height = 6)

    p_bubble <- ggplot(norm_long, aes(x = "ARI", y = MethodRef, size = ARI_Norm, color = ARI_Norm)) +
      geom_point() +
      geom_text(aes(label = paste0("AdjP=", formatC(ARI_pval, digits = 3, format = "f"))),
                vjust = -1, size = 3) +
      facet_wrap(~IterPair, scales = "free_y") +
      scale_size_continuous(range = c(1, 10)) +
      theme_minimal() +
      labs(title = paste("Bubble Plot of Normalized Metrics -", dataset_name),
          x = "Metric", y = "Unified Method", size = "Norm Value", color = "Norm Value")
    print(p_bubble)
    ggsave(file.path(results_folder, paste0("Bubble_Plot_", dataset_name, ".svg")),
           p_bubble, device = "svg", width = 9, height = 6)

    obs_rand <- all_results %>% 
      select(IterPair, UnifiedMethod, ARI_Real, ARI_RandMean, ARI_pval) %>% 
      pivot_longer(cols = c("ARI_Real", "ARI_RandMean"), names_to = "Type", values_to = "Value")
    p_obs_bar <- ggplot(obs_rand, aes(x = UnifiedMethod, y = Value, fill = Type)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(data = obs_rand %>% filter(Type == "ARI_Real"),
                aes(label = paste0("AdjP=", formatC(ARI_pval, digits = 3, format = "f")), vjust = -0.2),
                position = position_dodge(width = 0.9), size = 3) +
      facet_wrap(~IterPair, scales = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
      coord_cartesian(clip = "off") +
      theme_minimal() +
      labs(title = paste("Observed vs. Random Baseline -", dataset_name),
          x = "Unified Method", y = "Metric Value", fill = "Type")
    print(p_obs_bar)
    ggsave(file.path(results_folder, paste0("Observed_vs_Random_", dataset_name, ".svg")),
           p_obs_bar, device = "svg", width = 9, height = 6)

    dev.off()
    log_message("Completed all plotting and saved PDF and SVG files")

    ###############################
    # Process linked clusters for this unified method
    ###############################
    log_message("Extracting linked clusters from iterations for unified method (using global mapping).")
    linked_data_iterations <- process_linked_clusters(data_iterations_method, linked_cluster_col = "LinkedCluster", verbose = verbose)
    if (length(linked_data_iterations) < 2) {
      log_message("Fewer than two iterations have linked cluster assignments. Skipping linked cluster comparison for unified method.")
    } else {
      log_message("Performing cross-iteration comparison on linked clusters for unified method.")
      cross_iter_linked_results <- cross_iteration_comparison_with_betti(
        data_iterations = linked_data_iterations,
        group_by_col = "LinkedCluster",       # Clusters are now aligned according to the unified method mapping
        orig_ident_columns = "orig.ident",      # Adjust if needed
        tau_vals = seq(0, 1, length.out = 500),
        dimensions = c(0, 1),
        grid_points = 500,
        num_permutations = 10000,
        bootstrap_samples = 100,
        base_sigma = 1,
        output_folder = file.path(results_folder, "linked_cluster_comparisons", unified_method),
        num_cores = n_cores,
        verbose = verbose
      )
      linked_output_file <- file.path(results_folder, paste0(dataset_name, "_", unified_method, "_linked_cluster_comparison.rds"))
      saveRDS(cross_iter_linked_results, linked_output_file)
      log_message(sprintf("Saved linked cluster comparison results for unified method '%s' to: %s", unified_method, linked_output_file))
      rm(cross_iter_linked_results); gc()
    }
    rm(data_iterations_method, mapping_sub); gc()
  }

  ###############################
  # Return results
  ###############################
  list(
    raw_comparison = all_results,
    norm_metrics   = norm_df
    # Linked cluster comparisons for each unified method are saved separately.
  )
}


cross_iteration_comparison_with_betti(
  data_iterations = data_iterations,
  group_by_col = "Tissue",
  tau_vals = seq(0, 1, length.out = 500),
  dimensions = c(0, 1),
  grid_points = 500,
  num_permutations = 10000,
  bootstrap_samples = 100,
  base_sigma = 1,
  output_folder = file.path(results_folder, "cross_iteration_comparisons"),
  num_cores = 16
  # orig_ident_columns = c("SRA Number", "SRS Number")
)

cross_iteration_comparison_with_betti(
  data_iterations = data_iterations,
  group_by_col = "SRA Number",
  tau_vals = seq(0, 1, length.out = 500),
  dimensions = c(0, 1),
  grid_points = 500,
  num_permutations = 10000,
  bootstrap_samples = 100,
  base_sigma = 1,
  output_folder = file.path(results_folder, "cross_iteration_comparisons"),
  num_cores = 16
  # orig_ident_columns = c("SRA Number", "SRS Number")

)

cross_iteration_comparison_with_betti(
  data_iterations = data_iterations,
  group_by_col = "Approach",
  tau_vals = seq(0, 1, length.out = 500),
  dimensions = c(0, 1),
  grid_points = 500,
  num_permutations = 10000,
  bootstrap_samples = 100,
  base_sigma = 1,
  output_folder = file.path(results_folder, "cross_iteration_comparisons"),
  num_cores = 16
  # orig_ident_columns = c("SRA Number", "SRS Number")
)


cross_iteration_comparison_with_betti(
  data_iterations = data_iterations,
  group_by_col = "Random_Group_bootstrap_1",
  tau_vals = seq(0, 1, length.out = 500),
  dimensions = c(0, 1),
  grid_points = 500,
  num_permutations = 10000,
  bootstrap_samples = 100,
  base_sigma = 1,
  output_folder = file.path(results_folder, "cross_iteration_comparisons"),
  num_cores = 16
  # orig_ident_columns = c("SRA Number", "SRS Number")
)

# ###########################
# cross_iteration_comparison_with_betti(
#   data_iterations = data_iterations_bonemarrow,
#   group_by_col = "Tissue",
#   tau_vals = seq(0, 1, length.out = 500),
#   dimensions = c(0, 1),
#   grid_points = 500,
#   num_permutations = 10000,
#   bootstrap_samples = 100,
#   base_sigma = 1,
#   output_folder = file.path(results_folder, "cross_iteration_comparisons"),
#   num_cores = 16,
#   orig_ident_columns = c("orig.ident"),
#   metadata = "./data/GSE120221/metadata.csv"
# )

# cross_iteration_comparison_with_betti(
#   data_iterations = data_iterations_bonemarrow,
#   group_by_col = "SRA",
#   tau_vals = seq(0, 1, length.out = 500),
#   dimensions = c(0, 1),
#   grid_points = 500,
#   num_permutations = 10000,
#   bootstrap_samples = 100,
#   base_sigma = 1,
#   output_folder = file.path(results_folder, "cross_iteration_comparisons"),
#   num_cores = 16,
#   orig_ident_columns = c("orig.ident"),
#   metadata = "./data/GSE120221/metadata.csv"
# )

# cross_iteration_comparison_with_betti(
#   data_iterations = data_iterations_bonemarrow,
#   group_by_col = "Approach",
#   tau_vals = seq(0, 1, length.out = 500),
#   dimensions = c(0, 1),
#   grid_points = 500,
#   num_permutations = 10000,
#   bootstrap_samples = 100,
#   base_sigma = 1,
#   output_folder = file.path(results_folder, "cross_iteration_comparisons"),
#   num_cores = 16,
#   orig_ident_columns = c("orig.ident"),
#   metadata = "./data/GSE120221/metadata.csv"
# )

# cross_iteration_comparison_with_betti(
#   data_iterations = data_iterations_bonemarrow,
#   group_by_col = "orig.ident",
#   tau_vals = seq(0, 1, length.out = 500),
#   dimensions = c(0, 1),
#   grid_points = 500,
#   num_permutations = 10000,
#   bootstrap_samples = 100,
#   base_sigma = 1,
#   output_folder = file.path(results_folder, "cross_iteration_comparisons"),
#   num_cores = 16,
#   orig_ident_columns = c("orig.ident"),
#   metadata = "./data/GSE120221/metadata.csv"
# )

# cross_iteration_comparison_with_betti(
#   data_iterations = data_iterations_bonemarrow,
#   group_by_col = "Random_Group_bootstrap_1",
#   tau_vals = seq(0, 1, length.out = 500),
#   dimensions = c(0, 1),
#   grid_points = 500,
#   num_permutations = 10000,
#   bootstrap_samples = 100,
#   base_sigma = 1,
#   output_folder = file.path(results_folder, "cross_iteration_comparisons"),
#   num_cores = 16,
#   orig_ident_columns = c("orig.ident"),
#   metadata = "./data/GSE120221/metadata.csv"
# )

# # Suppose you have two iterations:
# # 1) Raw_bonemarrow: cluster column named "seurat_cluster_Raw_bonemarrow"
# # 2) Integrated_bonemarrow: cluster column named "seurat_cluster_Integrated_bonemarrow"

# method_mapping_list_bonemarrow <- list(
#   # 1) Seurat clustering
#   "seurat_cluster" = list(
#     "Raw_bonemarrow"         = "seurat_cluster_raw_bonemarrow",
#     "SCT_Individual_bonemarrow" = "seurat_cluster_sct_individual_bonemarrow",
#     "SCT_Whole_bonemarrow"   = "seurat_cluster_sct_whole_bonemarrow",
#     "Integrated_bonemarrow"  = "seurat_cluster_integrated_bonemarrow"
#   ),

#   # 2) K-means (tissue)
#   # "kmeans_cluster_tissue" = list(
#   #   "Raw_bonemarrow"         = "kmeans_cluster_raw_bonemarrow_tissue",
#   #   "SCT_Individual_bonemarrow" = "kmeans_cluster_sct_individual_bonemarrow_tissue",
#   #   "SCT_Whole_bonemarrow"   = "kmeans_cluster_sct_whole_bonemarrow_tissue",
#   #   "Integrated_bonemarrow"  = "kmeans_cluster_integrated_bonemarrow_tissue"
#   # ),

#   # 3) K-means (sra)
#   "kmeans_cluster_sra" = list(
#     "Raw_bonemarrow"         = "kmeans_cluster_raw_bonemarrow_sra",
#     "SCT_Individual_bonemarrow" = "kmeans_cluster_sct_individual_bonemarrow_sra",
#     "SCT_Whole_bonemarrow"   = "kmeans_cluster_sct_whole_bonemarrow_sra",
#     "Integrated_bonemarrow"  = "kmeans_cluster_integrated_bonemarrow_sra"
#   ),

#   # 4) K-means (approach)
#   # "kmeans_cluster_approach" = list(
#   #   "Raw_bonemarrow"         = "kmeans_cluster_raw_bonemarrow_approach",
#   #   "SCT_Individual_bonemarrow" = "kmeans_cluster_sct_individual_bonemarrow_approach",
#   #   "SCT_Whole_bonemarrow"   = "kmeans_cluster_sct_whole_bonemarrow_approach",
#   #   "Integrated_bonemarrow"  = "kmeans_cluster_integrated_bonemarrow_approach"
#   # ),

#   # 5) Hierarchical clustering (BDM + PH) (tissue)
#   # "hierarchical_cluster_bdm_ph_tissue" = list(
#   #   "Raw_bonemarrow"         = "hierarchical_cluster_bdm_ph_raw_bonemarrow_tissue",
#   #   "SCT_Individual_bonemarrow" = "hierarchical_cluster_bdm_ph_sct_individual_bonemarrow_tissue",
#   #   "SCT_Whole_bonemarrow"   = "hierarchical_cluster_bdm_ph_sct_whole_bonemarrow_tissue",
#   #   "Integrated_bonemarrow"  = "hierarchical_cluster_bdm_ph_integrated_bonemarrow_tissue"
#   # ),

#   # 6) Hierarchical clustering (BDM + PH) (sra)
#   "hierarchical_cluster_bdm_ph_sra" = list(
#     "Raw_bonemarrow"         = "hierarchical_cluster_bdm_ph_raw_bonemarrow_sra",
#     "SCT_Individual_bonemarrow" = "hierarchical_cluster_bdm_ph_sct_individual_bonemarrow_sra",
#     "SCT_Whole_bonemarrow"   = "hierarchical_cluster_bdm_ph_sct_whole_bonemarrow_sra",
#     "Integrated_bonemarrow"  = "hierarchical_cluster_bdm_ph_integrated_bonemarrow_sra"
#   ),

#   # 7) Hierarchical clustering (BDM + PH) (approach)
#   # "hierarchical_cluster_bdm_ph_approach" = list(
#   #   "Raw_bonemarrow"         = "hierarchical_cluster_bdm_ph_raw_bonemarrow_approach",
#   #   "SCT_Individual_bonemarrow" = "hierarchical_cluster_bdm_ph_sct_individual_bonemarrow_approach",
#   #   "SCT_Whole_bonemarrow"   = "hierarchical_cluster_bdm_ph_sct_whole_bonemarrow_approach",
#   #   "Integrated_bonemarrow"  = "hierarchical_cluster_bdm_ph_integrated_bonemarrow_approach"
#   # ),

#   # 8) Hierarchical clustering (SDM + PH) (tissue)
#   # "hierarchical_cluster_sdm_ph_tissue" = list(
#   #   "Raw_bonemarrow"         = "hierarchical_cluster_sdm_ph_raw_bonemarrow_tissue",
#   #   "SCT_Individual_bonemarrow" = "hierarchical_cluster_sdm_ph_sct_individual_bonemarrow_tissue",
#   #   "SCT_Whole_bonemarrow"   = "hierarchical_cluster_sdm_ph_sct_whole_bonemarrow_tissue",
#   #   "Integrated_bonemarrow"  = "hierarchical_cluster_sdm_ph_integrated_bonemarrow_tissue"
#   # ),

#   # 9) Hierarchical clustering (SDM + PH) (sra)
#   "hierarchical_cluster_sdm_ph_sra" = list(
#     "Raw_bonemarrow"         = "hierarchical_cluster_sdm_ph_raw_bonemarrow_sra",
#     "SCT_Individual_bonemarrow" = "hierarchical_cluster_sdm_ph_sct_individual_bonemarrow_sra",
#     "SCT_Whole_bonemarrow"   = "hierarchical_cluster_sdm_ph_sct_whole_bonemarrow_sra",
#     "Integrated_bonemarrow"  = "hierarchical_cluster_sdm_ph_integrated_bonemarrow_sra"
#   ),

#   # 10) Hierarchical clustering (SDM + PH) (approach)
#   # "hierarchical_cluster_sdm_ph_approach" = list(
#   #   "Raw_bonemarrow"         = "hierarchical_cluster_sdm_ph_raw_bonemarrow_approach",
#   #   "SCT_Individual_bonemarrow" = "hierarchical_cluster_sdm_ph_sct_individual_bonemarrow_approach",
#   #   "SCT_Whole_bonemarrow"   = "hierarchical_cluster_sdm_ph_sct_whole_bonemarrow_approach",
#   #   "Integrated_bonemarrow"  = "hierarchical_cluster_sdm_ph_integrated_bonemarrow_approach"
#   # ),

#   # 11) Spectral clustering (BDM) (tissue)
#   # "spectral_cluster_bdm_tissue" = list(
#   #   "Raw_bonemarrow"         = "spectral_cluster_bdm_raw_bonemarrow_tissue",
#   #   "SCT_Individual_bonemarrow" = "spectral_cluster_bdm_sct_individual_bonemarrow_tissue",
#   #   "SCT_Whole_bonemarrow"   = "spectral_cluster_bdm_sct_whole_bonemarrow_tissue",
#   #   "Integrated_bonemarrow"  = "spectral_cluster_bdm_integrated_bonemarrow_tissue"
#   # ),

#   # 12) Spectral clustering (BDM) (sra)
#   "spectral_cluster_bdm_sra" = list(
#     "Raw_bonemarrow"         = "spectral_cluster_bdm_raw_bonemarrow_sra",
#     "SCT_Individual_bonemarrow" = "spectral_cluster_bdm_sct_individual_bonemarrow_sra",
#     "SCT_Whole_bonemarrow"   = "spectral_cluster_bdm_sct_whole_bonemarrow_sra",
#     "Integrated_bonemarrow"  = "spectral_cluster_bdm_integrated_bonemarrow_sra"
#   ),

#     # 13) Spectral clustering (BDM) (approach)
#   #   "spectral_cluster_bdm_approach" = list(
#   #     "Raw_bonemarrow"         = "spectral_cluster_bdm_raw_bonemarrow_approach",
#   #     "SCT_Individual_bonemarrow" = "spectral_cluster_bdm_sct_individual_bonemarrow_approach",
#   #     "SCT_Whole_bonemarrow"   = "spectral_cluster_bdm_sct_whole_bonemarrow_approach",
#   #     "Integrated_bonemarrow"  = "spectral_cluster_bdm_integrated_bonemarrow_approach"
#   #  ),
#   # 14) Spectral clustering (SDM) (tissue)
#   # "spectral_cluster_sdm_tissue" = list(
#   #   "Raw_bonemarrow"         = "spectral_cluster_sdm_raw_bonemarrow_tissue",
#   #   "SCT_Individual_bonemarrow" = "spectral_cluster_sdm_sct_individual_bonemarrow_tissue",
#   #   "SCT_Whole_bonemarrow"   = "spectral_cluster_sdm_sct_whole_bonemarrow_tissue",
#   #   "Integrated_bonemarrow"  = "spectral_cluster_sdm_integrated_bonemarrow_tissue"
#   # ),

#   # 15) Spectral clustering (SDM) (sra)
#   "spectral_cluster_sdm_sra" = list(
#     "Raw_bonemarrow"         = "spectral_cluster_sdm_raw_bonemarrow_sra",
#     "SCT_Individual_bonemarrow" = "spectral_cluster_sdm_sct_individual_bonemarrow_sra",
#     "SCT_Whole_bonemarrow"   = "spectral_cluster_sdm_sct_whole_bonemarrow_sra",
#     "Integrated_bonemarrow"  = "spectral_cluster_sdm_integrated_bonemarrow_sra"
#   )
#   # # 16) Spectral clustering (SDM) (approach)
#   # "spectral_cluster_sdm_approach" = list(
#   #   "Raw_bonemarrow"         = "spectral_cluster_sdm_raw_bonemarrow_approach",
#   #   "SCT_Individual_bonemarrow" = "spectral_cluster_sdm_sct_individual_bonemarrow_approach",
#   #   "SCT_Whole_bonemarrow"   = "spectral_cluster_sdm_sct_whole_bonemarrow_approach",
#   #   "Integrated_bonemarrow"  = "spectral_cluster_sdm_integrated_bonemarrow_approach"
#   # ),
# )
# # all_cluster_comparisons <- compare_clusters_across_iterations_by_centroids_mapped(
# #   data_iterations    = data_iterations_bonemarrow,         # e.g. raw, integrated, etc.
# #   method_mapping_list = method_mapping_list_bonemarrow,
# #   dims_use          = 1:10,
# #   reduction_name    = "pca",
# #   results_folder    = paste0(results_folder, "/compare_clusters_across_iterations"),
# #   num_cores         = 4
# # )

# results <- compare_all_iterations_with_mapped_methods_random_baseline_including_purity_pvals (
#   data_iterations = data_iterations_bonemarrow,
#   method_mapping_list = method_mapping_list_bonemarrow,
#   dims_use         = 1:10,
#   reduction_name   = "pca",
#   n_random         = 100,
#   results_folder    = paste0(results_folder, "/compare_clusters_across_iterations"),
#   seed             = 42
# )
