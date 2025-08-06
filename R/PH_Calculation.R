# Internal helper to extract columns with constant values from metadata
.get_constant_metadata <- function(obj, count_field) {
  metadata <- obj@meta.data
  constant_cols <- sapply(metadata, function(col) length(unique(col)) == 1)
  constant_metadata <- metadata[1, constant_cols, drop = FALSE]  # Keep one row with constant columns
  constant_metadata$Sample <- obj@project.name  # Add the sample name
  constant_metadata[[count_field]] <- nrow(metadata)  # Add cell count with custom label
  return(constant_metadata)
}

#' PH Pipeline for Single-Cell RNA-Seq Data
#'
#' @description 
#' The `process_datasets_PH` function is designed to process single-cell RNA-seq datasets using 
#' Persistent Homology (PH), data integration techniques, and clustering methods. 
#' The pipeline handles both integrated and unintegrated data, calculates persistence diagrams, 
#' and performs statistical comparisons across different clustering methods such as UMAP 
#' and hierarchical clustering.
#' 
#' Batch effect correction is currently supported through Seurat integration for
#' multi-sample or multi-dataset workflows.
#' The final output includes persistence diagrams. Distance matrices
#' (BDM/SDM/LDM) are produced later during post-processing along with UMAP
#' embeddings and clustering visualizations.
#'
#' @section Workflow Overview:
#' \enumerate{
#'   \item \strong{Data Loading & Preprocessing:}
#'     \itemize{
#'       \item Loads single-cell RNA-seq data from `.RData` files based on paths in the input metadata.
#'       \item Creates `Seurat` objects for each sample and filters cells and features based on quality metrics.
#'       \item Calculates mitochondrial, ribosomal, and hemoglobin gene percentages for quality control.
#'       \item Normalizes data and identifies variable features for downstream analysis.
#'     }
#'   
#'   \item \strong{Persistent Homology (PH) Analysis:}
#'     \itemize{
#'       \item Calculates persistence diagrams for unintegrated data using `ripserr` and `TDA` packages.
#'       \item Integrates datasets using Seurat and recalculates persistence diagrams for integrated data.
#'     }
#'   
#'   \item \strong{Distance Matrix Generation (Post-Processing):}
#'     \itemize{
#'       \item Persistence diagrams created in this step are used by
#'         `process_iteration_calculate_matrices()` in
#'         `PH_PostProcessing_andAnalysis.R` to compute BDM, SDM and LDM once.
#'       \item Downsampling of persistence diagrams may be performed to optimise
#'         memory usage during these calculations.
#'     }
#'   
#'   \item \strong{Clustering & Visualization:}
#'     \itemize{
#'       \item Performs UMAP projections for visualizing data in two-dimensional space.
#'       \item Conducts hierarchical clustering on UMAP embeddings and compares results with PH-based clusters.
#'       \item Generates heatmaps, dendrograms, and annotations based on tissue types or experimental conditions.
#'     }
#'   
#'   \item \strong{Statistical Analysis:}
#'     \itemize{
#'       \item Performs statistical tests to compare clustering results (e.g., Adjusted Rand Index, Mutual Information, Chi-Square).
#'       \item Calculates cluster purity, Silhouette analysis, and Davies-Bouldin Index for validating clusters.
#'     }
#' }
#'
#' @section Dependencies:
#' The script requires the following R packages:
#' \code{tidyverse}, \code{Matrix}, \code{ripserr}, \code{TDA}, \code{foreach},
#' \code{doParallel}, \code{parallel}, \code{rliger}, \code{Seurat}, \code{stringr}, 
#' \code{ComplexHeatmap}, \code{mclust}, \code{aricode}, \code{clusterSim}, \code{Rtsne}.
#'
#' @section Usage:
#' Ensure all required packages are installed, and the input metadata CSV is correctly formatted with paths to `.RData` files. 
#' Execute the script in an environment with adequate computational resources, as the pipeline utilizes parallel processing.
#'
#' @section Example:
#' \dontrun{
#'   # Load metadata and process datasets with Persistent Homology
#'   metadata <- read.csv("./data/VastlyDifferentTissues/metadata.csv",
#'                        check.names = FALSE)
#'   process_datasets_PH(metadata, num_cores = 32)
#' }
#'
#' @section Author:
#' Jonah Daneshmand
#'
#' @section License:
#' MIT License

#' @param metadata Data frame with file paths and associated metadata.
#' @param integration_method Integration method to use ("seurat",
#'   "liger" or "mnn").
#' @param num_cores Number of cores for parallel processing.
#' @param MIN_CELLS Minimum number of cells to keep a dataset.
#' @param DIM Dimensionality for persistent homology calculation.
#' @param THRESHOLD Threshold for persistence diagram calculation.
#' @param dataset_tag Optional tag appended to output files.
#'
#' @return A list containing processed iterations as well as the detected
#'   column names for SRA, tissue and approach.
#'
#' @examples
#' \dontrun{
#' metadata <- read.csv("./data/VastlyDifferentTissues/metadata.csv",
#'                     check.names = FALSE)
#' results <- process_datasets_PH(metadata, num_cores = 32)
#' }
#' @export



process_datasets_PH <- function(metadata,
                                integration_method = "seurat",
                                num_cores = 16,
                                MIN_CELLS = 250,
                                DIM = 1,
                                THRESHOLD = -1,
                                dataset_tag = "dataset") {
  # Determine metadata column names and warn if missing
  sra_col <- intersect(c("SRA", "SRA_Number", "SRA Number"), colnames(metadata))[1]
  tissue_col <- intersect(c("Tissue", "tissue"), colnames(metadata))[1]
  approach_col <- intersect(c("Approach", "approach"), colnames(metadata))[1]
  if (all(is.na(c(sra_col, tissue_col, approach_col)))) {
    stop("Metadata must contain at least one of 'SRA', 'Tissue', or 'Approach'.")
  }
  missing <- c()
  if (is.na(sra_col)) missing <- c(missing, "SRA")
  if (is.na(tissue_col)) missing <- c(missing, "Tissue")
  if (is.na(approach_col)) missing <- c(missing, "Approach")
  if (length(missing) > 0) {
    message("Metadata missing columns: ", paste(missing, collapse = ", "))
  }

  dataset_suffix <- if (nzchar(dataset_tag)) paste0("_", dataset_tag) else ""
  
  # Start logging
  log_file <- paste0("PH_Pipeline_Log_", Sys.Date(), dataset_suffix, ".txt")
  sink(log_file, append = TRUE)
  on.exit(sink(), add = TRUE)
  
  # Helper functions are available once the package is loaded
  
  log_message("Processing series with CSV input")

  # Extract file paths from either 'File Path' or 'File.Path' column
  file_paths <- dplyr::coalesce(metadata[['File Path']], metadata[['File.Path']])
  if (all(is.na(file_paths))) {
    stop("Metadata must contain a 'File Path' or 'File.Path' column.")
  }

  loaded <- tryCatch(
    load_sparse_matrices(file_paths),
    error = function(e) {
      log_message(paste("Error in processing .RData files from CSV:", e$message))
      NULL
    }
  )

  if (is.null(loaded) || length(loaded$sample_names) == 0) {
    return(NULL)
  }

  my_sparse_matrices <- loaded$matrices
  sample_names <- loaded$sample_names
  
  # Create Seurat objects for each sample
  my_seurat_list <- tryCatch(
    {
      mclapply(seq_along(my_sparse_matrices), function(i) {
        spm <- my_sparse_matrices[[i]]
        name <- names(my_sparse_matrices)[i]
        CreateSeuratObject(counts = spm, project = name)
      }, mc.cores = num_cores)
    },
    error = function(e) {
      log_message(paste("Error in creating Seurat objects:", e$message))
      NULL
    }
  )
  
  if (is.null(my_seurat_list)) {
    return(NULL)
  }

  
  # Calculate percentage of mitochondrial, ribosomal, and hemoglobin genes
  my_seurat_list <- tryCatch(
    {
      mclapply(my_seurat_list, function(obj) {
        obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent_mito")
        obj <- PercentageFeatureSet(obj, pattern = "^RP[SL]", col.name = "percent_ribo")
        # Hemoglobin genes are typically annotated as 'HB' followed by a subunit
        # letter (e.g., HBA1, HBB). Exclude non-hemoglobin genes such as 'HBP'
        # by using a negative lookahead.
        obj <- PercentageFeatureSet(obj, pattern = "^HB(?!P)", perl = TRUE,
                                    col.name = "percent_hb")
        return(obj)
      }, mc.cores = num_cores)
    },
    error = function(e) {
      log_message(paste("Error in calculating percentage metrics:", e$message))
      NULL
    }
  )
  
  if (is.null(my_seurat_list)) {
    return(NULL)
  }
  
  
  # Add SRA, tissue, and approach metadata without automatic approach detection
  my_seurat_list <- tryCatch(
    {
      mclapply(seq_along(my_seurat_list), function(i) {
        seurat_obj <- my_seurat_list[[i]]

        if (!is.na(sra_col)) {
          seurat_obj@meta.data$SRA <- metadata[[sra_col]][i]
        }
        if (!is.na(tissue_col)) {
          seurat_obj@meta.data$Tissue <- metadata[[tissue_col]][i]
        }
        if (!is.na(approach_col)) {
          seurat_obj@meta.data$Approach <- metadata[[approach_col]][i]
        }

        return(seurat_obj)
      }, mc.cores = num_cores)
    },
    error = function(e) {
      log_message(paste("Error in adding SRA, tissue, and approach metadata:", e$message))
      NULL
    }
  )
  
  if (is.null(my_seurat_list)) {
    return(NULL)
  }
  

  
  prefiltered_cells_file <- paste0("prefiltered_cells", dataset_suffix, ".csv")

  # Extract constant metadata for each sample before filtering
  prefiltered_cells <- tryCatch(
    {
      # Combine metadata from all Seurat objects in the list
      do.call(
        rbind,
        lapply(
          my_seurat_list,
          .get_constant_metadata,
          count_field = "Number_of_Cells_Before_Filtering"
        )
      )
    },
    error = function(e) {
      log_message(paste("Error in creating prefiltered metadata table:", e$message))
      NULL
    }
  )

  # Save the prefiltering metadata to a CSV file
  if (!is.null(prefiltered_cells)) {
    tryCatch(
      {
        write.csv(prefiltered_cells, file = prefiltered_cells_file, row.names = FALSE)
        log_message(paste("Prefiltered metadata saved to", prefiltered_cells_file))
      },
      error = function(e) {
        log_message(paste("Error in saving prefiltered metadata to CSV:", e$message))
      }
    )
  }

  
  # Filter Seurat objects based on specific criteria
  my_seurat_list_filtered <- tryCatch(
    {
      mclapply(my_seurat_list, function(obj) {
        minGenesPerCell <- 500
        maxGenesPerCell <- 9000
        maxMitoPercent <- 20
        minRiboPercent <- 5
        minCellsPerGene <- 3
        
        # Filtering based on cell metrics
        cells_to_keep <- WhichCells(obj, expression = nFeature_RNA >= minGenesPerCell &
                                      nFeature_RNA <= maxGenesPerCell &
                                      percent_mito <= maxMitoPercent &
                                      percent_ribo > minRiboPercent)
        obj_filtered <- subset(obj, cells = cells_to_keep)
        
        # Filtering based on gene metrics
        counts_matrix <- GetAssayData(obj_filtered, slot = "counts")
        genes_to_keep <- rowSums(counts_matrix > 0) > minCellsPerGene
        obj_filtered <- subset(obj_filtered, features = names(genes_to_keep[genes_to_keep]))
        
        return(obj_filtered)
      }, mc.cores = num_cores)
    },
    error = function(e) {
      log_message(paste("Error in filtering Seurat objects for series:", e$message))
      NULL
    }
  )
  
  if (is.null(my_seurat_list_filtered)) {
    return(NULL)
  }
  
  # Remove datasets with insufficient cells after filtering
  my_seurat_list_filtered <- tryCatch(
    {
      discard(my_seurat_list_filtered, function(obj) {
        num_cells <- nrow(obj@meta.data)
        log_message(paste("Sample", obj@project.name, "has", num_cells, "cells"))
        if (num_cells < MIN_CELLS) {
          log_message(paste("Skipping sample", obj@project.name, "due to insufficient cells after filtering (", num_cells, "cells)"))
          return(TRUE)  # Discard this object
        }
        return(FALSE)  # Keep this object
      })
    },
    error = function(e) {
      log_message(paste("Error in removing insufficient cell datasets for series:", e$message))
      NULL
    }
  )


  filtered_cells_file <- paste0("filtered_cells", dataset_suffix, ".csv")

  # Extract constant metadata for each sample after filtering
  filtered_cells <- tryCatch(
    {
      # Combine metadata from all filtered Seurat objects in the list
      do.call(
        rbind,
        lapply(
          my_seurat_list_filtered,
          .get_constant_metadata,
          count_field = "Number_of_Cells_After_Filtering"
        )
      )
    },
    error = function(e) {
      log_message(paste("Error in creating filtered metadata table:", e$message))
      NULL
    }
  )

  # Save the filtered metadata to a CSV file
  if (!is.null(filtered_cells)) {
    tryCatch(
      {
        write.csv(filtered_cells, file = filtered_cells_file, row.names = FALSE)
        log_message(paste("Filtered metadata saved to", filtered_cells_file))
      },
      error = function(e) {
        log_message(paste("Error in saving filtered metadata to CSV:", e$message))
      }
    )
  }

  # Return NULL if the Seurat list is empty
  if (is.null(my_seurat_list_filtered)) {
    return(NULL)
  }

  
  # Normalize data according to the integration method
  if (integration_method == "seurat") {
    # Perform SCTransform normalization on each individually
    log_message("Performing SCTransform normalization for Seurat integration")
    my_seurat_list_filtered <- tryCatch(
      {
        mclapply(my_seurat_list_filtered, function(x) {
          # Ensure RNA assay has a data slot before SCTransform
          x <- NormalizeData(x, verbose = TRUE) # Creates a data slot if missing
          
          # Store the RNA assay's data slot
          rna_data <- GetAssayData(x, assay = "RNA", slot = "data")

          # Perform SCTransform
          x <- SCTransform(x, variable.features.n = 6000, verbose = TRUE)

          # Restore RNA assay's data slot
          x@assays$RNA$data <- rna_data
          DefaultAssay(x) <- "SCT"

          return(x)
        }, mc.cores = num_cores)
      },
      error = function(e) {
        log_message(paste("Error in SCTransform normalization for series:", e$message))
        NULL
      }
    )
  
    # Check if my_seurat_list_filtered exists and merged_seurat_unintegrated doesn't exist
    if (exists("my_seurat_list_filtered") && !exists("merged_seurat_unintegrated")) {
      
    ### 1. Rename cells in each original object so that cell names are unique
    names <- vapply(
      my_seurat_list_filtered,
      function(seurat_obj) seurat_obj@meta.data$orig.ident[1],
      FUN.VALUE = character(1)
    )
    my_seurat_list_filtered <- lapply(my_seurat_list_filtered, function(seurat_obj) {
      # Remove any existing prefixes
      RenameCells(seurat_obj, add.cell.id = NULL)
    })

    ### 2. Merge the objects (merge.data = TRUE keeps all per-sample data)
    merged_seurat_unintegrated <- merge(
      x = my_seurat_list_filtered[[1]],
      y = my_seurat_list_filtered[-1],
      add.cell.ids = as.character(names),  # This will prefix cell names (e.g. "Sample1_cellbarcode")
      project = "MergedUnintegrated",
      merge.data = TRUE                    # Keeps all raw data (though in a multi-layer Assay5 object)
    )

    ### 3. Create a combined “normal” RNA assay with one counts matrix and one data slot
    # Instead of trying to extract from the merged object's multi-layer RNA assay,
    # we’ll get the counts from each original object.

    # (a) For each original object, get its RNA counts and update cell names with the prefix.
    rna_counts_list <- lapply(my_seurat_list_filtered, function(x) {
      counts <- GetAssayData(x, assay = "RNA", slot = "counts")
      prefix <- x@meta.data$orig.ident[1]
      colnames(counts) <- paste0(prefix, "_", colnames(counts))
      return(counts)
    })

    # (b) Get the union of all gene names across objects
    all_genes <- unique(unlist(lapply(rna_counts_list, rownames)))

    # (c) For each counts matrix, add missing genes (filled with 0) so all matrices share the same gene set,
    #     and reorder rows to match the full gene set.
    rna_counts_list_aligned <- lapply(rna_counts_list, function(mat) {
      missing <- setdiff(all_genes, rownames(mat))
      if (length(missing) > 0) {
        add_mat <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                          dimnames = list(missing, colnames(mat)))
        mat <- rbind(mat, add_mat)
      }
      mat <- mat[all_genes, , drop = FALSE]
      return(mat)
    })

    # (d) Combine the counts matrices column-wise (each column = one cell)
    combined_counts <- do.call(cbind, rna_counts_list_aligned)

    # (e) Create a new standard RNA assay from the combined counts matrix and replace the merged object's RNA assay
    new_RNA_assay <- CreateAssayObject(counts = combined_counts)
        
    merged_seurat_unintegrated[["RNA"]] <- new_RNA_assay

    # (f) Normalize the new RNA assay to populate the "data" slot
    merged_seurat_unintegrated <- NormalizeData(merged_seurat_unintegrated, assay = "RNA")

    ### 4. Run SCTransform on the merged object (using the new RNA assay) to create an SCT assay
    merged_seurat_unintegrated <- SCTransform(
      merged_seurat_unintegrated,
      assay = "RNA",
      new.assay.name = "SCT",
      verbose = TRUE
    )

    ### 5. Create the SCT_Ind assay from the individually SCTransformed objects
    # (a) For each original object, extract its SCT "data" (normalized values) and update cell names.
    sct_data_list <- lapply(my_seurat_list_filtered, function(x) {
      prefix <- x@meta.data$orig.ident[1]
      mat <- GetAssayData(x, assay = "SCT", slot = "data")
      colnames(mat) <- paste0(prefix, "_", colnames(mat))
      return(mat)
    })

    # (b) Align features (genes) across these matrices
    all_features_sct <- unique(unlist(lapply(sct_data_list, rownames)))
    sct_data_list_aligned <- lapply(sct_data_list, function(mat) {
      missing <- setdiff(all_features_sct, rownames(mat))
      if (length(missing) > 0) {
        add_mat <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                          dimnames = list(missing, colnames(mat)))
        mat <- rbind(mat, add_mat)
      }
      mat <- mat[all_features_sct, , drop = FALSE]
      return(mat)
    })

    # (c) Combine the individually SCT-transformed matrices column-wise
    combined_sct_data <- do.call(cbind, sct_data_list_aligned)
    # Ensure the cell (column) order matches the merged object
    combined_sct_data <- combined_sct_data[, colnames(merged_seurat_unintegrated)]

    # (d) Create a new assay "SCT_Ind" from the combined SCT data
    merged_seurat_unintegrated[["SCT_Ind"]] <- CreateAssayObject(counts = combined_sct_data)

    ### (Optional) Set the default assay as needed
    DefaultAssay(merged_seurat_unintegrated) <- "RNA"

    ### Final Merged Object Contains:
    # - "RNA": A standard RNA assay with a combined counts matrix and normalized data.
    # - "SCT": An assay created from running SCTransform on the combined RNA assay.
    # - "SCT_Ind": An assay containing combined individually SCTransformed data.

    saveRDS(merged_seurat_unintegrated, file="merged_seurat_unintegrated.Rds")
    # merged_seurat_unintegrated <- readRDS("merged_seurat_unintegrated.Rds")

    } else {
      # If the conditions are not met, print appropriate messages
      if (!exists("my_seurat_list_filtered")) {
        message("The list 'my_seurat_list_filtered' does not exist. Please ensure it is created before running this code.")
      }

      if (exists("merged_seurat_unintegrated")) {
        message("The object 'merged_seurat_unintegrated' already exists. Skipping merging and SCTransform steps.")
      }
    }
  }

  #### RUN PH on UNINTEGRATED DATA ####
  if (integration_method == "seurat") {
    # Extract SCTransformed data for unintegrated PH
    log_message("Extracting SCTransformed data for unintegrated PH")
    
    expr_list_sctInd <- tryCatch(
      {
        # Extract expression matrices from each Seurat object using SCT slot
        expr_matrices <- lapply(my_seurat_list_filtered, function(obj) {
          as.matrix(GetAssayData(obj, assay = "SCT", slot = "data"))
        })
        
        # Extract the first level of the orig.ident factor from each Seurat object
        orig_idents <- sapply(my_seurat_list_filtered, function(obj) {
          levels(obj@meta.data$orig.ident)[1]
        })
        
        # Assign the extracted orig.ident values as names to the expression matrices
        names(expr_matrices) <- orig_idents
        
        # Return the named list of expression matrices
        expr_matrices
      },
      error = function(e) {
        log_message(paste("Error in extracting SCTransformed data:", e$message))
        return(NULL)
      }
    )
    
    # Verify the naming of expr_list_sctransformed
    if (!is.null(expr_list_sctInd)) {
      log_message(paste("Successfully named expr_list_sctransformed entries by orig.ident:",
                        paste(names(expr_list_sctInd), collapse = ", ")))
      saveRDS(expr_list_sctInd, file = paste0("expr_list_sctInd", dataset_suffix, ".Rds"))
    } else {
      log_message("expr_list_sctransformed is NULL due to an error in extraction.")
    }
    
    # Extract Raw Counts Data for Unintegrated PH
    log_message("Extracting raw count data for unintegrated PH")
    
    expr_list_raw <- tryCatch(
      {
        # Extract expression matrices from each Seurat object using counts slot
        expr_matrices <- lapply(my_seurat_list_filtered, function(obj) {
          as.matrix(GetAssayData(obj, assay = "RNA", slot = "counts"))
        })
        
        # Extract the first level of the orig.ident factor from each Seurat object
        orig_idents <- sapply(my_seurat_list_filtered, function(obj) {
          levels(obj@meta.data$orig.ident)[1]
        })
        
        # Assign the extracted orig.ident values as names to the expression matrices
        names(expr_matrices) <- orig_idents
        
        # Return the named list of expression matrices
        expr_matrices
      },
      error = function(e) {
        log_message(paste("Error in extracting raw count data:", e$message))
        return(NULL)
      }
    )
    
    # Verify the naming of expr_list_raw
    if (!is.null(expr_list_raw)) {
      log_message(paste("Successfully named expr_list_raw entries by orig.ident:",
                        paste(names(expr_list_raw), collapse = ", ")))
      saveRDS(expr_list_raw, file = paste0("expr_list_raw", dataset_suffix, ".Rds"))
    } else {
      log_message("expr_list_raw is NULL due to an error in extraction.")
    }
    
    # Ensure that SCT is the default assay after SCTransform
    DefaultAssay(merged_seurat_unintegrated) <- "SCT"
    
    # Verify that the 'orig.ident' column exists in the metadata of the merged object
    if (!"orig.ident" %in% colnames(merged_seurat_unintegrated@meta.data)) {
      stop("Error: The 'orig.ident' column is not found in the metadata. Ensure that the orig.ident information is present.")
    }
    
    # Split the merged Seurat object by orig.ident using the metadata column
    origIdent_list <- SplitObject(merged_seurat_unintegrated, split.by = "orig.ident")
    
    # Create an empty list to store the SCT-normalized expression data for each orig.ident
    expr_list_sctWhole <- list()
    
    # Loop through each orig.ident-specific Seurat object and extract the SCT-normalized data
    for (sra_name in names(origIdent_list)) {
      # Extract the SCT-normalized data from the 'data' slot of the SCT assay
      sra_data <- GetAssayData(origIdent_list[[sra_name]], slot = "data", assay = "SCT")
      
      # Store the extracted data in expr_list using the orig.ident name as the key
      expr_list_sctWhole[[sra_name]] <- sra_data
    }
    
    # Check the structure of expr_list to verify the extraction
    str(expr_list_sctWhole)
    saveRDS(expr_list_sctWhole, file = paste0("expr_list_sctWhole", dataset_suffix, ".Rds"))

    timeout_datasets <- NULL

    PD_result_unintegrated <- compute_ph_batch(
      expr_list_sctInd, DIM, log_message, dataset_suffix,
      "_unintegrated", max_cores = 6
    )
    save_ph_results(
      PD_result_unintegrated, expr_list_sctInd, DIM, THRESHOLD,
      dataset_suffix, "_unintegrated", log_message
    )

    PD_result_unintegrated_RAW <- compute_ph_batch(
      expr_list_raw, DIM, log_message, dataset_suffix,
      "_unintegrated_RAW", max_cores = 6
    )
    save_ph_results(
      PD_result_unintegrated_RAW, expr_list_raw, DIM, THRESHOLD,
      dataset_suffix, "_unintegrated_RAW", log_message
    )

    PD_result_unintegrated_sctWhole <- compute_ph_batch(
      expr_list_sctWhole, DIM, log_message, dataset_suffix,
      "_unintegrated_sctWhole", max_cores = 8
    )
    save_ph_results(
      PD_result_unintegrated_sctWhole, expr_list_sctWhole, DIM, THRESHOLD,
      dataset_suffix, "_unintegrated_sctWhole", log_message
    )

    # Integration routines are available in the package

    # Call the integration function
    result <- perform_integration(
      seurat_list = my_seurat_list_filtered,
      integration_features_n = 3000,
      dims = 50,
      min_cells_threshold = 200,
      num_cores = num_cores
      #checkpoint_dir = "./integration_checkpoints_bonemarrow",
      #output_path = "final_integrated_seurat_bonemarrow.rds",
      #anchor_batches_dir = "./intermediate_seuratIntegration_results/anchor_batches_bonemarrow"
    )

    # result <- readRDS("final_integrated_seurat_bonemarrow.rds")
    # result <- readRDS("final_integrated_seurat.rds")

    DefaultAssay(result) <- "integrated"

    # # Save the final integrated Seurat object
    # final_save_path <- "final_integrated_seurat_v5_multiStudySeriesOnly.rds"
    # saveRDS(merged_integrated, file = final_save_path)
    # log_message(paste("Final integrated Seurat object saved successfully at:", final_save_path))

    # Verify that the 'orig.ident' column exists in the metadata of the merged object
    if (!"orig.ident" %in% colnames(result@meta.data)) {
      stop("Error: The 'orig.ident' column is not found in the metadata. Ensure that the orig.ident information is present.")
    }

    # Split the merged Seurat object by orig.ident using the metadata column
    origIdent_list <- SplitObject(result, split.by = "orig.ident")

    # Create an empty list to store the SCT-normalized expression data for each orig.ident
    expr_list_integrated <- list()

    # Loop through each orig.ident-specific Seurat object and extract the SCT-normalized data
    for (sra_name in names(origIdent_list)) {
      # Extract the SCT-normalized data from the 'data' slot of the SCT assay
      sra_data <- GetAssayData(origIdent_list[[sra_name]], slot = "data", assay = "integrated")

      # Store the extracted data in expr_list using the orig.ident name as the key
      expr_list_integrated[[sra_name]] <- sra_data
    }

    # Check the structure of expr_list to verify the extraction
    str(expr_list_integrated)

    saveRDS(expr_list_integrated, file = paste0("expr_list_integrated", dataset_suffix, ".Rds"))

    PD_result_integrated <- compute_ph_batch(
      expr_list_integrated, DIM, log_message, dataset_suffix,
      "_integrated", max_cores = 8
    )
    save_ph_results(
      PD_result_integrated, expr_list_integrated, DIM, THRESHOLD,
      dataset_suffix, "_integrated", log_message
    )
  }

  ##
  ## Assemble output structure for post processing
  ## -------------------------------------------------------
  data_iterations <- assemble_ph_results(
    merged_seurat_unintegrated, result,
    expr_list_raw, expr_list_sctInd,
    expr_list_sctWhole, expr_list_integrated
  )

  ph_results <- list(
    data_iterations = data_iterations,
    SRA_col = sra_col,
    Tissue_col = tissue_col,
    Approach_col = approach_col
  )
  return(ph_results)
}

