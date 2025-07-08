#' PH Pipeline for Single-Cell RNA-Seq Data
#'
#' @description 
#' The `process_datasets_PH` function is designed to process single-cell RNA-seq datasets using 
#' Persistent Homology (PH), data integration techniques, and clustering methods. 
#' The pipeline handles both integrated and unintegrated data, calculates persistence diagrams, 
#' and performs statistical comparisons across different clustering methods such as UMAP 
#' and hierarchical clustering.
#' 
#' It supports batch integration methods including LIGER, Seurat, and MNN (Mutual Nearest Neighbors),
#' allowing for batch effect correction in multi-sample or multi-dataset workflows. 
#' The final output includes persistence diagrams, bottleneck distance matrices (BDMs), 
#' UMAP embeddings, and clustering visualizations with statistical tests.
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
#'       \item Converts `Seurat` objects into `LIGER` objects for batch correction.
#'       \item Calculates persistence diagrams for unintegrated data using `ripserr` and `TDA` packages.
#'       \item Integrates datasets using methods such as Seurat, LIGER, or MNN, and recalculates persistence diagrams for integrated data.
#'     }
#'   
#'   \item \strong{Bottleneck Distance Matrix (BDM) Calculation:}
#'     \itemize{
#'       \item Computes bottleneck distances between persistence diagrams for both unintegrated and integrated data.
#'       \item Downsamples persistence diagrams to optimize memory usage and speed up calculations.
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
#'   metadata <- read.csv("./data/VastlyDifferentTissues/metadata.csv")
#'   process_datasets_PH(metadata, num_cores = 32)
#' }
#'
#' @section Author:
#' Jonah Daneshmand
#'
#' @section License:
#' [Specify the license under which the script is distributed, e.g., GPL-3]

source("PH_Functions.R")
metadata <- read_csv("./data/VastlyDifferentTissues/metadata.csv")
metadata <- read_csv("./data/GSE120221/metadata.csv")

filtered_cells <- read_csv("filtered_cells.csv")

# datasets_to_drop <- c(12,15,18) #going to come back to these in a future run

process_datasets_PH <- function(metadata, integration_method = "seurat", num_cores = 16, MIN_CELLS = 250, DIM = 1, THRESHOLD = -1, datasets_to_drop = datasets_to_drop) {
  packages <- c(
    "tidyverse",
    "Matrix",
    "ripserr",
    "TDA",
    "foreach",
    "doParallel",
    "parallel",
    "rliger",
    "Seurat",
    "SeuratDisk",
    "stringr",
    "pheatmap",
    "mclust",
    "aricode",
    "clusterSim",
    "Rtsne",
    "batchelor",  
    "BiocSingular",
    "scCustomize",
    "kernlab",
    "igraph",
    "progressr"
  )
  
  # Check and load required packages
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      library(pkg, character.only = TRUE)
    }
  }
  
  # Start logging
  log_file <- paste0("PH_Pipeline_Log_", Sys.Date(), ".txt")
  sink(log_file, append = TRUE)
  
  source("PH_Functions.R")
  
  log_message("Processing series with CSV input")
  
  # Extract file paths directly from the CSV
  file_paths <- metadata$'File Path'
  
  # file_paths <- file_paths[-datasets_to_drop]
  
  # Read sparse matrices into a list
  my_sparse_matrices <- list()
  
  # Process .RData files directly
  if (length(file_paths) > 0) {
    tryCatch(
      {
        sparse_matrices <- lapply(file_paths, loadRData)
        my_sparse_matrices <- sparse_matrices
      },
      error = function(e) {
        log_message(paste("Error in processing .RData files from CSV:", e$message))
      }
    )
    
    sample_names <- tryCatch(
      {
        gsub(".*/|\\.sparse.RData$", "", file_paths)
      },
      error = function(e) {
        log_message(paste("Error in defining sample names from CSV:", e$message))
        NULL
      }
    )
  }
  
  if (is.null(sample_names)) {
    return(NULL)
  }
  
  names(my_sparse_matrices) <- sample_names
  
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
        obj <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]", col.name = "percent_hb")
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
  
  
  # Add SRA, tissue, and approach metadata with detection based on data characteristics
  my_seurat_list <- tryCatch(
    {
      mclapply(seq_along(my_seurat_list), function(i) {
        seurat_obj <- my_seurat_list[[i]]
        file_path <- file_paths[i]
        
        # Add SRA and tissue metadata
        seurat_obj@meta.data$SRA <- metadata$'SRA Number'[i]
        seurat_obj@meta.data$Tissue <- metadata$Tissue[i]
        
        #seurat_obj@meta.data$Tissue <- metadata$source_name[i]
        #seurat_obj@meta.data$sample <- metadata$Run[i]
        #seurat_obj@meta.data$SRA <- metadata$`SRA Study`[i]

        # Detect scRNA-seq or snRNA-seq based on data characteristics
        median_genes_per_cell <- median(seurat_obj$nFeature_RNA)  # Median number of genes per cell
        mt_gene_percent <- mean(seurat_obj$percent_mito)  # Percentage of mitochondrial genes (if calculated)
        median_umi_counts <- median(seurat_obj$nCount_RNA)  # Median UMI counts per cell
        
        # Apply multiple heuristics to detect scRNA-seq or snRNA-seq
        if (median_genes_per_cell < 1000 && mt_gene_percent < 5 && median_umi_counts < 5000) {
          seurat_obj@meta.data$Approach <- "snRNA-seq"
        } else {
          seurat_obj@meta.data$Approach <- "scRNA-seq"
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
  

  
  prefiltered_cells_file <- "prefiltered_cells_bonemarrow.csv"

  # Function to extract columns with constant values from metadata
  get_constant_metadata <- function(obj) {
    metadata <- obj@meta.data
    constant_cols <- sapply(metadata, function(col) length(unique(col)) == 1)
    constant_metadata <- metadata[1, constant_cols, drop = FALSE]  # Keep one row with constant columns
    constant_metadata$Sample <- obj@project.name  # Add the sample name
    constant_metadata$Number_of_Cells_Before_Filtering <- nrow(metadata)  # Add cell count
    return(constant_metadata)
  }

  # After filtering, extract the updated metadata for each remaining sample
  prefiltered_cells <- tryCatch(
    {
      # Combine metadata from all Seurat objects in the list
      do.call(
        rbind,
        lapply(my_seurat_list, get_constant_metadata)
      )
    },
    error = function(e) {
      log_message(paste("Error in creating filtered metadata table:", e$message))
      NULL
    }
  )

  # Save the filtered metadata to a CSV file
  if (!is.null(prefiltered_cells)) {
    tryCatch(
      {
        write.csv(prefiltered_cells, file = prefiltered_cells_file, row.names = FALSE)
        log_message(paste("Filtered metadata saved to", prefiltered_cells_file))
      },
      error = function(e) {
        log_message(paste("Error in saving filtered metadata to CSV:", e$message))
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
  
  
  filtered_cells_file <- "filtered_cells.csv"

  # Function to extract columns with constant values from metadata
  get_constant_metadata <- function(obj) {
    metadata <- obj@meta.data
    constant_cols <- sapply(metadata, function(col) length(unique(col)) == 1)
    constant_metadata <- metadata[1, constant_cols, drop = FALSE]  # Keep one row with constant columns
    constant_metadata$Sample <- obj@project.name  # Add the sample name
    constant_metadata$Number_of_Cells_After_Filtering <- nrow(metadata)  # Add cell count
    return(constant_metadata)
  }

  # After filtering, extract the updated metadata for each remaining sample
  filtered_cells <- tryCatch(
    {
      # Combine metadata from all Seurat objects in the list
      do.call(
        rbind,
        lapply(my_seurat_list_filtered, get_constant_metadata)
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
    names <- list()
    my_seurat_list_filtered <- lapply(my_seurat_list_filtered, function(seurat_obj) {
      # Remove any existing prefixes
      seurat_obj <- RenameCells(seurat_obj, add.cell.id = NULL)
      # Save the sample identifier (from orig.ident) for later use
      names <<- append(x = names, values = seurat_obj@meta.data$orig.ident[1])
      return(seurat_obj)
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
        cat("The list 'my_seurat_list_filtered' does not exist. Please ensure it is created before running this code.\n")
      }
      
      if (exists("merged_seurat_unintegrated")) {
        cat("The object 'merged_seurat_unintegrated' already exists. Skipping merging and SCTransform steps.\n")
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
      log_message("Successfully named expr_list_sctransformed entries by orig.ident:")
      print(names(expr_list_sctInd))
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
      log_message("Successfully named expr_list_raw entries by orig.ident:")
      print(names(expr_list_raw))
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
    saveRDS(expr_list_sctWhole, file="expr_list_sctWhole_bonemarrow.Rds")
}
  
  # saveRDS(expr_list_sctInd, file = 'expr_list_scInd_bonemarrow.Rds')
  # saveRDS(expr_list_sctWhole, file = 'expr_list_sctWhole_bonemarrow.Rds')
  # saveRDS(expr_list_raw, file = 'expr_list_raw_bonemarrow.Rds')

  # expr_list_sctInd <- readRDS(file = 'expr_list_scInd_bonemarrow.Rds')
  # expr_list_sctWhole <- readRDS(file = 'expr_list_sctWhole_bonemarrow.Rds')
  # expr_list_raw <- readRDS(file = 'expr_list_raw_bonemarrow.Rds')
  
  # expr_list_sctInd <- readRDS("expr_list_sctInd.Rds")
  # expr_list_raw <- readRDS("expr_list_raw.Rds")
  # expr_list_sctWhole <- readRDS("expr_list_sctWhole.Rds") 
  # expr_list_integrated <- readRDS("expr_list_integrated.Rds")

  
  # timeout_datasets <- c(62, 68, 69, 72, 67, 64) 
  
  # Step 1: Process unintegrated data with the new batching, memory checking, and preschedule approach
  PD_result_unintegrated <- tryCatch(
    {
      # Call the process_expression_list_with_monitoring function with timeout_datasets
      process_expression_list_with_monitoring(
        expr_list = expr_list_sctInd, 
        DIM = 1, 
        log_message = log_message, 
        max_cores = 6, 
        memory_threshold = 0.25,  # Example: 25% memory threshold
        log_file = "progress_log_bonemarrow.csv", 
        results_file = "intermediate_results_bonemarrow.rds",
        timeout_datasets = NULL  # Pass the timeout_datasets here
      )    
    },
    error = function(e) {
      log_message(paste("Error in calculating persistence diagrams for unintegrated data:", e$message))
      NULL
    }
  )
  
  PD_result_unintegrated_RAW <- tryCatch(
    {
      # Call the process_expression_list_with_monitoring function with timeout_datasets
      process_expression_list_with_monitoring(
        expr_list = expr_list_raw, 
        DIM = 1, 
        log_message = log_message, 
        max_cores = 6, 
        memory_threshold = 0.25,  # Example: 25% memory threshold
        log_file = "progress_log_RAW_bonemarrow.csv", 
        results_file = "intermediate_results_RAW_bonemarrow.rds",
        timeout_datasets = NULL  # Pass the timeout_datasets here
      )    
    },
    error = function(e) {
      log_message(paste("Error in calculating persistence diagrams for unintegrated data:", e$message))
      NULL
    }
  )
  
  PD_result_unintegrated_sctWhole <- tryCatch(
    {
      # Call the process_expression_list_with_monitoring function with timeout_datasets
      process_expression_list_with_monitoring(
        expr_list = expr_list_sctWhole,
        DIM = 1,
        log_message = log_message,
        max_cores = 8,
        memory_threshold = 0.25, # Example: 25% memory threshold
        log_file = "progress_log_sctWhole.csv",
        results_file = "intermediate_results_sctWhole.rds",
        timeout_datasets = 1:124 # Pass the timeout_datasets here
      )
    },
    error = function(e) {
      log_message(paste("Error in calculating persistence diagrams for unintegrated data:", e$message))
      NULL
    }
  )
  #PD_result_unintegrated_sctWhole <- readRDS("intermediate_results_sctWhole_bonemarrow.rds")

  # Extract PD_list and thresholds from unintegrated PH results
  if (!is.null(PD_result_unintegrated)) {
    thresholds_unintegrated <- PD_result_unintegrated$thresholds
    PD_list_unintegrated <- PD_result_unintegrated$PD_list
    
    names(PD_list_unintegrated) <- names(expr_list_sctInd)
    names(thresholds_unintegrated) <- names(expr_list_sctInd)
    
    # Save the unintegrated PD results
    tryCatch(
      {
        saveRDS(object = PD_list_unintegrated, file = paste0("PD_list_dim", DIM, "_th", THRESHOLD, "_unintegrated_bonemarrow", ".Rds"))
        saveRDS(object = thresholds_unintegrated, file = paste0("thresholds_dim", DIM, "_unintegrated_bonemarrow", ".Rds"))
      },
      error = function(e) {
        log_message(paste("Error in saving persistence diagrams for unintegrated data:", e$message))
      }
    )
  }
  
  # Extract PD_list and thresholds from unintegrated_RAW PH results
  if (!is.null(PD_result_unintegrated_RAW)) {
    thresholds_unintegrated_RAW <- PD_result_unintegrated_RAW$thresholds
    PD_list_unintegrated_RAW <- PD_result_unintegrated_RAW$PD_list
    
    names(PD_list_unintegrated_RAW) <- names(expr_list_raw)
    names(thresholds_unintegrated_RAW) <- names(expr_list_raw)
    
    # Save the unintegrated PD results
    tryCatch(
      {
        saveRDS(object = PD_list_unintegrated_RAW, file = paste0("PD_list_dim", DIM, "_th", THRESHOLD, "_unintegrated_RAW_bonemarrow", ".Rds"))
        saveRDS(object = thresholds_unintegrated_RAW, file = paste0("thresholds_dim", DIM, "_unintegrated_RAW_bonemarrow", ".Rds"))
      },
      error = function(e) {
        log_message(paste("Error in saving persistence diagrams for unintegrated data RAW:", e$message))
      }
    )
  }
  
  # Extract PD_list and thresholds from unintegrated_sctWhole PH results
  if (!is.null(PD_result_unintegrated_sctWhole)) {
    thresholds_unintegrated_sctWhole <- PD_result_unintegrated_sctWhole$thresholds
    PD_list_unintegrated_sctWhole <- PD_result_unintegrated_sctWhole$PD_list

    names(PD_list_unintegrated_sctWhole) <- names(expr_list_sctWhole)
    names(thresholds_unintegrated_sctWhole) <- names(expr_list_sctWhole)

    # Save the unintegrated PD results
    tryCatch(
      {
        saveRDS(object = PD_list_unintegrated_sctWhole, file = paste0("PD_list_dim", DIM, "_th", THRESHOLD, "_unintegrated_sctWhole_bonemarrow", ".Rds"))
        saveRDS(object = thresholds_unintegrated_sctWhole, file = paste0("thresholds_dim", DIM, "_unintegrated_sctWhole_bonemarrow", ".Rds"))
      },
      error = function(e) {
        log_message(paste("Error in saving persistence diagrams for unintegrated data sctWhole:", e$message))
      }
    )
  }
  
  #load integration function
  source("Integration_flexible.R")

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
    
    saveRDS(expr_list_integrated, file = "expr_list_integrated.rds")
  }
  
  # Step 8: Process integrated data with the new batching, memory checking, and preschedule approach
  PD_result_integrated <- tryCatch(
    {
      process_expression_list_with_monitoring(
        expr_list = expr_list_integrated, 
        DIM = DIM, 
        log_message = log_message, 
        max_cores = 8, 
        memory_threshold = 0.25,  # Example: 25% memory threshold
        log_file = "progress_log_integrated.csv", 
        results_file = "intermediate_results_integrated.rds"
      )
    },
    error = function(e) {
      log_message(paste("Error in calculating persistence diagrams for integrated data:", e$message))
      NULL
    }
  )
  
  # Extract PD_list and thresholds from integrated PH results
  if (!is.null(PD_result_integrated)) {
    thresholds_integrated <- PD_result_integrated$thresholds
    PD_list_integrated <- PD_result_integrated$PD_list
    
    # Save the integrated PD results
    tryCatch(
      {
        saveRDS(object = PD_list_integrated, file = paste0("PD_list_dim", DIM, "_th", THRESHOLD, "_integrated", ".Rds"))
        saveRDS(object = thresholds_integrated, file = paste0("thresholds_dim", DIM, "_integrated", ".Rds"))
      },
      error = function(e) {
        log_message(paste("Error in saving persistence diagrams for integrated data:", e$message))
      }
    )
  }
  
  # Retry for Unintegrated Data
  PD_list_after_retries_unintegrated_sctInd <- tryCatch(
    {
      retry_pd_calculation(
        progress_log = "progress_log_unintegrated_sctInd.csv",  # Original progress log for unintegrated data
        results_file = "intermediate_results_unintegrated_sctInd.rds",  # Original results file for unintegrated data
        expr_list = expr_list_sctInd,  # Expression list for unintegrated data
        DIM = 1,  # Dimension for PH calculation
        doubling_limit = 5,  # Retry limit for threshold doubling
        time_limit = 16 * 3600,  # 12 hours for each retry
        num_cores = 6,  # Number of cores for parallel processing
        log_message = log_message,  # Logging function
        retry_progress_log = "retry_progress_log_unintegrated_sctInd.csv",  # Separate retry progress log
        retry_results_file = "retry_results_unintegrated_sctInd.rds",  # Separate retry results file
        overwrite_original = FALSE  # Avoid overwriting original logs and results
      )
    },
    error = function(e) {
      log_message(paste("Error during retry of failed persistence diagrams for unintegrated data:", e$message))
      NULL
    }
  )

    # Retry for Unintegrated Data
  PD_list_after_retries_unintegrated_sctWhole <- tryCatch(
    {
      retry_pd_calculation(
        progress_log = "progress_log_unintegrated_sctWhole.csv",  # Original progress log for unintegrated data
        results_file = "intermediate_results_unintegrated_sctWhole.rds",  # Original results file for unintegrated data
        expr_list = expr_list_sctWhole,  # Expression list for unintegrated data
        DIM = 1,  # Dimension for PH calculation
        doubling_limit = 5,  # Retry limit for threshold doubling
        time_limit = 16 * 3600,  # 12 hours for each retry
        num_cores = 6,  # Number of cores for parallel processing
        log_message = log_message,  # Logging function
        retry_progress_log = "retry_progress_log_unintegrated_sctWhole.csv",  # Separate retry progress log
        retry_results_file = "retry_results_unintegrated_sctWhole.rds",  # Separate retry results file
        overwrite_original = FALSE  # Avoid overwriting original logs and results
      )
    },
    error = function(e) {
      log_message(paste("Error during retry of failed persistence diagrams for unintegrated data:", e$message))
      NULL
    }
  )
  
  # Retry for Unintegrated Data
  PD_list_after_retries_unintegrated_Raw <- tryCatch(
    {
      retry_pd_calculation(
        progress_log = "progress_log_RAW.csv",  # Original progress log for unintegrated data
        results_file = "intermediate_results_RAW.rds",  # Original results file for unintegrated data
        expr_list = expr_list_raw,  # Expression list for unintegrated data
        DIM = 1,  # Dimension for PH calculation
        doubling_limit = 5,  # Retry limit for threshold doubling
        time_limit = 16 * 3600,  # 12 hours for each retry
        num_cores = 6,  # Number of cores for parallel processing
        log_message = log_message,  # Logging function
        retry_progress_log = "retry_progress_log_unintegrated_RAW.csv",  # Separate retry progress log
        retry_results_file = "retry_results_unintegrated_RAW.rds",  # Separate retry results file
        overwrite_original = FALSE  # Avoid overwriting original logs and results
      )
    },
    error = function(e) {
      log_message(paste("Error during retry of failed persistence diagrams for unintegrated data:", e$message))
      NULL
    }
  )
  
  # Retry for Integrated Data
  PD_list_after_retries_integrated <- tryCatch(
    {
      retry_pd_calculation(
        progress_log = "progress_log_integrated.csv",  # Original progress log for integrated data
        results_file = "intermediate_results_integrated.rds",  # Original results file for integrated data
        expr_list = expr_list_integrated,  # Expression list for integrated data
        DIM = 1,  # Dimension for PH calculation
        doubling_limit = 5,  # Retry limit for threshold doubling
        time_limit = 12 * 3600,  # 12 hours for each retry
        num_cores = 32,  # Number of cores for parallel processing
        log_message = log_message,  # Logging function
        retry_progress_log = "retry_progress_log_integrated.csv",  # Separate retry progress log
        retry_results_file = "retry_results_integrated.rds",  # Separate retry results file
        overwrite_original = FALSE  # Avoid overwriting original logs and results
      )
    },
    error = function(e) {
      log_message(paste("Error during retry of failed persistence diagrams for integrated data:", e$message))
      NULL
    }
  )
  
    # Helper function to update names in PD list
  update_pd_list_names <- function(PD_list, expr_list_names, log_message) {
    tryCatch({
      if (length(PD_list) != length(expr_list_names)) {
        stop("Length of PD list does not match length of expression list names.")
      }
      names(PD_list) <- expr_list_names
      log_message("PD list names successfully updated.")
      return(PD_list)
    }, error = function(e) {
      log_message(paste("Error in updating PD list names:", e$message))
      return(PD_list)
    })
  }

  # Save results after retries for unintegrated SCT Individual
  if (!is.null(PD_list_after_retries_unintegrated_sctInd)) {
    tryCatch(
      {
        expr_list_names <- names(expr_list_sctInd) # Replace with the correct expression list
        PD_list_after_retries_unintegrated_sctInd <- update_pd_list_names(
          PD_list_after_retries_unintegrated_sctInd$PD_list,
          expr_list_names,
          log_message
        )
        saveRDS(object = PD_list_after_retries_unintegrated_sctInd, file = "PD_list_after_retries_unintegrated_sctInd.rds")
        log_message("Final PD list for unintegrated SCT Individual saved successfully.")
      },
      error = function(e) {
        log_message(paste("Error in saving final PD list for unintegrated SCT Individual:", e$message))
      }
    )
  }

  # Save results after retries for unintegrated SCT Whole
  if (!is.null(PD_list_after_retries_unintegrated_sctWhole)) {
    tryCatch(
      {
        expr_list_names <- names(expr_list_sctWhole) # Replace with the correct expression list
        PD_list_after_retries_unintegrated_sctWhole <- update_pd_list_names(
          PD_list_after_retries_unintegrated_sctWhole$PD_list,
          expr_list_names,
          log_message
        )
        saveRDS(object = PD_list_after_retries_unintegrated_sctWhole, file = "PD_list_after_retries_unintegrated_sctWhole.rds")
        log_message("Final PD list for unintegrated SCT Whole saved successfully.")
      },
      error = function(e) {
        log_message(paste("Error in saving final PD list for unintegrated SCT Whole:", e$message))
      }
    )
  }

  # Save results after retries for unintegrated Raw
  if (!is.null(PD_list_after_retries_unintegrated_Raw)) {
    tryCatch(
      {
        expr_list_names <- names(expr_list_raw) # Replace with the correct expression list
        PD_list_after_retries_unintegrated_Raw <- update_pd_list_names(
          PD_list_after_retries_unintegrated_Raw$PD_list,
          expr_list_names,
          log_message
        )
        saveRDS(object = PD_list_after_retries_unintegrated_Raw, file = "PD_list_after_retries_unintegrated_Raw.rds")
        log_message("Final PD list for unintegrated Raw saved successfully.")
      },
      error = function(e) {
        log_message(paste("Error in saving final PD list for unintegrated Raw:", e$message))
      }
    )
  }

  # Save results after retries for integrated
  if (!is.null(PD_list_after_retries_integrated)) {
    tryCatch(
      {
        expr_list_names <- names(expr_list_integrated) # Replace with the correct expression list
        PD_list_after_retries_integrated <- update_pd_list_names(
          PD_list_after_retries_integrated$PD_list,
          expr_list_names,
          log_message
        )
        saveRDS(object = PD_list_after_retries_integrated, file = "PD_list_after_retries_integrated.rds")
        log_message("Final PD list for integrated data saved successfully.")
      },
      error = function(e) {
        log_message(paste("Error in saving final PD list for integrated data:", e$message))
      }
    )
  }


  #   # Save results after retries for integrated
  # if (!is.null(PD_list_integrated)) {
  #   tryCatch(
  #     {
  #       expr_list_names <- names(expr_list_integrated) # Replace with the correct expression list
  #       PD_list_integrated <- update_pd_list_names(
  #         PD_list_integrated,
  #         expr_list_names,
  #         log_message
  #       )
  #       saveRDS(object = PD_list_integrated, file = "PD_list_dim1_th-1_integrated_bonemarrow.Rds")
  #       log_message("Final PD list for integrated data saved successfully.")
  #     },
  #     error = function(e) {
  #       log_message(paste("Error in saving final PD list for integrated data:", e$message))
  #     }
  #   )
  # }

#   # Load required libraries
# library(purrr)

# # Function to add missing names to PD list
# ensure_pd_list_names <- function(pd_list, expr_list_names, log_message) {
#   if (is.null(names(pd_list)) || length(names(pd_list)) != length(pd_list)) {
#     log_message("Names missing or mismatched in PD list. Updating names.")
#     names(pd_list) <- expr_list_names
#   } else {
#     log_message("PD list names are already correct.")
#   }
#   return(pd_list)
# }

# # Load PD_list_after_retries_integrated
# tryCatch(
#   {
#     PD_list_after_retries_integrated <- readRDS("PD_list_after_retries_integrated.rds")
#     log_message("Loaded PD list successfully.")

#     # Get the correct names from the expression list
#     expr_list_names <- names(expr_list_integrated) # Replace this with your actual expression list names source

#     # Ensure the PD list has the correct names
#     PD_list_after_retries_integrated <- ensure_pd_list_names(
#       PD_list_after_retries_integrated,
#       expr_list_names,
#       log_message
#     )

#     # Save the updated PD list
#     saveRDS(PD_list_after_retries_integrated, "PD_list_after_retries_integrated.rds")
#     log_message("Updated and saved PD list with correct names.")
#   },
#   error = function(e) {
#     log_message(paste("Error in loading or updating PD list:", e$message))
#   }
# )



  
  #### BDM CALCULATIONS ####
  ## Function to compute and save distance matrices (BDM and SDM)

compute_and_save_distance_matrices <- function(
    pd_list, expr_list_names, dataset_name, num_cores, log_message, 
    dimension = 0, save_sdm = TRUE, save_bdm = TRUE) {
  
  # Create Bottleneck Distance Matrix
  if (save_bdm) {
    log_message(paste("Creating Bottleneck Distance Matrix for", dataset_name))
    BDM <- tryCatch(
      {
        CreateBottleneckDistanceMatrixParallel(pd_list, max_cores = num_cores, log_message = log_message, dimension = dimension)
      },
      error = function(e) {
        log_message(paste("Error in creating BDM for", dataset_name, ":", e$message))
        NULL
      }
    )
    
    if (!is.null(BDM)) {
      tryCatch(
        {
          BDM_df <- as.data.frame(BDM)
          rownames(BDM_df) <- expr_list_names
          colnames(BDM_df) <- expr_list_names
          
          saveRDS(BDM_df, file = paste0("BDM_", dataset_name, ".Rds"))
          write.csv(BDM_df, file = paste0("BDM_", dataset_name, ".csv"))
          log_message(paste("BDM for", dataset_name, "saved successfully."))
        },
        error = function(e) {
          log_message(paste("Error in saving BDM for", dataset_name, ":", e$message))
        }
      )
    }
  }
  
  # Create Spectral Distance Matrix
  if (save_sdm) {
    log_message(paste("Creating Spectral Distance Matrix for", dataset_name))
    SDM <- tryCatch(
      {
        CreateSpectralDistanceMatrixFromPD(pd_list, num_eigen = 50, log_message = log_message, dimension = dimension)
      },
      error = function(e) {
        log_message(paste("Error in creating SDM for", dataset_name, ":", e$message))
        NULL
      }
    )
    
    if (!is.null(SDM)) {
      tryCatch(
        {
          SDM_df <- as.data.frame(SDM)
          rownames(SDM_df) <- expr_list_names
          colnames(SDM_df) <- expr_list_names
          
          saveRDS(SDM_df, file = paste0("SDM_", dataset_name, ".Rds"))
          write.csv(SDM_df, file = paste0("SDM_", dataset_name, ".csv"))
          log_message(paste("SDM for", dataset_name, "saved successfully."))
        },
        error = function(e) {
          log_message(paste("Error in saving SDM for", dataset_name, ":", e$message))
        }
      )
    }
  }
}

# Apply to different datasets
log_message("Starting distance matrix computation...")

compute_and_save_distance_matrices(
  pd_list = PD_list_after_retries_unintegrated_scTransformed,
  expr_list_names = names(expr_list),
  dataset_name = "unintegrated_sctInd",
  num_cores = num_cores,
  log_message = log_message
)

compute_and_save_distance_matrices(
  pd_list = PD_list_unintegrated_sctWhole,
  expr_list_names = names(expr_list_sctWhole),
  dataset_name = "unintegrated_scTransformed_whole",
  num_cores = num_cores,
  log_message = log_message
)

compute_and_save_distance_matrices(
  pd_list = PD_list_unintegrated_RAW,
  expr_list_names = names(expr_list_raw),
  dataset_name = "unintegrated_Raw",
  num_cores = num_cores,
  log_message = log_message
)

compute_and_save_distance_matrices(
  pd_list = PD_list_integrated,
  expr_list_names = names(expr_list_integrated),
  dataset_name = "integrated",
  num_cores = num_cores,
  log_message = log_message
)

log_message("Distance matrix computation completed.")