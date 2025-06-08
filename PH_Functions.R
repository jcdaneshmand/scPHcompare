# --- 0. Setup: Load Libraries and Source External Functions ---

# Recommended: Use pacman for easier package management
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(tidyverse, Matrix, ripserr, TDA, foreach, doParallel, parallel,
#                Seurat, SeuratObject, SeuratDisk, stringr, pheatmap, mclust, aricode,
#                clusterSim, Rtsne, batchelor, BiocSingular, scCustomize, kernlab,
#                igraph, progressr, plyr)

packages <- c(
  "tidyverse", "Matrix", "ripserr", "TDA", "foreach", "doParallel", "parallel",
  "Seurat", "SeuratObject", "SeuratDisk", "stringr", "pheatmap", "mclust", "aricode",
  "clusterSim", "Rtsne", "batchelor", "BiocSingular", "scCustomize", "kernlab",
  "igraph", "progressr", "plyr" # Added plyr for rbind.fill
)

# Check and load required packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Optionally install missing packages
    # install.packages(pkg)
    # BiocManager::install(pkg) # For Bioconductor packages if needed
    stop(paste("Package", pkg, "is required but not installed."))
  }
  library(pkg, character.only = TRUE)
}

# Source helper functions defined elsewhere
# Ensure these files exist and contain the necessary functions:
# - PH_Functions.R: Defines log_message, loadRData,
#                   process_expression_list_with_monitoring, retry_pd_calculation
# - Integration_flexible.R: Defines perform_integration (now Seurat-only)
source("PH_Functions.R")
source("Integration_flexible.R") # Source the updated Seurat-only version

# Built-in Seurat cell cycle genes (Human) - Define alternatives if needed
# s.genes <- cc.genes.updated.2019$s.genes
# g2m.genes <- cc.genes.updated.2019$g2m.genes


# --- Helper Function Definitions ---

#' Setup Environment and Logging
#' @param output_dir Directory for output.
#' @param log_file_base Base name for the log file.
#' @return Path to the log file.
#' @noRd
setup_environment <- function(output_dir, log_file_base = "PH_Pipeline_Log") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("Created output directory:", output_dir))
  } else {
    message(paste("Using existing output directory:", output_dir))
  }
  log_file <- file.path(output_dir, paste0(log_file_base, "_", Sys.Date(), ".txt"))
  # Ensure log_message is defined globally or sourced from PH_Functions.R
  # sink() will direct message() output to the log file as well if log_message uses message()
  sink(log_file, append = TRUE, split = TRUE) # split=TRUE shows output on console too
  # Define log_message or ensure it's sourced before first use
  if (!exists("log_message", mode = "function")) {
    stop("log_message function not found. Please define it or source PH_Functions.R")
  }
  log_message(paste("Logging started:", Sys.time()))
  log_message(paste("Output directory:", output_dir))
  log_message(paste("Log file:", log_file))
  return(log_file)
}

#' Load Metadata
#' @param metadata_path Path to the metadata CSV.
#' @return A data frame or tibble, or NULL on error.
#' @noRd
load_metadata <- function(metadata_path) {
  log_message(paste("Loading metadata from:", metadata_path))
  tryCatch({
    if (!file.exists(metadata_path)) {
      stop(paste("Metadata file not found at:", metadata_path))
    }
    metadata <- read_csv(metadata_path, show_col_types = FALSE)
    # Basic validation - Adjust required columns as needed
    required_cols <- c('File Path', 'Sample Name', 'SRA Number', 'Tissue') # Add 'Approach' if strictly required
    missing_cols <- setdiff(required_cols, names(metadata))
    if (length(missing_cols) > 0) {
      stop(paste("Metadata CSV is missing required column(s):", paste(missing_cols, collapse = ", ")))
    }
    log_message("Metadata loaded successfully.")
    return(metadata)
  }, error = function(e) {
    log_message(paste("Error loading metadata:", e$message), type = "ERROR")
    return(NULL)
  })
}

#' Load Data and Create Seurat Objects
#' @param metadata Metadata dataframe containing 'File Path' and 'Sample Name'.
#' @param num_cores Number of cores.
#' @return A list of Seurat objects, or NULL on error.
#' @noRd
load_and_create_seurat <- function(metadata, num_cores) {
  log_message("Loading sparse matrices and creating Seurat objects...")
  file_paths <- metadata$'File Path'
  sample_names <- metadata$'Sample Name' # Use the dedicated sample name column

  if (length(file_paths) == 0) {
    log_message("No file paths found in metadata.", type = "WARN")
    return(list()) # Return empty list
  }
  if (length(file_paths) != length(sample_names)) {
    log_message("Mismatch between number of file paths and sample names in metadata.", type = "ERROR")
    return(NULL)
  }

  # Ensure loadRData is defined (e.g., in PH_Functions.R)
  if (!exists("loadRData", mode = "function")) stop("loadRData function not found.")

  sparse_matrices <- tryCatch({
    log_message(paste("Loading", length(file_paths), ".RData files..."))
    mclapply(file_paths, function(fpath) {
      if (!file.exists(fpath)) {
        log_message(paste("File not found:", fpath), type = "WARN")
        return(NULL)
      }
      loadRData(fpath)
    }, mc.cores = num_cores)
  }, error = function(e) {
    log_message(paste("Error loading .RData files:", e$message), type = "ERROR")
    return(NULL)
  })

  if (is.null(sparse_matrices)) return(NULL)

  # Assign sample names from metadata
  names(sparse_matrices) <- sample_names

  seurat_list <- tryCatch({
    log_message("Creating Seurat objects...")
    lapply(seq_along(sparse_matrices), function(i) {
      spm <- sparse_matrices[[i]]
      name <- names(sparse_matrices)[i]
      if (is.null(spm)) {
        log_message(paste("Warning: No matrix loaded for sample:", name, "- Skipping."), type = "WARN")
        return(NULL) # Return NULL for this element
      }
      if (!inherits(spm, "dgCMatrix")) {
        log_message(paste("Warning: Loaded object for sample", name, "is not a dgCMatrix. Attempting conversion."), type = "WARN")
        spm <- tryCatch(as(spm, "dgCMatrix"), error = function(e) NULL)
        if (is.null(spm)) {
          log_message(paste("Error: Could not convert object for sample", name, "to dgCMatrix. Skipping."), type = "ERROR")
          return(NULL)
        }
      }
      if (nrow(spm) == 0 || ncol(spm) == 0) {
        log_message(paste("Warning: Matrix for sample", name, "is empty. Skipping."), type = "WARN")
        return(NULL)
      }
      CreateSeuratObject(counts = spm, project = name)
    }) # Consider mc.cores=num_cores for large lists
  }, error = function(e) {
    log_message(paste("Error creating Seurat objects:", e$message), type = "ERROR")
    return(NULL)
  })

  if (is.null(seurat_list)) return(NULL)

  # Filter out NULLs
  original_count <- length(seurat_list)
  seurat_list <- seurat_list[!sapply(seurat_list, is.null)]
  final_count <- length(seurat_list)

  if (final_count > 0) {
    log_message(paste("Successfully created", final_count, "Seurat objects (out of", original_count, "initial attempts)."))
  } else {
    log_message("No Seurat objects were successfully created.", type = "ERROR")
    return(NULL)
  }
  return(seurat_list)
}


#' Add QC Metrics and Metadata to Seurat Objects
#' @param seurat_list List of Seurat objects.
#' @param metadata Metadata dataframe.
#' @param metadata_sample_col Column name in 'metadata' matching Seurat Project names.
#' @param num_cores Number of cores.
#' @return Updated list of Seurat objects, or NULL on error.
#' @noRd
add_qc_and_metadata <- function(seurat_list, metadata, metadata_sample_col = "Sample Name", num_cores) {
  log_message("Calculating QC metrics (Mito, Ribo, HB)...")
  seurat_list_updated <- tryCatch({
    mclapply(seurat_list, function(obj) {
      DefaultAssay(obj) <- "RNA"
      obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent_mito", assay = "RNA")
      obj <- PercentageFeatureSet(obj, pattern = "^RP[SL]", col.name = "percent_ribo", assay = "RNA")
      obj <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]", col.name = "percent_hb", assay = "RNA")
      return(obj)
    }, mc.cores = num_cores)
  }, error = function(e) {
    log_message(paste("Error calculating QC metrics:", e$message), type = "ERROR");
    return(NULL)
  })
  if (is.null(seurat_list_updated)) return(NULL)

  log_message("Adding metadata (SRA, Tissue, Approach)...")
  if (!metadata_sample_col %in% names(metadata)) {
    log_message(paste("Error: Metadata matching column '", metadata_sample_col, "' not found."), type = "ERROR");
    return(NULL)
  }
  metadata_lookup <- setNames(as.list(1:nrow(metadata)), metadata[[metadata_sample_col]])

  seurat_list_final <- tryCatch({
    lapply(seq_along(seurat_list_updated), function(i) {
      seurat_obj <- seurat_list_updated[[i]]
      obj_project_name <- Project(seurat_obj)
      metadata_row_idx <- metadata_lookup[[obj_project_name]]

      if (is.null(metadata_row_idx)) {
        log_message(paste("Warning: No metadata for sample:", obj_project_name), type = "WARN")
        seurat_obj@meta.data$SRA <- NA_character_;
        seurat_obj@meta.data$Tissue <- NA_character_;
        seurat_obj@meta.data$Approach <- NA_character_
      } else {
        metadata_row <- metadata[metadata_row_idx,]
        seurat_obj@meta.data$SRA <- .subset2(metadata_row, 'SRA Number')[1]
        seurat_obj@meta.data$Tissue <- .subset2(metadata_row, 'Tissue')[1]
        # Add other relevant columns if they exist
        for (col_name in c("Series", "Platform", "Instrument")) {
          # Add others as needed
          if (col_name %in% names(metadata_row)) {
            seurat_obj@meta.data[[col_name]] <- .subset2(metadata_row, col_name)[1]
          }
        }

        approach_metadata <- if ("Approach" %in% names(metadata_row)) .subset2(metadata_row, 'Approach')[1] else NA_character_
        if (!is.na(approach_metadata) && nzchar(trimws(approach_metadata))) {
          seurat_obj@meta.data$Approach <- approach_metadata
        } else {
          seurat_obj@meta.data$Approach <- NA_character_
          log_message(paste("Approach missing/NA for sample:", obj_project_name, ". Using heuristic."), type = "INFO")
        }
      }

      if (is.na(seurat_obj@meta.data$Approach) || !nzchar(trimws(seurat_obj@meta.data$Approach))) {
        median_genes <- median(seurat_obj$nFeature_RNA, na.rm = TRUE)
        mean_mito <- mean(seurat_obj$percent_mito, na.rm = TRUE)
        median_umi <- median(seurat_obj$nCount_RNA, na.rm = TRUE)
        if (any(is.na(c(median_genes, mean_mito, median_umi)))) {
          log_message(paste("WARN: Cannot run heuristic for", obj_project_name, "missing QC."), type = "WARN")
          detected_approach <- "Unknown"
        } else {
          detected_approach <- ifelse(median_genes < 1000 && mean_mito < 5 && median_umi < 5000, "snRNA-seq", "scRNA-seq")
          log_message(paste("Heuristic detected:", detected_approach, "for sample:", obj_project_name))
        }
        if (is.na(seurat_obj@meta.data$Approach) || !nzchar(trimws(seurat_obj@meta.data$Approach))) {
          seurat_obj@meta.data$Approach <- detected_approach
        }
      }
      return(seurat_obj)
    })
  }, error = function(e) {
    log_message(paste("Error adding metadata:", e$message), type = "ERROR");
    return(NULL)
  })

  if (is.null(seurat_list_final)) { log_message("Failed during metadata add step.", type = "ERROR"); return(NULL) }
  log_message("QC metrics and metadata added.")
  return(seurat_list_final)
}

#' Save Metadata Summary
#' @param seurat_list List of Seurat objects.
#' @param file_path Path to save the CSV summary.
#' @param before_filtering Boolean indicating if this is before or after filtering.
#' @return Invisibly returns the summary dataframe, or NULL on error.
#' @noRd
save_metadata_summary <- function(seurat_list, file_path, before_filtering = TRUE) {
  stage <- ifelse(before_filtering, "Before", "After")
  cell_count_col <- paste0("Number_of_Cells_", stage, "_Filtering")
  log_message(paste("Creating metadata summary", stage, "filtering..."))
  if (is.null(seurat_list) || length(seurat_list) == 0) { log_message("List empty.", type = "WARN"); return(invisible(NULL)) }

  get_constant_metadata <- function(obj) {
    metadata <- obj@meta.data
    constant_cols <- sapply(names(metadata), function(col_name) {
      if (col_name %in% c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb", "S.Score", "G2M.Score")) return(FALSE)
      length(unique(na.omit(metadata[[col_name]]))) <= 1
    })
    constant_col_names <- names(constant_cols[constant_cols])
    essential_cols <- c("orig.ident", "SRA", "Tissue", "Approach")
    cols_to_keep <- unique(c(constant_col_names, intersect(essential_cols, names(metadata))))
    summary_data <- metadata[1, cols_to_keep, drop = FALSE]
    summary_data$Sample <- Project(obj)
    summary_data[[cell_count_col]] <- ncol(obj)
    return(summary_data)
  }

  summary_table <- tryCatch({
    plyr::rbind.fill(lapply(seurat_list, get_constant_metadata))
  }, error = function(e) { log_message(paste("Error creating summary:", e$message), type = "ERROR"); return(NULL) })

  if (!is.null(summary_table)) {
    tryCatch({ write.csv(summary_table, file = file_path, row.names = FALSE); log_message(paste("Metadata summary saved:", file_path)) },
                 error = function(e) { log_message(paste("Error saving summary CSV:", e$message), type = "ERROR") })
  }
  invisible(summary_table)
}

#' Filter Seurat Objects Based on QC Metrics
#' @param seurat_list List of Seurat objects.
#' @param params List of filtering parameters.
#' @param min_cells Minimum cells required per dataset after filtering.
#' @param num_cores Number of cores.
#' @return Filtered list of Seurat objects, or NULL on error.
#' @noRd
filter_seurat_objects <- function(seurat_list, params, min_cells, num_cores) {
  log_message("Filtering Seurat objects based on cell and gene metrics...")
  if (is.null(seurat_list) || length(seurat_list) == 0) return(list())

  seurat_list_filtered <- tryCatch({
    mclapply(seurat_list, function(obj) {
      sample_name <- Project(obj)
      # log_message(paste("Filtering sample:", sample_name)) # Reduce verbosity
      cells_before <- ncol(obj);
      genes_before <- nrow(obj)
      filter_expr <- expression(nFeature_RNA >= params$minGenesPerCell & nFeature_RNA <= params$maxGenesPerCell & percent_mito <= params$maxMitoPercent & percent_ribo >= params$minRiboPercent)
      qc_cols_present <- c("nFeature_RNA", "percent_mito", "percent_ribo") %in% colnames(obj@meta.data)
      if (!all(qc_cols_present)) log_message(paste("WARN: Missing QC columns for", sample_name), type = "WARN")
      cells_to_keep <- WhichCells(obj, expression = filter_expr)
      if (length(cells_to_keep) == 0) { log_message(paste("WARN: No cells passed cell filtering:", sample_name), type = "WARN"); return(NULL) }
      obj_filtered_cells <- subset(obj, cells = cells_to_keep)
      cells_after_cellfilt <- ncol(obj_filtered_cells)

      if (!"RNA" %in% Assays(obj_filtered_cells) || !"counts" %in% slotNames(GetAssay(obj_filtered_cells, "RNA"))) {
        log_message(paste("WARN: RNA/counts missing for", sample_name, "Cannot filter genes."), type = "WARN");
        obj_filtered_genes <- obj_filtered_cells
      } else {
        counts_matrix <- GetAssayData(obj_filtered_cells, assay = "RNA", slot = "counts")
        if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
          log_message(paste("WARN: Counts matrix empty for", sample_name), type = "WARN");
          obj_filtered_genes <- obj_filtered_cells
        } else {
          genes_to_keep_logical <- Matrix::rowSums(counts_matrix > 0) >= params$minCellsPerGene
          genes_to_keep_names <- rownames(counts_matrix)[genes_to_keep_logical]
          if (length(genes_to_keep_names) == 0) { log_message(paste("WARN: No genes passed filtering:", sample_name), type = "WARN"); return(NULL) }
          obj_filtered_genes <- subset(obj_filtered_cells, features = genes_to_keep_names)
        }
        rm(counts_matrix);
        gc()
      }
      log_message(sprintf("Filtered sample: %s | Cells: %d->%d->%d | Genes: %d->%d", sample_name, cells_before, cells_after_cellfilt, ncol(obj_filtered_genes), genes_before, nrow(obj_filtered_genes)))
      return(obj_filtered_genes)
    }, mc.cores = num_cores)
  }, error = function(e) { log_message(paste("Error parallel filtering:", e$message), type = "ERROR"); return(NULL) })

  if (is.null(seurat_list_filtered)) { log_message("Filtering failed.", type = "ERROR"); return(NULL) }
  seurat_list_filtered <- seurat_list_filtered[!sapply(seurat_list_filtered, is.null)]
  if (length(seurat_list_filtered) == 0) { log_message("All datasets removed during QC filtering.", type = "WARN"); return(NULL) }
  log_message(paste("QC filtering complete.", length(seurat_list_filtered), "datasets remaining."))

  log_message(paste("Removing datasets with <", min_cells, "cells..."))
  initial_count <- length(seurat_list_filtered)
  seurat_list_final <- Filter(function(obj) { keep <- ncol(obj) >= min_cells; if (!keep) log_message(paste("Discarding", Project(obj), ":", ncol(obj), "cells."), type = "INFO"); return(keep) }, seurat_list_filtered)
  final_count <- length(seurat_list_final)
  log_message(paste("Removed", initial_count - final_count, "datasets.", final_count, "datasets remain."))
  if (final_count == 0) { log_message("No datasets meet min cell requirement.", type = "WARN"); return(NULL) }
  return(seurat_list_final)
}


#' Normalize (SCTransform) and Prepare Seurat Data
#' @param seurat_list Filtered list of Seurat objects.
#' @param sct_params Parameters for SCTransform.
#' @param regress_params List controlling regression (mito, ribo, cell_cycle).
#' @param species Species for cell cycle genes.
#' @param num_cores Number of cores.
#' @param output_dir Output directory.
#' @return List containing 'seurat_list_normalized' and 'merged_seurat_unintegrated', or NULL.
#' @noRd
normalize_and_prepare_seurat <- function(seurat_list, sct_params, regress_params, species, num_cores, output_dir) {
  log_message("Running individual SCTransform normalization...")
  if (is.null(seurat_list) || length(seurat_list) == 0) return(NULL)

  seurat_list_normalized <- tryCatch({
    mclapply(seurat_list, function(x) {
      sample_name <- Project(x);
      log_message(paste("SCTransform:", sample_name))
      DefaultAssay(x) <- "RNA"
      if (!"data" %in% slotNames(GetAssay(x, "RNA"))) x <- NormalizeData(x, verbose = FALSE)
      x <- SCTransform(x, assay = "RNA", new.assay.name = "SCT", variable.features.n = sct_params$variable.features.n, verbose = FALSE) # Verbose FALSE
      DefaultAssay(x) <- "SCT"
      return(x)
    }, mc.cores = num_cores)
  }, error = function(e) { log_message(paste("Error individual SCT:", e$message), type = "ERROR"); return(NULL) })
  if (is.null(seurat_list_normalized) || length(seurat_list_normalized) == 0) { log_message("Individual SCT failed.", type = "ERROR"); return(NULL) }
  log_message("Individual SCTransform complete.")

  merged_unintegrated_path <- file.path(output_dir, "merged_seurat_unintegrated.rds")
  merged_seurat_unintegrated <- NULL
  # Add logic to load existing if desired here...

  if (is.null(merged_seurat_unintegrated)) {
    log_message("Creating merged Seurat object with RNA, SCT(regressed), SCT_Ind assays...")
    project_names <- sapply(seurat_list_normalized, Project)
    seurat_list_renamed <- lapply(seq_along(seurat_list_normalized), function(i) RenameCells(seurat_list_normalized[[i]], new.names = paste0(project_names[i], "_", colnames(seurat_list_normalized[[i]]))))
    names(seurat_list_renamed) <- project_names
    if (length(seurat_list_renamed) == 0) { log_message("No objects for merging.", type = "ERROR"); return(NULL) }

    if (length(seurat_list_renamed) == 1) {
      log_message("Single dataset: creating merged structure...")
      merged_seurat_unintegrated <- seurat_list_renamed[[1]]
      # Ensure RNA data slot exists
      if (!"data" %in% slotNames(GetAssay(merged_seurat_unintegrated, "RNA"))) merged_seurat_unintegrated <- NormalizeData(merged_seurat_unintegrated, assay = "RNA", verbose = FALSE)
      # Cell Cycle (Single)
      log_message("Cell cycle scoring (single)...")
      tryCatch({
        if (tolower(species) == "human") { s.genes <- cc.genes.updated.2019$s.genes; g2m.genes <- cc.genes.updated.2019$g2m.genes } else { stop("Custom species CC genes needed.") }
        merged_seurat_unintegrated <- CellCycleScoring(merged_seurat_unintegrated, s.features = s.genes, g2m.features = g2m.genes, assay = 'RNA', set.ident = FALSE)
      }, error = function(e) { log_message(paste("WARN: CC scores failed:", e$message), type = "WARN") })
      # Assays (Single)
      merged_seurat_unintegrated[["SCT_Ind"]] <- merged_seurat_unintegrated[["SCT"]]
      log_message("Running SCT on RNA assay (single)...")
      vars_sct <- c();
      if (isTRUE(regress_params$regress_mito) && "percent_mito" %in% colnames(merged_seurat_unintegrated@meta.data)) vars_sct <- c(vars_sct, "percent_mito");
      # ... add ribo, cc checks ...
      if (isTRUE(regress_params$regress_ribo) && "percent_ribo" %in% colnames(merged_seurat_unintegrated@meta.data)) vars_sct <- c(vars_sct, "percent_ribo")
      if (isTRUE(regress_params$regress_cell_cycle) && all(c("S.Score", "G2M.Score") %in% colnames(merged_seurat_unintegrated@meta.data))) vars_sct <- c(vars_sct, "S.Score", "G2M.Score")
      vars_sct <- unique(vars_sct);
      if (length(vars_sct) == 0) vars_sct <- NULL else log_message(paste("Regressing:", paste(vars_sct, collapse = ", ")))
      merged_seurat_unintegrated <- SCTransform(merged_seurat_unintegrated, assay = "RNA", new.assay.name = "SCT_temp", vars.to.regress = vars_sct, verbose = FALSE)
      merged_seurat_unintegrated[["SCT"]] <- merged_seurat_unintegrated[["SCT_temp"]];
      merged_seurat_unintegrated[["SCT_temp"]] <- NULL;
      DefaultAssay(merged_seurat_unintegrated) <- "RNA"

    } else {
      # Multiple datasets
      log_message("Merging multiple objects...")
      merged_seurat_unintegrated <- tryCatch({ merge(x = seurat_list_renamed[[1]], y = seurat_list_renamed[-1], project = "MergedUnintegrated", merge.data = FALSE) }, error = function(e) { log_message(paste("ERROR merging:", e$message), type = "ERROR"); return(NULL) })
      if (is.null(merged_seurat_unintegrated)) return(NULL)

      # Reconstruct RNA
      log_message("Reconstructing RNA assay...");
      all_rna_counts <- list()
      for (i in seq_along(seurat_list_normalized)) { counts <- GetAssayData(seurat_list_normalized[[i]], "RNA", "counts"); colnames(counts) <- colnames(seurat_list_renamed[[i]]); all_rna_counts[[project_names[i]]] <- counts }
      combined_counts <- Seurat:::RowMergeSparseMatrices(all_rna_counts);
      rm(all_rna_counts);
      gc()
      if (!all(colnames(combined_counts) %in% colnames(merged_seurat_unintegrated))) combined_counts <- combined_counts[, colnames(merged_seurat_unintegrated)]
      merged_seurat_unintegrated[["RNA"]] <- CreateAssayObject(counts = combined_counts);
      merged_seurat_unintegrated <- NormalizeData(merged_seurat_unintegrated, assay = "RNA", verbose = FALSE);
      rm(combined_counts);
      gc()

      # Cell Cycle (Merged)
      log_message("Cell cycle scoring (merged)...")
      tryCatch({ DefaultAssay(merged_seurat_unintegrated) <- "RNA"; if (tolower(species) == "human") { s.genes <- cc.genes.updated.2019$s.genes; g2m.genes <- cc.genes.updated.2019$g2m.genes } else { stop("Custom species CC genes needed.") }; merged_seurat_unintegrated <- CellCycleScoring(merged_seurat_unintegrated, s.features = s.genes, g2m.features = g2m.genes, assay = 'RNA', set.ident = FALSE) }, error = function(e) { log_message(paste("WARN: CC scores failed:", e$message), type = "WARN") })

      # Create SCT_Ind
      log_message("Creating SCT_Ind assay...");
      all_sct_data <- list()
      for (i in seq_along(seurat_list_renamed)) { sct_assay <- GetAssay(seurat_list_renamed[[i]], "SCT"); all_sct_data[[project_names[i]]] <- GetAssayData(sct_assay, "data") }
      combined_sct_ind_data <- Seurat:::RowMergeSparseMatrices(all_sct_data);
      rm(all_sct_data);
      gc()
      if (!all(colnames(combined_sct_ind_data) %in% colnames(merged_seurat_unintegrated))) combined_sct_ind_data <- combined_sct_ind_data[, colnames(merged_seurat_unintegrated)]
      merged_seurat_unintegrated[["SCT_Ind"]] <- CreateAssayObject(data = combined_sct_ind_data);
      rm(combined_sct_ind_data);
      gc()

      # Run SCT on Merged RNA
      log_message("Running SCTransform on merged RNA...")
      DefaultAssay(merged_seurat_unintegrated) <- "RNA"
      vars_to_regress <- c();
      if (isTRUE(regress_params$regress_mito) && "percent_mito" %in% colnames(merged_seurat_unintegrated@meta.data)) vars_to_regress <- c(vars_to_regress, "percent_mito");
      # Add ribo/cc checks...
      if (isTRUE(regress_params$regress_ribo) && "percent_ribo" %in% colnames(merged_seurat_unintegrated@meta.data)) vars_to_regress <- c(vars_to_regress, "percent_ribo")
      if (isTRUE(regress_params$regress_cell_cycle) && all(c("S.Score", "G2M.Score") %in% colnames(merged_seurat_unintegrated@meta.data))) vars_to_regress <- c(vars_to_regress, "S.Score", "G2M.Score")
      vars_to_regress <- unique(vars_to_regress);
      if (length(vars_to_regress) == 0) vars_to_regress <- NULL else log_message(paste("Regressing:", paste(vars_to_regress, collapse = ", ")))
      merged_seurat_unintegrated <- SCTransform(merged_seurat_unintegrated, assay = "RNA", new.assay.name = "SCT", vars.to.regress = vars_to_regress, verbose = FALSE) # Verbose FALSE
      DefaultAssay(merged_seurat_unintegrated) <- "RNA"
    }
    # End multi-object block

    log_message("Saving merged unintegrated object...")
    tryCatch({ saveRDS(merged_seurat_unintegrated, file = merged_unintegrated_path) }, error = function(e) { log_message(paste("Error saving merged obj:", e$message), type = "ERROR") })
  }
  # End merged object creation block

  log_message("Normalization and preparation complete.")
  return(list(seurat_list_normalized = seurat_list_normalized, merged_seurat_unintegrated = merged_seurat_unintegrated))
}

#' Extract Expression Data Matrices for PH Analysis
#' @param seurat_list_normalized List of individual SCT objects (for sctInd, raw).
#' @param merged_seurat_object Merged Seurat object (for sctWhole or integrated).
#' @param data_types Character vector (e.g., "raw", "sctInd", "sctWhole", "integrated").
#' @param output_dir Directory to save lists.
#' @return List named by data_types containing expression lists, or NULL.
#' @noRd
extract_expression_data <- function(seurat_list_normalized = NULL, merged_seurat_object = NULL, data_types = c("raw", "sctInd", "sctWhole", "integrated"), output_dir) {
  log_message(paste("Extracting expression data types:", paste(data_types, collapse = ", ")))
  results <- list();
  overall_success <- TRUE
  extract_and_save <- function(expr_list, label) {
    # Internal helper
    success_flag <- FALSE
    if (!is.null(expr_list) && length(expr_list) > 0 && !all(sapply(expr_list, is.null))) {
      expr_list <- expr_list[!sapply(expr_list, is.null)]
      if (length(expr_list) == 0) { log_message(paste("WARN: List empty after NULL removal:", label), type = "WARN"); return(FALSE) }
      file_path <- file.path(output_dir, paste0("expr_list_", label, ".rds"))
      tryCatch({
        log_message(paste("Saving", label, "list (", length(expr_list), "datasets):", basename(file_path)));
        saveRDS(expr_list, file = file_path);
        results[[label]] <<- expr_list;
        success_flag <- TRUE
      }, error = function(e) { log_message(paste("ERROR saving", label, "list:", e$message), type = "ERROR"); success_flag <- FALSE })
    } else { log_message(paste("WARN: Skipping empty/NULL list for", label), type = "WARN"); success_flag <- FALSE }
    return(success_flag)
  }

  if ("sctInd" %in% data_types) {
    if (!is.null(seurat_list_normalized)) {
      log_message("Extracting sctInd...");
      expr_list <- tryCatch({ lapply(seurat_list_normalized, function(obj) if ("SCT" %in% Assays(obj)) as.matrix(GetAssayData(obj, "SCT", "data")) else NULL) }, error = function(e) { log_message(paste("ERROR extracting sctInd:", e$message), type = "ERROR"); NULL });
      if (!extract_and_save(expr_list, "sctInd")) overall_success <- FALSE
    } else { log_message("WARN: Cannot extract 'sctInd': list missing.", type = "WARN"); if ("sctInd" %in% data_types) overall_success <- FALSE }
  }
  if ("raw" %in% data_types) {
    if (!is.null(seurat_list_normalized)) {
      log_message("Extracting raw...");
      expr_list <- tryCatch({ lapply(seurat_list_normalized, function(obj) if ("RNA" %in% Assays(obj)) as.matrix(GetAssayData(obj, "RNA", "counts")) else NULL) }, error = function(e) { log_message(paste("ERROR extracting raw:", e$message), type = "ERROR"); NULL });
      if (!extract_and_save(expr_list, "raw")) overall_success <- FALSE
    } else { log_message("WARN: Cannot extract 'raw': list missing.", type = "WARN"); if ("raw" %in% data_types) overall_success <- FALSE }
  }

  split_targets <- intersect(data_types, c("sctWhole", "integrated"))
  if (length(split_targets) > 0) {
    if (!is.null(merged_seurat_object)) {
      split_col <- if ("orig.ident" %in% colnames(merged_seurat_object@meta.data)) "orig.ident" else "Sample";
      if (!split_col %in% colnames(merged_seurat_object@meta.data)) { log_message(paste("ERROR: Split column", split_col, "missing."), type = "ERROR"); overall_success <- FALSE } else {
        merged_object_list <- tryCatch({ log_message(paste("Splitting merged object by", split_col, "...")); SplitObject(merged_seurat_object, split.by = split_col) }, error = function(e) { log_message(paste("ERROR splitting:", e$message), type = "ERROR"); NULL })
        if (!is.null(merged_object_list)) {
          for (target_type in split_targets) {
            assay_name <- ifelse(target_type == "integrated", "integrated", "SCT");
            data_label <- target_type
            if (!assay_name %in% Assays(merged_seurat_object)) { log_message(paste("WARN: Cannot extract", data_label, "- Assay", assay_name, "missing."), type = "WARN"); overall_success <- FALSE; next }
            log_message(paste("Extracting", assay_name, "data for", data_label, "..."))
            expr_list <- tryCatch({ lapply(merged_object_list, function(obj) as.matrix(GetAssayData(obj, assay = assay_name, slot = "data"))) }, error = function(e) { log_message(paste("ERROR extracting split data for", data_label, ":", e$message), type = "ERROR"); NULL })
            if (!extract_and_save(expr_list, data_label)) overall_success <- FALSE
          }
          rm(merged_object_list);
          gc()
        } else { overall_success <- FALSE }
      }
    } else { log_message("WARN: Cannot extract 'sctWhole' or 'integrated': merged object missing.", type = "WARN"); overall_success <- FALSE }
  }

  if (length(results) == 0 && length(data_types) > 0) { log_message("ERROR: Failed to extract any data.", type = "ERROR"); return(NULL) }
  if (!overall_success && length(results) > 0) { log_message("WARN: Extraction finished, some types failed.", type = "WARN") } else if (overall_success) { log_message("Expression data extraction finished.") }
  return(results)
}

#' Run Persistent Homology Calculation with Monitoring
#' @param expr_list Named list of expression matrices.
#' @param data_label Label for dataset (e.g., "raw", "sctInd", "integrated").
#' @param ph_dim Dimension for PH.
#' @param ph_threshold Threshold for PH (-1 for auto).
#' @param num_cores Max cores for PH calculation.
#' @param output_dir Directory for logs and results.
#' @param memory_threshold Memory threshold for monitoring.
#' @param timeout_datasets Datasets to potentially skip initially.
#' @return List containing PD_list, thresholds, status, log_file, results_file, or NULL.
#' @noRd
run_ph_analysis <- function(expr_list, data_label, ph_dim, ph_threshold, num_cores, output_dir, memory_threshold = 0.25, timeout_datasets = NULL) {
  log_message(paste("Starting PH calculation for:", data_label))
  if (is.null(expr_list) || length(expr_list) == 0) { log_message("WARN: List empty.", type = "WARN"); return(NULL) }
  if (!requireNamespace("TDA", quietly = TRUE) || !requireNamespace("ripserr", quietly = TRUE)) { log_message("ERROR: TDA/ripserr missing.", type = "ERROR"); return(NULL) }
  if (!exists("process_expression_list_with_monitoring", mode = "function")) { log_message("ERROR: process_expression_list_with_monitoring not found.", type = "ERROR"); return(NULL) }

  file_suffix <- paste0("_dim", ph_dim, "_th", ph_threshold, "_", data_label)
  log_file <- file.path(output_dir, paste0("progress_log", file_suffix, ".csv"))
  results_file <- file.path(output_dir, paste0("intermediate_results", file_suffix, ".rds"))
  final_pd_file <- file.path(output_dir, paste0("PD_list", file_suffix, ".rds"))
  final_thresh_file <- file.path(output_dir, paste0("thresholds", file_suffix, ".rds"))

  if (file.exists(final_pd_file) && file.exists(final_thresh_file)) {
    log_message(paste("Loading existing PH results:", data_label));
    pd_list <- readRDS(final_pd_file);
    thresholds <- readRDS(final_thresh_file)
    status_log <- NULL;
    if (file.exists(log_file)) tryCatch({ status_log <- read.csv(log_file) }, error = function(e) { })
    if (is.null(status_log)) status_log <- data.frame(dataset = names(pd_list), status = "completed (loaded)", stringsAsFactors = FALSE)
    return(list(PD_list = pd_list, thresholds = thresholds, status = status_log, log_file = log_file, results_file = results_file, data_label = data_label)) # Add data_label here
  }

  log_message(paste("Running PH monitoring function for", data_label))
  ph_result <- tryCatch({
    process_expression_list_with_monitoring(expr_list = expr_list, DIM = ph_dim, THRESHOLD = ph_threshold, log_message = log_message, max_cores = num_cores, memory_threshold = memory_threshold, log_file = log_file, results_file = results_file, timeout_datasets = timeout_datasets)
  }, error = function(e) { log_message(paste("ERROR during PH run:", e$message), type = "ERROR"); progress_log <- NULL; if (file.exists(log_file)) tryCatch(progress_log <- read.csv(log_file), error = function(e2) { }); return(list(PD_list = NULL, thresholds = NULL, status = progress_log, log_file = log_file, results_file = results_file, data_label = data_label)) }) # Add data_label here

  pd_list <- ph_result$PD_list;
  thresholds <- ph_result$thresholds;
  progress_log <- NULL
  if (file.exists(log_file)) tryCatch(progress_log <- read.csv(log_file), error = function(e2) { log_message(paste("WARN: Could not read progress log:", log_file), type = "WARN") }) else { log_message(paste("WARN: Log file not found:", log_file), type = "WARN") }

  if (is.null(pd_list) || length(pd_list) == 0) { log_message(paste("WARN: PH did not produce PD list for", data_label), type = "WARN"); return(list(PD_list = NULL, thresholds = NULL, status = progress_log, log_file = log_file, results_file = results_file, data_label = data_label)) }
  # Add data_label here

  if (is.null(names(pd_list)) || length(names(pd_list)) != length(names(expr_list))) {
    log_message(paste("WARN: PD list name/length mismatch for", data_label), type = "WARN")
    if (length(pd_list) == length(names(expr_list))) { names(pd_list) <- names(expr_list); if (!is.null(thresholds) && length(thresholds) == length(names(expr_list))) names(thresholds) <- names(expr_list); log_message("Assigned names based on order.") } else { log_message("ERROR: Cannot assign names.", type = "ERROR") }
  } else if (!is.null(thresholds) && (is.null(names(thresholds)) || length(names(thresholds)) != length(names(expr_list)))) { if (!is.null(thresholds) && length(thresholds) == length(names(expr_list))) names(thresholds) <- names(expr_list) }

  log_message(paste("Saving final PH results for", data_label))
  tryCatch({ saveRDS(object = pd_list, file = final_pd_file); saveRDS(object = thresholds, file = final_thresh_file) }, error = function(e) { log_message(paste("ERROR saving PH results:", e$message), type = "ERROR") })

  log_message(paste("PH calculation finished for:", data_label))
  return(list(PD_list = pd_list, thresholds = thresholds, status = progress_log, log_file = log_file, results_file = results_file, data_label = data_label)) # Add data_label here
}

#' Retry Failed PH Calculations and Update Results
#' @param expr_list Original expression list.
#' @param initial_ph_result List returned by run_ph_analysis.
#' @param data_label Label for the dataset.
#' @param ph_dim Dimension for PH.
#' @param ph_threshold Original PH Threshold.
#' @param num_cores Cores for retry.
#' @param output_dir Output directory.
#' @param retry_params List with time_limit, doubling_limit.
#' @param overwrite_original Boolean.
#' @return Updated PD list or original list on failure.
#' @noRd
retry_and_update_ph <- function(expr_list, initial_ph_result, data_label, ph_dim, ph_threshold, num_cores, output_dir, retry_params = list(time_limit = 12 * 3600, doubling_limit = 5), overwrite_original = FALSE) {
  log_message(paste("Checking retries for:", data_label))
  if (is.null(initial_ph_result) || is.null(initial_ph_result$log_file) || is.null(initial_ph_result$results_file)) { log_message("ERROR: Initial result invalid.", type = "ERROR"); return(initial_ph_result$PD_list %||% NULL) }
  progress_log_path <- initial_ph_result$log_file;
  results_file_path <- initial_ph_result$results_file;
  original_pd_list <- initial_ph_result$PD_list %||% list()
  if (!file.exists(progress_log_path)) { log_message(paste("WARN: Log not found:", basename(progress_log_path)), type = "WARN"); return(original_pd_list) }
  progress_data <- tryCatch(read.csv(progress_log_path), error = function(e) { log_message(paste("ERROR reading log", basename(progress_log_path)), type = "ERROR"); NULL })
  if (is.null(progress_data)) return(original_pd_list)
  failed_datasets <- progress_data$dataset[progress_data$status != "completed"]
  if (length(failed_datasets) == 0) { log_message("No failures found."); return(original_pd_list) }
  log_message(paste("Retrying", length(failed_datasets), "failed datasets:", paste(failed_datasets, collapse = ", ")))
  expr_list_failed <- expr_list[names(expr_list) %in% failed_datasets]
  if (length(expr_list_failed) == 0) { log_message("ERROR: Name mismatch, cannot subset failed list.", type = "ERROR"); return(original_pd_list) }
  if (!exists("retry_pd_calculation", mode = "function")) { log_message("ERROR: retry_pd_calculation not found.", type = "ERROR"); return(original_pd_list) }

  file_suffix <- paste0("_dim", ph_dim, "_th", ph_threshold, "_", data_label)
  retry_log_file <- file.path(output_dir, paste0("retry_progress_log", file_suffix, ".csv"))
  retry_results_file <- file.path(output_dir, paste0("retry_results", file_suffix, ".rds"))
  log_message("Calling retry_pd_calculation...")
  retry_result <- tryCatch({
    retry_pd_calculation(progress_log = progress_log_path, results_file = results_file_path, expr_list = expr_list_failed, DIM = ph_dim, doubling_limit = retry_params$doubling_limit, time_limit = retry_params$time_limit, num_cores = num_cores, log_message = log_message, retry_progress_log = retry_log_file, retry_results_file = retry_results_file, overwrite_original = overwrite_original)
  }, error = function(e) { log_message(paste("ERROR during retry execution:", e$message), type = "ERROR"); return(NULL) })

  if (is.null(retry_result) || is.null(retry_result$PD_list) || length(retry_result$PD_list) == 0) { log_message("Retry failed or returned no results."); final_pd_file_final <- file.path(output_dir, paste0("PD_list", file_suffix, "_final.rds")); if (!file.exists(final_pd_file_final)) tryCatch({ saveRDS(original_pd_list, final_pd_file_final) }, error = function(e) { }); return(original_pd_list) }

  log_message("Merging retry results...");
  updated_pd_list <- original_pd_list;
  retried_pd_list <- retry_result$PD_list
  successful_retries <- 0
  for (dataset_name in names(retried_pd_list)) { if (!is.null(retried_pd_list[[dataset_name]])) { if (dataset_name %in% failed_datasets) { log_message(paste("Successfully retried:", dataset_name)); updated_pd_list[[dataset_name]] <- retried_pd_list[[dataset_name]]; successful_retries <- successful_retries + 1 } else { log_message(paste("WARN: Retry returned non-failed:", dataset_name), type = "WARN") }} else { log_message(paste("WARN: Retry yielded NULL for", dataset_name), type = "WARN") }}
  still_failed <- setdiff(failed_datasets, names(updated_pd_list[sapply(updated_pd_list[failed_datasets], function(x)!is.null(x))]))
  if (length(still_failed) > 0) log_message(paste("WARN:", length(still_failed), "still failed after retry:", paste(still_failed, collapse = ", ")), type = "WARN")
  log_message(paste("Retries finished.", successful_retries, "updated."))

  final_pd_file_final <- file.path(output_dir, paste0("PD_list", file_suffix, "_final.rds"))
  log_message(paste("Saving final updated PD list:", basename(final_pd_file_final)))
  tryCatch({ saveRDS(object = updated_pd_list, file = final_pd_file_final) }, error = function(e) { log_message(paste("ERROR saving final PD list:", e$message), type = "ERROR") })
  # Add threshold merging/saving here if needed
  return(updated_pd_list)
}