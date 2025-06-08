# --- 0. Setup: Load Libraries and Source External Functions ---

# Recommended: Use pacman for easier package management
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(tidyverse, Matrix, ripserr, TDA, foreach, doParallel, parallel,
#                Seurat, SeuratObject, SeuratDisk, stringr, pheatmap, mclust, aricode,
#                clusterSim, Rtsne, batchelor, BiocSingular, scCustomize, kernlab,
#                igraph, progressr, plyr, digest, transport, kernlab, svglite,
#                ComplexHeatmap, circlize, viridis, dendextend) # Added more deps potentially used by moved funcs

packages <- c(
  "tidyverse", "Matrix", "ripserr", "TDA", "foreach", "doParallel", "parallel",
  "Seurat", "SeuratObject", "SeuratDisk", "stringr", "pheatmap", "mclust", "aricode",
  "clusterSim", "Rtsne", "batchelor", "BiocSingular", "scCustomize", "kernlab",
  "igraph", "progressr", "plyr", "digest", "transport", "kernlab", "svglite",
  "ComplexHeatmap", "circlize", "viridis", "dendextend", "ps", "processx", "TDAstats" # Added all potential deps
)

# Check and load required packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
  library(pkg, character.only = TRUE)
}

# Source helper functions defined elsewhere
source("PH_Functions.R") # Assumes cleaned version
source("Integration_flexible.R") # Seurat-only version

# Define %||% operator if not using rlang
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Helper Function Definitions ---

# (Keep setup_environment, load_metadata, load_and_create_seurat,
#  add_qc_and_metadata, save_metadata_summary, filter_seurat_objects,
#  normalize_and_prepare_seurat, extract_expression_data, run_ph_analysis,
#  retry_and_update_ph from previous version)
# ...

# --- NEW: Add Matrix Calculation Functions Here (Moved from Post-Processing) ---

#' Create Bottleneck Distance Matrix using Parallel External Processes
#' (Copied and cleaned from previous post-processing script version)
#' @noRd # Mark as internal helper
CreateBottleneckDistanceMatrixParallel <- function(
    PD, DIM = 1, dataset_name, safety_factor = 0.8, max_cores = 12,
    log_message, save_progress = TRUE, progress_file = "BDM_progress.rds") {
  # Use specific name

  n <- length(PD)
  if (n < 2) { log_message("Need >= 2 PDs for BDM.", type = "WARN"); return(matrix(0, n, n)) }
  log_message(sprintf("Creating BDM (%dx%d) Dim %d for %s", n, n, DIM, dataset_name))

  # --- Initialize or Load Progress ---
  distanceMatrix <- matrix(NA_real_, nrow = n, ncol = n);
  diag(distanceMatrix) <- 0
  if (save_progress && file.exists(progress_file)) {
    log_message(paste("Loading existing BDM progress:", basename(progress_file)))
    loaded_matrix <- tryCatch(readRDS(progress_file), error = function(e) { log_message("ERROR loading BDM progress.", type = "ERROR"); NULL })
    if (!is.null(loaded_matrix) && is.matrix(loaded_matrix) && all(dim(loaded_matrix) == c(n, n))) {
      distanceMatrix <- loaded_matrix
    } else { log_message(paste("WARN: BDM progress file invalid.", basename(progress_file)), type = "WARN") }
  }

  # --- Setup Temporary Directory ---
  session_id <- paste0(dataset_name, "_BDM_", Sys.getpid())
  temp_dir <- file.path(tempdir(), "PD_temp_files", session_id)
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
  on.exit({ if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE) }, add = TRUE)

  # --- Save PDs to Disk ---
  log_message("Saving PDs to temp disk storage...")
  PD_file_paths <- character(n)
  pd_names <- names(PD) %||% as.character(1:n) # Use names if available
  names(PD_file_paths) <- pd_names
  for (i in seq_along(PD)) {
    pd_id <- pd_names[i]
    file_path <- file.path(temp_dir, paste0("PD_", gsub("[^A-Za-z0-9_.-]", "_", pd_id), ".rds"))
    tryCatch(saveRDS(PD[[i]], file_path), error = function(e) { stop("Failed to save temp PD file: ", e$message) })
    PD_file_paths[i] <- file_path
  }
  log_message("Temporary PDs saved.")

  # --- Determine Batching/Cores ---
  # Ensure calculate_bdm_batch_size_and_cores is available (should be in PH_Functions.R)
  if (!exists("calculate_bdm_batch_size_and_cores", mode = "function")) stop("Helper calculate_bdm_batch_size_and_cores not found.")
  batch_info <- calculate_bdm_batch_size_and_cores(PD, safety_factor = safety_factor, max_cores = max_cores, log_message = log_message)
  num_cores_to_use <- batch_info$num_cores

  # --- Identify Pairs to Compute ---
  all_pairs_indices <- which(upper.tri(distanceMatrix), arr.ind = TRUE)
  na_indices_linear <- which(is.na(distanceMatrix[upper.tri(distanceMatrix)]))
  if (length(na_indices_linear) == 0) { log_message("BDM already complete."); return(distanceMatrix) }
  pairs_to_compute <- all_pairs_indices[na_indices_linear,, drop = FALSE]
  total_pairs_to_compute <- nrow(pairs_to_compute)
  log_message(sprintf("Need to compute BDM distances for %d pairs.", total_pairs_to_compute))

  # --- Process Pairs in Chunks ---
  chunk_size <- num_cores_to_use * 2
  pair_chunks <- split(seq_len(total_pairs_to_compute), ceiling(seq_len(total_pairs_to_compute) / chunk_size))
  total_chunks <- length(pair_chunks)
  log_message(paste("Processing BDM in", total_chunks, "chunks using up to", num_cores_to_use, "cores."))

  results_list <- list()

  for (chunk_idx in seq_along(pair_chunks)) {
    current_pair_indices_in_chunk <- pair_chunks[[chunk_idx]]
    chunk_pairs <- pairs_to_compute[current_pair_indices_in_chunk,, drop = FALSE]
    n_pairs_in_chunk <- nrow(chunk_pairs)
    log_message(sprintf("BDM Chunk %d/%d (%d pairs)...", chunk_idx, total_chunks, n_pairs_in_chunk))

    job_list <- vector("list", n_pairs_in_chunk)
    procs_launched_map <- list() # Use list to track PIDs

    for (k in seq_len(n_pairs_in_chunk)) {
      # Simple concurrency limit: wait if max cores busy
      while (length(procs_launched_map) >= num_cores_to_use) {
        Sys.sleep(1)
        pids_alive <- names(procs_launched_map)[sapply(procs_launched_map, function(p) p$is_alive())]
        procs_launched_map <- procs_launched_map[pids_alive] # Keep only alive ones
      }

      i <- chunk_pairs[k, 1];
      j <- chunk_pairs[k, 2]
      result_file <- file.path(temp_dir, paste0("bdm_result_", i, "_", j, ".rds"))
      r_command <- sprintf("tryCatch({library(TDA); pd1<-readRDS('%s'); pd2<-readRDS('%s'); d<-TDA::bottleneck(pd1,pd2,dimension=%d); saveRDS(d,'%s'); quit(save='no',status=0);}, error=function(e){cat('ERROR:',conditionMessage(e),'\\n'); saveRDS(NA_real_,'%s'); quit(save='no',status=1);})", PD_file_paths[i], PD_file_paths[j], DIM, result_file, result_file)

      px <- tryCatch(processx::process$new("Rscript", c("-e", r_command), stdout = "|", stderr = "|"), error = function(e) { log_message(paste("ERROR launching processx:", e$message), type = "ERROR"); NULL })
      if (!is.null(px)) {
        job_list[[k]] <- list(process = px, i = i, j = j, result_file = result_file)
        procs_launched_map[[as.character(px$get_pid())]] <- px # Track launched proc
      } else { job_list[[k]] <- NULL }
    }

    # --- Wait for and Collect Results ---
    log_message("Collecting results for BDM chunk...")
    for (k in seq_along(job_list)) {
      job_info <- job_list[[k]];
      if (is.null(job_info)) next
      px <- job_info$process;
      i <- job_info$i;
      j <- job_info$j;
      result_file <- job_info$result_file
      job_label <- sprintf("BDM Pair (%d, %d)", i, j)
      tryCatch(px$wait(timeout = 3 * 3600 * 1000), error = function(e) { log_message(paste("WARN:", job_label, "wait error/timeout:", e$message), type = "WARN"); if (px$is_alive()) try(px$kill(), silent = TRUE) })
      exit_status <- px$get_exit_status();
      dist_val <- NA_real_
      if (!is.na(exit_status) && exit_status != 0) log_message(paste("WARN:", job_label, "non-zero exit:", exit_status), type = "WARN")
      if (file.exists(result_file)) { res <- tryCatch(readRDS(result_file), error = function(e) NA_real_); if (is.numeric(res) && length(res) == 1) dist_val <- res; try(file.remove(result_file), silent = TRUE) } else { log_message(paste("WARN:", job_label, "result file missing."), type = "WARN") }
      results_list[[paste0(i, "_", j)]] <- dist_val
    }

    # --- Update Distance Matrix ---
    log_message("Updating BDM matrix...")
    updated_count = 0
    for (name in names(results_list)) {
      coords <- as.integer(strsplit(name, "_")[[1]]);
      i <- coords[1];
      j <- coords[2]
      dist_val <- results_list[[name]]
      if (!is.na(dist_val) && is.na(distanceMatrix[i, j])) { distanceMatrix[i, j] <- distanceMatrix[j, i] <- dist_val; updated_count <- updated_count + 1 }
      else if (is.na(dist_val)) { log_message(paste("WARN: NA BDM result for pair", i, j), type = "WARN") }
    }
    results_list <- list();
    log_message(paste("Updated", updated_count, "BDM entries."))

    # --- Save Progress ---
    if (save_progress) { tryCatch(saveRDS(distanceMatrix, progress_file), error = function(e) { log_message("ERROR saving BDM progress.", type = "ERROR") }); log_message(paste("BDM Progress saved after chunk", chunk_idx)) }
    clean_memory()
  }
  # End loop over chunks

  num_na <- sum(is.na(distanceMatrix[upper.tri(distanceMatrix)]));
  if (num_na > 0) log_message(paste("WARN: BDM finished with", num_na, "NA pairs."), type = "WARN") else log_message("BDM computation complete.")
  if (save_progress) { tryCatch(saveRDS(distanceMatrix, progress_file), error = function(e) { }); log_message(paste("Final BDM saved:", basename(progress_file))) }
  return(distanceMatrix)
}


#' Create Spectral Distance Matrix from Persistence Diagrams
#' (Copied and cleaned from previous post-processing script version)
#' @noRd # Mark as internal helper
CreateSpectralDistanceMatrixFromPD <- function(
    pd_list, bdm_matrix, num_eigen = 50, scale_factor = 0.5, log_message,
    dimension = 1, save_progress = TRUE, progress_dir = "spectral_progress",
    dataset_name = "default", use_normalized = TRUE) {
  # Default Dim 1

  n <- length(pd_list)
  if (n < 2) { log_message("Need >= 2 PDs.", type = "WARN"); return(matrix(0, n, n)) }
  if (!is.matrix(bdm_matrix) || !all(dim(bdm_matrix) == c(n, n))) { log_message("Invalid BDM matrix.", type = "ERROR"); return(NULL) }
  log_message(sprintf("Creating SDM (%dx%d) from BDM (Dim %d) for: %s", n, n, dimension, dataset_name))

  progress_file <- file.path(progress_dir, paste0("spectralMatrix_", dataset_name, "_dim", dimension, "_progress.rds"))
  if (save_progress) {
    ensure_dir(progress_dir)
    if (file.exists(progress_file)) {
      spectralMatrix <- tryCatch(readRDS(progress_file), error = function(e) NULL)
      if (!is.null(spectralMatrix) && is.matrix(spectralMatrix) && all(dim(spectralMatrix) == c(n, n))) { log_message("Loaded existing SDM."); return(spectralMatrix) } else { log_message("WARN: Existing SDM file invalid.", type = "WARN") }
    }
  }

  log_message("Validating BDM & computing similarity...")
  bdm_matrix_finite <- bdm_matrix;
  is_infinite <- !is.finite(bdm_matrix_finite)
  if (any(is_infinite)) { log_message("WARN: Replacing non-finite BDM values.", type = "WARN"); max_finite_val <- max(bdm_matrix_finite[is.finite(bdm_matrix_finite)], na.rm = TRUE) %||% 1; bdm_matrix_finite[is_infinite] <- max_finite_val * 1.5 }
  bdm_range <- diff(range(bdm_matrix_finite[is.finite(bdm_matrix_finite)], na.rm = TRUE));
  sigma <- if (bdm_range > 1e-9) bdm_range * scale_factor else 1.0
  log_message(sprintf("Using sigma = %.4f", sigma))
  similarity_matrix <- exp(-bdm_matrix_finite ^ 2 / (2 * sigma ^ 2));
  diag(similarity_matrix) <- 1;
  similarity_matrix <- (similarity_matrix + t(similarity_matrix)) / 2

  log_message(paste("Computing Laplacian:", ifelse(use_normalized, "normalized", "unnormalized")))
  D_vec <- Matrix::rowSums(similarity_matrix);
  if (any(D_vec <= 1e-9)) { D_vec[D_vec <= 1e-9] <- 1e-9 }
  if (use_normalized) { L <- Matrix::Diagonal(n) - Matrix::Diagonal(n, x = 1 / sqrt(D_vec)) %*% similarity_matrix %*% Matrix::Diagonal(n, x = 1 / sqrt(D_vec)) } else { L <- Matrix::Diagonal(n, x = D_vec) - similarity_matrix }
  L <- (L + t(L)) / 2

  log_message("Performing spectral decomposition...")
  num_eigen_requested <- min(num_eigen, n - 1);
  num_eigen_requested <- max(2, num_eigen_requested)
  eigen_decomp <- NULL
  if (requireNamespace("RSpectra", quietly = TRUE)) { eigen_decomp <- tryCatch(RSpectra::eigs_sym(L, k = num_eigen_requested + 1, which = "SM"), error = function(e) { log_message(paste("WARN: RSpectra failed:", e$message), type = "WARN"); NULL }) }
  if (is.null(eigen_decomp)) { eigen_decomp <- tryCatch(eigen(as.matrix(L), symmetric = TRUE), error = function(e) { log_message(paste("ERROR: Eigen decomp failed:", e$message), type = "ERROR"); NULL }) }
  if (is.null(eigen_decomp)) return(NULL)

  eigen_order <- order(eigen_decomp$values);
  eigenvalues_sorted <- eigen_decomp$values[eigen_order];
  eigenvectors_sorted <- eigen_decomp$vectors[, eigen_order, drop = FALSE]
  if (ncol(eigenvectors_sorted) <= 1) { log_message("ERROR: Not enough eigenvectors.", type = "ERROR"); return(NULL) }
  num_eigen_to_use <- min(num_eigen_requested, ncol(eigenvectors_sorted) - 1)
  spectral_embedding <- eigenvectors_sorted[, 2:(num_eigen_to_use + 1), drop = FALSE]

  log_message(paste("Computing SDM using", num_eigen_to_use, "eigenvectors."))
  spectralMatrix <- tryCatch(as.matrix(dist(spectral_embedding)), error = function(e) { log_message("ERROR computing dist:", type = "ERROR"); NULL })
  if (is.null(spectralMatrix)) return(NULL)

  if (save_progress) { tryCatch(saveRDS(spectralMatrix, progress_file), error = function(e) { log_message("ERROR saving SDM:", type = "ERROR") }); log_message(paste("SDM saved:", basename(progress_file))) }
  return(spectralMatrix)
}

#' Compute Persistence Landscapes
#' @param pd Persistence diagram matrix.
#' @param grid Sequence of tau values. Default seq(0, 1, length.out = 100).
#' @return List with dim0 and dim1 landscapes, or NULL on error.
#' @noRd # Mark as internal helper
ComputePersistenceLandscapes <- function(pd) {
  # The 'grid' argument is removed
  library("TDA")

  res <- tryCatch({
    # Let TDA choose the optimal grid by NOT providing the tseq argument.
    # We now compute the first 5 landscape functions by setting KK = 1:5.
    # The output will be a 5-row matrix.
    landscape0 <- TDA::landscape(Diag = pd, dimension = 0, KK = 1:5)
    landscape1 <- TDA::landscape(Diag = pd, dimension = 1, KK = 1:5)

    list(dim0 = landscape0, dim1 = landscape1)
  },
    error = function(e) {
      log_message(paste("Error computing persistence landscape for a PD:", e$message))
      return(NULL)
    }
  )
  return(res)
}

#' Compute and Save Landscape List and Combined L2 Distance Matrix
#' @param pd_list Named list of persistence diagrams.
#' @param landscape_list_path Path to save the list of landscapes (RDS).
#' @param landscape_distance_matrix_path Path to save the L2 distance matrix (RDS).
#' @param log_message Logging function.
#' @param num_cores Number of cores for parallel landscape computation.
#' @return Invisible NULL. Saves results to files.
#' @noRd # Mark as internal helper
compute_and_save_landscape_matrices <- function(pd_list, landscape_list_path, landscape_distance_matrix_path, log_message, num_cores = 6) {
  log_message("Computing persistence landscapes...")
  if (is.null(pd_list) || length(pd_list) == 0) { log_message("WARN: PD list empty, skipping landscape.", type = "WARN"); return(invisible(NULL)) }

  # Use mclapply for parallel computation
  landscape_list <- mclapply(pd_list, ComputePersistenceLandscapes, mc.cores = num_cores)

  # Basic validation
  if (length(landscape_list) != length(pd_list) || any(sapply(landscape_list, is.null))) {
    log_message("WARN: Some landscapes failed to compute.", type = "WARN")
    # Optional: filter out NULLs? Or proceed with NAs in distance matrix?
    # Let's proceed and handle NAs in distance calculation.
  }
  names(landscape_list) <- names(pd_list) # Ensure names are preserved

  # Save landscape list
  tryCatch(saveRDS(landscape_list, file = landscape_list_path), error = function(e) log_message("ERROR saving landscape list.", type = "ERROR"))
  log_message(paste("Saved landscape list to:", basename(landscape_list_path)))

  # Compute L2 distance matrix
  n <- length(landscape_list)
  combined_distance_matrix <- matrix(NA_real_, n, n) # Initialize with NA
  rownames(combined_distance_matrix) <- names(landscape_list)
  colnames(combined_distance_matrix) <- names(landscape_list)

  grid_len <- 100 # Assuming default grid length used in ComputePersistenceLandscapes

  for (i in 1:n) {
    combined_distance_matrix[i, i] <- 0 # Distance to self is 0
    if (i < n) {
      for (j in (i + 1):n) {
        land1 <- landscape_list[[i]]
        land2 <- landscape_list[[j]]
        # Check if landscapes are valid before computing distance
        if (!is.null(land1) && !is.null(land2) && !is.null(land1$dim0) && !is.null(land1$dim1) && !is.null(land2$dim0) && !is.null(land2$dim1)) {
          # Compute L2 distance for dim0 (use first level k=1)
          d0_1 <- if (nrow(land1$dim0) > 0) land1$dim0[1,] else rep(0, grid_len)
          d0_2 <- if (nrow(land2$dim0) > 0) land2$dim0[1,] else rep(0, grid_len)
          dist_dim0_sq <- sum((d0_1 - d0_2) ^ 2)
          # Compute L2 distance for dim1 (use first level k=1)
          d1_1 <- if (nrow(land1$dim1) > 0) land1$dim1[1,] else rep(0, grid_len)
          d1_2 <- if (nrow(land2$dim1) > 0) land2$dim1[1,] else rep(0, grid_len)
          dist_dim1_sq <- sum((d1_1 - d1_2) ^ 2)
          # Combined Euclidean distance in (L2(dim0), L2(dim1)) space
          combined_dist <- sqrt(dist_dim0_sq + dist_dim1_sq)
        } else {
          combined_dist <- NA_real_ # Set NA if any landscape is invalid
        }
        combined_distance_matrix[i, j] <- combined_distance_matrix[j, i] <- combined_dist
      }
    }
  }

  log_message("Saving combined landscape L2 distance matrix...")
  tryCatch({ saveRDS(combined_distance_matrix, file = landscape_distance_matrix_path); write.csv(combined_distance_matrix, file = sub(".Rds", ".csv", landscape_distance_matrix_path)) }, error = function(e) log_message("ERROR saving landscape matrix.", type = "ERROR"))
  log_message(paste("Saved landscape distance matrix:", basename(landscape_distance_matrix_path)))
  invisible(NULL)
}

#' Compute and Save Distance Matrices (BDM, SDM) for an Iteration
#' Wrapper function called after PH analysis for an iteration.
#' @noRd # Mark as internal helper
compute_and_save_distance_matrices <- function(
    ph_results, # Output from run_ph_analysis (contains PD_list, data_label etc)
    expr_list_names, # Names corresponding to PD list order
    output_dir,
    num_cores,
    log_message,
    recompute_bdm = TRUE, # Flags to force recomputation
    recompute_sdm = TRUE) {

  if (is.null(ph_results) || is.null(ph_results$PD_list) || length(ph_results$PD_list) == 0) {
    log_message("WARN: No valid PD list found, skipping distance matrix calculation.", type = "WARN")
    return(list(BDM = NULL, SDM = NULL))
  }

  pd_list <- ph_results$PD_list
  dataset_name_suffix <- ph_results$data_label # e.g., "raw", "sctInd", "integrated_seurat"
  ph_dim <- 1 # Assuming Dim 1 for now, parameterize if needed

  # Ensure PD list has correct names
  if (length(pd_list) != length(expr_list_names) || is.null(names(pd_list)) || !all(names(pd_list) == expr_list_names)) {
    if (length(pd_list) == length(expr_list_names)) {
      log_message("WARN: Assigning names to PD list based on expression list order.", type = "WARN")
      names(pd_list) <- expr_list_names
    } else {
      log_message("ERROR: PD list length mismatch, cannot assign names reliably. Skipping matrix calculation.", type = "ERROR")
      return(list(BDM = NULL, SDM = NULL))
    }
  }

  # Define file paths
  bdm_filename <- paste0("BDM_", dataset_name_suffix, ".rds")
  sdm_filename <- paste0("SDM_", dataset_name_suffix, ".rds")
  bdm_path <- file.path(output_dir, bdm_filename)
  sdm_path <- file.path(output_dir, sdm_filename)
  bdm_progress_file <- file.path(output_dir, paste0("BDM_progress_", dataset_name_suffix, ".rds"))
  sdm_progress_dir <- file.path(output_dir, "spectral_progress") # Directory for SDM progress

  # --- Calculate BDM ---
  BDM <- NULL
  bdm_exists <- file.exists(bdm_path)
  if (bdm_exists && !recompute_bdm) {
    log_message(paste("Loading existing BDM:", basename(bdm_path)))
    BDM <- tryCatch(readRDS(bdm_path), error = function(e) NULL)
    if (!is.null(BDM) && is.matrix(BDM) && all(dim(BDM) == length(expr_list_names))) {
      # Ensure names match current list
      if (!identical(rownames(BDM), expr_list_names)) {
        log_message("Updating names on loaded BDM.", type = "INFO")
        dimnames(BDM) <- list(expr_list_names, expr_list_names)
        # Optionally re-save with correct names
        tryCatch({ saveRDS(BDM, bdm_path); write.csv(as.data.frame(BDM), sub(".rds", ".csv", bdm_path)) }, error = function(e) { })
      }
    } else { BDM <- NULL; log_message("Existing BDM invalid, recomputing.", type = "WARN") }
  }
  if (is.null(BDM)) {
    log_message(paste("Computing BDM for", dataset_name_suffix))
    BDM <- CreateBottleneckDistanceMatrixParallel(
                 PD = pd_list, DIM = ph_dim, dataset_name = dataset_name_suffix,
                 max_cores = num_cores, log_message = log_message,
                 progress_file = bdm_progress_file, save_progress = TRUE # Save progress internally
             )
    if (!is.null(BDM)) {
      dimnames(BDM) <- list(expr_list_names, expr_list_names)
      tryCatch({ saveRDS(BDM, file = bdm_path); write.csv(as.data.frame(BDM), file = sub(".rds", ".csv", bdm_path)) }, error = function(e) { log_message(paste("ERROR saving BDM:", e$message), type = "ERROR") })
    } else { log_message(paste("ERROR: BDM computation failed for", dataset_name_suffix), type = "ERROR") }
  }

  # --- Calculate SDM ---
  SDM <- NULL
  sdm_exists <- file.exists(sdm_path)
  if (sdm_exists && !recompute_sdm) {
    log_message(paste("Loading existing SDM:", basename(sdm_path)))
    SDM <- tryCatch(readRDS(sdm_path), error = function(e) NULL)
    if (!is.null(SDM) && is.matrix(SDM) && all(dim(SDM) == length(expr_list_names))) {
      if (!identical(rownames(SDM), expr_list_names)) {
        log_message("Updating names on loaded SDM.", type = "INFO")
        dimnames(SDM) <- list(expr_list_names, expr_list_names)
        tryCatch({ saveRDS(SDM, sdm_path); write.csv(as.data.frame(SDM), sub(".rds", ".csv", sdm_path)) }, error = function(e) { })
      }
    } else { SDM <- NULL; log_message("Existing SDM invalid, recomputing.", type = "WARN") }
  }

  if (is.null(SDM)) {
    if (is.null(BDM)) {
      log_message("WARN: Skipping SDM calculation because BDM is missing.", type = "WARN")
    } else {
      log_message(paste("Computing SDM for", dataset_name_suffix))
      SDM <- CreateSpectralDistanceMatrixFromPD(
                      pd_list = pd_list, # Pass pd_list for count N
                      bdm_matrix = BDM,
                      num_eigen = 50, # Parameterize?
                      log_message = log_message,
                      dataset_name = dataset_name_suffix,
                      dimension = ph_dim, # Pass dimension
                      progress_dir = sdm_progress_dir, # Separate progress dir
                      save_progress = TRUE
                  )
      if (!is.null(SDM)) {
        dimnames(SDM) <- list(expr_list_names, expr_list_names)
        tryCatch({ saveRDS(SDM, file = sdm_path); write.csv(as.data.frame(SDM), file = sub(".rds", ".csv", sdm_path)) }, error = function(e) { log_message(paste("ERROR saving SDM:", e$message), type = "ERROR") })
      } else { log_message(paste("ERROR: SDM computation failed for", dataset_name_suffix), type = "ERROR") }
    }
  }

  return(list(BDM = BDM, SDM = SDM)) # Return computed/loaded matrices
}

# --- Main Pipeline Function Definition ---
#' Process Single-Cell Datasets with Persistent Homology (Seurat Integration Only)
#' @param metadata_path Path to metadata CSV.
#' @param output_dir Output directory.
#' @param metadata_sample_col Column name in metadata for sample IDs.
#' @param integration_method Integration method ("seurat" or "none").
#' @param num_cores Number of CPU cores.
#' @param min_cells Minimum cells per dataset after filtering.
#' @param ph_dim Dimension for PH.
#' @param ph_threshold Threshold for PH (-1 for auto).
#' @param filter_params List of filtering parameters.
#' @param sct_params List of SCTransform parameters.
#' @param regress_params List controlling regression variables.
#' @param species Species for cell cycle genes.
#' @param integration_params List of Seurat integration parameters.
#' @param ph_memory_threshold Memory threshold for PH monitoring.
#' @param perform_retries Whether to retry failed PH jobs.
#' @param ph_retry_params List of retry parameters.
#' @param compute_matrices Boolean: Whether to compute BDM/SDM/Landscapes in this run. Default TRUE.
#' @param run_downstream_analysis Whether to run downstream steps (placeholder).
#' @return List summarizing results and output paths.
#' @export
process_datasets_PH <- function(
    metadata_path, output_dir, metadata_sample_col = "Sample Name",
    integration_method = c("seurat", "none"),
    num_cores = parallel::detectCores() %/% 2,
    min_cells = 250,
    filter_params = list(minGenesPerCell = 500, maxGenesPerCell = 9000, maxMitoPercent = 20, minRiboPercent = 5, minCellsPerGene = 3),
    sct_params = list(variable.features.n = 3000),
    regress_params = list(regress_mito = TRUE, regress_ribo = FALSE, regress_cell_cycle = TRUE),
    species = "human",
    integration_params = list(dims = 30, k.anchor = 5, integration_features_n = 3000),
    ph_dim = 1, ph_threshold = -1, ph_memory_threshold = 0.25,
    perform_retries = TRUE, ph_retry_params = list(time_limit = 12 * 3600, doubling_limit = 5),
    compute_matrices = TRUE, # New flag to control matrix computation
    run_downstream_analysis = TRUE
    ) {

  # --- 1. Setup ---
  integration_method <- match.arg(integration_method)
  log_file <- setup_environment(output_dir)
  log_message("Starting PH Pipeline (Seurat Integration Only)")
  log_message(paste("Time:", Sys.time()))
  log_message(paste("Integration Method:", integration_method))

  # --- 2. Load Metadata ---
  metadata <- load_metadata(metadata_path);
  if (is.null(metadata)) { sink(); return(NULL) }

  # --- 3. Load/Create Seurat Objects ---
  seurat_list_initial <- load_and_create_seurat(metadata, num_cores)
  if (is.null(seurat_list_initial) || length(seurat_list_initial) == 0) { log_message("ERROR: Load/Create failed.", type = "ERROR"); sink(); return(NULL) }

  # --- 4. QC & Metadata ---
  seurat_list_qc <- add_qc_and_metadata(seurat_list_initial, metadata, metadata_sample_col, num_cores)
  if (is.null(seurat_list_qc)) { sink(); return(NULL) }
  rm(seurat_list_initial);
  gc()

  # --- 5. Pre-filter Summary ---
  save_metadata_summary(seurat_list_qc, file.path(output_dir, "prefiltering_summary.csv"), TRUE)

  # --- 6. Filtering ---
  seurat_list_filtered <- filter_seurat_objects(seurat_list_qc, filter_params, min_cells, num_cores)
  if (is.null(seurat_list_filtered) || length(seurat_list_filtered) == 0) { log_message("WARN: No datasets remain after filtering.", type = "WARN"); sink(); return(NULL) }
  rm(seurat_list_qc);
  gc()

  # --- 7. Post-filter Summary ---
  save_metadata_summary(seurat_list_filtered, file.path(output_dir, "postfiltering_summary.csv"), FALSE)

  # --- 8. Normalization & Preparation ---
  log_message("Running normalization & preparation...")
  normalization_result <- normalize_and_prepare_seurat(seurat_list = seurat_list_filtered, sct_params = sct_params, regress_params = regress_params, species = species, num_cores = num_cores, output_dir = output_dir)
  if (is.null(normalization_result)) { log_message("ERROR: Normalization failed.", type = "ERROR"); sink(); return(NULL) }
  seurat_list_normalized <- normalization_result$seurat_list_normalized # Keep this list needed for extraction
  merged_seurat_unintegrated <- normalization_result$merged_seurat_unintegrated
  rm(seurat_list_filtered);
  gc()

  # --- 9. Data Integration ---
  integrated_seurat_object <- NULL;
  integration_ran_successfully <- FALSE
  if (integration_method == "seurat") {
    log_message("Performing Seurat integration...")
    if (!is.null(seurat_list_normalized) && length(seurat_list_normalized) >= 2) {
      if (exists("perform_integration", mode = "function")) {
        # Define params for Seurat integration call
        seurat_output_path = file.path(output_dir, "final_integrated_seurat.rds")
        seurat_checkpoint_dir = file.path(output_dir, "integration_checkpoints")
        seurat_anchor_dir = file.path(output_dir, "anchor_batches")
        # Call perform_integration with specific params
        integrated_seurat_object <- perform_integration(
                  seurat_list = seurat_list_normalized,
                  integration_method = "seurat",
                  integration_params = integration_params,
                  output_path = seurat_output_path,
                  checkpoint_dir = seurat_checkpoint_dir,
                  anchor_batches_dir = seurat_anchor_dir,
                  num_cores = num_cores # Pass other relevant params...
              )
        if (!is.null(integrated_seurat_object)) { integration_ran_successfully <- TRUE; log_message("Seurat integration successful.") }
        else { log_message("WARN: Seurat integration failed.", type = "WARN") }
      } else { log_message("ERROR: perform_integration function missing.", type = "ERROR") }
    } else { log_message("WARN: Need >= 2 datasets for Seurat integration.", type = "WARN") }
  } else { log_message("Skipping integration (method='none').") }
  # Keep seurat_list_normalized until after extraction

  # --- 10. Extract Expression Data ---
  data_to_extract <- c("raw", "sctInd", "sctWhole")
  if (integration_ran_successfully) data_to_extract <- c(data_to_extract, "integrated")
  log_message(paste("Extracting expression data:", paste(data_to_extract, collapse = ", ")))
  extraction_results <- extract_expression_data(
      seurat_list_normalized = seurat_list_normalized, # NEEDED for raw, sctInd
      merged_seurat_object = if (integration_ran_successfully) integrated_seurat_object else merged_seurat_unintegrated,
      data_types = data_to_extract,
      output_dir = output_dir
  )
  if (is.null(extraction_results) || length(extraction_results) == 0) { log_message("ERROR: Failed extracting data.", type = "ERROR"); sink(); return(NULL) }
  rm(seurat_list_normalized);
  gc() # NOW safe to remove

  # --- 11. Run PH Analysis & Matrix Computations ---
  ph_results <- list()
  dist_matrix_results <- list() # To store computed BDM/SDM

  for (dtype in names(extraction_results)) {
    current_expr_list <- extraction_results[[dtype]]
    data_label <- dtype
    if (!is.null(current_expr_list) && length(current_expr_list) > 0) {
      log_message(paste("--- Processing PH & Matrices for:", data_label, "---"))
      ph_run_result <- run_ph_analysis(expr_list = current_expr_list, data_label = data_label, ph_dim = ph_dim, ph_threshold = ph_threshold, num_cores = num_cores, output_dir = output_dir, memory_threshold = ph_memory_threshold)
      ph_results[[data_label]] <- ph_run_result

      if (perform_retries && !is.null(ph_run_result) && !is.null(ph_run_result$status) && any(ph_run_result$status$status != "completed")) {
        updated_pd_list <- retry_and_update_ph(expr_list = current_expr_list, initial_ph_result = ph_run_result, data_label = data_label, ph_dim = ph_dim, ph_threshold = ph_threshold, num_cores = num_cores, output_dir = output_dir, retry_params = ph_retry_params)
        if (!is.null(updated_pd_list)) { ph_results[[data_label]]$PD_list <- updated_pd_list }
        # Update PD list in results
        rm(updated_pd_list);
        gc()
      }

      # --- Compute Distance Matrices if requested and PD list is valid ---
      if (compute_matrices && !is.null(ph_results[[data_label]]$PD_list) && length(ph_results[[data_label]]$PD_list) > 0) {
        log_message(paste("Computing distance matrices for:", data_label))
        # Get names in correct order
        current_expr_names <- names(current_expr_list) %||% as.character(seq_along(current_expr_list))

        # Compute BDM and SDM
        matrices <- compute_and_save_distance_matrices(
                   ph_results = ph_results[[data_label]], # Pass the result list containing PD_list etc.
                   expr_list_names = current_expr_names,
                   output_dir = output_dir,
                   num_cores = num_cores,
                   log_message = log_message,
                   recompute_bdm = TRUE, # Control recomputation if needed
                   recompute_sdm = TRUE
               )
        dist_matrix_results[[paste0(data_label, "_BDM")]] <- matrices$BDM
        dist_matrix_results[[paste0(data_label, "_SDM")]] <- matrices$SDM

        # Compute Landscapes and L2 Matrix
        landscape_list_path <- file.path(output_dir, paste0("landscape_list_", data_label, ".Rds"))
        landscape_matrix_path <- file.path(output_dir, paste0("landscape_l2_distance_matrix_", data_label, ".Rds"))
        pd_list_path_final <- file.path(output_dir, paste0("PD_list_dim", ph_dim, "_th", ph_threshold, "_", data_label, "_final.rds"))
        if (!file.exists(pd_list_path_final)) pd_list_path_final <- file.path(output_dir, paste0("PD_list_dim", ph_dim, "_th", ph_threshold, "_", data_label, ".rds")) # Fallback

        if (file.exists(pd_list_path_final)) {
          pd_list_final <- tryCatch(readRDS(pd_list_path_final),
                                    error = function(e) {
                                      log_message(paste("ERROR reading PD list:", e$message),
                                                  type = "ERROR");
                                      NULL
                                    })
          if (!is.null(pd_list_final) && length(pd_list_final) > 0) {
            compute_and_save_landscape_matrices(
              pd_list = pd_list_final,
              landscape_list_path = landscape_list_path,
              landscape_distance_matrix_path = landscape_matrix_path,
              log_message = log_message,
              num_cores = num_cores
            )
          } else {
            log_message(paste("WARN: Failed to load PD list for landscapes:",
                              basename(pd_list_path_final)), type = "WARN")
          }
        } else {
          log_message(paste("WARN: Final PD list not found for landscape calculation:",
                           basename(pd_list_path_final)), type = "WARN")
        }

      } else if (compute_matrices) {
        log_message(paste("Skipping matrix computation for", data_label, "- PD list invalid or empty."), type = "WARN")
      }
      rm(ph_run_result);
      gc()
    } else { log_message(paste("WARN: Skipping PH & Matrices for empty list:", data_label), type = "WARN"); ph_results[[data_label]] <- NULL }
  }
  rm(extraction_results);
  gc()


  # --- 12. Run Downstream Analysis ---
  if (run_downstream_analysis) {
    log_message("Proceeding to downstream analysis...")
    # Ensure objects are loaded if needed by downstream functions
    if (!exists("merged_seurat_unintegrated") || is.null(merged_seurat_unintegrated)) merged_seurat_unintegrated <- tryCatch(readRDS(file.path(output_dir, "merged_seurat_unintegrated.rds")), error = function(e) NULL)
    if (!exists("integrated_seurat_object") || is.null(integrated_seurat_object)) integrated_seurat_object <- if (integration_ran_successfully) tryCatch(readRDS(file.path(output_dir, "final_integrated_seurat.rds")), error = function(e) NULL) else NULL

    # Placeholder calls - these functions need the results (PH lists, matrices, Seurat objects)
    bdm_results_paths <- list.files(output_dir, pattern = "^BDM_.*\\.rds$", full.names = TRUE) # Example: get paths
    # Modify downstream function signatures as needed
    # clustering_results <- perform_clustering_and_visualization(merged_seurat_unintegrated, integrated_seurat_object, bdm_results_paths, metadata, output_dir, integration_method)
    # statistical_analysis_results <- run_statistical_tests(clustering_results, metadata, output_dir)
  } else { log_message("Skipping downstream analysis.") }


  # --- 13. Finalize and Return ---
  log_message("Pipeline finished.");
  log_message(paste("Time:", Sys.time()));
  sink()
  results_summary <- list(
      output_directory = output_dir, log_file = log_file,
      merged_unintegrated_object_path = file.path(output_dir, "merged_seurat_unintegrated.rds"),
      final_integrated_object_path = if (integration_ran_successfully) file.path(output_dir, "final_integrated_seurat.rds") else NULL,
      ph_results_summary = lapply(ph_results, function(res) { if (!is.null(res) && !is.null(res$data_label)) { file_suffix <- paste0("_dim", ph_dim, "_th", ph_threshold, "_", res$data_label); final_pd_path <- file.path(output_dir, paste0("PD_list", file_suffix, "_final.rds")); if (!file.exists(final_pd_path)) final_pd_path <- file.path(output_dir, paste0("PD_list", file_suffix, ".rds")); return(list(data_type = res$data_label, pd_list_path = if (file.exists(final_pd_path)) final_pd_path else NA_character_, status_log_path = res$log_file)) } else { return(NULL) }}),
      distance_matrix_paths = list.files(output_dir, pattern = "^(BDM|SDM|landscape_l2_distance_matrix)_.*\\.rds$", full.names = TRUE) # Add matrix paths
  )
  results_summary$ph_results_summary <- results_summary$ph_results_summary[!sapply(results_summary$ph_results_summary, is.null)]
  return(results_summary)

}

# --- End of process_datasets_PH function ---

# --- Example Usage ---
# Define parameters
# metadata_file <- "./path/to/your/metadata.csv"
# output_directory <- "./PH_Pipeline_Seurat_Output"
# cores <- 8

# Run the pipeline
# pipeline_results <- process_datasets_PH(
#     metadata_path = metadata_file,
#     output_dir = output_directory,
#     num_cores = cores,
#     integration_method = "seurat", # Or "none"
#     min_cells = 100 # Adjust as needed
# )
# print(pipeline_results)
