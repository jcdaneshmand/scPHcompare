# PH Functions Script
# This script is intended to be sourced and used in the PH pipeline. It includes functions for processing
# persistence homology calculations, managing memory, handling parallel processing, and dynamically adjusting thresholds.

# Load necessary libraries
# General Utility Functions
#----------------------------------------------

# Function to load an RData file
# Input: fileName (string) - Path to the RData file
# Output: The loaded data from the RData file
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"]) # Return the loaded data, excluding the fileName variable itself
}

# Function to log messages to the console with flushing
log_message <- function(message) {
  cat(message, "\n")
  flush.console()
}

# Function to update the progress log file
# Now includes the best threshold used for the dataset
# Input: log_file - Path to the log file
#        job_id - The ID of the job being processed
#        status - Status of the job (e.g., "completed", "deferred", "error")
#        threshold - The best threshold used for the dataset
# Modified update_progress_log function to create log file if it does not exist
update_progress_log <- function(log_file, job_id, status, threshold = NA) {
  log_entry <- data.frame(
    job_id = job_id,
    status = status,
    best_threshold = threshold,
    timestamp = Sys.time(),
    stringsAsFactors = FALSE
  )

  tryCatch({
    write.table(log_entry, file = log_file, sep = ",", append = TRUE, col.names = !file.exists(log_file), row.names = FALSE)
    log_message(paste("Logged progress for job", job_id, "with threshold:", threshold, "Status:", status))
  },
    error = function(e) {
      log_message(paste("Failed to log progress for job", job_id, ":", e$message))
    }
  )
}

# Function to get available memory using ps package
get_available_memory <- function() {
  mem_info <- ps::ps_system_memory()
  # Available memory in bytes
  memfree_bytes <- mem_info$free
  return(memfree_bytes)
}

# Function to save intermediate results for persistence diagrams
# Input: results_file - Path to the file where results should be saved
#        PD_list - List of persistence diagrams
#        job_id - The ID of the current job
#        PD - The persistence diagram for the current job
# Function to save only the best intermediate results for persistence diagrams
# Function to save intermediate results for persistence diagrams
# Input: results_file - Path to the file where results should be saved
#        job_id - The ID of the current job
#        best_PD - The persistence diagram for the current job
# Modified save_intermediate_results function with explicit results_file parameter
save_intermediate_results <- function(results_file, job_id, best_PD) {
  # Load existing intermediate results
  if (file.exists(results_file)) {
    PD_list <- readRDS(results_file)
  } else {
    PD_list <- list()
  }

  # Ensure PD_list is a named list
  if (is.null(names(PD_list))) {
    names(PD_list) <- as.character(seq_along(PD_list))
  }

  # Assign the PD to the correct job_id slot
  PD_list[[as.character(job_id)]] <- best_PD

  # Save the updated list back to the results file
  saveRDS(PD_list, results_file)
}


# Modified load_intermediate_results_and_identify_unfinished_jobs function
load_intermediate_results_and_identify_unfinished_jobs <- function(results_file, log_file, expr_list, log_message) {
  PD_list <- list() # Initialize the PD list
  job_indices <- seq_along(expr_list) # Get indices for all jobs

  log_message("Starting to load intermediate results and identify unfinished jobs.")

  # Step 1: Load the progress log and filter completed jobs
  if (file.exists(log_file)) {
    log_message(paste("Loading progress log from:", log_file))
    log_data <- tryCatch({
      read.csv(log_file, stringsAsFactors = FALSE)
    }, error = function(e) {
      log_message(paste("Error reading log file:", log_file, "-", e$message))
      data.frame(job_id = integer(0), status = character(0), best_threshold = numeric(0), timestamp = as.POSIXct(character(0)))
    })
  } else {
    log_message(paste("Progress log file", log_file, "does not exist. Initializing for a fresh run."))
    log_data <- data.frame(job_id = integer(0), status = character(0), best_threshold = numeric(0), timestamp = as.POSIXct(character(0)))
    write.csv(log_data, log_file, row.names = FALSE)
  }

  # Ensure that completed_jobs is properly initialized
  if (nrow(log_data) > 0 && "job_id" %in% colnames(log_data) && "status" %in% colnames(log_data)) {
    log_message(paste(nrow(log_data), "entries found in progress log."))
    completed_jobs <- log_data$job_id[log_data$status == "completed"] # Filter jobs marked as completed
    log_message(paste(length(completed_jobs), "jobs are marked as completed in the progress log."))
  } else {
    log_message("No valid entries found in the progress log. This is normal for a fresh run or when no previous jobs were completed.")
    completed_jobs <- integer(0) # Initialize as empty vector if no valid entries found
  }

  job_indices <- setdiff(job_indices, completed_jobs) # Only keep jobs that are not marked as completed

  # Step 2: Load intermediate results if needed, but only for jobs that are incomplete
  if (file.exists(results_file)) {
    log_message(paste("Loading intermediate results from", results_file))
    PD_list <- tryCatch({
      readRDS(results_file) # Load intermediate results
    }, error = function(e) {
      log_message(paste("Error reading intermediate results file:", results_file, "-", e$message))
      vector("list", length(expr_list)) # Return an empty list of the same length as expr_list
    })
  } else {
    log_message("No intermediate results file found. Initializing an empty list for persistence diagrams.")
    PD_list <- vector("list", length(expr_list)) # Initialize an empty list if no results file is found
    saveRDS(PD_list, results_file) # Save the initialized empty list
  }

  # Ensure PD_list has names matching expr_list
  if (length(expr_list) > 0) {
    names(PD_list) <- names(expr_list)
  } else if (is.null(names(PD_list))) {
    names(PD_list) <- as.character(seq_along(PD_list))
  }

  # Identify jobs with incomplete or invalid PD results
  incomplete_jobs <- which(sapply(PD_list, is.null) | sapply(PD_list, function(x) is.matrix(x) && nrow(x) == 0))
  job_indices <- intersect(job_indices, incomplete_jobs) # Focus only on jobs that still need processing

  log_message(paste(length(job_indices), "jobs need further processing."))

  return(list(PD_list = PD_list, job_indices = job_indices))
}


#----------------------------------------------
# Persistence Diagram (PD) Functions
#----------------------------------------------
# Estimate the memory usage of a persistence diagram (PD)
# Input: pd - Persistence diagram (data frame)
# Output: Memory usage in MB
estimate_data_size_pd <- function(pd) {
  object.size(pd) / (1024 ^ 2) # Convert bytes to MB
}

#----------------------------------------------
# Batch Size & Core Allocation for Memory Management
#----------------------------------------------

# Function to calculate batch size and allocate cores based on available memory
# Input: pd_list - List of persistence diagrams
#        safety_factor - Fraction of available memory to use (default: 0.8)
#        max_cores - Maximum number of cores to use (default: 12)
# Output: A list containing the calculated batch size and the number of cores to use
calculate_bdm_batch_size_and_cores <- function(pd_list, safety_factor = 0.8, max_cores = 12) {
  available_memory <- get_available_memory() # Get available memory in MB
  log_message(paste("Available memory for BDM:", available_memory, "MB"))

  # Estimate memory size for each PD and find the largest
  pd_sizes <- sapply(pd_list, estimate_data_size_pd)
  max_pd_size <- max(pd_sizes)

  # Estimate memory required for bottleneck distance calculations
  total_memory_per_pd <- max_pd_size * 2
  max_pd_in_memory <- floor((available_memory * safety_factor) / total_memory_per_pd)

  # Set the number of cores and batch size based on memory constraints
  num_cores <- min(max_pd_in_memory, max_cores)
  batch_size <- num_cores

  log_message(paste("Max PD size:", max_pd_size, "MB"))
  log_message(paste("Using", num_cores, "cores and batch size of", batch_size, "for BDM calculations."))

  return(list(batch_size = batch_size, num_cores = num_cores))
}

#----------------------------------------------
# Bottleneck Distance Matrix (BDM) Creation in Parallel
#----------------------------------------------

# Function to create a Bottleneck Distance Matrix (BDM) in parallel with memory-aware batching
# Input: PD - List of persistence diagrams
#        DIM - Dimension for bottleneck distance (default: 0)
#        dataset_name - Name of the dataset to make temporary directories unique
#        safety_factor - Fraction of available memory to use (default: 0.8)
#        max_cores - Maximum number of cores to use (default: 12)
#        total_memory_threshold - Total memory threshold in MB (default: 1.3TB)
#        log_message - Logging function
# Output: Symmetric Bottleneck Distance Matrix
CreateBottleneckDistanceMatrixParallel <- function(
    PD, DIM = 1, dataset_name, safety_factor = 0.8, max_cores = 12,
    total_memory_threshold = 1300000, max_time_per_calculation = Inf,
    log_message, save_progress = TRUE, progress_file = "distanceMatrix_progress.rds") {
  n <- length(PD) # Number of persistence diagrams
  log_message(paste0("Creating Output Matrix of Size ", n, "x", n))

  # Initialize or load the distance matrix
  if (save_progress && file.exists(progress_file)) {
    # Load existing distance matrix
    distanceMatrix <- readRDS(progress_file)
    log_message(paste("Loaded existing distance matrix from", progress_file))
  } else {
    # Initialize an empty distance matrix with NA values
    distanceMatrix <- matrix(NA, nrow = n, ncol = n)
  }

  # Create a unique temporary directory for storing persistence diagrams (PDs) on disk
  temp_dir <- file.path("temp_PD_files", dataset_name)
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
    log_message(paste("Created temporary directory for PD files:", temp_dir))
  }

  # Save each PD to disk as an RDS file (if not already saved)
  PD_file_paths <- sapply(seq_along(PD), function(i) {
    file_path <- file.path(temp_dir, paste0("PD_", i, ".rds"))
    if (!file.exists(file_path)) {
      saveRDS(PD[[i]], file_path)
    }
    file_path
  })
  log_message("All persistence diagrams saved to disk for processx access.")

  # Dynamically calculate batch size and cores
  batch_info <- calculate_bdm_batch_size_and_cores(
    PD,
    safety_factor = safety_factor, max_cores = max_cores
  )
  batch_size <- batch_info$batch_size
  num_cores <- batch_info$num_cores

  log_message(paste("Using batch size:", batch_size, "and number of cores:", num_cores))

  # Prepare a list of all (i, j) pairs to compute (upper triangle)
  pair_indices <- which(upper.tri(distanceMatrix), arr.ind = TRUE)

  # Filter out pairs that have already been computed
  incomplete_pairs <- which(is.na(distanceMatrix[upper.tri(distanceMatrix)]))
  if (length(incomplete_pairs) == 0) {
    log_message("All pairs have been computed. Nothing to do.")
    return(distanceMatrix)
  }

  pair_indices <- pair_indices[incomplete_pairs,, drop = FALSE]
  total_pairs <- nrow(pair_indices)

  # Process pairs in chunks to limit the number of concurrent processes
  pair_chunks <- split(seq_len(nrow(pair_indices)), ceiling(seq_along(seq_len(nrow(pair_indices))) / num_cores))

  total_chunks <- length(pair_chunks)
  log_message(paste("Total number of chunks to process:", total_chunks))

  for (chunk_idx in seq_along(pair_chunks)) {
    chunk <- pair_indices[pair_chunks[[chunk_idx]],, drop = FALSE]
    log_message(paste("Processing chunk", chunk_idx, "of", total_chunks, "with", nrow(chunk), "pairs..."))

    # Start processes for the current chunk
    job_list <- list()
    for (k in seq_len(nrow(chunk))) {
      i <- chunk[k, 1]
      j <- chunk[k, 2]

      if (!is.na(distanceMatrix[i, j])) {
        log_message(paste("Skipping pair (", i, ",", j, ") as it has already been calculated."))
        next # Skip already calculated pairs
      }

      # Start the bottleneck distance calculation for each pair using `processx`
      PD1_path <- PD_file_paths[i]
      PD2_path <- PD_file_paths[j]

      # Create a temporary result file path
      result_file <- file.path(temp_dir, paste0("BDM_result_", i, "_", j, ".txt"))

      # Create the Rscript command using the saved PD file paths
      ph_job <- processx::process$new(
        "Rscript",
        c("-e", paste0("\n          PD1 <- readRDS('", PD1_path, "');\n          PD2 <- readRDS('", PD2_path, "');\n          bottleneck_distance <- TDA::bottleneck(PD1, PD2, dimension = ", DIM, ");\n          writeLines(as.character(bottleneck_distance), '", result_file, "')\n        ")),
        stdout = "|", stderr = "|"
      )

      # Log process ID and store the job
      log_message(paste("Started bottleneck distance calculation for pair (", i, ",", j, ") with PID:", ph_job$get_pid()))
      job_list[[k]] <- list(process = ph_job, i = i, j = j, result_file = result_file)
    }

    # Wait for all processes in the chunk to complete
    for (job in job_list) {
      ph_job <- job$process
      i <- job$i
      j <- job$j
      result_file <- job$result_file

      start_time <- Sys.time()

      # Wait for the process to finish
      ph_job$wait()

      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      log_message(paste("Process for pair (", i, ",", j, ") with PID:", ph_job$get_pid(), "completed in", round(elapsed_time), "seconds"))

      # Collect the result from the temporary result file
      if (file.exists(result_file)) {
        distance <- as.numeric(readLines(result_file))
        distanceMatrix[i, j] <- distance
        distanceMatrix[j, i] <- distance # Mirror the result
        log_message(paste("Collected bottleneck distance for pair (", i, ",", j, "):", distance))

        # Remove the temporary result file
        file.remove(result_file)
      } else {
        log_message(paste("Result file for pair (", i, ",", j, ") not found."))
      }
    }

    # Save progress after each chunk
    if (save_progress) {
      saveRDS(distanceMatrix, progress_file)
      log_message(paste("Progress saved to", progress_file))
    }

    # Trigger garbage collection to free up memory
    gc()
  }

  # Clean up temporary PD files
  sapply(PD_file_paths, file.remove)
  log_message(paste("Temporary PD files removed from directory:", temp_dir))

  # Optionally remove the temporary directory
  unlink(temp_dir, recursive = TRUE)
  log_message(paste("Temporary directory removed:", temp_dir))

  # Final save of the completed distance matrix
  if (save_progress) {
    saveRDS(distanceMatrix, progress_file)
    log_message(paste("Final distance matrix saved to", progress_file))
  }

  return(distanceMatrix)
}
CreateSpectralDistanceMatrixFromPD <- function(
    pd_list,
    bdm_matrix,
    num_eigen = 50,
    scale_factor = 0.5,
    log_message = print,
    dimension = 0,
    save_progress = TRUE,
    progress_dir = "spectral_progress",
    dataset_name = "default",
    use_normalized = TRUE) {
  n <- length(pd_list)
  log_message(paste(
    "Creating Spectral Distance Matrix for", n,
    "persistence diagrams for dataset:", dataset_name
  ))

  # Ensure progress directory exists.
  if (save_progress && !dir.exists(progress_dir)) {
    dir.create(progress_dir, recursive = TRUE)
    log_message(paste("Created progress directory:", progress_dir))
  }

  # Unique progress file for the dataset.
  progress_file <- file.path(progress_dir, paste0("spectralMatrix_", dataset_name, "_progress.rds"))

  # Check for existing progress.
  if (save_progress && file.exists(progress_file)) {
    spectralMatrix <- readRDS(progress_file)
    log_message(paste("Loaded existing spectral distance matrix from", progress_file))
    return(spectralMatrix)
  } else {
    spectralMatrix <- matrix(NA, nrow = n, ncol = n)
  }

  # Validate the Bottleneck Distance Matrix.
  log_message("Validating Bottleneck Distance Matrix.")
  if (any(!is.finite(as.numeric(bdm_matrix)))) {
    log_message("WARNING: Found infinite or missing values in the Bottleneck Distance Matrix.")
    bdm_matrix[!is.finite(bdm_matrix)] <- max(bdm_matrix[is.finite(bdm_matrix)], na.rm = TRUE) * 1.5
  }

  # Compute sigma based on the range of the finite BDM values.
  bdm_finite <- bdm_matrix[is.finite(bdm_matrix)]
  bdm_range <- max(bdm_finite) - min(bdm_finite)
  if (bdm_range == 0) {
    log_message("WARNING: BDM range is zero. Using default sigma value of 1.")
    sigma <- 1
  } else {
    sigma <- bdm_range * scale_factor
  }
  log_message(paste(
    "Using sigma =", sigma,
    "based on BDM range =", round(bdm_range, 4),
    "with scale factor", scale_factor
  ))

  # Convert BDM to a similarity matrix using an RBF (Gaussian) kernel.
  log_message("Converting Bottleneck Distance Matrix to Similarity Matrix.")
  similarity_matrix <- exp(-bdm_matrix ^ 2 / (2 * sigma ^ 2))
  diag(similarity_matrix) <- 1 # Ensure self-similarity is exactly 1.

  # Compute the Laplacian.
  if (use_normalized) {
    log_message("Using normalized Laplacian.")
    # Compute the normalized Laplacian: L_sym = I - D^{-1/2} * W * D^{-1/2}
    D <- diag(rowSums(similarity_matrix))
    D_inv_sqrt <- diag(1 / sqrt(diag(D)))
    L <- diag(n) - D_inv_sqrt %*% similarity_matrix %*% D_inv_sqrt
  } else {
    log_message("Using unnormalized Laplacian.")
    D <- diag(rowSums(similarity_matrix))
    L <- D - similarity_matrix
  }

  # Perform spectral decomposition.
  log_message("Performing spectral decomposition.")
  eigen_decomp <- tryCatch(
    eigen(L, symmetric = TRUE),
    error = function(e) {
      log_message(paste("ERROR: Spectral decomposition failed:", e$message))
      return(NULL)
    }
  )
  if (is.null(eigen_decomp)) {
    return(NULL)
  }

  # Order eigenvalues (and corresponding eigenvectors) in increasing order.
  order_index <- order(eigen_decomp$values)
  eigenvalues <- eigen_decomp$values[order_index]
  eigenvectors <- eigen_decomp$vectors[, order_index]

  # Skip the trivial eigenvector (typically the constant one corresponding to eigenvalue 0).
  if (ncol(eigenvectors) < 2) {
    log_message("Not enough eigenvectors to skip the trivial eigenvector. Using available eigenvectors.")
    spectral_embedding <- eigenvectors
  } else {
    num_eigen_used <- min(num_eigen, ncol(eigenvectors) - 1)
    spectral_embedding <- eigenvectors[, 2:(num_eigen_used + 1)]
  }

  # Compute the pairwise Euclidean distances in the spectral embedding space.
  log_message("Computing Spectral Distance Matrix.")
  spectralMatrix <- as.matrix(dist(spectral_embedding))

  # Save progress if requested.
  if (save_progress) {
    saveRDS(spectralMatrix, progress_file)
    log_message(paste("Spectral distance matrix saved to", progress_file))
  }

  return(spectralMatrix)
}


#----------------------------------------------
# Processing Datasets with Tiered Thresholds for Complex Datasets
#----------------------------------------------

# Function to process a list of datasets with tiered thresholds for complex datasets
# Input: expr_list - List of expression datasets
#        DIM - Dimension for persistence homology
#        complex_datasets - List of datasets considered complex
#        log_message - Logging function
#        max_cores - Maximum number of cores for parallel processing
#        total_memory_threshold - Total memory threshold in MB
#        log_file - File path for logging progress
#        results_file - File path for saving intermediate results
#        THRESHOLD - Default threshold to use for non-complex datasets
# Output: List of persistence diagrams
# Modified function to process the expression list and track the best threshold for each dataset
# Function to process a list of datasets with tiered thresholds for complex datasets
process_expression_list_with_monitoring <- function(expr_list, DIM, log_message, max_cores, memory_threshold,
                                                    log_file, results_file, batch_size = 64, timeout_datasets = NULL) {
  # Load intermediate results and identify unfinished jobs initially
  load_res <- load_intermediate_results_and_identify_unfinished_jobs(results_file, log_file, expr_list, log_message)
  PD_list <- load_res$PD_list
  job_indices <- load_res$job_indices
  threshold_tracker <- list()

  # Ensure PD_list is a named list matching expr_list
  if (length(expr_list) > 0) {
    names(PD_list) <- names(expr_list)
  } else if (is.null(names(PD_list))) {
    names(PD_list) <- as.character(seq_along(PD_list))
  }

  log_message(paste("Processing datasets in batches of", batch_size, "with", max_cores, "cores."))
  batches <- split(job_indices, ceiling(seq_along(job_indices) / batch_size)) # Split job indices into batches

  # Process each batch of datasets
  for (batch in batches) {
    log_message(paste("Processing batch:", paste(batch, collapse = ", ")))
    batch_results <- process_batch_of_datasets(batch, expr_list, DIM, log_message, max_cores, memory_threshold,
                                               results_file, PD_list, log_file, timeout_datasets
    )

    # Update PD_list and threshold_tracker based on the results from process_batch_of_datasets
    for (i in seq_along(batch)) {
      job_id <- batch[i]
      if (!is.null(batch_results[[i]]$PD)) {
        PD_list[[as.character(job_id)]] <- batch_results[[i]]$PD
        threshold_tracker[[as.character(job_id)]] <- batch_results[[i]]$threshold
      }
      # No need to save or log here since it's handled within process_and_monitor
    }
  }

  # Reload intermediate results to capture the final state after processing all batches
  PD_list <- readRDS(results_file) # Ensure the latest version of PD_list is loaded

  # Extract final thresholds from the progress log file
  log_message("Extracting thresholds from the progress log file.")
  threshold_tracker <- extract_thresholds_from_log(log_file)

  return(list(PD_list = PD_list, thresholds = threshold_tracker))
}


# Function to extract thresholds from the progress log
extract_thresholds_from_log <- function(log_file) {
  # Check if the log file exists
  if (!file.exists(log_file)) {
    log_message(paste("Log file", log_file, "does not exist. Returning an empty threshold tracker."))
    return(list()) # Return an empty list if the log file is not found
  }

  # Read the log file into a data frame
  log_data <- tryCatch({
    read.csv(log_file, stringsAsFactors = FALSE)
  },
    error = function(e) {
      log_message(paste("Error reading log file:", log_file, "-", e$message))
      return(data.frame(job_id = integer(0), best_threshold = numeric(0), status = character(0), timestamp = as.POSIXct(character(0))))
    }
  )

  # Check if the expected columns are present in the log
  if (!all(c("job_id", "best_threshold") %in% colnames(log_data))) {
    stop("Log file does not contain expected columns: 'job_id', 'best_threshold'")
  }

  # Filter to only completed datasets and extract thresholds
  completed_entries <- log_data[log_data$status == "completed", c("job_id", "best_threshold")]

  # Create a named list to store thresholds
  threshold_tracker <- as.list(setNames(completed_entries$best_threshold, completed_entries$job_id))
  return(threshold_tracker)
}

# Update the process_batch_of_datasets function to include a check for completed jobs before reprocessing
process_batch_of_datasets <- function(batch_indices, expr_list, DIM, log_message, max_cores, memory_threshold,
                                      results_file, PD_list, log_file, timeout_datasets = NULL) {

  log_message(paste("Starting parallel processing for batch with", length(batch_indices), "datasets."))

  # Parallel processing for each dataset in the batch
  results <- mclapply(seq_along(batch_indices), function(i) {
    global_job_id <- batch_indices[i] # Get the global job index

    # Check if the job is already completed using the new is_job_completed function
    if (is_job_completed(global_job_id, PD_list, log_file)) {
      log_message(paste("Skipping dataset", global_job_id, "as it is already completed."))
      return(list(PD = PD_list[[as.character(global_job_id)]], threshold = NA)) # Skip and return already completed PD
    }

    log_message(paste("Processing dataset", global_job_id))
    # Process the dataset and monitor
    result <- tryCatch({
      process_and_monitor(
          expr_matrix = expr_list[[global_job_id]],
          i = global_job_id,
          DIM = DIM,
          log_message = log_message,
          memory_threshold = memory_threshold,
          timeout_datasets = timeout_datasets,
          results_file = results_file, # Ensure results_file is passed here
          log_file = log_file
        )
    },
      error = function(e) {
        log_message(paste("Error in processing dataset", global_job_id, ":", e$message))
        return(list(PD = NULL, threshold = NULL)) # Return NULL if there is an error
      }
    )

    # Return the result from process_and_monitor
    return(list(PD = result$PD, threshold = result$threshold))
  }, mc.cores = min(length(batch_indices), max_cores)) # Parallel processing

  # Return the processed results for the current batch
  return(results)
}

# Function to verify if a job is completed based on both the log file and the intermediate results
# Input: job_id - Unique identifier for the job/dataset
#        PD_list - List of persistence diagrams loaded from intermediate results
#        log_file - Path to the progress log file
# Output: TRUE if the job is completed, FALSE otherwise
is_job_completed <- function(job_id, PD_list, log_file) {
  # Check if the job is completed in the PD list
  if (!is.null(PD_list[[as.character(job_id)]]) && nrow(PD_list[[as.character(job_id)]]) > 0) {
    return(TRUE)
  }

  # Check if the job is completed in the progress log
  if (file.exists(log_file)) {
    log_data <- read.csv(log_file, stringsAsFactors = FALSE)
    completed_jobs <- log_data$job_id[log_data$status == "completed"]
    if (job_id %in% completed_jobs) {
      return(TRUE)
    }
  }

  return(FALSE)
}

#----------------------------------------------
# Dataset Processing Helper Functions
#----------------------------------------------

#' process_and_monitor
#'
#' Run persistent homology (PH) on a single expression matrix using
#' `ripserr::vietoris_rips` from the **ripserr** package. The function chooses an
#' initial threshold of `-1` unless the dataset previously timed out, in
#' which case a data-driven estimate is used. The PH computation runs in a
#' separate process and is killed if it exceeds `max_time_per_iteration`.
#'
#' @param expr_matrix A numeric matrix or `dgCMatrix` of expression values.
#' @param i Identifier of the dataset.
#' @param DIM Dimension for PH computation.
#' @param log_message Function for logging messages.
#' @param memory_threshold Reserved for future memory checks.
#' @param max_time_per_iteration Maximum runtime in seconds (default 12 hours).
#' @param timeout_datasets Vector of dataset indices that previously timed out.
#'
#' @return A list with
#' \describe{
#'   \item{PD}{The persistence diagram matrix or `NULL` if no valid result was produced.}
#'   \item{threshold}{The threshold used for the PH run, or `NULL` on failure.}
#' }
# Function to process and monitor a single dataset for Persistent Homology (PH) calculation
process_and_monitor <- function(expr_matrix, i, DIM, log_message, memory_threshold,
                                max_time_per_iteration = 20 * 24 * 3600, # 10*24 hours
                                timeout_datasets = NULL,
                                results_file = NULL,
                                log_file = NULL,
                                temp_dataset_dir = "temp_datasets_sct_whole") {
  tryCatch({
    is_timeout <- i %in% timeout_datasets
    if (is_timeout) {
      # log_message(paste("Dataset", i, "is marked as timed out. Setting initial threshold based on maximum of non-zero values."))
      # max_val <- calculate_max_threshold(expr_matrix)
      # threshold_increment <- determine_proportional_increment(max_val, proportion = 0.1)
      # current_threshold <- max_val / 4
      # log_message(paste("Initial threshold for dataset", i, "set to:", current_threshold))


      # -------------------------------------------------------------
      # 1. run a quick PCA on expr_matrix
      pcs <- prcomp(expr_matrix, center = TRUE, scale. = TRUE)$x[, seq_len(DIM)]

      # 2. choose k for the kNN (e.g. 10)
      k_nn <- 100

      # 3. compute the k‐th nearest neighbor distances
      knn_info <- FNN::get.knn(pcs, k = k_nn)

      # 4. extract each point’s distance to its k-th neighbor
      kth_dists <- knn_info$nn.dist[, k_nn]

      # after getting your k-th neighbor distances
      kth_dists_nonzero <- kth_dists[kth_dists > 0]
      if (length(kth_dists_nonzero) == 0) {
        stop("All k-NN distances are zero—check for duplicates or increase DIM/k.")
      }
      # 5. pick your threshold as the median (or another percentile)
      ideal_thresh <- quantile(kth_dists_nonzero, 0.9)

      log_message(paste("Auto‐selected threshold (k =", k_nn, "90th Percentile):", round(ideal_thresh, 3)))

      # 6. decide on your increment strategy
      threshold_increment <- ideal_thresh * 0.1 # e.g. 10% of that local scale

      # and feed that into your PH routine:
      current_threshold <- ideal_thresh
      # -------------------------------------------------------------

    } else {
      current_threshold <- -1 # Infinite threshold
      threshold_increment <- 50
      log_message(paste("Dataset", i, "has an initial threshold of -1 (infinite)."))
    }

    best_PD <- NULL

    log_message(paste("Starting process for dataset", i, "with threshold", current_threshold))

    if (!dir.exists(temp_dataset_dir)) {
      dir.create(temp_dataset_dir, recursive = TRUE)
      log_message(paste("Created temporary directory for datasets:", temp_dataset_dir))
    }
    dataset_file <- file.path(temp_dataset_dir, paste0("dataset_", i, ".rds"))
    saveRDS(expr_matrix, dataset_file)
    log_message(paste("Dataset", i, "saved to", dataset_file))

    pd_file <- file.path(temp_dataset_dir, paste0("PD_", i, "_", current_threshold, ".rds"))

    ph_job <- processx::process$new(
      "Rscript",
      c("-e", paste0(
        "tryCatch({",
        "dataset <- readRDS('", dataset_file, "');",
        "dataset <- as.matrix(dataset);",
        "PD <- ripserr::vietoris_rips(dataset = dataset, max_dim = ", DIM, ", threshold = ", current_threshold, ", return_format = 'mat');",
        "saveRDS(PD, '", pd_file, "')",
        "}, error = function(e) { cat('Error:', e$message, '\n') })",
        ";quit(save='no')"
      )),
      stdout = "|", stderr = "|"
    )

    start_time <- Sys.time()
    log_message(paste("PH process for dataset", i, "started with PID:", ph_job$get_pid()))

    while (ph_job$is_alive()) {
      Sys.sleep(60)
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      log_message(paste("Elapsed time for dataset", i, "PH calculation:", round(elapsed_time), "seconds"))

      if (elapsed_time > max_time_per_iteration) {
        log_message(paste("Timeout exceeded for dataset", i, "at threshold:", current_threshold, "- killing process."))
        ph_job$kill()
        break
      }
    }

    if (file.exists(pd_file)) {
      best_PD <- readRDS(pd_file)
    }

    if (!is.null(best_PD) && nrow(best_PD) > 0) {
      save_intermediate_results(results_file = results_file, job_id = i, best_PD = best_PD)
      update_progress_log(log_file, i, "completed", current_threshold)
      return(list(PD = best_PD, threshold = current_threshold))
    } else {
      update_progress_log(log_file, i, "failed", current_threshold)
      return(list(PD = NULL, threshold = NULL))
    }

  }, error = function(e) {
    log_message(paste("Unexpected error in processing dataset", i, ":", e$message))
    update_progress_log(log_file, i, "error", threshold = NA)
    return(list(PD = NULL, threshold = NULL))
  })
}

# Helper function to calculate max threshold
calculate_max_threshold <- function(expr_matrix) {
  # Return the maximum of the non-zero values in the expression matrix
  max_val <- max(expr_matrix[expr_matrix > 0], na.rm = TRUE)
  return(max_val)
}



# Modified get_anchors function to handle missing anchor files
get_anchors <- function(seurat_list, integration_features, dims, anchor_file) {
  if (file.exists(anchor_file)) {
    log_message(paste("Loading precomputed anchor batch from:", anchor_file))
    anchors <- readRDS(anchor_file)
  } else {
    log_message(paste("Anchor batch file", anchor_file, "does not exist. Initializing a new anchor batch."))
    anchors <- FindIntegrationAnchors(
      object.list = seurat_list,
      normalization.method = "SCT",
      anchor.features = integration_features,
      dims = 1:dims,
      verbose = FALSE
    )
    saveRDS(anchors, file = anchor_file)
    log_message(paste("Saved new anchor batch to:", anchor_file))
  }
  return(anchors)
}
