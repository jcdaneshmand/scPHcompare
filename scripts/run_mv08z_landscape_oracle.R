#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(paste(
  "usage: run_mv08z_landscape_oracle.R <prefreeze> <private-sentinel>",
  "<mv08s-private> <mv08v-private> <rust-library> <output-csv>",
  "<execution-head> <run-id>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
selection_path <- normalizePath(args[[2L]], mustWork = TRUE)
s_root <- normalizePath(args[[3L]], mustWork = TRUE)
v_root <- normalizePath(args[[4L]], mustWork = TRUE)
rust_library <- normalizePath(args[[5L]], mustWork = TRUE)
output <- normalizePath(args[[6L]], mustWork = FALSE)
execution_head <- tolower(args[[7L]])
run_id <- args[[8L]]
if (file.exists(output) || file.exists(paste0(output, ".partial")) ||
    !grepl("^[0-9a-f]{40}$", execution_head) || !nzchar(run_id)) {
  stop("MV8-Z oracle requires a fresh valid output", call. = FALSE)
}

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(prefreeze, "mv08z-artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(prefreeze, "mv08z-contract.csv"))
inputs <- .mv08z_read_csv(file.path(prefreeze, "mv08z-input-manifest.csv"))
selection <- .mv08z_read_csv(selection_path)
if (nrow(contract) != 1L || nrow(selection) != 1L ||
    .mv08z_sha256_file(selection_path) !=
      inputs$sha256[inputs$role == "private_sentinel_selection"] ||
    .mv08z_sha256_file(rust_library) != contract$rust_library_sha256) {
  stop("MV8-Z oracle binding drift", call. = FALSE)
}
root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("MV8-Z oracle unknown PH source role", call. = FALSE)
}
load_one <- function(prefix) {
  path <- file.path(root_for(selection[[paste0(prefix, "_source_role")]]),
                    selection[[paste0(prefix, "_output_file")]])
  if (!file.exists(path) ||
      as.numeric(file.info(path)$size) !=
        as.numeric(selection[[paste0(prefix, "_output_bytes")]]) ||
      .mv08z_sha256_file(path) !=
        selection[[paste0(prefix, "_output_sha256")]]) {
    stop("MV8-Z oracle PH artifact drift", call. = FALSE)
  }
  record <- readRDS(path)
  mv08s_validate_ph_record_v1(record)
  if (record$topology_result$provenance$diagram_sha256 !=
      selection[[paste0(prefix, "_diagram_sha256")]]) {
    stop("MV8-Z oracle diagram drift", call. = FALSE)
  }
  record
}
first_record <- load_one("first")
second_record <- load_one("second")
dimension <- as.integer(sub("H", "", selection$homology_dimension, fixed = TRUE))
first <- .mv08z_finite_intervals(first_record, selection$homology_dimension)
second <- .mv08z_finite_intervals(second_record, selection$homology_dimension)
reference_route <- if (max(nrow(first), nrow(second)) <= 500L) {
  "r_exact_breakpoint"
} else "r_adaptive_certified"
reference <- if (reference_route == "r_exact_breakpoint") {
  landscape_reference_exact_dimension(
    first_record$topology_result$diagram,
    second_record$topology_result$diagram,
    dimension, exact_max_intervals = 500L
  )
} else {
  landscape_reference_adaptive_dimension(
    first_record$topology_result$diagram,
    second_record$topology_result$diagram,
    dimension, abs_tol = 1e-8, rel_tol = 1e-8, subdivisions = 200L
  )
}
candidate <- landscape_rust_prototype_dimension(
  first, second, dimension, library = rust_library
)
threshold <- 1e-10 * max(1, abs(reference$squared_distance))
absolute_error <- abs(candidate$squared_distance - reference$squared_distance)
passed <- isTRUE(reference$within_requested_tolerance) &&
  isTRUE(candidate$rust_used) && candidate$status == 0L &&
  absolute_error <= threshold &&
  as.integer(candidate$active_levels) ==
    max(.mv08z_active_depth(first), .mv08z_active_depth(second))
result <- data.frame(
  contract_id = "mv08z_sentinel_oracle_v1", execution_head = execution_head,
  run_id = run_id, pair_identity_sha256 = selection$pair_identity_sha256,
  group_order = selection$group_order,
  homology_dimension = selection$homology_dimension,
  first_finite_intervals = nrow(first), second_finite_intervals = nrow(second),
  reference_route = reference_route, reference_method = reference$method,
  reference_exact = reference$exact,
  reference_error_estimate = reference$achieved_absolute_error_estimate,
  reference_squared_distance = reference$squared_distance,
  candidate_squared_distance = candidate$squared_distance,
  absolute_error = absolute_error, acceptance_threshold = threshold,
  active_levels = as.integer(candidate$active_levels),
  expected_active_levels = max(.mv08z_active_depth(first),
                               .mv08z_active_depth(second)),
  rust_status = candidate$status, passed = passed,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
)
if (!passed) stop("MV8-Z sentinel canonical-R oracle failed", call. = FALSE)
.mv08z_atomic_csv(result, output)
message("MV8-Z sentinel oracle passed")
