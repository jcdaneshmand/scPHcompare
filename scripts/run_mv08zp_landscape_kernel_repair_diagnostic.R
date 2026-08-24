#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv08zp_landscape_kernel_repair_diagnostic.R <prefreeze>",
  "<private-bindings> <mv08s-private> <mv08v-private> <candidate-library>",
  "<output-dir> <execution-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
bindings_path <- normalizePath(args[[2L]], mustWork = TRUE)
s_root <- normalizePath(args[[3L]], mustWork = TRUE)
v_root <- normalizePath(args[[4L]], mustWork = TRUE)
candidate <- normalizePath(args[[5L]], mustWork = TRUE)
output <- normalizePath(args[[6L]], mustWork = FALSE)
execution_head <- tolower(args[[7L]])
if (dir.exists(output) || !grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV8-ZP requires fresh output and exact execution head", call. = FALSE)
}

source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv08z_landscape_production.R")
.mv08z_verify_manifest(prefreeze, "mv08zp-artifact-manifest.csv")
contract <- .mv08z_read_csv(file.path(prefreeze, "mv08zp-contract.csv"))
inputs <- .mv08z_read_csv(file.path(prefreeze, "mv08zp-input-bindings.csv"))
decision <- .mv08z_read_csv(file.path(prefreeze, "mv08zp-decision.csv"))
if (nrow(contract) != 1L || nrow(decision) != 1L ||
    !.mv08z_truth(decision$diagnostic_execution_authorized) ||
    decision$private_failed_pair_jobs_authorized != 1L ||
    .mv08z_sha256_file(candidate) != contract$candidate_sha256 ||
    .mv08z_sha256_file(bindings_path) !=
      inputs$sha256[inputs$role == "private_unit_bindings"] ||
    execution_head != contract$execution_head) {
  stop("MV8-ZP diagnostic binding drift", call. = FALSE)
}

z_root <- normalizePath(contract$mv08z_root, mustWork = TRUE)
zo_root <- normalizePath(contract$mv08zo_root, mustWork = TRUE)
.mv08z_verify_manifest(z_root, "mv08z-artifact-manifest.csv")
.mv08z_verify_manifest(zo_root, "mv08zo-artifact-manifest.csv")
failure <- .mv08z_read_csv(file.path(zo_root, "mv08zo-failure-evidence.csv"))
groups <- .mv08z_read_csv(file.path(z_root, "mv08z-group-queue.csv"))
bindings <- .mv08z_read_csv(bindings_path)
group <- groups[as.integer(groups$group_order) == failure$failed_group_order,
                , drop = FALSE]
group_bindings <- bindings[
  as.integer(bindings$group_order) == failure$failed_group_order,
  , drop = FALSE
]
pairs <- .mv08z_group_pairs(group_bindings)
pair <- pairs[as.integer(pairs$pair_ordinal) == failure$failed_pair_ordinal,
              , drop = FALSE]

root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("MV8-ZP unknown PH source role", call. = FALSE)
}
intervals_for <- function(axis) {
  row <- group_bindings[
    as.integer(group_bindings$axis_order) == as.integer(axis), , drop = FALSE
  ]
  path <- file.path(root_for(row$source_role), row$output_file)
  if (nrow(row) != 1L || !file.exists(path) ||
      as.numeric(file.info(path)$size) != as.numeric(row$output_bytes) ||
      .mv08z_sha256_file(path) != row$output_sha256) {
    stop("MV8-ZP private PH input drift", call. = FALSE)
  }
  record <- readRDS(path)
  .canonical_rust_intervals(.mv08z_finite_intervals(
    record, group$homology_dimension
  ))
}
first <- intervals_for(pair$first_axis_order)
second <- intervals_for(pair$second_axis_order)
empty <- matrix(numeric(), nrow = 0L, ncol = 2L,
                dimnames = list(NULL, c("birth", "death")))
synthetic <- matrix(c(1, 4, 3, 4, 0, 2, 1, 2), ncol = 2L, byrow = TRUE,
                    dimnames = list(NULL, c("birth", "death")))
started <- proc.time()[["elapsed"]]
pair_result <- landscape_rust_prototype_dimension(first, second, 1L, candidate)
first_norm <- landscape_rust_prototype_dimension(first, empty, 1L, candidate)
second_norm <- landscape_rust_prototype_dimension(second, empty, 1L, candidate)
synthetic_norm <- landscape_rust_prototype_dimension(synthetic, empty, 1L, candidate)
canonical_norm <- function(x) sum((x[, "death"] - x[, "birth"]) ^ 3 / 12)
elapsed <- proc.time()[["elapsed"]] - started

result <- data.frame(
  contract_id = "mv08zp_private_diagnostic_v1",
  check_id = c("failed_pair", "first_norm", "second_norm", "synthetic_norm"),
  rust_status = c(pair_result$status, first_norm$status, second_norm$status,
                  synthetic_norm$status),
  engine_version = c(pair_result$engine_version, first_norm$engine_version,
                     second_norm$engine_version, synthetic_norm$engine_version),
  active_levels = c(pair_result$active_levels, first_norm$active_levels,
                    second_norm$active_levels, synthetic_norm$active_levels),
  expected_active_levels = c(
    max(.mv08z_active_depth(first), .mv08z_active_depth(second)),
    .mv08z_active_depth(first), .mv08z_active_depth(second),
    .mv08z_active_depth(synthetic)
  ),
  squared_distance = c(pair_result$squared_distance,
                       first_norm$squared_distance, second_norm$squared_distance,
                       synthetic_norm$squared_distance),
  exact_norm = c(NA_real_, canonical_norm(first), canonical_norm(second),
                 canonical_norm(synthetic)),
  absolute_norm_error = c(
    NA_real_, abs(first_norm$squared_distance - canonical_norm(first)),
    abs(second_norm$squared_distance - canonical_norm(second)),
    abs(synthetic_norm$squared_distance - canonical_norm(synthetic))
  ),
  passed = c(
    pair_result$status == 0L && pair_result$engine_version == 2L &&
      pair_result$active_levels ==
        max(.mv08z_active_depth(first), .mv08z_active_depth(second)),
    first_norm$status == 0L && first_norm$engine_version == 2L &&
      abs(first_norm$squared_distance - canonical_norm(first)) <= 1e-15,
    second_norm$status == 0L && second_norm$engine_version == 2L &&
      abs(second_norm$squared_distance - canonical_norm(second)) <= 1e-15,
    synthetic_norm$status == 0L && synthetic_norm$engine_version == 2L &&
      abs(synthetic_norm$squared_distance - canonical_norm(synthetic)) <= 1e-12
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (!all(result$passed) || elapsed > contract$diagnostic_elapsed_cap_seconds) {
  stop("MV8-ZP private diagnostic failed", call. = FALSE)
}
resource <- data.frame(
  contract_id = "mv08zp_private_diagnostic_resource_v1",
  execution_head = execution_head, elapsed_seconds = elapsed,
  elapsed_cap_seconds = contract$diagnostic_elapsed_cap_seconds,
  workers = 1L, retries = 0L, production_landscape_jobs = 0L,
  comparison_jobs = 0L, clustering_jobs = 0L, fusion_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
dir.create(output, recursive = TRUE)
.mv08z_atomic_csv(result, file.path(output, "diagnostic-results.csv"))
.mv08z_atomic_csv(resource, file.path(output, "resource.csv"))
manifest <- data.frame(
  artifact = c("diagnostic-results.csv", "resource.csv"),
  bytes = as.numeric(file.info(file.path(output,
    c("diagnostic-results.csv", "resource.csv")))$size),
  sha256 = vapply(file.path(output, c("diagnostic-results.csv", "resource.csv")),
                  .mv08z_sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
.mv08z_atomic_csv(manifest, file.path(output, "artifact-manifest.csv"))
cat("MV8-ZP private diagnostic passed 4/4; engine=2; production=0\n")
