#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    "usage: validate_mv05j_admission.R FIRST_RDS SECOND_RDS FOLD_RDS ",
    "MEAN_RDS G_RECORD_RDS I_COMPONENT_GZ OUTPUT_CSV SAMPLE_COUNT",
    call. = FALSE
  )
}
first_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
second_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
fold_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
mean_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
g_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
i_path <- normalizePath(args[[6L]], winslash = "/", mustWork = TRUE)
output_path <- args[[7L]]
sample_count <- as.integer(args[[8L]])

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05d5_retrieval_inputs.R")
source("R/mv05f_integration_gate.R")
source("R/mv05h_integrated_ph_production.R")
source("R/mv05j_integrated_retrieval_inputs.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
first <- readRDS(first_path)
second <- readRDS(second_path)
mv05j_validate_group_bundle_v1(first)
mv05j_validate_group_bundle_v1(second)
if (!identical(first, second) || !identical(file_sha(first_path),
                                            file_sha(second_path))) {
  stop("Independent MV5-J admission bundles are not byte-identical.",
       call. = FALSE)
}
fold <- readRDS(fold_path)
mean_bundle <- readRDS(mean_path)
mv05d1_validate_cell_fold_record_v1(fold)
mv05d5_validate_mean_profile_bundle_v1(mean_bundle)
g_record <- readRDS(g_path)
mv05f_validate_group_record_v1(g_record)
i_components <- utils::read.csv(
  i_path, stringsAsFactors = FALSE, check.names = FALSE
)
i_components <- i_components[i_components$fold_id == first$identity$fold_id &
           i_components$seed == first$identity$seed, , drop = FALSE]
pairs <- first$payload$pairs
key <- function(data) paste(
  data$query_sample_id, data$training_sample_id, sep = "\r"
)
i_components <- i_components[
  match(unique(key(pairs)), key(i_components)), , drop = FALSE
]
topology_max <- numeric()
for (method_id in c("integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1",
                    "integrated_cell_landscape_h0_h1_raw_euclidean_v1")) {
  observed <- pairs[pairs$method_id == method_id, , drop = FALSE]
  observed <- observed[match(key(i_components), key(observed)), , drop = FALSE]
  expected <- switch(
    method_id,
    integrated_cell_landscape_h0_v1 = i_components$h0_distance,
    integrated_cell_landscape_h1_v1 = i_components$h1_distance,
    integrated_cell_landscape_h0_h1_raw_euclidean_v1 =
      sqrt(i_components$h0_distance^2 + i_components$h1_distance^2)
  )
  topology_max[[method_id]] <- max(abs(observed$distance - expected))
}

energy <- pairs[
  pairs$method_id == "integrated_cell_distribution_energy_v1",
  , drop = FALSE
]
energy <- energy[order(energy$query_sample_id, energy$training_sample_id,
                       method = "radix"), , drop = FALSE]
indices <- unique(round(seq(1, nrow(energy), length.out = sample_count)))
energy_expected <- vapply(indices, function(index) {
  .mv05_empirical_energy_distance(
    g_record$payload$coordinates[[energy$query_sample_id[[index]]]],
    g_record$payload$coordinates[[energy$training_sample_id[[index]]]]
  )
}, numeric(1L))
energy_max <- max(abs(energy$distance[indices] - energy_expected))

profile_vector <- function(sample_id) {
  profile <- mean_bundle$payload$profiles[[sample_id]]
  panel <- fold$payload$panel
  present <- panel$feature_id %in% names(profile)
  result <- numeric(nrow(panel))
  result[present] <- (profile[panel$feature_id[present]] -
                        fold$payload$center[present]) /
    fold$payload$scale[present]
  result
}
pseudobulk <- pairs[
  pairs$method_id == "pseudobulk_training_standardized_panel_v1",
  , drop = FALSE
]
pseudobulk <- pseudobulk[order(
  pseudobulk$query_sample_id, pseudobulk$training_sample_id,
  method = "radix"
), , drop = FALSE]
pseudobulk_expected <- vapply(indices, function(index) {
  sqrt(sum((profile_vector(pseudobulk$query_sample_id[[index]]) -
              profile_vector(pseudobulk$training_sample_id[[index]]))^2))
}, numeric(1L))
pseudobulk_max <- max(abs(
  pseudobulk$distance[indices] - pseudobulk_expected
))

rank_recomputed <- mv05j_rank_pairs_v1(pairs)
result <- data.frame(
  contract_id = "mv05j_admission_validation_v1",
  group_id = first$identity$group_id,
  retrieval_rows = nrow(pairs), methods = length(unique(pairs$method_id)),
  biological_pairs = nrow(pairs) / length(unique(pairs$method_id)),
  exact_repeat_file_sha256 = file_sha(first_path),
  repeat_byte_identical = identical(file_sha(first_path), file_sha(second_path)),
  h0_max_absolute_difference = topology_max[["integrated_cell_landscape_h0_v1"]],
  h1_max_absolute_difference = topology_max[["integrated_cell_landscape_h1_v1"]],
  combined_max_absolute_difference =
    topology_max[["integrated_cell_landscape_h0_h1_raw_euclidean_v1"]],
  energy_oracle_pairs = length(indices),
  energy_max_absolute_difference = energy_max,
  pseudobulk_oracle_pairs = length(indices),
  pseudobulk_max_absolute_difference = pseudobulk_max,
  rankings_canonical = identical(rank_recomputed, pairs),
  outcome_columns_absent = !any(c("tissue", "approach") %in% names(pairs)),
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
if (!result$repeat_byte_identical ||
    any(unlist(result[c(
      "h0_max_absolute_difference", "h1_max_absolute_difference",
      "combined_max_absolute_difference", "energy_max_absolute_difference",
      "pseudobulk_max_absolute_difference"
    )], use.names = FALSE) > 1e-12) || !result$rankings_canonical ||
    !result$outcome_columns_absent) {
  stop("MV5-J admission validation failed.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("MV5-J admission validation passed.")
