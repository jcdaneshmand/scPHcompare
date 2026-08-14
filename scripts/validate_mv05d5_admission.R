#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: validate_mv05d5_admission.R FIRST_RDS SECOND_RDS FOLD_RDS ",
    "MEAN_RDS D4_COMPONENT_GZ OUTPUT_CSV SAMPLE_COUNT", call. = FALSE
  )
}
first_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
second_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
fold_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
mean_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
d4_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output_path <- args[[6L]]
sample_count <- as.integer(args[[7L]])

source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05d5_retrieval_inputs.R")
file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
first <- readRDS(first_path)
second <- readRDS(second_path)
mv05d5_validate_group_bundle_v1(first)
mv05d5_validate_group_bundle_v1(second)
if (!identical(first, second) || !identical(file_sha(first_path),
                                            file_sha(second_path))) {
  stop("Independent MV5-D5 admission bundles are not byte-identical.",
       call. = FALSE)
}
fold <- readRDS(fold_path)
mean_bundle <- readRDS(mean_path)
mv05d1_validate_cell_fold_record_v1(fold)
mv05d5_validate_mean_profile_bundle_v1(mean_bundle)
d4 <- utils::read.csv(
  d4_path, stringsAsFactors = FALSE, check.names = FALSE
)
d4 <- d4[d4$fold_id == first$identity$fold_id &
           d4$seed == first$identity$seed, , drop = FALSE]
pairs <- first$payload$pairs
key <- function(data) paste(
  data$query_sample_id, data$training_sample_id, sep = "\r"
)
d4 <- d4[match(unique(key(pairs)), key(d4)), , drop = FALSE]
topology_max <- numeric()
for (method_id in c("cell_landscape_h0_v1", "cell_landscape_h1_v1",
                    "cell_landscape_h0_h1_raw_euclidean_v1")) {
  observed <- pairs[pairs$method_id == method_id, , drop = FALSE]
  observed <- observed[match(key(d4), key(observed)), , drop = FALSE]
  expected <- switch(
    method_id, cell_landscape_h0_v1 = d4$h0_distance,
    cell_landscape_h1_v1 = d4$h1_distance,
    cell_landscape_h0_h1_raw_euclidean_v1 =
      sqrt(d4$h0_distance^2 + d4$h1_distance^2)
  )
  topology_max[[method_id]] <- max(abs(observed$distance - expected))
}

energy <- pairs[
  pairs$method_id == "cell_distribution_energy_shared_pca_v1",
  , drop = FALSE
]
energy <- energy[order(energy$query_sample_id, energy$training_sample_id,
                       method = "radix"), , drop = FALSE]
indices <- unique(round(seq(1, nrow(energy), length.out = sample_count)))
energy_expected <- vapply(indices, function(index) {
  .mv05_empirical_energy_distance(
    fold$payload$cell_views[[energy$query_sample_id[[index]]]]$payload,
    fold$payload$cell_views[[energy$training_sample_id[[index]]]]$payload
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
  pairs$method_id == "pseudobulk_shared_panel_euclidean_v1",
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

rank_recomputed <- mv05d5_rank_pairs_v1(pairs)
result <- data.frame(
  contract_id = "mv05d5_admission_validation_v1",
  group_id = first$identity$group_id,
  retrieval_rows = nrow(pairs), methods = length(unique(pairs$method_id)),
  biological_pairs = nrow(pairs) / length(unique(pairs$method_id)),
  exact_repeat_file_sha256 = file_sha(first_path),
  repeat_byte_identical = identical(file_sha(first_path), file_sha(second_path)),
  h0_max_absolute_difference = topology_max[["cell_landscape_h0_v1"]],
  h1_max_absolute_difference = topology_max[["cell_landscape_h1_v1"]],
  combined_max_absolute_difference =
    topology_max[["cell_landscape_h0_h1_raw_euclidean_v1"]],
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
  stop("MV5-D5 admission validation failed.", call. = FALSE)
}
write_provenance_csv(result, output_path)
message("MV5-D5 admission validation passed.")
