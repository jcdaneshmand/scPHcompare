#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(paste(
    "usage: validate_mv05n_baseline_admission.R PAIRS_A RESOURCES_A PAIRS_B",
    "RESOURCES_B LANDSCAPE_PROJECTION VALIDATION_OUTPUT RESOURCE_OUTPUT",
    "COMBINED_PROJECTION_OUTPUT"
  ), call. = FALSE)
}
source("R/provenance_utils.R")
file_sha <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
read_public <- function(path) utils::read.csv(
  normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE,
  check.names = FALSE
)
pairs_a <- read_public(args[[1L]])
resources_a <- read_public(args[[2L]])
pairs_b <- read_public(args[[3L]])
resources_b <- read_public(args[[4L]])
landscape <- read_public(args[[5L]])

required_methods <- c("cell_distribution_energy_v1",
                      "pseudobulk_training_standardized_panel_v1")
if (nrow(pairs_a) != 384L || nrow(pairs_b) != 384L ||
    nrow(resources_a) != 12L || nrow(resources_b) != 12L ||
    !setequal(pairs_a$method_id, required_methods) ||
    any(!is.finite(pairs_a$distance)) || any(pairs_a$distance < 0) ||
    any(pairs_a$outcome_label_state != "closed") ||
    any(as.logical(pairs_a$biological_outcomes_computed)) ||
    any(pairs_a$clustering_jobs_executed != 0L) ||
    any(!resources_a$elapsed_cap_passed)) {
  stop("MV5-N baseline admission violates its frozen boundary.", call. = FALSE)
}

pair_key <- paste(pairs_a$profile, pairs_a$pair_ordinal,
                  pairs_a$first_sample_id, pairs_a$second_sample_id, sep = "\r")
pseudobulk <- pairs_a[
  pairs_a$method_id == "pseudobulk_training_standardized_panel_v1", ]
pseudo_key <- paste(pseudobulk$profile, pseudobulk$pair_ordinal,
                    pseudobulk$first_sample_id, pseudobulk$second_sample_id,
                    sep = "\r")
pseudo_split <- split(pseudobulk$distance, pseudo_key)
pseudobulk_identity <- all(vapply(pseudo_split, function(values) {
  length(values) == 2L && identical(values[[1L]], values[[2L]])
}, logical(1L)))
byte_repeat <- identical(file_sha(args[[1L]]), file_sha(args[[3L]]))
if (!pseudobulk_identity || !byte_repeat) {
  stop("MV5-N baseline identity or byte-repeat validation failed.", call. = FALSE)
}

validation <- data.frame(
  contract_id = "mv05n_baseline_admission_validation_v1",
  baseline_rows = nrow(pairs_a), unique_training_pairs = length(unique(pair_key)),
  profiles = length(unique(pairs_a$profile)),
  representations = length(unique(pairs_a$representation)),
  methods = length(unique(pairs_a$method_id)),
  finite_nonnegative_distances = all(is.finite(pairs_a$distance) &
                                       pairs_a$distance >= 0),
  pseudobulk_representation_identity_pairs = length(pseudo_split),
  pseudobulk_representation_identity_passed = pseudobulk_identity,
  pair_artifact_primary_sha256 = file_sha(args[[1L]]),
  pair_artifact_repeat_sha256 = file_sha(args[[3L]]),
  pair_artifact_byte_identical = byte_repeat,
  elapsed_caps_passed = all(resources_a$elapsed_cap_passed),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  clustering_jobs_executed = 0L, stringsAsFactors = FALSE
)

resources <- resources_a
resources$contract_id <- "mv05n_baseline_admission_resource_v1"

projection_rows <- list()
cursor <- 0L
for (representation in c("sct_whole", "inductive_integrated")) {
  selected <- resources[resources$representation == representation &
                          resources$method_id ==
                            "cell_distribution_energy_v1", ]
  cursor <- cursor + 1L
  projection_rows[[cursor]] <- data.frame(
    contract_id = "mv05n_combined_full_matrix_resource_projection_v1",
    component = paste0(representation, "_energy"),
    projected_pairs = 262675L,
    projected_worker_hours = max(selected$operation_seconds_per_pair) *
      262675 / 3600,
    projection_rule = "maximum_minimum_representative_maximum_seconds_per_pair_v1",
    stringsAsFactors = FALSE
  )
}
pseudo <- resources[resources$method_id ==
                      "pseudobulk_training_standardized_panel_v1", ]
cursor <- cursor + 1L
projection_rows[[cursor]] <- data.frame(
  contract_id = "mv05n_combined_full_matrix_resource_projection_v1",
  component = "shared_pseudobulk_reused_across_representations",
  projected_pairs = 262675L,
  projected_worker_hours = max(pseudo$operation_seconds_per_pair) *
    262675 / 3600,
  projection_rule = "maximum_profile_seconds_per_pair_compute_once_reuse_v1",
  stringsAsFactors = FALSE
)
topology <- landscape[landscape$representation == "all_representations", ]
cursor <- cursor + 1L
projection_rows[[cursor]] <- data.frame(
  contract_id = "mv05n_combined_full_matrix_resource_projection_v1",
  component = "sct_and_integrated_h0_h1_landscapes",
  projected_pairs = 525350L,
  projected_worker_hours = topology$projected_worker_hours,
  projection_rule = topology$projection_rule,
  stringsAsFactors = FALSE
)
projection <- do.call(rbind, projection_rows)
core_hours <- sum(projection$projected_worker_hours)
reserve_hours <- core_hours * 0.10
projection <- rbind(projection, data.frame(
  contract_id = "mv05n_combined_full_matrix_resource_projection_v1",
  component = "validation_io_and_slow_tail_reserve_10_percent",
  projected_pairs = 0L, projected_worker_hours = reserve_hours,
  projection_rule = "ten_percent_of_measured_core_projection_v1",
  stringsAsFactors = FALSE
))
total_hours <- core_hours + reserve_hours
projection <- rbind(projection, data.frame(
  contract_id = "mv05n_combined_full_matrix_resource_projection_v1",
  component = "all_required_distance_matrices_with_reserve",
  projected_pairs = 1050700L + 262675L * 3L,
  projected_worker_hours = total_hours,
  projection_rule = "sum_landscape_two_energy_shared_pseudobulk_plus_reserve_v1",
  stringsAsFactors = FALSE
))
projection$planning_cap_worker_hours <- 21.6
projection$planning_cap_passed <-
  projection$projected_worker_hours <= projection$planning_cap_worker_hours
projection$full_production_authorized <- FALSE
projection$outcome_label_state <- "closed"
projection$biological_outcomes_computed <- FALSE

write_provenance_csv(validation, args[[6L]])
write_provenance_csv(resources, args[[7L]])
write_provenance_csv(projection, args[[8L]])
message("Validated MV5-N matched baselines and combined resource projection.")
