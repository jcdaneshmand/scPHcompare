#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "usage: validate_mv05l_representation_comparison.R RESULT_DIR OUTPUT_CSV",
    call. = FALSE
  )
}
result_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
source("R/provenance_utils.R")

file_sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
expected <- c(
  mv05e_query_endpoints = "b65dfd0bae4889fd6297372e7ea180ec8f8465487f68e3c6a5bafd4b9920565f",
  mv05e_prediction_lock = "ed2ff4e8d5c52ddabb527853ba5ff60f81dc6017f79bbea068bdcd07939cd053",
  mv05e_independent_validation = "0f8173d57e1bef27311d8e655de867e2b694c6e39c05e7253a52361c3eb2561b",
  mv05e_deterministic_repeat = "8974196f2e28f712b3b726181b7ebdcc056b7ef8cf779bc4c53ac362e9e7bd07",
  mv05e_artifact_manifest = "2a2401ac30a52067530139fcfb4c6d642d681da3e4ada5a8270db13edf97329a",
  mv05k_query_endpoints = "8fa9f70a64f35b87c5e408d32bd7fc4440d5df611c87cd201ff68968b959df2b",
  mv05k_prediction_lock = "e19334a5adfc2fed09dfd0297958f8a837ddbf5ff6753d97ba8e4ea60696ae56",
  mv05k_independent_validation = "bd80254875753b4488b7488545f348ac7ccf478a1b9f2b89a4cc87aefeffbfd2",
  mv05k_deterministic_repeat = "296ed828e30436782b6053c8ca5a35b71e844d1bd6387ac9ef02c0b10d205f49",
  mv05k_artifact_manifest = "983b6be20b1a5858d14b4d08959f382fb09453b64f08c70ed3931e8d10c51a06"
)
audit_dir <- "docs/audits"
paths <- c(
  mv05e_query_endpoints = file.path(audit_dir, "mv05e-query-endpoints-2026-08-08.csv"),
  mv05e_prediction_lock = file.path(audit_dir, "mv05e-prediction-lock-2026-08-08.csv"),
  mv05e_independent_validation = file.path(audit_dir, "mv05e-independent-validation-2026-08-08.csv"),
  mv05e_deterministic_repeat = file.path(audit_dir, "mv05e-deterministic-repeat-2026-08-08.csv"),
  mv05e_artifact_manifest = file.path(audit_dir, "mv05e-artifact-manifest-2026-08-08.csv"),
  mv05k_query_endpoints = file.path(audit_dir, "mv05k-query-endpoints-2026-08-10.csv"),
  mv05k_prediction_lock = file.path(audit_dir, "mv05k-prediction-lock-2026-08-10.csv"),
  mv05k_independent_validation = file.path(audit_dir, "mv05k-independent-validation-2026-08-10.csv"),
  mv05k_deterministic_repeat = file.path(audit_dir, "mv05k-deterministic-repeat-2026-08-10.csv"),
  mv05k_artifact_manifest = file.path(audit_dir, "mv05k-artifact-manifest-2026-08-10.csv")
)
observed_hash <- vapply(paths, file_sha, character(1L))
if (!identical(unname(observed_hash), unname(expected))) {
  stop("Independent MV5-L input hash validation failed.", call. = FALSE)
}

read_result <- function(name) {
  utils::read.csv(file.path(result_dir, name), stringsAsFactors = FALSE,
                  check.names = FALSE)
}
lock <- read_result("mv05l-input-lock-2026-08-10.csv")
method_map_public <- read_result("mv05l-method-map-2026-08-10.csv")
compatibility_public <- read_result("mv05l-endpoint-compatibility-2026-08-10.csv")
paired_public <- read_result("mv05l-paired-query-endpoints-2026-08-10.csv")
identity_public <- read_result("mv05l-pseudobulk-identity-control-2026-08-10.csv")
sample_public <- read_result("mv05l-sample-estimands-2026-08-10.csv")
tissue_public <- read_result("mv05l-tissue-estimands-2026-08-10.csv")
macro_public <- read_result("mv05l-macro-estimands-2026-08-10.csv")
interval_public <- read_result("mv05l-estimand-intervals-2026-08-10.csv")
primary_public <- read_result("mv05l-primary-contrasts-2026-08-10.csv")
bootstrap_public <- read_result("mv05l-bootstrap-audit-2026-08-10.csv")
random_public <- read_result("mv05l-randomization-audit-2026-08-10.csv")
production_public <- read_result("mv05l-production-summary-2026-08-10.csv")

if (nrow(lock) != 10L ||
    !identical(lock$input_id, names(expected)) ||
    !identical(lock$expected_sha256, unname(expected)) ||
    !all(lock$lock_passed_before_endpoint_read) ||
    !all(lock$marginal_aggregate_outcomes_known_at_specification) ||
    any(lock$joint_sample_contrasts_known_at_specification)) {
  stop("Independent MV5-L lock-record validation failed.", call. = FALSE)
}

map <- data.frame(
  family_id = c("h0", "h1", "raw_composite", "energy", "pseudobulk"),
  sct_method_id = c(
    "cell_landscape_h0_v1", "cell_landscape_h1_v1",
    "cell_landscape_h0_h1_raw_euclidean_v1",
    "cell_distribution_energy_shared_pca_v1",
    "pseudobulk_shared_panel_euclidean_v1"),
  integrated_method_id = c(
    "integrated_cell_landscape_h0_v1", "integrated_cell_landscape_h1_v1",
    "integrated_cell_landscape_h0_h1_raw_euclidean_v1",
    "integrated_cell_distribution_energy_v1",
    "pseudobulk_training_standardized_panel_v1"),
  stringsAsFactors = FALSE
)
if (!identical(method_map_public$family_id, map$family_id) ||
    !identical(method_map_public$sct_method_id, map$sct_method_id) ||
    !identical(method_map_public$integrated_method_id,
               map$integrated_method_id) ||
    !all(method_map_public$mapping_frozen_before_cross_representation_join)) {
  stop("Independent MV5-L method-map validation failed.", call. = FALSE)
}

sct <- utils::read.csv(paths[["mv05e_query_endpoints"]],
                       stringsAsFactors = FALSE, check.names = FALSE)
integrated <- utils::read.csv(paths[["mv05k_query_endpoints"]],
                              stringsAsFactors = FALSE, check.names = FALSE)
keys <- c("fold_id", "held_out_study", "seed", "query_sample_id",
          "query_tissue")
common <- c("training_samples", "training_studies", "endpoint_status",
            "labels_opened_after_prediction_lock", "upstream_refit",
            "reranked_after_label_open")
parts <- vector("list", nrow(map))
for (index in seq_len(nrow(map))) {
  left <- sct[sct$method_id == map$sct_method_id[[index]], , drop = FALSE]
  right <- integrated[
    integrated$method_id == map$integrated_method_id[[index]], , drop = FALSE]
  left <- left[order(left$fold_id, left$seed, left$query_sample_id,
                     method = "radix"), , drop = FALSE]
  right <- right[order(right$fold_id, right$seed, right$query_sample_id,
                       method = "radix"), , drop = FALSE]
  if (nrow(left) != 450L || nrow(right) != 450L ||
      !all(vapply(c(keys, common), function(field) {
        identical(left[[field]], right[[field]])
      }, logical(1L)))) {
    stop("Independent MV5-L endpoint pairing failed for ",
         map$family_id[[index]], ".", call. = FALSE)
  }
  parts[[index]] <- data.frame(
    family_id = map$family_id[[index]], left[keys],
    training_samples = left$training_samples,
    training_studies = left$training_studies,
    sct_first_same_tissue_sample_id = left$first_same_tissue_sample_id,
    integrated_first_same_tissue_sample_id = right$first_same_tissue_sample_id,
    sct_first_same_tissue_rank = left$first_same_tissue_rank,
    integrated_first_same_tissue_rank = right$first_same_tissue_rank,
    sct_reciprocal_rank = left$reciprocal_rank,
    integrated_reciprocal_rank = right$reciprocal_rank,
    direct_reciprocal_rank_difference = right$reciprocal_rank - left$reciprocal_rank,
    sct_nearest_sample_id = left$nearest_sample_id,
    integrated_nearest_sample_id = right$nearest_sample_id,
    sct_nearest_tissue = left$nearest_tissue,
    integrated_nearest_tissue = right$nearest_tissue,
    sct_one_nn_correct = left$one_nn_correct,
    integrated_one_nn_correct = right$one_nn_correct,
    direct_one_nn_difference = as.numeric(right$one_nn_correct) -
      as.numeric(left$one_nn_correct),
    sct_nearest_distance_tied = left$nearest_distance_tied,
    integrated_nearest_distance_tied = right$nearest_distance_tied,
    endpoint_status = left$endpoint_status,
    stringsAsFactors = FALSE
  )
}
paired <- do.call(rbind, parts)
paired$family_id <- factor(paired$family_id, levels = map$family_id)
paired <- paired[order(paired$family_id, paired$fold_id, paired$seed,
                       paired$query_sample_id, method = "radix"), , drop = FALSE]
paired$family_id <- as.character(paired$family_id)
rownames(paired) <- NULL

pair_character <- c(
  "family_id", keys, "sct_first_same_tissue_sample_id",
  "integrated_first_same_tissue_sample_id", "sct_nearest_sample_id",
  "integrated_nearest_sample_id", "sct_nearest_tissue",
  "integrated_nearest_tissue", "endpoint_status"
)
if (!all(vapply(pair_character, function(field) {
  identical(paired[[field]], paired_public[[field]])
}, logical(1L)))) {
  stop("Independent MV5-L paired identity reconstruction failed.", call. = FALSE)
}
pair_numeric <- c(
  "training_samples", "training_studies", "sct_first_same_tissue_rank",
  "integrated_first_same_tissue_rank", "sct_reciprocal_rank",
  "integrated_reciprocal_rank", "direct_reciprocal_rank_difference",
  "direct_one_nn_difference"
)
pair_diff <- max(vapply(pair_numeric, function(field) {
  max(abs(paired[[field]] - paired_public[[field]]))
}, numeric(1L)))
if (!is.finite(pair_diff) || pair_diff > 5e-15) {
  stop("Independent MV5-L paired numerical reconstruction failed.", call. = FALSE)
}

pseudo <- paired[paired$family_id == "pseudobulk", , drop = FALSE]
pseudo_checks <- c(
  all(pseudo$sct_first_same_tissue_sample_id ==
        pseudo$integrated_first_same_tissue_sample_id),
  all(pseudo$sct_first_same_tissue_rank ==
        pseudo$integrated_first_same_tissue_rank),
  all(pseudo$sct_reciprocal_rank == pseudo$integrated_reciprocal_rank),
  all(pseudo$sct_nearest_sample_id == pseudo$integrated_nearest_sample_id),
  all(pseudo$sct_nearest_tissue == pseudo$integrated_nearest_tissue),
  all(pseudo$sct_one_nn_correct == pseudo$integrated_one_nn_correct),
  all(pseudo$sct_nearest_distance_tied ==
        pseudo$integrated_nearest_distance_tied)
)
if (!all(pseudo_checks) || !all(identity_public$exact_identity_passed) ||
    any(identity_public$mismatches != 0L)) {
  stop("Independent MV5-L pseudobulk identity control failed.", call. = FALSE)
}
if (nrow(compatibility_public) != 5L ||
    any(compatibility_public$paired_query_seed_rows != 450L) ||
    any(compatibility_public$identity_mismatches != 0L) ||
    any(compatibility_public$denominator_mismatches != 0L) ||
    any(compatibility_public$nonestimable_pairs != 0L)) {
  stop("Independent MV5-L compatibility validation failed.", call. = FALSE)
}

endpoints <- c("cross_study_tissue_mrr_v1",
               "cross_study_tissue_1nn_balanced_accuracy_v1")
estimands <- c(
  "did_h0_topology_minus_energy", "did_h1_topology_minus_energy",
  "direct_h0_integrated_minus_sct", "direct_h1_integrated_minus_sct",
  "direct_raw_composite_integrated_minus_sct",
  "direct_energy_integrated_minus_sct",
  "direct_pseudobulk_integrated_minus_sct"
)
direct_names <- c(
  h0 = estimands[[3L]], h1 = estimands[[4L]],
  raw_composite = estimands[[5L]], energy = estimands[[6L]],
  pseudobulk = estimands[[7L]]
)
query_rows <- list()
counter <- 0L
for (family in names(direct_names)) {
  part <- paired[paired$family_id == family, , drop = FALSE]
  for (endpoint_index in seq_along(endpoints)) {
    counter <- counter + 1L
    query_rows[[counter]] <- data.frame(
      part[keys], endpoint_id = endpoints[[endpoint_index]],
      estimand_id = direct_names[[family]],
      value = if (endpoint_index == 1L)
        part$direct_reciprocal_rank_difference else part$direct_one_nn_difference,
      stringsAsFactors = FALSE
    )
  }
}
energy <- paired[paired$family_id == "energy", , drop = FALSE]
energy <- energy[order(energy$fold_id, energy$seed, energy$query_sample_id,
                       method = "radix"), , drop = FALSE]
for (family in c("h0", "h1")) {
  topology <- paired[paired$family_id == family, , drop = FALSE]
  topology <- topology[order(topology$fold_id, topology$seed,
                             topology$query_sample_id, method = "radix"), ]
  for (endpoint_index in seq_along(endpoints)) {
    counter <- counter + 1L
    topology_direct <- if (endpoint_index == 1L)
      topology$direct_reciprocal_rank_difference else
        topology$direct_one_nn_difference
    energy_direct <- if (endpoint_index == 1L)
      energy$direct_reciprocal_rank_difference else energy$direct_one_nn_difference
    query_rows[[counter]] <- data.frame(
      topology[keys], endpoint_id = endpoints[[endpoint_index]],
      estimand_id = if (family == "h0") estimands[[1L]] else estimands[[2L]],
      value = topology_direct - energy_direct, stringsAsFactors = FALSE
    )
  }
}
query <- do.call(rbind, query_rows)
query$endpoint_id <- factor(query$endpoint_id, levels = endpoints)
query$estimand_id <- factor(query$estimand_id, levels = estimands)
query <- query[order(query$endpoint_id, query$estimand_id,
                     query$query_tissue, query$held_out_study,
                     query$query_sample_id, query$seed,
                     method = "radix"), , drop = FALSE]
query$endpoint_id <- as.character(query$endpoint_id)
query$estimand_id <- as.character(query$estimand_id)

sample_calc <- stats::aggregate(
  value ~ endpoint_id + estimand_id + query_sample_id + query_tissue +
    held_out_study, data = query, FUN = mean
)
names(sample_calc)[names(sample_calc) == "value"] <- "estimate"
sample_calc <- sample_calc[order(
  match(sample_calc$endpoint_id, endpoints),
  match(sample_calc$estimand_id, estimands), sample_calc$query_tissue,
  sample_calc$held_out_study, sample_calc$query_sample_id, method = "radix"
), , drop = FALSE]
sample_public <- sample_public[order(
  match(sample_public$endpoint_id, endpoints),
  match(sample_public$estimand_id, estimands), sample_public$query_tissue,
  sample_public$held_out_study, sample_public$query_sample_id, method = "radix"
), , drop = FALSE]
sample_diff <- max(abs(sample_calc$estimate - sample_public$estimate))
if (sample_diff > 5e-15 || nrow(sample_calc) != 1260L ||
    any(sample_public$completed_seeds != 5L)) {
  stop("Independent MV5-L sample aggregation failed.", call. = FALSE)
}

tissue_calc <- stats::aggregate(
  estimate ~ endpoint_id + estimand_id + query_tissue,
  data = sample_calc, FUN = mean
)
tissue_calc <- tissue_calc[order(
  match(tissue_calc$endpoint_id, endpoints),
  match(tissue_calc$estimand_id, estimands), tissue_calc$query_tissue,
  method = "radix"), , drop = FALSE]
tissue_public <- tissue_public[order(
  match(tissue_public$endpoint_id, endpoints),
  match(tissue_public$estimand_id, estimands), tissue_public$query_tissue,
  method = "radix"), , drop = FALSE]
tissue_diff <- max(abs(tissue_calc$estimate - tissue_public$estimate))
macro_calc <- stats::aggregate(
  estimate ~ endpoint_id + estimand_id, data = tissue_calc, FUN = mean
)
macro_calc <- macro_calc[order(match(macro_calc$endpoint_id, endpoints),
                               match(macro_calc$estimand_id, estimands)), ]
macro_public <- macro_public[order(match(macro_public$endpoint_id, endpoints),
                                   match(macro_public$estimand_id, estimands)), ]
macro_diff <- max(abs(macro_calc$estimate - macro_public$estimate))
if (tissue_diff > 5e-15 || macro_diff > 5e-15 ||
    nrow(tissue_calc) != 70L || nrow(macro_calc) != 14L) {
  stop("Independent MV5-L tissue/macro aggregation failed.", call. = FALSE)
}

sample_ids <- sort(unique(sample_calc$query_sample_id), method = "radix")
base <- unique(sample_calc[, c("query_sample_id", "query_tissue",
                               "held_out_study")])
base <- base[match(sample_ids, base$query_sample_id), , drop = FALSE]
arrays <- lapply(endpoints, function(endpoint) {
  matrix(vapply(estimands, function(estimand) {
    part <- sample_calc[sample_calc$endpoint_id == endpoint &
                          sample_calc$estimand_id == estimand, , drop = FALSE]
    part$estimate[match(sample_ids, part$query_sample_id)]
  }, numeric(length(sample_ids))), nrow = length(sample_ids),
  dimnames = list(sample_ids, estimands))
})
names(arrays) <- endpoints
strata_base <- unique(base[, c("query_tissue", "held_out_study")])
strata <- split(strata_base, strata_base$query_tissue)
strata <- strata[sort(names(strata), method = "radix")]
set.seed(20260810L)
bootstrap <- vector("list", length(endpoints))
for (endpoint_index in seq_along(endpoints)) {
  values <- arrays[[endpoint_index]]
  estimates <- matrix(NA_real_, 2000L, length(estimands),
                      dimnames = list(NULL, estimands))
  for (replicate_id in seq_len(2000L)) {
    tissue_values <- matrix(NA_real_, length(strata), length(estimands))
    for (tissue_index in seq_along(strata)) {
      studies <- sort(strata[[tissue_index]]$held_out_study, method = "radix")
      drawn <- sample(studies, length(studies), replace = TRUE)
      rows <- unlist(lapply(drawn, function(study) {
        which(base$query_tissue == names(strata)[[tissue_index]] &
                base$held_out_study == study)
      }), use.names = FALSE)
      tissue_values[tissue_index, ] <- colMeans(values[rows, , drop = FALSE])
    }
    estimates[replicate_id, ] <- colMeans(tissue_values)
  }
  bootstrap[[endpoint_index]] <- estimates
}
names(bootstrap) <- endpoints
bootstrap_hash <- vapply(bootstrap, digest::digest, character(1L),
                         algo = "sha256", serialize = TRUE)
if (!identical(unname(bootstrap_hash),
               bootstrap_public$replicate_matrix_sha256)) {
  stop("Independent MV5-L bootstrap hash reconstruction failed.", call. = FALSE)
}
interval_calc <- do.call(rbind, lapply(endpoints, function(endpoint) {
  do.call(rbind, lapply(seq_along(estimands), function(index) {
    interval <- stats::quantile(bootstrap[[endpoint]][, index],
                                c(0.025, 0.975), names = FALSE, type = 7)
    data.frame(endpoint_id = endpoint, estimand_id = estimands[[index]],
               ci_lower = interval[[1L]], ci_upper = interval[[2L]])
  }))
}))
interval_public <- interval_public[order(
  match(interval_public$endpoint_id, endpoints),
  match(interval_public$estimand_id, estimands)), , drop = FALSE]
interval_diff <- max(abs(c(
  interval_calc$ci_lower - interval_public$ci_lower,
  interval_calc$ci_upper - interval_public$ci_upper
)))
if (interval_diff > 5e-15) {
  stop("Independent MV5-L interval reconstruction failed.", call. = FALSE)
}

study_ids <- sort(unique(base$held_out_study), method = "radix")
set.seed(20260811L)
signs <- matrix(sample(c(-1, 1), 9999L * length(study_ids), replace = TRUE),
                nrow = 9999L, ncol = length(study_ids),
                dimnames = list(NULL, study_ids))
if (!all(random_public$sign_matrix_sha256 == digest::digest(
  signs, algo = "sha256", serialize = TRUE))) {
  stop("Independent MV5-L sign-matrix reconstruction failed.", call. = FALSE)
}
primary_mrr <- primary_public[
  primary_public$family_id == "F1_representation_did_mrr", , drop = FALSE]
null_hash <- character(2L)
counts <- integer(2L)
p_values <- numeric(2L)
for (index in seq_len(2L)) {
  values <- arrays[[endpoints[[1L]]]][, primary_mrr$estimand_id[[index]]]
  null <- vapply(seq_len(9999L), function(replicate_id) {
    signed <- values * signs[replicate_id, match(base$held_out_study, study_ids)]
    by_tissue <- tapply(signed, base$query_tissue, mean)
    mean(by_tissue[sort(names(by_tissue), method = "radix")])
  }, numeric(1L))
  tolerance <- 64 * .Machine$double.eps *
    pmax(1, abs(null), abs(primary_mrr$estimate[[index]]))
  counts[[index]] <- sum(abs(null) + tolerance >=
                           abs(primary_mrr$estimate[[index]]))
  p_values[[index]] <- (counts[[index]] + 1) / 10000
  null_hash[[index]] <- digest::digest(null, algo = "sha256", serialize = TRUE)
}
random_diff <- max(abs(c(
  counts - random_public$exceedance_count,
  p_values - primary_mrr$raw_p_value,
  stats::p.adjust(p_values, method = "holm") - primary_mrr$holm_p_value
)))
if (!identical(null_hash, random_public$null_distribution_sha256) ||
    random_diff > 5e-15 ||
    any(random_public$boundary_policy !=
          "absolute_two_sided_64eps_ties_count_as_exceedance_v1")) {
  stop("Independent MV5-L randomization/Holm reconstruction failed.",
       call. = FALSE)
}

counter_fields <- c(
  "endpoint_recomputations", "reranking_operations",
  "method_tuning_operations", "method_selection_operations",
  "tissue_selection_operations", "clustering_jobs_executed",
  "gene_view_jobs_executed", "fusion_jobs_executed",
  "new_data_jobs_executed", "optimization_jobs_executed",
  "source_evaluations_modified"
)
if (any(unlist(production_public[counter_fields]) != 0L) ||
    !isTRUE(production_public$pseudobulk_identity_passed) ||
    !isTRUE(production_public$marginal_aggregate_outcomes_known_at_specification) ||
    isTRUE(production_public$joint_sample_contrasts_known_at_specification)) {
  stop("Independent MV5-L boundary-counter validation failed.", call. = FALSE)
}
if (!identical(vapply(paths, file_sha, character(1L)), observed_hash)) {
  stop("An immutable MV5-L source changed during validation.", call. = FALSE)
}

validation <- data.frame(
  contract_id = "mv05l_independent_validation_v1",
  validation = c(
    "input_lock_and_timing_disclosure", "method_map_and_endpoint_pairing",
    "common_identity_and_denominators", "pseudobulk_identity_control",
    "query_direct_and_did_formulas", "seed_within_sample_aggregation",
    "tissue_and_macro_aggregation", "blocked_bootstrap_and_intervals",
    "randomization_boundary_and_holm", "forbidden_operation_counters",
    "source_hashes_unchanged_after_comparison"
  ),
  status = "passed",
  rows_checked = c(10L, 2250L, 2250L, 450L, 6300L, 1260L, 84L,
                   28014L, 19998L, length(counter_fields), 10L),
  max_absolute_difference = c(
    0, pair_diff, 0, 0, pair_diff, sample_diff,
    max(tissue_diff, macro_diff), interval_diff, random_diff, 0, 0
  ),
  evidence = c(
    "ten immutable hashes; pre-join freeze and known-marginal disclosure",
    "five fixed role families; 450 exact pairs per family",
    "fold/study/seed/sample/tissue and training denominators exact",
    "450 shared pseudobulk endpoints exact across seven fields",
    "five direct and two topology-minus-energy DID estimands reconstructed",
    "five technical seeds averaged within each of 90 biological samples",
    "sample means within five tissues then equal-tissue macro means",
    "two 2000-replicate paired study-block matrices and 14 intervals exact",
    "one 9999-row sign matrix, two nulls, 64-epsilon counts, and Holm exact",
    "zero recomputation/tuning/selection/downstream operation counters",
    "all ten accepted source artifact hashes unchanged"
  ),
  stringsAsFactors = FALSE
)
write_provenance_csv(validation, output_path)
message("MV5-L independent representation comparison validation passed.")
