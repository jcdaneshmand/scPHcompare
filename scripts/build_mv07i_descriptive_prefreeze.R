#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: build_mv07i_descriptive_prefreeze.R MV7H_VALIDATION ",
       "MV7E_PREFREEZE MV7D_PREFREEZE OUTPUT", call. = FALSE)
}
source("R/mv07h_full_topology.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
mv07h <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07e <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07d <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output <- args[[4L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-I prefreeze output must be empty.", call. = FALSE)
}
input_paths <- c(
  file.path(mv07h, "mv07h-landscape-complete-validation-decision.csv"),
  file.path(mv07h, "mv07h-landscape-complete-independent-validation.csv"),
  file.path(mv07h, "mv07h-landscape-complete-group-inventory.csv"),
  file.path(mv07e, "mv07e-sample-seed-axis.csv"),
  file.path(mv07e, "mv07e-canonical-approach.csv"),
  file.path(mv07e, "mv07e-approach-field-lineage.csv"),
  file.path(mv07d, "mv07d-sample-reconciliation.csv"),
  file.path(mv07d, "mv07d-estimand-populations.csv"),
  file.path(mv07d, "mv07d-tissue-study-summary.csv")
)
if (!all(file.exists(input_paths))) stop("MV7-I prefreeze input is missing.")
decision_h <- read_csv(input_paths[[1L]])
validation_h <- read_csv(input_paths[[2L]])
inventory <- read_csv(input_paths[[3L]])
axis <- read_csv(input_paths[[4L]])
approach <- read_csv(input_paths[[5L]])
lineage <- read_csv(input_paths[[6L]])
reconciliation <- read_csv(input_paths[[7L]])
estimands <- read_csv(input_paths[[8L]])
tissues <- read_csv(input_paths[[9L]])
input_manifest <- data.frame(
  contract_id = "mv07i_input_manifest_v1", input_order = seq_along(input_paths),
  path = input_paths, sha256 = vapply(input_paths, .mv07h_sha256, character(1L)),
  bytes = as.numeric(file.info(input_paths)$size), access_role = c(
    rep("label_closed_scientific_admission", 4L),
    "metadata_lineage_only_no_outcome_execution",
    "metadata_lineage_only_no_outcome_execution",
    "population_structure_only_no_outcome_execution",
    "population_structure_only_no_outcome_execution",
    "population_structure_only_no_outcome_execution"),
  stringsAsFactors = FALSE
)

representation <- data.frame(
  contract_id = "mv07i_representation_registry_v1",
  representation_order = 1:6,
  representation_id = c("cell_H0", "cell_H1", "cell_H0_H1_secondary",
                        "gene_H0", "gene_H1", "gene_H0_H1_secondary"),
  view_id = rep(c("cell_topology_v1", "gene_topology_v1"), each = 3L),
  component = rep(c("H0", "H1", "H0_H1_unweighted"), 2L),
  scientific_role = rep(c("primary_component", "primary_component",
                          "secondary_descriptive_composite"), 2L),
  seed_specific_matrices = 5L, samples = 124L,
  construction = rep(c("direct_validated_distance", "direct_validated_distance",
                       "sqrt(H0_squared_distance + H1_squared_distance)"), 2L),
  dimension_replacement_allowed = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
aggregation <- data.frame(
  contract_id = "mv07i_seed_aggregation_contract_v1",
  quantity = c("distance", "combined_distance", "h1_squared_fraction",
               "cluster_stability", "later_label_metric"),
  computation_order = c("direct_per_seed", "within_seed_then_aggregate",
                        "within_seed_then_aggregate", "across_seed_partitions",
                        "per_seed_then_jackknife_after_separate_authorization"),
  center = c(rep("median_across_five_seeds", 3L),
             "mean_of_ten_pairwise_seed_ARIs",
             "mean_across_five_seeds"),
  uncertainty = c(rep("min;max;IQR;raw_MAD_constant_1", 3L),
                  "delete_one_seed_jackknife_SE",
                  "delete_one_seed_jackknife_SE"),
  zero_rule = c(NA, NA, "fraction_zero_only_when_H0_and_H1_are_both_zero",
                NA, NA),
  select_favorable_seed = FALSE, stringsAsFactors = FALSE
)
clustering <- data.frame(
  contract_id = "mv07i_clustering_contract_v1",
  representation_id = representation$representation_id,
  candidate_k = "2:10", seeds = "20260805;20260806;20260807;20260808;20260809",
  primary_algorithm = "PAM_dissimilarity_v1",
  stability_statistic = "mean_10_pairwise_seed_adjusted_Rand_indices",
  stability_uncertainty = "five_seed_delete_one_seed_jackknife_SE",
  k_selection = "smallest_k_within_one_SE_of_maximum_mean_stability",
  failure_state = "no_stable_k",
  sensitivity_algorithm = "average_linkage_at_PAM_selected_k",
  canonical_cluster_ids = "sorted_member_signature_v1",
  spectral_status = "deferred_pending_affinity_eigengap_stability_prefreeze",
  ward_status = "prohibited_for_arbitrary_dissimilarity",
  kmeans_status = "prohibited_on_distance_matrix",
  label_values_used = FALSE, stringsAsFactors = FALSE
)
populations <- estimands[estimands$population_id %in% c(
  "corrected_primary_cross_study", "corrected_full_corpus_descriptive",
  "below_threshold_sensitivity"), , drop = FALSE]
populations$mv07i_role <- c(
  "cross_study_context_and_prior_primary_comparison",
  "all_eight_tissue_descriptive_topology_and_clustering",
  "excluded_unless_separate_threshold_sensitivity_authorized")
populations$may_support_new_cross_study_tissue_claim <-
  populations$population_id == "corrected_primary_cross_study"
populations$outcome_label_state <- "closed"
populations$biological_outcomes_computed <- FALSE

metadata_registry <- data.frame(
  contract_id = "mv07i_metadata_registry_v1",
  metadata_id = c("tissue_and_study", "canonical_approach"),
  source_path = c(input_paths[[7L]], input_paths[[5L]]),
  source_sha256 = vapply(c(input_paths[[7L]], input_paths[[5L]]),
                         .mv07h_sha256, character(1L)),
  join_key = "sample_id", expected_rows = 124L,
  authoritative_fields = c("tissue;study", "canonical_approach"),
  prohibited_fields = c("approach;approach_historical_retained",
                        "historical_heuristic_approach"),
  access_stage = "after_immutable_clustering_and_separate_outcome_prefreeze",
  may_select_method_k_view_dimension_seed = FALSE,
  stringsAsFactors = FALSE
)
interpretation <- data.frame(
  contract_id = "mv07i_interpretation_boundary_v1",
  scope = c("all_124", "primary_90", "added_34", "approach",
            "H1", "external_generalization"),
  allowed = c("descriptive_association_and_topology_location",
              "contextualize_against_prior_blocked_cross_study_results",
              "descriptive_only_complete_reporting",
              "descriptive_association_using_canonical_field_only",
              "distance_contribution_and_stability_only",
              "none_without_external_data_gate"),
  prohibited = c("confirmatory_or_causal_claim",
                 "retroactive_method_selection",
                 "cross_study_tissue_generalization",
                 "causal_technology_effect",
                 "biological_cycle_or_mechanism_claim",
                 "generalization_claim_from_current_corpus"),
  stringsAsFactors = FALSE
)
stages <- data.frame(
  contract_id = "mv07i_stage_authorization_v1",
  stage_order = 1:4,
  stage = c("matrix_reconstruction", "label_closed_clustering",
            "outcome_prefreeze", "descriptive_outcome_execution"),
  authorized_now = c(TRUE, TRUE, FALSE, FALSE),
  label_access = c(FALSE, FALSE, "structural_only", TRUE),
  prerequisite = c("MV7H_13_of_13", "immutable_matrices",
                   "immutable_partitions", "separate_outcome_prefreeze_pass"),
  output_role = c("private_matrices_public_manifest",
                  "private_partitions_public_label_free_stability",
                  "fixed_metadata_endpoints_and_reporting_order",
                  "complete_descriptive_results"),
  stringsAsFactors = FALSE
)
resources <- data.frame(
  contract_id = "mv07i_resource_resume_contract_v1",
  maximum_workers = 1L, elapsed_cap_seconds = 3600,
  process_tree_rss_cap_bytes = 4 * 1024^3,
  candidate_pam_fits = 6L * 5L * 9L,
  selected_average_linkage_fits = 6L * 5L,
  atomic_private_artifacts = TRUE, deterministic_repeat_required = TRUE,
  immutable_resume_required = TRUE, public_results_allowed = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
checks <- data.frame(
  contract_id = "mv07i_descriptive_prefreeze_check_v1",
  category = c("mv07h_admission", "landscape_inventory", "sample_seed_axis",
               "representation_scope", "dimension_separation",
               "seed_aggregation", "label_free_k", "algorithm_scope",
               "population_scope", "metadata_lineage", "single_study_limit",
               "label_firewall", "resource_resume"),
  passed = c(
    decision_h$decision == "authorize_MV7I_descriptive_prefreeze_only" &&
      all(validation_h$passed),
    nrow(inventory) == 20L && sum(inventory$component_rows) == 152520L,
    nrow(axis) == 620L && length(unique(axis$sample_id)) == 124L &&
      length(unique(axis$seed)) == 5L,
    nrow(representation) == 6L,
    sum(representation$scientific_role == "primary_component") == 4L &&
      !any(representation$dimension_replacement_allowed),
    all(aggregation$select_favorable_seed == FALSE),
    all(clustering$candidate_k == "2:10") &&
      !any(clustering$label_values_used),
    all(clustering$primary_algorithm == "PAM_dissimilarity_v1") &&
      all(clustering$sensitivity_algorithm ==
            "average_linkage_at_PAM_selected_k"),
    setequal(populations$population_id, c("corrected_primary_cross_study",
      "corrected_full_corpus_descriptive", "below_threshold_sensitivity")),
    nrow(approach) == 124L && nrow(lineage) >= 3L &&
      all(metadata_registry$expected_rows == 124L),
    nrow(tissues) == 8L && sum(!tissues$primary_cross_study_eligible) == 3L,
    !any(stages$label_access[stages$authorized_now] != "FALSE") &&
      !any(stages$authorized_now[stages$stage %in%
        c("outcome_prefreeze", "descriptive_outcome_execution")]),
    resources$maximum_workers == 1L &&
      resources$candidate_pam_fits == 270L &&
      resources$deterministic_repeat_required &&
      resources$immutable_resume_required
  ),
  detail = c("MV7-H full landscape decision and 13/13 validation",
             "20 groups and 152,520 rows", "124 samples x five seeds",
             "four primary components plus two secondary composites",
             "H0/H1 primary components never replaced",
             "within-seed derivation and complete five-seed summaries",
             "k=2:10 five-seed one-SE without labels",
             "PAM primary and average-linkage sensitivity only",
             "90 primary context, 124 descriptive, three excluded sensitivity",
             "canonical sample/tissue/study/approach lineage",
             "three single-study tissues descriptive only",
             "labels and outcomes closed for authorized stages",
             "one worker, caps, repeat and immutable resume"),
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) {
  stop("MV7-I descriptive prefreeze failed: ",
       paste(checks$category[!checks$passed], collapse = ", "),
       call. = FALSE)
}
decision <- data.frame(
  contract_id = "mv07i_descriptive_prefreeze_decision_v1",
  decision = "authorize_label_closed_matrix_and_clustering_production_only",
  representations = 6L, seed_specific_matrices = 30L,
  candidate_pam_fits = 270L, average_linkage_fits_after_k_selection = 30L,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  claims_authorized = FALSE, external_data_authorized = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
outputs <- list(
  "mv07i-input-manifest.csv" = input_manifest,
  "mv07i-representation-registry.csv" = representation,
  "mv07i-seed-aggregation-contract.csv" = aggregation,
  "mv07i-clustering-contract.csv" = clustering,
  "mv07i-populations.csv" = populations,
  "mv07i-metadata-registry.csv" = metadata_registry,
  "mv07i-interpretation-boundaries.csv" = interpretation,
  "mv07i-stage-authorization.csv" = stages,
  "mv07i-resource-resume-contract.csv" = resources,
  "mv07i-prefreeze-checks.csv" = checks,
  "mv07i-decision.csv" = decision
)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (name in names(outputs)) write.csv(
  outputs[[name]], file.path(output, name), row.names = FALSE, na = "")
message("MV7-I descriptive prefreeze: 13/13 pass")
