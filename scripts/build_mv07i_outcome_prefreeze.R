#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv07i_outcome_prefreeze.R LABEL_CLOSED_VALIDATION ",
       "MV7I_PREFREEZE MV7D_PREFREEZE MV7E_PREFREEZE SELECTED_PARTITIONS ",
       "OUTPUT", call. = FALSE)
}
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
sha256 <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
lc_validation <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07i <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07d <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
mv07e <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
selected_path <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
output <- args[[6L]]
if (dir.exists(output) && length(list.files(
    output, all.files = TRUE, no.. = TRUE))) {
  stop("MV7-I outcome prefreeze output must be empty.", call. = FALSE)
}

input_paths <- c(
  file.path(lc_validation, "mv07i-label-closed-decision.csv"),
  file.path(lc_validation, "mv07i-label-closed-validation.csv"),
  file.path(lc_validation, "mv07i-label-closed-artifact-manifest.csv"),
  file.path(mv07i, "mv07i-metadata-registry.csv"),
  file.path(mv07i, "mv07i-representation-registry.csv"),
  file.path(mv07d, "mv07d-sample-reconciliation.csv"),
  file.path(mv07e, "mv07e-canonical-approach.csv"), selected_path)
if (!all(file.exists(input_paths))) stop("MV7-I outcome input is missing.")
lc_decision <- read_csv(input_paths[[1L]])
lc_checks <- read_csv(input_paths[[2L]])
lc_manifest <- read_csv(input_paths[[3L]])
metadata_registry <- read_csv(input_paths[[4L]])
representation_registry <- read_csv(input_paths[[5L]])
reconciliation <- read_csv(input_paths[[6L]])
approach <- read_csv(input_paths[[7L]])
selected <- read_csv(input_paths[[8L]])

selected_manifest <- lc_manifest[
  lc_manifest$artifact == "selected_partitions", , drop = FALSE]
descriptive <- reconciliation[as.logical(
  reconciliation$corrected_descriptive_124), , drop = FALSE]
primary <- descriptive[as.logical(descriptive$corrected_primary_90), , drop = FALSE]
approach <- approach[order(approach$sample_id, method = "radix"), ]
descriptive <- descriptive[order(descriptive$sample_id, method = "radix"), ]
primary <- primary[order(primary$sample_id, method = "radix"), ]
if (nrow(lc_decision) != 1L || lc_decision$decision !=
      "authorize_MV7I_outcome_prefreeze_only" || !all(lc_checks$passed) ||
    nrow(selected_manifest) != 1L ||
    !identical(tolower(sha256(selected_path)),
               tolower(selected_manifest$production_sha256)) ||
    nrow(metadata_registry) != 2L || nrow(representation_registry) != 6L ||
    !identical(tolower(sha256(input_paths[[6L]])), tolower(
      metadata_registry$source_sha256[
        metadata_registry$metadata_id == "tissue_and_study"])) ||
    !identical(tolower(sha256(input_paths[[7L]])), tolower(
      metadata_registry$source_sha256[
        metadata_registry$metadata_id == "canonical_approach"])) ||
    nrow(descriptive) != 124L || nrow(primary) != 90L || nrow(approach) != 124L ||
    anyDuplicated(descriptive$sample_id) || anyDuplicated(approach$sample_id) ||
    !identical(descriptive$sample_id, approach$sample_id) ||
    !identical(descriptive$tissue, approach$tissue) ||
    !identical(descriptive$study, approach$study) || nrow(selected) != 7440L ||
    any(selected$outcome_label_state != "closed") ||
    any(as.logical(selected$biological_outcomes_computed)) ||
    anyDuplicated(selected[c("representation_id", "seed", "algorithm_id",
                             "sample_id")]) ||
    !setequal(unique(selected$sample_id), descriptive$sample_id)) {
  stop("MV7-I outcome-prefreeze admission is stale.", call. = FALSE)
}

input_manifest <- data.frame(
  contract_id = "mv07i_outcome_input_manifest_v1",
  input_order = seq_along(input_paths),
  input_id = c("label_closed_decision", "label_closed_checks",
               "label_closed_artifact_manifest", "metadata_registry",
               "representation_registry", "sample_reconciliation",
               "canonical_approach", "selected_partitions_private"),
  path = input_paths, sha256 = vapply(input_paths, sha256, character(1L)),
  bytes = as.numeric(file.info(input_paths)$size),
  access_role = c(rep("public_admission", 5L),
                  "metadata_structure_only_no_cluster_join",
                  "metadata_structure_only_no_cluster_join",
                  "label_closed_partition_structure_only_no_metadata_join"),
  cluster_metadata_join_executed = FALSE,
  association_computed = FALSE, stringsAsFactors = FALSE)

populations <- data.frame(
  contract_id = "mv07i_outcome_population_registry_v1",
  population_order = 1:2,
  population_id = c("full124_descriptive",
                    "primary90_context_restriction"),
  samples = c(124L, 90L),
  membership_source = c("corrected_descriptive_124_TRUE",
                        "corrected_primary_90_TRUE"),
  partition_policy = c("immutable_full124_partition",
                       "immutable_full124_partition_restricted_no_refit"),
  scientific_role = c("all_eight_tissue_description",
                      "five_tissue_cross_study_context_only"),
  confirmatory = FALSE, external_generalization = FALSE,
  stringsAsFactors = FALSE)

endpoint_rows <- list(); cursor <- 0L
for (population in populations$population_id) for (axis in
    c("tissue", "study", "canonical_approach")) {
  cursor <- cursor + 1L
  estimable <- !(population == "primary90_context_restriction" &&
                 axis == "canonical_approach")
  endpoint_rows[[cursor]] <- data.frame(
    contract_id = "mv07i_outcome_endpoint_registry_v1",
    endpoint_order = cursor, population_id = population, label_axis = axis,
    endpoint_id = paste(population, axis, sep = "__"),
    execution_status = if (estimable) "scheduled" else
      "structurally_not_estimable_single_class",
    scientific_role = if (axis == "tissue")
      "descriptive_biological_alignment" else if (axis == "study")
        "descriptive_technical_cohort_alignment" else
          "descriptive_confounded_technology_proxy",
    causal_interpretation_allowed = FALSE,
    cross_study_generalization_allowed = FALSE,
    association_computed = FALSE, stringsAsFactors = FALSE)
}
endpoints <- do.call(rbind, endpoint_rows)

metadata_counts <- list(); cursor <- 0L
for (population in populations$population_id) {
  value <- if (population == "full124_descriptive") approach else
    approach[as.logical(approach$corrected_primary_90), , drop = FALSE]
  for (axis in c("tissue", "study", "canonical_approach")) {
    counts <- sort(table(value[[axis]]))
    for (index in seq_along(counts)) {
      cursor <- cursor + 1L
      metadata_counts[[cursor]] <- data.frame(
        contract_id = "mv07i_outcome_metadata_count_v1",
        population_id = population, label_axis = axis,
        label_value = names(counts)[[index]], samples = as.integer(counts[[index]]),
        structural_only = TRUE, cluster_metadata_join_executed = FALSE,
        association_computed = FALSE, stringsAsFactors = FALSE)
    }
  }
}
metadata_counts <- do.call(rbind, metadata_counts)

metrics <- data.frame(
  contract_id = "mv07i_outcome_metric_registry_v1",
  metric_order = 1:2,
  metric_id = c("adjusted_rand_index", "normalized_mutual_information_max"),
  definition = c("Hubert_Arabie_complete_contingency_table",
                 "mutual_information_divided_by_maximum_partition_entropy"),
  implementation_contract = c("mv05s_adjusted_rand_index_v1",
                              "mv05s_nmi_max_v1"),
  p_value_authorized = FALSE, multiplicity_family = "none_complete_reporting",
  degenerate_policy = "retain_not_identifiable_never_replace",
  association_computed = FALSE, stringsAsFactors = FALSE)

partition_groups <- unique(selected[c(
  "representation_id", "algorithm_id", "selected_k")])
partition_groups <- partition_groups[order(
  match(partition_groups$representation_id,
        representation_registry$representation_id),
  match(partition_groups$algorithm_id,
        c("pam_stability_k_v1", "hclust_average_v1"))), ]
partition_registry <- data.frame(
  contract_id = "mv07i_outcome_partition_registry_v1",
  partition_order = seq_len(nrow(partition_groups)),
  representation_id = partition_groups$representation_id,
  algorithm_id = partition_groups$algorithm_id,
  algorithm_role = ifelse(partition_groups$algorithm_id ==
    "pam_stability_k_v1", "primary", "sensitivity"),
  selected_k = partition_groups$selected_k,
  seeds = "20260805;20260806;20260807;20260808;20260809",
  samples_per_seed = 124L,
  selected_partition_sha256 = sha256(selected_path), refit_authorized = FALSE,
  outcome_driven_selection_authorized = FALSE, association_computed = FALSE,
  stringsAsFactors = FALSE)

scheduled <- endpoints[endpoints$execution_status == "scheduled", ]
queue <- merge(scheduled[c("endpoint_order", "endpoint_id", "population_id",
                           "label_axis")],
               partition_registry[c("partition_order", "representation_id",
                                    "algorithm_id", "algorithm_role",
                                    "selected_k")], all = TRUE)
queue <- merge(queue, metrics[c("metric_order", "metric_id")], all = TRUE)
queue <- queue[order(queue$endpoint_order, queue$partition_order,
                     queue$metric_order), ]
queue$contract_id <- "mv07i_outcome_execution_queue_v1"
queue$execution_order <- seq_len(nrow(queue))
queue$evaluation_unit_id <- vapply(seq_len(nrow(queue)), function(index) {
  paste0("mv07i_outcome_v1:", digest::digest(list(
    endpoint_id = queue$endpoint_id[[index]],
    representation_id = queue$representation_id[[index]],
    algorithm_id = queue$algorithm_id[[index]],
    metric_id = queue$metric_id[[index]],
    selected_partition_sha256 = sha256(selected_path)),
    algo = "sha256", serialize = TRUE))
}, character(1L))
queue$expected_seeds <- 5L
queue$expected_samples_per_seed <- ifelse(
  queue$population_id == "full124_descriptive", 124L, 90L)
queue$cluster_metadata_join_executed <- FALSE
queue$association_computed <- FALSE
queue$method_selection_executed <- FALSE
queue <- queue[c("contract_id", "evaluation_unit_id", "execution_order",
                 "endpoint_id", "population_id", "label_axis",
                 "representation_id", "algorithm_id", "algorithm_role",
                 "selected_k", "metric_id", "expected_seeds",
                 "expected_samples_per_seed", "cluster_metadata_join_executed",
                 "association_computed", "method_selection_executed")]

aggregation <- data.frame(
  contract_id = "mv07i_outcome_seed_aggregation_v1",
  values_reported = "all_five_seed_estimates",
  summaries = "mean;median;minimum;maximum;delete_one_seed_jackknife_SE",
  seed_interpretation = "technical_realizations_not_biological_replicates",
  favorable_seed_selection = FALSE, algorithm_pooling = FALSE,
  representation_pooling = FALSE, p_value_computed = FALSE,
  association_computed = FALSE, stringsAsFactors = FALSE)

confounding <- data.frame(
  contract_id = "mv07i_outcome_confounding_boundary_v1",
  boundary_id = c("full124_tissue", "full124_study", "full124_approach",
                  "primary90_tissue", "primary90_study",
                  "primary90_approach", "H1_interpretation"),
  required_statement = c(
    "three_added_tissues_are_single_study_no_cross_study_generalization",
    "study_alignment_is_technical_or_cohort_association_not_batch_mechanism",
    "six_snRNA_samples_are_nested_in_substantia_nigra_and_SRA850958",
    "full124_fit_restricted_to_90_context_not_a_90_sample_refit",
    "tissue_homogeneous_studies_preclude_independent_tissue_study_attribution",
    "all_90_are_scRNA_seq_endpoint_not_estimable",
    "H1_alignment_is_not_a_biological_cycle_or_mechanism"),
  causal_claim_allowed = FALSE, association_computed = FALSE,
  stringsAsFactors = FALSE)

publication <- data.frame(
  contract_id = "mv07i_outcome_publication_contract_v1",
  artifact = c("seed_metrics", "unit_summaries", "contingency_long",
               "sample_metadata_join", "artifact_manifest", "resource_summary"),
  publication_state = c("public_after_validation", "public_after_validation",
                        "private", "private", "public_after_validation",
                        "public_after_validation"),
  may_contain_label_values = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  may_contain_sample_ids = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
  complete_reporting_required = TRUE,
  association_computed = FALSE, stringsAsFactors = FALSE)

resources <- data.frame(
  contract_id = "mv07i_outcome_resource_resume_v1",
  maximum_workers = 1L, elapsed_cap_seconds = 900L,
  process_tree_rss_cap_bytes = 2 * 1024^3,
  evaluation_units = nrow(queue), expected_seed_metric_rows = nrow(queue) * 5L,
  private_atomic_artifacts = TRUE, deterministic_repeat_required = TRUE,
  immutable_resume_required = TRUE, independent_recomputation_required = TRUE,
  association_computed = FALSE, stringsAsFactors = FALSE)

checks <- data.frame(
  contract_id = "mv07i_outcome_prefreeze_check_v1",
  check = c("label_closed_admission", "selected_artifact_hash",
            "selected_partition_axis", "metadata_authority",
            "metadata_sample_axis", "full124_category_support",
            "primary90_category_support", "representation_algorithm_scope",
            "endpoint_scope", "metric_scope", "queue_scope",
            "seed_aggregation", "confounding_boundaries", "publication_firewall",
            "resource_repeat_resume"),
  passed = c(
    lc_decision$decision == "authorize_MV7I_outcome_prefreeze_only" &&
      all(lc_checks$passed),
    identical(tolower(sha256(selected_path)),
              tolower(selected_manifest$production_sha256)),
    nrow(selected) == 7440L && length(unique(selected$representation_id)) == 6L &&
      length(unique(selected$seed)) == 5L &&
      length(unique(selected$algorithm_id)) == 2L,
    identical(descriptive$tissue, approach$tissue) &&
      identical(descriptive$study, approach$study) &&
      identical(tolower(sha256(input_paths[[6L]])), tolower(
        metadata_registry$source_sha256[
          metadata_registry$metadata_id == "tissue_and_study"])) &&
      identical(tolower(sha256(input_paths[[7L]])), tolower(
        metadata_registry$source_sha256[
          metadata_registry$metadata_id == "canonical_approach"])) &&
      all(approach$historical_heuristic_use ==
            "prohibited_for_scientific_approach_labels"),
    nrow(descriptive) == 124L && nrow(primary) == 90L &&
      setequal(selected$sample_id, descriptive$sample_id),
    length(unique(approach$tissue)) == 8L &&
      length(unique(approach$study)) == 18L &&
      identical(as.integer(table(approach$canonical_approach)[
        c("scRNA-seq", "snRNA-seq")]), c(118L, 6L)),
    length(unique(primary$tissue)) == 5L &&
      length(unique(primary$study)) == 15L &&
      identical(unique(approach$canonical_approach[
        as.logical(approach$corrected_primary_90)]), "scRNA-seq"),
    nrow(partition_registry) == 12L &&
      setequal(partition_registry$representation_id,
               representation_registry$representation_id) &&
      setequal(partition_registry$algorithm_id,
               c("pam_stability_k_v1", "hclust_average_v1")),
    nrow(endpoints) == 6L && sum(endpoints$execution_status == "scheduled") == 5L,
    nrow(metrics) == 2L && !any(metrics$p_value_authorized),
    nrow(queue) == 120L && anyDuplicated(queue$evaluation_unit_id) == 0L &&
      all(queue$expected_seeds == 5L),
    !aggregation$favorable_seed_selection && !aggregation$algorithm_pooling &&
      !aggregation$representation_pooling,
    nrow(confounding) == 7L && !any(confounding$causal_claim_allowed),
    !any(publication$may_contain_label_values[
      publication$publication_state == "public_after_validation"]) &&
      !any(publication$may_contain_sample_ids[
        publication$publication_state == "public_after_validation"]),
    resources$maximum_workers == 1L && resources$evaluation_units == 120L &&
      resources$expected_seed_metric_rows == 600L &&
      resources$deterministic_repeat_required &&
      resources$immutable_resume_required),
  detail = c(
    "13/13 label-closed validation authorizes outcome prefreeze only",
    "private selected partitions match the public production SHA-256",
    "6 representations x 2 algorithms x 5 seeds x 124 samples",
    "MV7-D tissue/study and MV7-E canonical approach only",
    "124 full and 90 contextual sample axes match immutable partitions",
    "8 tissues, 18 studies, 118 scRNA-seq and 6 nested snRNA-seq",
    "5 tissues, 15 studies, and one approach class",
    "all six representations; PAM primary and average sensitivity",
    "five scheduled endpoints plus explicit non-estimable primary approach",
    "ARI and max-NMI; no p-values",
    "120 fixed units and 600 expected seed rows",
    "all seeds plus mean/median/range/jackknife SE; no pooling",
    "seven mandatory interpretation limits",
    "public summaries contain neither label values nor sample IDs",
    "one worker, 900 seconds, 2 GiB, repeat and immutable resume"),
  stringsAsFactors = FALSE)
if (any(!checks$passed)) {
  stop("MV7-I outcome prefreeze failed: ",
       paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE)
}
decision <- data.frame(
  contract_id = "mv07i_outcome_prefreeze_decision_v1",
  decision = "authorize_MV7I_descriptive_outcome_execution_only",
  checks_passed = sum(checks$passed), checks_total = nrow(checks),
  evaluation_units = nrow(queue), expected_seed_metric_rows = nrow(queue) * 5L,
  labels_joined_to_clusters_now = FALSE, associations_computed_now = FALSE,
  labels_authorized_for_next_stage = TRUE,
  outcomes_authorized_for_next_stage = TRUE,
  claims_authorized = FALSE, method_selection_authorized = FALSE,
  external_data_authorized = FALSE, stringsAsFactors = FALSE)

outputs <- list(
  "mv07i-outcome-input-manifest.csv" = input_manifest,
  "mv07i-outcome-populations.csv" = populations,
  "mv07i-outcome-endpoints.csv" = endpoints,
  "mv07i-outcome-metadata-counts.csv" = metadata_counts,
  "mv07i-outcome-metrics.csv" = metrics,
  "mv07i-outcome-partitions.csv" = partition_registry,
  "mv07i-outcome-queue.csv" = queue,
  "mv07i-outcome-seed-aggregation.csv" = aggregation,
  "mv07i-outcome-confounding-boundaries.csv" = confounding,
  "mv07i-outcome-publication-contract.csv" = publication,
  "mv07i-outcome-resource-resume.csv" = resources,
  "mv07i-outcome-prefreeze-checks.csv" = checks,
  "mv07i-outcome-decision.csv" = decision)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
for (name in names(outputs)) write.csv(
  outputs[[name]], file.path(output, name), row.names = FALSE, na = "")
message("MV7-I descriptive outcome prefreeze passed 15/15.")
