#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: build_mv10p_clustering_outcome_prefreeze.R <private-partitions>",
  "<mv10e-public> <mv10o-closure> <mv07d> <mv07e> <output-dir>"
), call. = FALSE)
partitions_path <- normalizePath(args[[1L]], mustWork = TRUE)
mv10e <- normalizePath(args[[2L]], mustWork = TRUE)
mv10o <- normalizePath(args[[3L]], mustWork = TRUE)
mv07d <- normalizePath(args[[4L]], mustWork = TRUE)
mv07e <- normalizePath(args[[5L]], mustWork = TRUE)
output <- args[[6L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-P prefreeze")
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
source("R/mv05s_outcome_execution.R")
source("R/mv10p_clustering_outcomes.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
execution_head <- tolower(trimws(Sys.getenv("MV10P_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("MV10P_GIT_HEAD must bind one 40-character commit")
}
.mv08z_verify_manifest(mv10e, "mv10e-artifact-manifest.csv")
.mv08z_verify_manifest(mv10o, "mv10o-artifact-manifest.csv")
partitions <- readc(partitions_path)
selected <- mv10p_select_outcome_partitions_v1(partitions)
reconciliation_path <- file.path(mv07d, "mv07d-sample-reconciliation.csv")
approach_path <- file.path(mv07e, "mv07e-canonical-approach.csv")
reconciliation <- readc(reconciliation_path)
approach <- readc(approach_path)
descriptive <- reconciliation[as.logical(reconciliation$corrected_descriptive_124),
                               , drop = FALSE]
primary <- descriptive[as.logical(descriptive$corrected_primary_90), ,
                       drop = FALSE]
descriptive <- descriptive[order(descriptive$sample_id), , drop = FALSE]
primary <- primary[order(primary$sample_id), , drop = FALSE]
approach <- approach[order(approach$sample_id), , drop = FALSE]
if (nrow(descriptive) != 124L || nrow(primary) != 90L || nrow(approach) != 124L ||
    !identical(descriptive$sample_id, approach$sample_id) ||
    !identical(descriptive$tissue, approach$tissue) ||
    !identical(descriptive$study, approach$study) ||
    !setequal(selected$sample_id, descriptive$sample_id)) {
  stop("MV10-P metadata/sample axis drifted.")
}
methods <- mv10p_method_registry_v1()
groups <- unique(selected[c("stack_id", "homology_dimension", "method_id", "k")])
groups <- merge(groups, methods, by = "method_id", all.x = TRUE)
stacks <- c("allqc_data_exact500", "allqc_residual_exact500",
            "existing_selectedfit_data_exact500")
groups <- groups[order(match(groups$stack_id, stacks),
                       groups$homology_dimension, groups$method_order), ]
group_hash <- vapply(seq_len(nrow(groups)), function(index) {
  value <- selected[
    selected$stack_id == groups$stack_id[[index]] &
      selected$homology_dimension == groups$homology_dimension[[index]] &
      selected$method_id == groups$method_id[[index]],
    c("stack_id", "homology_dimension", "method_id", "k", "seed",
      "sample_id", "cluster"), drop = FALSE
  ]
  digest::digest(value, algo = "sha256", serialize = TRUE)
}, character(1L))
partition_registry <- data.frame(
  contract_id = "mv10p_partition_registry_v1",
  partition_order = seq_len(nrow(groups)), stack_id = groups$stack_id,
  homology_dimension = groups$homology_dimension,
  method_id = groups$method_id, method_role = groups$method_role,
  selected_k = groups$k, seeds = "20260805;20260806;20260807;20260808;20260809",
  samples_per_seed = 124L, partition_group_sha256 = group_hash,
  source_file_sha256 = sha(partitions_path), refit_authorized = FALSE,
  outcome_driven_selection_authorized = FALSE, association_computed = FALSE,
  stringsAsFactors = FALSE
)
populations <- data.frame(
  contract_id = "mv10p_population_registry_v1", population_order = 1:2,
  population_id = c("full124_descriptive", "primary90_context_restriction"),
  samples = c(124L, 90L),
  partition_policy = c("immutable_full124_partition",
                       "immutable_full124_partition_restricted_no_refit"),
  scientific_role = c("all_eight_tissue_description",
                      "five_tissue_cross_study_context_only"),
  confirmatory = FALSE, external_generalization = FALSE,
  stringsAsFactors = FALSE
)
endpoint_rows <- list(); cursor <- 0L
for (population in populations$population_id) for (axis in
    c("tissue", "study", "canonical_approach")) {
  cursor <- cursor + 1L
  estimable <- !(population == "primary90_context_restriction" &&
                 axis == "canonical_approach")
  endpoint_rows[[cursor]] <- data.frame(
    contract_id = "mv10p_endpoint_registry_v1", endpoint_order = cursor,
    population_id = population, label_axis = axis,
    endpoint_id = paste(population, axis, sep = "__"),
    execution_status = if (estimable) "scheduled" else
      "structurally_not_estimable_single_class",
    scientific_role = if (axis == "tissue")
      "descriptive_biological_alignment" else if (axis == "study")
        "descriptive_technical_cohort_alignment" else
          "descriptive_confounded_technology_proxy",
    causal_interpretation_allowed = FALSE,
    cross_study_generalization_allowed = FALSE,
    association_computed = FALSE, stringsAsFactors = FALSE
  )
}
endpoints <- do.call(rbind, endpoint_rows)
metrics <- data.frame(
  contract_id = "mv10p_metric_registry_v1", metric_order = 1:2,
  metric_id = c("adjusted_rand_index", "normalized_mutual_information_max"),
  definition = c("Hubert_Arabie_complete_contingency_table",
                 "mutual_information_divided_by_maximum_partition_entropy"),
  p_value_authorized = FALSE, multiplicity_family = "none_complete_reporting",
  association_computed = FALSE, stringsAsFactors = FALSE
)
queue <- mv10p_build_queue_v1(partition_registry, endpoints, metrics)
aggregation <- data.frame(
  contract_id = "mv10p_seed_aggregation_v1",
  values_reported = "all_five_seed_estimates",
  summaries = "mean;median;minimum;maximum;delete_one_seed_jackknife_SE",
  seed_interpretation = "technical_realizations_not_biological_replicates",
  favorable_seed_selection = FALSE, method_pooling = FALSE,
  representation_pooling = FALSE, homology_pooling = FALSE,
  p_value_computed = FALSE, stringsAsFactors = FALSE
)
confounding <- readc(
  "docs/audits/mv07i-outcome-prefreeze-evidence/mv07i-outcome-confounding-boundaries.csv"
)
confounding$contract_id <- "mv10p_confounding_boundary_v1"
confounding$association_computed <- FALSE
publication <- data.frame(
  contract_id = "mv10p_publication_contract_v1",
  artifact = c("seed_metrics", "unit_summaries", "structural_status",
               "provenance", "metadata_join", "contingency_long"),
  publication_state = c(rep("public_after_validation", 4L), "private", "private"),
  may_contain_label_values = c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE),
  may_contain_sample_ids = c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE),
  complete_reporting_required = TRUE, stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv10p_clustering_outcomes.R",
  "scripts/build_mv10p_clustering_outcome_prefreeze.R",
  "scripts/run_mv10q_clustering_outcomes.R",
  "scripts/build_mv10r_clustering_outcome_closure.R"
)
implementation <- data.frame(
  contract_id = "mv10p_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  partitions_path, file.path(mv10e, "mv10e-artifact-manifest.csv"),
  file.path(mv10o, "mv10o-artifact-manifest.csv"),
  file.path(mv10o, "mv10o-validation.csv"), reconciliation_path, approach_path,
  "docs/audits/mv07i-outcome-prefreeze-evidence/mv07i-outcome-metadata-counts.csv",
  "docs/audits/mv07i-outcome-prefreeze-evidence/mv07i-outcome-metrics.csv",
  "docs/audits/mv07i-outcome-prefreeze-evidence/mv07i-outcome-confounding-boundaries.csv"
)
source_freeze <- data.frame(
  contract_id = "mv10p_source_freeze_v1", source_order = seq_along(source_files),
  artifact = source_files, bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)), stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv10p_clustering_outcome_prefreeze_v1",
  execution_head = execution_head, source_assignment_rows = nrow(partitions),
  selected_assignment_rows = nrow(selected), partition_families = 30L,
  evaluation_units = 300L, expected_seed_metric_rows = 1500L,
  selected_H0_k = 2L, selected_H1_k = 3L,
  labels_joined_to_clusters_now = FALSE, associations_computed_now = FALSE,
  p_values_authorized = FALSE, method_selection_authorized = FALSE,
  biological_claims_authorized = FALSE, manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
resources <- data.frame(
  contract_id = "mv10p_resource_contract_v1", maximum_workers = 1L,
  automatic_retries = 0L, elapsed_cap_seconds = 900L,
  process_tree_rss_cap_bytes = 2 * 1024^3,
  evaluation_units = nrow(queue), expected_seed_metric_rows = nrow(queue) * 5L,
  public_storage_cap_bytes = 16 * 1024^2,
  private_storage_cap_bytes = 16 * 1024^2,
  independent_recomputation_required = TRUE,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10p_validation_v1",
  check_id = c("mv10e_manifest", "mv10o_manifest", "private_source_hash",
               "source_rows", "selected_rows", "sample_axis_124",
               "primary_axis_90", "metadata_join_axis", "three_stacks",
               "two_dimensions", "five_methods", "five_seeds",
               "thirty_partition_groups", "selected_K_2_3", "six_endpoints",
               "five_scheduled", "two_metrics", "three_hundred_units",
               "fifteen_hundred_seed_rows", "complete_reporting",
               "no_method_specific_retuning", "no_method_ranking",
               "no_representation_ranking", "no_H0_H1_pooling",
               "confounding_boundaries", "public_privacy", "implementation_bound",
               "sources_bound", "resource_contract", "no_association_now"),
  passed = c(TRUE, TRUE, grepl("^[0-9a-f]{64}$", sha(partitions_path)),
             nrow(partitions) == 167400L, nrow(selected) == 18600L,
             nrow(descriptive) == 124L, nrow(primary) == 90L,
             identical(descriptive$sample_id, approach$sample_id),
             length(unique(selected$stack_id)) == 3L,
             length(unique(selected$homology_dimension)) == 2L,
             length(unique(selected$method_id)) == 5L,
             length(unique(selected$seed)) == 5L,
             nrow(partition_registry) == 30L,
             all(partition_registry$selected_k[
               partition_registry$homology_dimension == "H0"] == 2L) &&
               all(partition_registry$selected_k[
                 partition_registry$homology_dimension == "H1"] == 3L),
             nrow(endpoints) == 6L,
             sum(endpoints$execution_status == "scheduled") == 5L,
             nrow(metrics) == 2L && !any(metrics$p_value_authorized),
             nrow(queue) == 300L && !anyDuplicated(queue$evaluation_unit_id),
             nrow(queue) * 5L == 1500L,
             !aggregation$favorable_seed_selection,
             all(!partition_registry$outcome_driven_selection_authorized),
             !aggregation$method_pooling, !aggregation$representation_pooling,
             !aggregation$homology_pooling,
             nrow(confounding) == 7L && !any(confounding$causal_claim_allowed),
             !any(publication$may_contain_label_values[
               publication$publication_state == "public_after_validation"]) &&
               !any(publication$may_contain_sample_ids[
                 publication$publication_state == "public_after_validation"]),
             nrow(implementation) == 4L && all(file.exists(implementation$file)),
             length(source_files) == 9L && all(file.exists(source_files)),
             resources$maximum_workers == 1L && resources$automatic_retries == 0L,
             !any(queue$association_computed)),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-P validation failed")
decision <- data.frame(
  contract_id = "mv10p_decision_v1",
  decision = "authorize_fixed_descriptive_outcome_execution_after_commit",
  execution_authorized_after_commit = TRUE, evaluation_units = 300L,
  expected_seed_rows = 1500L, p_values_authorized = FALSE,
  method_selection_authorized = FALSE, biological_claims_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  next_if_closed = "separate_complete_outcome_review_prefreeze",
  stringsAsFactors = FALSE
)
artifacts <- list(
  "mv10p-contract.csv" = contract,
  "mv10p-partition-registry.csv" = partition_registry,
  "mv10p-populations.csv" = populations, "mv10p-endpoints.csv" = endpoints,
  "mv10p-metrics.csv" = metrics, "mv10p-queue.csv" = queue,
  "mv10p-seed-aggregation.csv" = aggregation,
  "mv10p-confounding-boundaries.csv" = confounding,
  "mv10p-publication-contract.csv" = publication,
  "mv10p-implementation-bindings.csv" = implementation,
  "mv10p-source-freeze.csv" = source_freeze,
  "mv10p-resource-contract.csv" = resources,
  "mv10p-validation.csv" = validation, "mv10p-decision.csv" = decision
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-P clustering-outcome prefreeze", "",
  "This extends the already validated MV7-I descriptive label contract to the",
  "immutable MV10 partitions. All five fixed method roles are reported at the",
  "PAM-selected H0=2 and H1=3 K values. No method-specific retuning, ranking,",
  "p-value, homology pooling, biological claim, or manuscript claim is allowed."
), file.path(output, "MV10P_CLUSTERING_OUTCOME_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10p-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10p_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10p-artifact-manifest.csv"))
message("Built MV10-P clustering-outcome prefreeze; checks=30")
