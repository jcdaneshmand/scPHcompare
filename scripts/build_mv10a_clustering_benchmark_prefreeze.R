#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: build_mv10a_clustering_benchmark_prefreeze.R <mv07h-root>",
  "<mv08zu-private-root> <mv08zx-private-root> <output-dir>",
  "<implementation-head>"
), call. = FALSE)
roots <- vapply(args[1:3], normalizePath, character(1L), mustWork = TRUE)
output <- args[[4L]]
implementation_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", implementation_head)) stop("invalid head")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("refusing to overwrite MV10-A prefreeze")
}
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv08z_landscape_production.R")
source("R/mv08zy_distance_comparison.R")
source("R/mv10_clustering_benchmark.R")
readc <- .mv08z_read_csv; sha <- .mv08z_sha256_file; atomic <- .mv08z_atomic_csv
mv08zy <- "docs/audits/mv08zy-distance-comparison-execution-prefreeze-v1"
mv08zz <- "docs/audits/mv08zz-distance-comparison-closure-v1"
mv09m <- "docs/audits/mv09m-corrected-review-visual-inspection-v1"
.mv08z_verify_manifest(mv08zy, "mv08zy-artifact-manifest.csv")
.mv08z_verify_manifest(mv08zz, "mv08zz-artifact-manifest.csv")
.mv08z_verify_manifest(mv09m, "mv09m-artifact-manifest.csv")
catalog_all <- readc(file.path(mv08zy, "mv08zy-stack-bindings.csv"))
catalog <- catalog_all[catalog_all$dataset_scope == "internal124", , drop = FALSE]
catalog <- catalog[order(catalog$catalog_order), , drop = FALSE]
mv10_validate_stack_catalog_v1(catalog)
verification <- lapply(seq_len(nrow(catalog)), function(i) {
  binding <- catalog[i, , drop = FALSE]
  loaded <- mv08zy_read_distance_stack_v1(binding, roots[[1L]], roots[[2L]],
                                           roots[[3L]])
  data.frame(
    contract_id = "mv10a_stack_verification_v1",
    verification_order = i, catalog_order = binding$catalog_order,
    stack_id = binding$stack_id, seed = binding$seed,
    homology_dimension = binding$homology_dimension,
    units = length(unique(c(loaded$pairs$first_unit_id,
                            loaded$pairs$second_unit_id))),
    unordered_pairs = nrow(loaded$pairs),
    payload_set_sha256 = loaded$payload_set_sha256,
    pair_axis_sha256 = loaded$pair_axis_sha256,
    payload_binding_passed = loaded$payload_set_sha256 ==
      binding$payload_set_sha256,
    pair_axis_binding_passed = loaded$pair_axis_sha256 ==
      binding$pair_axis_sha256,
    labels_loaded = FALSE, outcomes_loaded = FALSE,
    stringsAsFactors = FALSE
  )
})
verification <- do.call(rbind, verification)
methods <- mv10_method_registry_v1()
distances <- mv10_distance_registry_v1()
literature <- mv10_literature_registry_v1()
analyses <- mv10_analysis_registry_v1(catalog)
owner_acceptance <- data.frame(
  contract_id = "mv10a_owner_acceptance_v1",
  accepted_on = "2026-08-25",
  accepted_artifact = "MV9 corrected four-figure review set",
  owner_statement_verbatim =
    "Yes they look good as long as they are comprehensive enough",
  acceptance_interpretation =
    "figures_look_good_and_are_comprehensive_for_current_scope",
  presentation_gate = "accepted",
  final_manuscript_comprehensiveness_implied = FALSE,
  biological_interpretation_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
view_contract <- data.frame(
  contract_id = "mv10a_view_contract_v1",
  view = c("gene_topology", "cell_topology", "cell_gene_fusion", "external8"),
  mv10b_role = c(
    "current_full124_benchmark", "separate_prior_evidence_no_recompute",
    "prohibited_pending_separate_multiview_gate",
    "descriptive_validation_only_no_clustering_claim"
  ),
  H0_H1_separate = TRUE,
  labels_allowed = FALSE, outcomes_allowed = FALSE,
  stringsAsFactors = FALSE
)
k_contract <- data.frame(
  contract_id = "mv10a_k_contract_v1", candidate_k = .mv10_k_grid,
  all_methods_report_complete_grid = TRUE,
  primary_selection_method = "pam_dissimilarity_v1",
  selection_rule = "smallest_k_within_one_SE_of_maximum_five_seed_mean_ARI",
  silhouette_role = "descriptive_internal_validation_not_selection",
  result_dependent_thresholds = FALSE, labels_allowed = FALSE,
  stringsAsFactors = FALSE
)
output_contract <- data.frame(
  contract_id = "mv10a_output_contract_v1",
  artifact = c(
    "private_sample_partitions", "public_partition_quality",
    "public_seed_stability", "public_primary_k_selection",
    "public_method_agreement", "public_resource_receipt"
  ),
  expected_rows = c(167400L, 1350L, 270L, 2L, 2700L, 1L),
  contains_sample_ids = c(TRUE, rep(FALSE, 5L)),
  public = c(FALSE, rep(TRUE, 5L)),
  labels_allowed = FALSE, outcomes_allowed = FALSE,
  stringsAsFactors = FALSE
)
resource <- data.frame(
  contract_id = "mv10a_resource_contract_v1",
  matrices = 30L, units_per_matrix = 124L, pairs_per_matrix = 7626L,
  authorized_methods = sum(methods$authorized_for_mv10b), k_values = 9L,
  partition_fits = 30L * 5L * 9L,
  private_assignment_rows = 30L * 5L * 9L * 124L,
  workers = 1L, automatic_retries = 0L,
  child_elapsed_cap_seconds = 1800L, process_tree_rss_cap_bytes = 4 * 1024^3,
  private_storage_cap_bytes = 512 * 1024^2,
  execution_state = "closed_until_runner_closure_implementation_committed",
  stringsAsFactors = FALSE
)
implementation_files <- c(
  "R/mv05_benchmark_contract.R",
  "R/mv05n_clustering_gate.R",
  "R/mv08zy_distance_comparison.R",
  "R/mv10_clustering_benchmark.R",
  "scripts/build_mv10a_clustering_benchmark_prefreeze.R"
)
implementation <- data.frame(
  contract_id = "mv10a_implementation_binding_v1",
  implementation_order = seq_along(implementation_files),
  file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
source_files <- c(
  file.path(mv08zy, "mv08zy-artifact-manifest.csv"),
  file.path(mv08zz, "mv08zz-artifact-manifest.csv"),
  file.path(mv09m, "mv09m-artifact-manifest.csv"),
  "docs/specifications/MV05N_LABEL_CLOSED_CLUSTERING_RESOURCE_GATE_SPECIFICATION_V1.md",
  "docs/specifications/MV05Q_LABEL_CLOSED_CLUSTERING_ARTIFACT_PRODUCTION_SPECIFICATION_V1.md",
  "docs/architecture/ADR-024-MV05N-LABEL-CLOSED-CLUSTERING-GATE.md"
)
source_freeze <- data.frame(
  contract_id = "mv10a_source_freeze_v1",
  source_order = seq_along(source_files), artifact = source_files,
  bytes = as.numeric(file.info(source_files)$size),
  sha256 = vapply(source_files, sha, character(1L)),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv10a_clustering_benchmark_prefreeze_v1",
  implementation_head = implementation_head,
  internal_units = 124L, distance_matrices = 30L,
  representations = 3L, homology_dimensions = 2L, seeds = 5L,
  authorized_clustering_methods = 5L, k_grid = "2:10",
  primary_distance = "exact_all_level_landscape_l2_H0_H1_separate",
  primary_representation = "allqc_residual_exact500",
  external_clustering = FALSE, cell_gene_fusion = FALSE,
  labels_used = FALSE, outcomes_used = FALSE, inference_performed = FALSE,
  biological_claims = FALSE, manuscript_claims = FALSE,
  result_dependent_selection = FALSE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv10a_decision_v1",
  decision = "freeze_design_and_inputs_before_execution_implementation",
  mv10b_execution_authorized = FALSE,
  next_required_stage = "commit_runner_independent_closure_and_resource_sentinel",
  label_opening_authorized = FALSE, outcome_evaluation_authorized = FALSE,
  biological_interpretation_authorized = FALSE,
  manuscript_claims_authorized = FALSE,
  stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv10a_validation_v1",
  check_id = c(
    "mv08zy_manifest", "mv08zz_manifest", "mv09m_manifest",
    "owner_presentation_acceptance", "thirty_internal_matrices",
    "three_representations", "five_seeds", "H0_H1_separate",
    "exact500_gene_view", "all_124_units", "all_7626_pairs",
    "all_payload_bindings", "one_pair_axis", "five_authorized_methods",
    "pam_primary", "average_complete_diana_sensitivities",
    "single_linkage_diagnostic", "ward_kmeans_excluded",
    "spectral_hdbscan_deferred", "complete_k2_k10_grid",
    "primary_k_rule_frozen", "silhouette_not_selection",
    "six_distance_dispositions", "landscape_H0_H1_only_authorized",
    "no_H0_H1_combination", "cell_view_separate", "no_fusion",
    "external_no_clustering_claim", "private_sample_ids_only",
    "one_worker_zero_retry", "resource_caps", "implementation_bound",
    "source_evidence_bound", "label_outcome_firewall", "claim_firewall",
    "execution_remains_closed"
  ),
  passed = c(
    TRUE, TRUE, TRUE, owner_acceptance$presentation_gate == "accepted",
    nrow(catalog) == 30L, length(unique(catalog$stack_id)) == 3L,
    identical(sort(unique(as.integer(catalog$seed))), .mv10_required_seeds),
    setequal(catalog$homology_dimension, c("H0", "H1")),
    all(catalog$panel_id == "exact500") &&
      all(catalog$view_kind == "gene_topology_v1"),
    all(verification$units == 124L), all(verification$unordered_pairs == 7626L),
    all(verification$payload_binding_passed &
          verification$pair_axis_binding_passed),
    length(unique(verification$pair_axis_sha256)) == 1L,
    sum(methods$authorized_for_mv10b) == 5L,
    methods$role[methods$method_id == "pam_dissimilarity_v1"] == "primary",
    all(methods$authorized_for_mv10b[methods$method_id %in% c(
      "hclust_average_v1", "hclust_complete_v1", "diana_dissimilarity_v1"
    )]),
    methods$role[methods$method_id == "hclust_single_diagnostic_v1"] ==
      "diagnostic_chaining",
    all(!methods$authorized_for_mv10b[methods$method_id %in% c(
      "ward_arbitrary_distance_v0", "kmeans_distance_matrix_v0"
    )]),
    all(!methods$authorized_for_mv10b[methods$method_id %in% c(
      "spectral_affinity_v0", "hdbscan_dissimilarity_v0"
    )]),
    identical(k_contract$candidate_k, .mv10_k_grid),
    all(k_contract$selection_rule ==
          "smallest_k_within_one_SE_of_maximum_five_seed_mean_ARI"),
    all(k_contract$silhouette_role ==
          "descriptive_internal_validation_not_selection"),
    nrow(distances) == 6L,
    sum(distances$authorized_for_mv10b) == 2L,
    all(!analyses$H0_H1_combined),
    view_contract$mv10b_role[view_contract$view == "cell_topology"] ==
      "separate_prior_evidence_no_recompute",
    !contract$cell_gene_fusion,
    !contract$external_clustering,
    all(!output_contract$public[output_contract$contains_sample_ids]),
    resource$workers == 1L && resource$automatic_retries == 0L,
    resource$process_tree_rss_cap_bytes == 4 * 1024^3 &&
      resource$private_storage_cap_bytes == 512 * 1024^2,
    all(file.exists(implementation$file)) &&
      all(vapply(implementation$file, sha, character(1L)) ==
            implementation$sha256),
    all(file.exists(source_freeze$artifact)),
    !contract$labels_used && !contract$outcomes_used,
    !contract$biological_claims && !contract$manuscript_claims,
    !decision$mv10b_execution_authorized
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV10-A prefreeze validation failed")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
artifacts <- list(
  "mv10a-contract.csv" = contract,
  "mv10a-owner-acceptance.csv" = owner_acceptance,
  "mv10a-internal-stack-bindings.csv" = catalog,
  "mv10a-stack-verification.csv" = verification,
  "mv10a-analysis-registry.csv" = analyses,
  "mv10a-method-registry.csv" = methods,
  "mv10a-distance-registry.csv" = distances,
  "mv10a-literature-registry.csv" = literature,
  "mv10a-view-contract.csv" = view_contract,
  "mv10a-k-contract.csv" = k_contract,
  "mv10a-output-contract.csv" = output_contract,
  "mv10a-resource-contract.csv" = resource,
  "mv10a-implementation-bindings.csv" = implementation,
  "mv10a-source-freeze.csv" = source_freeze,
  "mv10a-decision.csv" = decision,
  "mv10a-validation.csv" = validation
)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c(
  "# MV10-A label-closed clustering benchmark prefreeze", "",
  "Jonah accepted the corrected MV9 presentation for its current robustness",
  "scope. MV10-A binds 30 internal exact-landscape matrices across three gene",
  "representations, five seeds, and separate H0/H1. PAM is primary; average,",
  "complete, single-linkage diagnostic, and DIANA form the authorized method",
  "grid at K=2:10. Labels, outcomes, execution, fusion, and claims remain closed."
), file.path(output, "MV10A_CLUSTERING_BENCHMARK_PREFREEZE_2026-08-25.md"))
files <- sort(setdiff(list.files(output), "mv10a-artifact-manifest.csv"))
manifest <- data.frame(
  contract_id = "mv10a_artifact_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(output, files))$size),
  sha256 = vapply(file.path(output, files), sha, character(1L)),
  stringsAsFactors = FALSE
)
atomic(manifest, file.path(output, "mv10a-artifact-manifest.csv"))
message("Built MV10-A clustering benchmark prefreeze; checks=36")
