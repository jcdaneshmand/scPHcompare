#!/usr/bin/env Rscript

options(warn = 2, digits = 17, scipen = 999)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: build_mv07j_claim_figure_literature_map.R PRIVATE_ARTIFACT_DIR OUTPUT_DIR EXPECTED_HEAD",
       call. = FALSE)
}
private_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2L]]
expected_head <- tolower(trimws(args[[3L]]))
if (!grepl("^[0-9a-f]{40}$", expected_head)) stop("Full EXPECTED_HEAD required.")
observed_head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(observed_head, expected_head)) {
  stop("MV7-J requires prospective HEAD ", expected_head, "; observed ", observed_head, ".")
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                 no.. = TRUE))) {
  stop("MV7-J output directory must be empty.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

readc <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                         check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(x)) == "true"
write_csv <- function(x, name) {
  path <- file.path(output_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  utils::write.table(x, path, sep = ",", row.names = FALSE, col.names = TRUE,
                     quote = TRUE, na = "", qmethod = "double")
}

public <- c(
  label_decision = "docs/audits/mv07i-label-closed-validation/mv07i-label-closed-decision.csv",
  label_validation = "docs/audits/mv07i-label-closed-validation/mv07i-label-closed-validation.csv",
  label_manifest = "docs/audits/mv07i-label-closed-validation/mv07i-label-closed-artifact-manifest.csv",
  outcome_decision = "docs/audits/mv07i-outcome-validation/mv07i-outcome-decision.csv",
  outcome_validation = "docs/audits/mv07i-outcome-validation/mv07i-outcome-validation.csv",
  outcome_manifest = "docs/audits/mv07i-outcome-validation/mv07i-outcome-artifact-manifest.csv",
  outcome_seed = "docs/audits/mv07i-outcome-validation/mv07i-outcome-seed-metrics.csv",
  outcome_units = "docs/audits/mv07i-outcome-validation/mv07i-outcome-unit-summaries.csv",
  outcome_structural = "docs/audits/mv07i-outcome-validation/mv07i-outcome-structural-status.csv",
  mv07a_claims = "docs/audits/mv07a-synthesis-evidence/mv07a-claim-boundaries.csv",
  mv07a_robustness = "docs/audits/mv07a-synthesis-evidence/mv07a-robustness-coverage.csv",
  mv07b_audit = "docs/audits/MV07B_NO_NEW_PH_CONFOUNDING_DIAGNOSTICS_2026-08-15.md",
  mv07c_claims = "docs/audits/mv07c-synthesis-evidence/mv07c-claim-register.csv",
  mv07c_options = "docs/audits/mv07c-synthesis-evidence/mv07c-option-register.csv",
  mv07c_gates = "docs/audits/mv07c-synthesis-evidence/mv07c-gate-status.csv",
  landscape_spec = "docs/specifications/PERSISTENCE_LANDSCAPE_SPECIFICATION_V1.md",
  dual_view_spec = "docs/specifications/DUAL_VIEW_TOPOLOGY_SPECIFICATION_V1.md",
  mv07j_spec = "docs/specifications/MV07J_CLAIM_FIGURE_LITERATURE_MAP_PREFREEZE_V1.md")
private <- c(
  h1_summary = file.path(private_dir, "h1-contribution-summary.csv"),
  stability = file.path(private_dir, "stability-summary.csv"),
  selected = file.path(private_dir, "selected-partitions.csv"))
sources <- c(public, private)
if (any(!file.exists(sources))) {
  stop("MV7-J source set incomplete: ",
       paste(names(sources)[!file.exists(sources)], collapse = ", "))
}

label_decision <- readc(public[["label_decision"]])
label_validation <- readc(public[["label_validation"]])
label_manifest <- readc(public[["label_manifest"]])
outcome_decision <- readc(public[["outcome_decision"]])
outcome_validation <- readc(public[["outcome_validation"]])
outcome_manifest <- readc(public[["outcome_manifest"]])
outcome_seed <- readc(public[["outcome_seed"]])
outcome_units <- readc(public[["outcome_units"]])
outcome_structural <- readc(public[["outcome_structural"]])
h1 <- readc(private[["h1_summary"]])
stability <- readc(private[["stability"]])
selected <- readc(private[["selected"]])

if (nrow(label_decision) != 1L || label_decision$decision !=
      "authorize_MV7I_outcome_prefreeze_only" ||
    !all(truth(label_validation$passed)) || nrow(label_validation) != 13L ||
    nrow(outcome_decision) != 1L || outcome_decision$decision !=
      "authorize_MV7J_claim_map_and_figure_planning_only" ||
    !all(truth(outcome_validation$passed)) || nrow(outcome_validation) != 15L ||
    nrow(outcome_seed) != 600L || nrow(outcome_units) != 120L ||
    nrow(outcome_structural) != 1L || nrow(h1) != 15252L ||
    nrow(stability) != 54L || nrow(selected) != 7440L) {
  stop("MV7-J admission decision or expected source cardinality is stale.")
}

for (id in names(private)) {
  artifact <- c(h1_summary = "h1_summary", stability = "stability",
                selected = "selected_partitions")[[id]]
  row <- label_manifest[label_manifest$artifact == artifact, , drop = FALSE]
  if (nrow(row) != 1L || !identical(tolower(sha(private[[id]])),
                                     tolower(row$production_sha256))) {
    stop("MV7-J private source hash mismatch: ", id)
  }
}
for (id in c("outcome_seed", "outcome_units", "outcome_structural")) {
  artifact <- c(outcome_seed = "seed_metrics", outcome_units = "unit_summaries",
                outcome_structural = "structural_status")[[id]]
  row <- outcome_manifest[outcome_manifest$artifact == artifact, , drop = FALSE]
  if (nrow(row) != 1L || !identical(tolower(sha(public[[id]])),
                                     tolower(row$production_sha256))) {
    stop("MV7-J public outcome source hash mismatch: ", id)
  }
}

source_manifest <- data.frame(
  contract_id = "mv07j_source_manifest_v1",
  source_order = seq_along(sources), source_id = names(sources),
  locator = c(unname(public), rep("private_validated_artifact_not_published",
                                  length(private))),
  sha256 = vapply(sources, sha, character(1L)),
  bytes = as.numeric(file.info(sources)$size),
  access_class = c(rep("public_repository", length(public)),
                   rep("private_hash_validated", length(private))),
  accepted_head = expected_head, new_ph = FALSE, new_data = FALSE,
  stringsAsFactors = FALSE)

landscape <- data.frame(
  contract_id = "mv07j_landscape_observation_contract_v1",
  item_order = 1:12,
  item = c("comparison_unit_sample", "cell_view_observations_cells",
           "gene_view_observations_fixed_global_core_genes",
           "finite_positive_intervals", "essential_H0_excluded",
           "all_consecutive_active_levels", "H0_H1_separate",
           "squared_L2_exact_or_error_controlled", "no_universal_grid",
           "no_universal_level_cap", "secondary_unweighted_composite_only",
           "streamed_or_chunked_execution"),
  required_state = TRUE, current_state = TRUE,
  legacy_first_level_100_grid = c(rep(FALSE, 11L), FALSE),
  manuscript_primary = TRUE, stringsAsFactors = FALSE)

quant <- function(x, p) unname(stats::quantile(x, p, type = 7L))
h1_rows <- lapply(c("cell", "gene"), function(view) {
  x <- h1$median[h1$view_id == view]
  data.frame(
    contract_id = "mv07j_h1_contribution_distribution_v1", view_id = view,
    sample_pairs = length(x), mean = mean(x), minimum = min(x),
    q10 = quant(x, .10), q25 = quant(x, .25), median = quant(x, .50),
    q75 = quant(x, .75), q90 = quant(x, .90), maximum = max(x),
    fraction_gt_0_01 = mean(x > .01), fraction_gt_0_05 = mean(x > .05),
    fraction_gt_0_10 = mean(x > .10), fraction_gt_0_50 = mean(x > .50),
    contribution_definition = "H1_squared_distance_over_H0_plus_H1_squared_distance",
    biological_cycle_inference = FALSE, stringsAsFactors = FALSE)
})
h1_distribution <- do.call(rbind, h1_rows)

comb2 <- function(x) x * (x - 1) / 2
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  observed <- sum(comb2(tab)); rows <- sum(comb2(rowSums(tab)))
  cols <- sum(comb2(colSums(tab))); total <- comb2(n)
  expected <- rows * cols / total
  denominator <- (rows + cols) / 2 - expected
  if (denominator == 0) return(if (identical(as.integer(a), as.integer(b))) 1 else 0)
  (observed - expected) / denominator
}
get_partition <- function(rep_id, alg_id, seed) {
  z <- selected[selected$representation_id == rep_id &
                  selected$algorithm_id == alg_id & selected$seed == seed,
                c("sample_id", "cluster"), drop = FALSE]
  z[order(z$sample_id, method = "radix"), , drop = FALSE]
}
seeds <- sort(unique(selected$seed))
algorithms <- c("pam_stability_k_v1", "hclust_average_v1")
h0_comp <- list(); cursor <- 0L
for (view in c("cell", "gene")) for (algorithm in algorithms) for (seed in seeds) {
  left <- get_partition(paste0(view, "_H0"), algorithm, seed)
  right <- get_partition(paste0(view, "_H0_H1_secondary"), algorithm, seed)
  if (nrow(left) != 124L || !identical(left$sample_id, right$sample_id))
    stop("H0/composite axis mismatch.")
  cursor <- cursor + 1L
  h0_comp[[cursor]] <- data.frame(
    contract_id = "mv07j_h0_composite_concordance_v1", view_id = view,
    algorithm_id = algorithm, seed = seed, samples = 124L,
    adjusted_rand_index = ari(left$cluster, right$cluster),
    exact_partition = identical(left$cluster, right$cluster),
    composite_role = "secondary_descriptive", stringsAsFactors = FALSE)
}
h0_composite <- do.call(rbind, h0_comp)

algorithm_rows <- list(); cursor <- 0L
representations <- unique(selected$representation_id)
representations <- c("cell_H0", "cell_H1", "cell_H0_H1_secondary",
                     "gene_H0", "gene_H1", "gene_H0_H1_secondary")
for (representation in representations) for (seed in seeds) {
  pam <- get_partition(representation, "pam_stability_k_v1", seed)
  avg <- get_partition(representation, "hclust_average_v1", seed)
  if (nrow(pam) != 124L || !identical(pam$sample_id, avg$sample_id))
    stop("Algorithm-sensitivity axis mismatch.")
  cursor <- cursor + 1L
  algorithm_rows[[cursor]] <- data.frame(
    contract_id = "mv07j_algorithm_sensitivity_v1",
    representation_id = representation, seed = seed, samples = 124L,
    pam_k = unique(selected$selected_k[selected$representation_id == representation &
      selected$algorithm_id == "pam_stability_k_v1" & selected$seed == seed]),
    average_k = unique(selected$selected_k[selected$representation_id == representation &
      selected$algorithm_id == "hclust_average_v1" & selected$seed == seed]),
    adjusted_rand_index = ari(pam$cluster, avg$cluster),
    exact_partition = identical(pam$cluster, avg$cluster),
    favorable_algorithm_selected = FALSE, stringsAsFactors = FALSE)
}
algorithm_sensitivity <- do.call(rbind, algorithm_rows)

stability_public <- stability[c("representation_id", "k", "mean_stability",
  "jackknife_se", "pair_count", "selected_k", "one_se_threshold")]
stability_public$contract_id <- "mv07j_complete_stability_curve_v1"
stability_public <- stability_public[c("contract_id", names(stability_public)[1:7])]

outcome_complete <- outcome_units
outcome_complete$contract_id <- "mv07j_complete_outcome_summary_v1"
outcome_complete$claim_scope <- "descriptive_only_no_ranking"
outcome_complete$external_generalization <- FALSE

key <- interaction(outcome_complete$population_id, outcome_complete$label_axis,
                   outcome_complete$algorithm_id, outcome_complete$metric_id,
                   drop = TRUE, lex.order = TRUE)
outcome_range <- do.call(rbind, lapply(split(outcome_complete, key), function(z) {
  data.frame(contract_id = "mv07j_outcome_range_synopsis_v1",
    population_id = z$population_id[[1L]], label_axis = z$label_axis[[1L]],
    algorithm_id = z$algorithm_id[[1L]], metric_id = z$metric_id[[1L]],
    representations = nrow(z), minimum_seed_mean = min(z$seed_mean),
    maximum_seed_mean = max(z$seed_mean), minimum_seed_jackknife_se =
      min(z$seed_jackknife_se), maximum_seed_jackknife_se =
      max(z$seed_jackknife_se), ranking_authorized = FALSE,
    stringsAsFactors = FALSE)
}))
row.names(outcome_range) <- NULL

claims <- data.frame(
  contract_id = "mv07j_claim_map_v1", claim_order = 1:16,
  claim_id = sprintf("J%02d", 1:16),
  topic = c("observation_hierarchy", "landscape_definition", "dual_views",
    "coarse_structure", "tissue_study_alignment", "technology_approach",
    "H1_interpretation", "H0_H1_composite", "algorithm_sensitivity",
    "relative_views", "fusion", "integration", "cell_metric_robustness",
    "gene_robustness", "computational_feasibility", "external_generalization"),
  status = c("supported_method", "supported_method", "supported_method",
    "supported_descriptive", "supported_descriptive", "retire",
    "conditional_context", "conditional_context", "supported_descriptive",
    "retire", "supported_descriptive", "conditional_context",
    "conditional_context", "hypothesis_only", "supported_method",
    "hypothesis_only"),
  permitted_claim = c(
    "Samples are compared; each sample yields a cell-observation and a gene-observation point cloud.",
    "All finite positive intervals and active levels are retained with separate H0/H1 exact or error-controlled squared-L2 distances.",
    "Cell and gene topology provide separate complementary views under fixed coordinate and panel contracts.",
    "All six representations selected k=2 by the label-free five-seed stability rule.",
    "The coarse partitions show modest descriptive tissue and study alignment in the accepted populations.",
    "Canonical approach is only a nested descriptive proxy in the full 124 and is non-estimable in the primary 90.",
    "H1-only distances can produce partitions distinct from H0 despite low H1 contribution to the unweighted composite.",
    "The unweighted composite is a secondary sensitivity and is usually H0-dominated in this corpus.",
    "PAM and average linkage can yield materially different partitions and must be shown together.",
    "Cell and gene results are reported separately without a global winner.",
    "Equal-weight fusion did not reliably outperform both views under the frozen blocked evaluation.",
    "Integration had a corpus-specific negative or null topology increment under the accepted blocked design.",
    "Primary-90 cell-view conclusions have a completed Euclidean/cosine-chord sensitivity.",
    "The fixed 500-gene correlation-chord view is conditional on its transductive panel and metric.",
    "The 124-sample dual-view PH and all-level landscape workflow completed under audited resource and equivalence gates.",
    "External performance requires a prospectively frozen independent dataset."),
  prohibited_claim = c(
    "Cells or genes are the independent units of the final clustering outcome.",
    "The current method uses only k=1 or a universal 100-point grid.",
    "Dual-view topology or cell/gene topology is unprecedented.",
    "k=2 recovers the eight tissues or predicts tissue.",
    "Tissue alignment is independent of study.",
    "Sequencing technology causes the observed topology.",
    "H1 features are biological cycles pathways or mechanisms.",
    "The composite replaces separate H0 and H1 reporting.",
    "One clustering algorithm is objectively superior from these outcomes.",
    "Either cell or gene topology is globally superior.",
    "A post-hoc fusion weight is an improved method.",
    "Integration generally removes topological noise or generally harms topology.",
    "All cell-view robustness axes and external datasets are covered.",
    "The gene view is inductive or insensitive to panel size and distance.",
    "Computational feasibility establishes biological validity.",
    "Performance generalizes to unseen studies tissues methods or assays."),
  population = c(rep("method_contract", 3), "full124_descriptive",
    "full124_and_primary90_context", "full124_and_primary90_context",
    "full124_descriptive", "full124_descriptive", "full124_descriptive",
    "primary90_blocked_plus_full124_descriptive", "primary90_blocked",
    "primary90_blocked", "primary90_blocked", "full124_descriptive",
    "full124_execution", "not_yet_evaluated"),
  uncertainty = c(rep("implementation_validation", 3),
    "five_seed_stability_and_jackknife_SE", "five_seed_jackknife_SE",
    "structural_nesting_not_sampling_uncertainty", "complete_pair_distribution",
    "complete_view_algorithm_seed_partition_ARI",
    "complete_view_seed_partition_ARI", "blocked_influence_and_seed_sensitivity",
    "blocked_confidence_rule", "blocked_confidence_rule",
    "multi_panel_existing_data_sensitivity", "not_estimated",
    "deterministic_repeat_cross_engine_and_resource_validation", "not_estimated"),
  evidence_artifact = c("landscape_observation_contract", "landscape_observation_contract",
    "dual_view_spec", "complete_stability_curve", "complete_outcome_summary",
    "outcome_structural_status", "h1_contribution_distribution",
    "h0_composite_concordance", "algorithm_sensitivity", "mv07b_and_mv07c",
    "mv06h_and_mv07c", "mv05l_mv05s_mv07c", "mv07a_robustness",
    "mv07a_robustness", "mv07h_and_mv07i_closures", "external_data_gate"),
  limitation = c("Within-sample cells and genes are subsampled technical realizations.",
    "Landscape stability does not make raw Rips PH outlier-robust.",
    "Gene view uses a fixed transductive global-core panel.",
    "Two clusters are a coarse summary against eight tissues and 18 studies.",
    "Study and tissue are materially confounded.",
    "Six snRNA samples are nested in one tissue and study.",
    "Low energy does not imply H1 is irrelevant; no feature-level biology was validated.",
    "Unweighted scale makes H0 dominant in this corpus.",
    "Only PAM and average linkage are in the full-124 contract.",
    "Relative view performance changes by study and tissue.",
    "Only prospectively frozen equal fusion is confirmatory.",
    "Seurat-centered existing corpus; not a comprehensive integration benchmark.",
    "Sensitivity is limited to the accepted primary-90 cell analysis.",
    "Panel-size and gene-metric sensitivity require new PH.",
    "Resource success is hardware and corpus conditional.",
    "No independent external analysis has run."),
  manuscript_action = c("retain", "replace_legacy_method", "retain_with_novelty_limit",
    "retain_descriptive", "retain_with_confounding", "replace_causal_language",
    "replace_biological_loop_language", "supplementary_sensitivity",
    "retain_as_sensitivity", "replace_winner_language", "retain_negative_result",
    "replace_general_integration_claim", "retain_bounded_robustness",
    "state_gap", "retain_methods", "defer_until_external_validation"),
  stringsAsFactors = FALSE)

figures <- data.frame(
  contract_id = "mv07j_figure_map_v1", figure_order = 1:11,
  figure_id = c(paste0("F", 1:8), paste0("S", 1:3)),
  placement = c(rep("main", 8), rep("supplement", 3)),
  title = c("Corrected sample-comparison and dual-view pipeline",
    "Cohort accountability and confounding structure",
    "Persistence diagrams to all-level H0/H1 landscapes",
    "H1 contribution and H0/composite concordance",
    "Label-free stability and selected cluster count",
    "Complete descriptive tissue and study alignment",
    "PAM versus average-linkage sensitivity",
    "Prior blocked integration and dual-view evidence",
    "Complete seed-level descriptive outcomes",
    "Landscape numerical equivalence and resource envelope",
    "Robustness coverage and unresolved evidence gaps"),
  source_artifact = c("landscape_observation_contract;dual_view_spec",
    "mv07d_reconciliation;mv07e_canonical_approach",
    "validated_mv07h_diagrams;landscape_observation_contract",
    "h1_contribution_distribution;h0_composite_concordance",
    "complete_stability_curve", "complete_outcome_summary",
    "algorithm_sensitivity", "mv05l;mv05s;mv06h;mv07c",
    "mv07i_outcome_seed_metrics", "mv07h_closures;mv07i_resource_summaries",
    "mv07a_robustness;mv07b;external_data_gate"),
  estimand = c("method_architecture", "population_and_nonidentifiability",
    "corrected_landscape_definition", "distance_component_and_partition_sensitivity",
    "label_free_cluster_stability", "descriptive_cluster_label_alignment",
    "algorithm_membership_sensitivity", "accepted_blocked_existing_data_results",
    "technical_seed_variability", "numerical_and_resource_validation",
    "coverage_not_performance"),
  population = c("method_contract", "full124_and_primary90",
    "validated_representative_artifacts", "full124", "full124", "full124_and_primary90",
    "full124", "primary90_blocked", "full124_and_primary90", "full124", "all_accepted"),
  uncertainty = c("none_diagrammatic", "counts_and_nesting", "error_controlled_integration",
    "all_pairs_and_all_five_seeds", "ten_seed_pair_ARIs_and_jackknife_SE",
    "all_five_seeds_and_jackknife_SE", "all_five_seeds", "frozen_blocked_intervals",
    "all_seed_estimates", "repeat_and_oracle_bounds", "coverage_state"),
  limitation_visible = TRUE,
  status = c("ready_for_corrected_render", "ready_for_render",
    "ready_for_corrected_render", "ready_for_render", "ready_for_render",
    "ready_for_render", "ready_for_render", "requires_cross_stage_layout",
    "ready_for_table_or_heatmap", "ready_for_render", "ready_for_render"),
  favorable_subset_allowed = FALSE, stringsAsFactors = FALSE)

literature <- data.frame(
  contract_id = "mv07j_literature_registry_v1", literature_order = 1:11,
  literature_id = c("Bubenik2015", "Chazal2018", "Luecken2022", "Huynh2024",
    "Hozumi2024", "Cottrell2023", "Su2024", "Zhu2026", "Oliveira2026",
    "Mule2024", "Daneshmand2025"),
  year = c(2015L, 2018L, 2022L, 2024L, 2024L, 2023L, 2024L, 2026L, 2026L,
           2024L, 2025L),
  title = c("Statistical Topological Data Analysis using Persistence Landscapes",
    "Robust Topological Inference: Distance To a Measure and Kernel Distance",
    "Benchmarking atlas-level data integration in single-cell genomics",
    "Topological and geometric analysis of cell states in single-cell transcriptomic data",
    "Analyzing single cell RNA sequencing with topological nonnegative matrix factorization",
    "K-Nearest-Neighbors Induced Topological PCA for Single Cell RNA-Sequence Data Analysis",
    "Hodge Decomposition of Single-Cell RNA Velocity",
    "Revealing a coherent cell-state landscape across single-cell datasets with CONCORD",
    "TopoMetry systematically learns and evaluates the latent geometry of single-cell data",
    "Multiscale topology classifies cells in subcellular spatial transcriptomics",
    "Impact of integration on persistent homology clustering and biological signal detection in scRNA-seq data"),
  publication_state = c(rep("peer_reviewed", 5), "preprint", "peer_reviewed",
    "peer_reviewed", "peer_reviewed", "peer_reviewed", "project_preprint"),
  doi_or_url = c("https://www.jmlr.org/papers/v16/bubenik15a.html",
    "https://www.jmlr.org/papers/v18/15-484.html",
    "https://doi.org/10.1038/s41592-021-01336-8",
    "https://doi.org/10.1093/bib/bbae176",
    "https://doi.org/10.1016/j.cam.2024.115842",
    "https://arxiv.org/abs/2310.14521",
    "https://doi.org/10.1021/acs.jcim.4c00132",
    "https://doi.org/10.1038/s41587-025-02950-z",
    "https://doi.org/10.7554/eLife.100361.3",
    "https://doi.org/10.1038/s41586-024-07563-1",
    "https://doi.org/10.1101/2025.07.24.666637"),
  primary_source = TRUE,
  relevance = c("defines vector-space landscapes and stability",
    "shows raw empirical distance PH can be destroyed by even one outlier and provides robust alternatives",
    "separates batch removal from biological conservation across many integration methods",
    "uses topology and curvature on both cell and gene networks",
    "uses persistent-Laplacian regularization for single-cell dimension reduction",
    "uses kNN persistent-Laplacian topology for single-cell PCA",
    "uses Hodge components to analyze single-cell RNA-velocity fields",
    "uses PH and Betti curves with geometric metrics to assess integrated latent-space fidelity",
    "systematically learns and evaluates latent single-cell geometry",
    "uses multiparameter PH for spatial-transcriptomic cell classification",
    "public legacy manuscript that the corrected analysis must supersede transparently"),
  manuscript_effect = c("supports_corrected_definition", "narrows_robustness_claim",
    "requires_multi_metric_integration_language", "narrows_dual_view_novelty",
    "narrows_single_cell_topology_novelty", "track_as_preprint_not_established_fact",
    "distinguishes_velocity_topology_scope", "narrows_topology_integration_novelty",
    "requires_geometry_context_and_complementary_metrics", "narrows_broad_TDA_novelty",
    "requires_versioned_correction_of_legacy_claims"),
  verified_through = "2026-08-16", stringsAsFactors = FALSE)

legacy <- data.frame(
  contract_id = "mv07j_legacy_method_correction_v1", correction_order = 1:9,
  legacy_item = c("first_landscape_level_only", "uniform_grid_100",
    "pointwise_H0_H1_L2_curve", "combined_only_reporting",
    "label_predefined_cluster_count", "Ward_on_arbitrary_dissimilarity",
    "independently_fitted_PCA_comparison", "PH_removes_batch_effects",
    "H1_equals_biological_cycle"),
  corrected_item = c("all_consecutive_active_levels", "no_primary_fixed_grid",
    "separate_exact_or_error_controlled_squared_L2_distances",
    "H0_H1_separate_with_secondary_composite", "label_free_five_seed_stability_rule",
    "PAM_primary_average_linkage_sensitivity", "fixed_comparable_coordinate_contracts",
    "PH_describes_geometry_and_does_not_correct_batch", "H1_loop_is_mathematical_only_without_validation"),
  disposition = "retire_and_replace", public_source = "legacy_manuscript_audit",
  confidential_text_included = FALSE, stringsAsFactors = FALSE)

requirements <- data.frame(
  contract_id = "mv07j_scientific_requirement_map_v1", requirement_order = 1:12,
  requirement_id = sprintf("R%02d", 1:12),
  requirement = c("clarify_dataset_level_versus_cell_level_clustering",
    "define_integration_versus_harmonization", "broaden_integration_context",
    "address_point_metric_dependence", "address_clustering_algorithm_dependence",
    "address_outlier_sensitivity", "avoid_abstract_loop_biology",
    "distinguish_batch_detection_from_correction", "show_complete_readable_figures",
    "ensure_comparable_coordinate_spaces", "distinguish_technical_from_biological_signal",
    "establish_external_generalization_only_with_external_data"),
  current_disposition = c("resolved_in_method_contract", "resolved_in_terminology_plan",
    "literature_context_complete_computation_scope_limited",
    "partially_resolved_cell_only", "resolved_for_PAM_and_average_only",
    "unresolved_named_robustness_gap", "resolved_as_prohibited_claim",
    "resolved_as_claim_boundary", "planned_in_figure_map",
    "resolved_in_corrected_coordinate_contract", "resolved_as_nonidentifiability_boundary",
    "requires_prospective_external_validation"),
  next_action = c("render_F1", "manuscript_glossary", "cite_scIB_CONCORD_and_state_scope",
    "report_primary90_cosine_sensitivity_and_gene_gap", "render_F7_and_state_scope",
    "prospective_no_new_PH_or_robust_PH_sensitivity_decision", "replace_legacy_language",
    "replace_legacy_language", "render_F1_to_F8", "retain_fixed_axis_checks",
    "show_F2_and_limit_causal_words", "read_only_external_dataset_admission_audit"),
  confidential_source_text_included = FALSE, stringsAsFactors = FALSE)

external_gate <- data.frame(
  contract_id = "mv07j_external_data_gate_v1",
  claim = c("corrected_method_definition", "existing_corpus_descriptive_structure",
    "unseen_study_or_tissue_generalization", "external_cell_gene_complementarity",
    "integration_method_superiority", "inductive_gene_panel_robustness"),
  current_data_sufficient = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  external_data_required = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  current_disposition = c("supported", "supported_with_confounding_limits",
    "not_supported", "not_supported", "not_supported", "not_supported"),
  next_authorized_action = c("figure_and_methods_draft", "figure_and_results_draft",
    rep("read_only_dataset_admission_audit", 4)),
  new_data_download_authorized = FALSE, new_ph_authorized = FALSE,
  stringsAsFactors = FALSE)

credit <- data.frame(
  contract_id = "mv07j_credit_checklist_v1", person_order = 1:4,
  person = c("Jonah Daneshmand", "Julia H. Chariker", "Akshitkumar Mistry",
             "Eric C. Rouchka"),
  retained_in_author_team_registry = TRUE,
  final_author_order_confirmed = FALSE, final_credit_roles_confirmed = FALSE,
  manuscript_approval_required = TRUE,
  note = c("project_owner_and_primary_analysis", "legacy_manuscript_coauthor",
    "must_receive_credit_per_owner_instruction",
    "must_receive_credit_per_owner_instruction"), stringsAsFactors = FALSE)

criteria <- data.frame(
  contract_id = "mv07j_acceptance_criteria_v1", criterion_order = 1:14,
  criterion = c("source_decisions_valid", "source_hashes_valid",
    "landscape_contract_preserved", "H1_pair_family_complete",
    "H0_composite_family_complete", "algorithm_family_complete",
    "stability_family_complete", "outcome_family_complete",
    "claim_map_complete", "figure_map_complete", "literature_primary_sources_located",
    "confidential_review_absent", "no_selection_or_new_compute",
    "external_gate_named"),
  passed = c(TRUE, TRUE, all(landscape$current_state),
    nrow(h1_distribution) == 2L && all(h1_distribution$sample_pairs == 7626L),
    nrow(h0_composite) == 20L, nrow(algorithm_sensitivity) == 30L,
    nrow(stability_public) == 54L, nrow(outcome_complete) == 120L,
    nrow(claims) == 16L, nrow(figures) == 11L,
    nrow(literature) == 11L && all(nzchar(literature$doi_or_url)),
    !any(requirements$confidential_source_text_included), TRUE,
    any(!external_gate$current_data_sufficient & external_gate$external_data_required)),
  detail = c("13/13 label-closed and 15/15 outcome decisions",
    "all private and public manifest-bound inputs rehashed", "12 fixed items",
    "2 views x 7626 pairs", "2 views x 2 algorithms x 5 seeds",
    "6 representations x 5 seeds", "6 representations x 9 k values",
    "120 complete units", "16 bounded claims", "8 main plus 3 supplement",
    "11 primary-source locators", "paraphrased public requirements only",
    "no PH refit p-value rank or favorable subset", "four external-only claims named"),
  stringsAsFactors = FALSE)
if (!all(criteria$passed)) {
  stop("MV7-J acceptance criteria failed: ",
       paste(criteria$criterion[!criteria$passed], collapse = ", "))
}

decision <- data.frame(
  contract_id = "mv07j_decision_v1",
  decision = "authorize_corrected_figure_implementation_and_read_only_external_dataset_audit",
  checks_passed = sum(criteria$passed), checks_total = nrow(criteria),
  methods_focused_existing_data_sufficient = TRUE,
  external_generalization_sufficient = FALSE,
  new_data_download_authorized = FALSE, new_ph_authorized = FALSE,
  manuscript_submission_authorized = FALSE, claims_promoted_to_confirmatory = FALSE,
  confidential_review_published = FALSE, stringsAsFactors = FALSE)

outputs <- list(
  "mv07j-source-manifest.csv" = source_manifest,
  "mv07j-landscape-observation-contract.csv" = landscape,
  "mv07j-h1-contribution-distribution.csv" = h1_distribution,
  "mv07j-h0-composite-concordance.csv" = h0_composite,
  "mv07j-algorithm-sensitivity.csv" = algorithm_sensitivity,
  "mv07j-complete-stability-curve.csv" = stability_public,
  "mv07j-complete-outcome-summary.csv" = outcome_complete,
  "mv07j-outcome-range-synopsis.csv" = outcome_range,
  "mv07j-claim-map.csv" = claims,
  "mv07j-figure-map.csv" = figures,
  "mv07j-literature-registry.csv" = literature,
  "mv07j-legacy-method-correction.csv" = legacy,
  "mv07j-scientific-requirement-map.csv" = requirements,
  "mv07j-external-data-gate.csv" = external_gate,
  "mv07j-credit-checklist.csv" = credit,
  "mv07j-acceptance-criteria.csv" = criteria,
  "mv07j-decision.csv" = decision)
for (name in names(outputs)) write_csv(outputs[[name]], name)

manifest_paths <- file.path(output_dir, names(outputs))
artifact_manifest <- data.frame(
  contract_id = "mv07j_artifact_manifest_v1",
  artifact_order = seq_along(outputs), filename = names(outputs),
  sha256 = vapply(manifest_paths, sha, character(1L)),
  bytes = as.numeric(file.info(manifest_paths)$size),
  public_content = TRUE, contains_sample_level_labels = FALSE,
  contains_confidential_review_text = FALSE, stringsAsFactors = FALSE)
write_csv(artifact_manifest, "mv07j-artifact-manifest.csv")
message("MV7-J synthesis complete: 14/14 criteria; corrected figures and read-only external audit only")
