# MV7-E metadata provenance and 124-sample descriptive-estimand helpers.

.mv07e_require_columns <- function(x, required, label) {
  if (!is.data.frame(x) || !all(required %in% names(x))) {
    stop(label, " is missing required columns: ",
         paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

.mv07e_truth <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

#' Reproduce Seurat's historical feature-name normalization for availability.
mv07e_seurat_feature_ids_v1 <- function(ids) {
  ids <- as.character(ids)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids))) {
    stop("Feature identifiers must be non-empty and non-missing.", call. = FALSE)
  }
  gsub("_", "-", ids, fixed = TRUE)
}

#' Validate the accession-level resolution table for all disputed approaches.
mv07e_validate_accession_evidence_v1 <- function(evidence) {
  required <- c("sample_id", "srs", "srr", "geo_sample", "geo_series",
    "public_approach", "historical_heuristic_approach", "canonical_approach",
    "official_method_indicator", "official_url", "accessed_date")
  .mv07e_require_columns(evidence, required, "accession evidence")
  if (nrow(evidence) != 16L || anyDuplicated(evidence$sample_id) ||
      any(evidence$canonical_approach != evidence$public_approach) ||
      any(evidence$canonical_approach == evidence$historical_heuristic_approach) ||
      sum(evidence$canonical_approach == "scRNA-seq") != 14L ||
      sum(evidence$canonical_approach == "snRNA-seq") != 2L ||
      any(!grepl("^https://www[.]ncbi[.]nlm[.]nih[.]gov/geo/", evidence$official_url)) ||
      any(!grepl("^GSM[0-9]+$", evidence$geo_sample))) {
    stop("Accession evidence does not resolve the fixed 16-row conflict set.",
         call. = FALSE)
  }
  invisible(evidence)
}

#' Resolve approach labels without using expression heuristics.
mv07e_resolve_approaches_v1 <- function(reconciliation, retained, evidence) {
  .mv07e_require_columns(reconciliation,
    c("sample_id", "study", "tissue", "approach_public",
      "approach_historical_retained", "approach_metadata_conflict",
      "historical_retained_124", "corrected_primary_90"), "reconciliation")
  .mv07e_require_columns(retained,
    c("orig.ident", "Approach.x", "Approach.y"), "retained metadata")
  mv07e_validate_accession_evidence_v1(evidence)
  retained_rows <- reconciliation[.mv07e_truth(
    reconciliation$historical_retained_124), , drop = FALSE]
  if (nrow(reconciliation) != 127L || nrow(retained_rows) != 124L ||
      nrow(retained) != 124L || anyDuplicated(retained$orig.ident)) {
    stop("MV7-E requires the fixed 127/124 reconciliation.", call. = FALSE)
  }
  index <- match(retained_rows$sample_id, retained$orig.ident)
  if (anyNA(index)) stop("Retained metadata join is incomplete.", call. = FALSE)
  x <- trimws(retained$Approach.x[index])
  y <- trimws(retained$Approach.y[index])
  public <- trimws(retained_rows$approach_public)
  conflict <- x != public
  if (sum(conflict) != 16L || any(y != public) ||
      !setequal(retained_rows$sample_id[conflict], evidence$sample_id)) {
    stop("Approach.x/Approach.y lineage does not match the accession audit.",
         call. = FALSE)
  }
  ev <- evidence[match(retained_rows$sample_id, evidence$sample_id), , drop = FALSE]
  resolved <- !is.na(ev$sample_id)
  if (any(conflict & (!resolved | ev$canonical_approach != public))) {
    stop("A disputed approach lacks matching official resolution.", call. = FALSE)
  }
  result <- data.frame(
    contract_id = "mv07e_canonical_approach_v1",
    sample_id = retained_rows$sample_id, study = retained_rows$study,
    tissue = retained_rows$tissue, corrected_primary_90 =
      .mv07e_truth(retained_rows$corrected_primary_90),
    historical_heuristic_approach = x, public_metadata_approach = y,
    canonical_approach = public, metadata_conflict = conflict,
    official_accession_resolution = conflict & resolved,
    canonical_source = ifelse(conflict,
      "public_metadata_confirmed_by_official_geo_methods", "public_metadata"),
    historical_heuristic_use = "prohibited_for_scientific_approach_labels",
    stringsAsFactors = FALSE)
  result[order(result$sample_id, method = "radix"), , drop = FALSE]
}

#' Summarize the correction's impact on approach-only diagnostics.
mv07e_approach_correction_v1 <- function(metadata) {
  .mv07e_require_columns(metadata,
    c("sample_id", "study", "canonical_approach", "corrected_primary_90"),
    "canonical approach metadata")
  primary <- metadata[.mv07e_truth(metadata$corrected_primary_90), , drop = FALSE]
  mixed <- names(which(vapply(split(primary$canonical_approach, primary$study),
    function(value) length(unique(value)) > 1L, logical(1L))))
  data.frame(
    contract_id = "mv07e_mv07b_approach_correction_v1",
    affected_result = "MV7-B mixed-study approach association only",
    primary_samples = nrow(primary),
    primary_scrna = sum(primary$canonical_approach == "scRNA-seq"),
    primary_snrna = sum(primary$canonical_approach == "snRNA-seq"),
    primary_mixed_approach_studies = length(mixed),
    corrected_disposition = if (length(mixed)) "estimable_descriptive_only" else
      "not_estimable_zero_mixed_approach_studies",
    original_approach_result_status = "superseded_metadata_field_error",
    topology_or_landscape_recalculation_required = FALSE,
    study_tissue_cell_count_results_affected = FALSE,
    causal_technology_claim_authorized = FALSE,
    stringsAsFactors = FALSE)
}

#' Audit one added sample against the accepted 500-feature panel.
mv07e_panel_availability_row_v1 <- function(sample_id, source_feature_ids,
                                             panel) {
  .mv07e_require_columns(panel, c("feature_id", "panel_order"), "panel")
  if (nrow(panel) != 500L || anyDuplicated(panel$feature_id) ||
      anyDuplicated(panel$panel_order)) stop("Frozen panel is invalid.", call. = FALSE)
  normalized <- mv07e_seurat_feature_ids_v1(source_feature_ids)
  missing <- panel$feature_id[!panel$feature_id %in% normalized]
  data.frame(
    contract_id = "mv07e_added_sample_panel_availability_v1",
    sample_id = sample_id, source_features = length(source_feature_ids),
    normalized_unique_features = length(unique(normalized)),
    panel_features = nrow(panel), panel_features_present = nrow(panel) - length(missing),
    missing_features = length(missing),
    missing_feature_ids = paste(missing, collapse = ";"),
    expression_values_used = FALSE, expression_values_exported = FALSE,
    topology_read = FALSE,
    biological_labels_read = FALSE, stringsAsFactors = FALSE)
}

#' Select the availability-only panel branch.
mv07e_panel_decision_v1 <- function(availability, panel_sha256) {
  .mv07e_require_columns(availability,
    c("sample_id", "panel_features", "missing_features"), "availability")
  if (nrow(availability) != 34L || anyDuplicated(availability$sample_id) ||
      any(availability$panel_features != 500L) ||
      !grepl("^[0-9a-f]{64}$", panel_sha256)) {
    stop("Panel availability decision inputs are invalid.", call. = FALSE)
  }
  complete <- all(availability$missing_features == 0L)
  data.frame(
    contract_id = "mv07e_panel_branch_decision_v1",
    accepted_90_panel_sha256 = panel_sha256,
    added_samples = 34L, complete_samples = sum(availability$missing_features == 0L),
    incomplete_samples = sum(availability$missing_features > 0L),
    missing_feature_occurrences = sum(availability$missing_features),
    branch = if (complete) "retain_accepted_90_derived_panel" else
      "fit_deterministic_global_core_over_124",
    selection_basis = "feature_availability_only_no_expression_or_outcomes",
    exact_final_panel_status = if (complete) "locked_existing_panel" else
      "pending_620_sct_cache_inventory_after_mv07f",
    fallback_algorithm = "mv06c_global_core_panel_v1_over_124_samples_five_seeds",
    ph_authorized = FALSE, stringsAsFactors = FALSE)
}

#' Return the full sample/seed axis without biological labels.
mv07e_sample_seed_axis_v1 <- function(sample_ids,
                                      seeds = 20260805:20260809) {
  sample_ids <- sort(unique(as.character(sample_ids)), method = "radix")
  seeds <- as.integer(seeds)
  if (length(sample_ids) != 124L || length(seeds) != 5L ||
      !identical(seeds, 20260805:20260809)) {
    stop("MV7-E requires 124 samples and the five accepted seeds.", call. = FALSE)
  }
  grid <- expand.grid(sample_id = sample_ids, seed = seeds,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- grid[order(grid$seed, grid$sample_id, method = "radix"), , drop = FALSE]
  grid$contract_id <- "mv07e_sample_seed_axis_v1"
  grid$axis_order <- seq_len(nrow(grid)); grid$selected_cells <- 384L
  grid$outcome_label_state <- "closed"
  grid$biological_outcomes_computed <- FALSE
  grid[c("contract_id", "axis_order", "sample_id", "seed", "selected_cells",
         "outcome_label_state", "biological_outcomes_computed")]
}

#' Freeze one transductive descriptive transform per seed.
mv07e_fit_scopes_v1 <- function(seeds = 20260805:20260809) {
  seeds <- as.integer(seeds)
  if (!identical(seeds, 20260805:20260809)) stop("Seed axis changed.", call. = FALSE)
  data.frame(
    contract_id = "mv07e_global_descriptive_fit_scope_v1",
    fit_scope_order = seq_along(seeds),
    fit_scope_id = paste0("mv07e_full124_seed_", seeds), seed = seeds,
    fit_samples = 124L, fit_cells = 124L * 384L, panel_genes = 500L,
    pca_components = 30L,
    standardization_scope = "all_124_equal_depth_cells_within_seed",
    pca_scope = "all_124_equal_depth_cells_within_seed",
    transductive = TRUE, comparison_to_loso_transform = "prohibited_as_same_estimand",
    labels_used_for_fit = FALSE, stringsAsFactors = FALSE)
}

#' Freeze typed views, PH, and pairwise landscape cardinalities.
mv07e_topology_contract_v1 <- function() {
  data.frame(
    contract_id = "mv07e_typed_topology_contract_v1",
    view = c("cell_topology_v1", "gene_topology_v1"),
    points = c(384L, 500L), ambient_dimensions = c(30L, 384L),
    point_identity = c("selected_cells", "frozen_panel_genes"),
    geometry = c("euclidean_in_global_descriptive_pcs",
                 "correlation_chord_across_matched_cells"),
    correlation_rule = c("not_applicable",
      "sqrt_2_times_1_minus_pearson_r_clamped_0_2"),
    filtration = "complete_vietoris_rips", homology_dimensions = "H0;H1",
    essential_h0 = "exclude_infinite_interval",
    zero_variance_policy = "stop_no_imputation", labels_used = FALSE,
    stringsAsFactors = FALSE)
}

mv07e_pair_scope_v1 <- function() {
  data.frame(
    contract_id = "mv07e_pair_scope_v1", samples = 124L, seeds = 5L,
    unordered_pairs_per_seed = choose(124L, 2L),
    biological_pairs_per_view = choose(124L, 2L) * 5L,
    views = 2L, homology_dimensions = 2L,
    component_distance_rows = choose(124L, 2L) * 5L * 2L * 2L,
    cross_seed_pairs = FALSE,
    pair_identity = "lexicographic_sample_a_less_than_sample_b_within_seed",
    stringsAsFactors = FALSE)
}

mv07e_landscape_contract_v1 <- function() {
  data.frame(
    contract_id = "mv07e_landscape_contract_v1",
    item = c("finite_intervals", "essential_h0", "level_policy", "integration",
      "dimension_policy", "grid_policy", "level_cap_policy", "streaming"),
    required_state = c("all_finite_positive_persistence_intervals",
      "exclude_infinite_interval", "all_consecutive_active_levels",
      "exact_or_error_controlled_squared_l2_on_dimension_support",
      "h0_h1_separate", "no_universal_fixed_grid", "no_universal_level_cap",
      "stream_or_chunk_without_dense_landscape_materialization"),
    changed_by_mv07e = FALSE, stringsAsFactors = FALSE)
}

mv07e_resource_contract_v1 <- function() {
  data.frame(
    contract_id = "mv07e_resource_and_resume_contract_v1",
    stage = c("mv07f_raw_child", "mv07f_sct_child", "mv07f_aggregate",
      "mv07g_sentinel", "mv07h_full_ph_landscape"),
    elapsed_cap_seconds = c(1800, 1800, 14400, NA, NA),
    rss_cap_bytes = c(8, 8, 8, 8, NA) * 1024^3,
    storage_cap_bytes = c(NA, NA, 4, NA, NA) * 1024^3,
    cap_policy = c("stop_child", "stop_child", "stop_before_publish",
      "set_from_measured_mv07f_then_stop_or_authorize",
      "set_only_after_mv07g_measured_projection"),
    atomic_write = TRUE, immutable_resume = TRUE,
    partial_state_publishable = FALSE, stringsAsFactors = FALSE)
}

mv07e_decision_v1 <- function(panel_decision, correction) {
  if (nrow(panel_decision) != 1L || nrow(correction) != 1L ||
      panel_decision$branch != "fit_deterministic_global_core_over_124" ||
      correction$corrected_disposition !=
        "not_estimable_zero_mixed_approach_studies") {
    stop("MV7-E continuation decision inputs changed.", call. = FALSE)
  }
  data.frame(
    contract_id = "mv07e_prefreeze_decision_v1",
    decision = "authorize_mv07f_upstream_only_then_finalize_124_panel",
    next_stage = "MV7-F",
    raw_shards_authorized = 34L, sct_caches_authorized = 170L,
    final_panel_selection_authorized_after_sct = TRUE,
    pca_authorized = FALSE, ph_authorized = FALSE,
    landscape_authorized = FALSE, clustering_authorized = FALSE,
    outcomes_authorized = FALSE, label_access_authorized = FALSE,
    primary_90_recalculation_authorized = FALSE,
    external_data_authorized = FALSE, owner_input_required = FALSE,
    stringsAsFactors = FALSE)
}
