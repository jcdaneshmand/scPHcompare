# MV7-D full-corpus reconciliation helpers.

.mv07d_require_columns <- function(x, required, label) {
  if (!is.data.frame(x) || !all(required %in% names(x))) {
    stop(label, " is missing required columns: ",
         paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
  }
  invisible(x)
}

.mv07d_truth <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) == "true"
}

#' Reconcile the 127 candidates, 124 retained samples, and primary 90.
mv07d_reconcile_samples_v1 <- function(candidates, retained, accepted_ids) {
  .mv07d_require_columns(candidates,
    c("Sample Name", "SRA", "Tissue", "Approach"), "candidate metadata")
  .mv07d_require_columns(retained,
    c("orig.ident", "SRA", "Tissue.x", "Approach.x",
      "Number_of_Cells_After_Filtering"), "retained metadata")
  if (nrow(candidates) != 127L || nrow(retained) != 124L) {
    stop("MV7-D requires the historical 127-candidate and 124-retained axes.",
         call. = FALSE)
  }
  candidate <- data.frame(
    candidate_order = seq_len(nrow(candidates)),
    sample_id = paste(candidates$SRA, candidates[["Sample Name"]], sep = "_"),
    study = trimws(candidates$SRA), tissue = trimws(candidates$Tissue),
    approach = trimws(candidates$Approach), stringsAsFactors = FALSE)
  retained_clean <- data.frame(
    sample_id = trimws(retained[["orig.ident"]]),
    study_retained = trimws(retained$SRA),
    tissue_retained = trimws(retained[["Tissue.x"]]),
    approach_retained = trimws(retained[["Approach.x"]]),
    post_qc_cells = as.integer(retained[["Number_of_Cells_After_Filtering"]]),
    stringsAsFactors = FALSE)
  if (anyDuplicated(candidate$sample_id) || anyDuplicated(retained_clean$sample_id) ||
      anyNA(retained_clean$post_qc_cells) ||
      !all(retained_clean$sample_id %in% candidate$sample_id)) {
    stop("Candidate or retained sample identity is invalid.", call. = FALSE)
  }
  result <- merge(candidate, retained_clean, by = "sample_id", all.x = TRUE,
                  sort = FALSE)
  result <- result[match(candidate$sample_id, result$sample_id), , drop = FALSE]
  retained_flag <- !is.na(result$post_qc_cells)
  if (any(result$study[retained_flag] != result$study_retained[retained_flag]) ||
      any(result$tissue[retained_flag] != result$tissue_retained[retained_flag]) ||
      any(result$approach[retained_flag] != result$approach_retained[retained_flag])) {
    stop("Retained metadata conflicts with the public candidate identity.",
         call. = FALSE)
  }
  excluded <- data.frame(
    sample_id = c("SRA850958_SRS4386092", "SRA850958_SRS4386091",
                  "SRA850958_SRS4386107"),
    post_qc_cells = c(169L, 197L, 166L), stringsAsFactors = FALSE)
  idx <- match(excluded$sample_id, result$sample_id)
  if (anyNA(idx) || any(retained_flag[idx])) {
    stop("Historical below-threshold sample identity changed.", call. = FALSE)
  }
  result$post_qc_cells[idx] <- excluded$post_qc_cells
  accepted_ids <- sort(unique(trimws(as.character(accepted_ids))), method = "radix")
  if (length(accepted_ids) != 90L || !all(accepted_ids %in% result$sample_id)) {
    stop("Accepted primary sample axis must contain exactly 90 candidates.",
         call. = FALSE)
  }
  result$historical_min_cells <- 250L
  result$historical_retained_124 <- retained_flag
  result$historical_ph_output_expected <- retained_flag
  result$corrected_primary_90 <- result$sample_id %in% accepted_ids
  result$corrected_descriptive_124 <- retained_flag
  result$threshold_sensitivity_only <- !retained_flag
  result$corpus_class <- ifelse(
    !retained_flag, "below_250_pre_ph_exclusion",
    ifelse(result$corrected_primary_90, "primary_cross_study_eligible",
           "single_study_tissue_descriptive_only"))
  result$reason_code <- ifelse(
    !retained_flag, "excluded_post_qc_min_cells",
    ifelse(result$corrected_primary_90, "included_primary_cross_study",
           "single_study_tissue_not_estimable"))
  result$outcome_use <- ifelse(
    !retained_flag, "separate_threshold_sensitivity_not_authorized",
    ifelse(result$corrected_primary_90, "primary_and_descriptive",
           "descriptive_only"))
  keep <- c("candidate_order", "sample_id", "study", "tissue", "approach",
            "post_qc_cells", "historical_min_cells", "historical_retained_124",
            "historical_ph_output_expected", "corrected_primary_90",
            "corrected_descriptive_124", "threshold_sensitivity_only",
            "corpus_class", "reason_code", "outcome_use")
  result <- result[order(result$candidate_order), keep, drop = FALSE]
  rownames(result) <- NULL
  if (sum(result$historical_retained_124) != 124L ||
      sum(result$corrected_primary_90) != 90L ||
      sum(result$corrected_descriptive_124) != 124L ||
      sum(result$threshold_sensitivity_only) != 3L) {
    stop("Reconciled corpus cardinalities are not 127 -> 124 -> 90.",
         call. = FALSE)
  }
  result
}

#' Summarize tissue/study coverage and endpoint eligibility.
mv07d_tissue_study_summary_v1 <- function(samples) {
  .mv07d_require_columns(samples,
    c("sample_id", "study", "tissue", "historical_retained_124",
      "corrected_primary_90"), "reconciled samples")
  rows <- lapply(sort(unique(samples$tissue), method = "radix"), function(tissue) {
    x <- samples[samples$tissue == tissue, , drop = FALSE]
    retained <- x[.mv07d_truth(x$historical_retained_124), , drop = FALSE]
    data.frame(
      contract_id = "mv07d_tissue_study_summary_v1", tissue = tissue,
      candidate_samples = nrow(x), retained_samples = nrow(retained),
      retained_studies = length(unique(retained$study)),
      primary_samples = sum(.mv07d_truth(x$corrected_primary_90)),
      primary_cross_study_eligible = length(unique(retained$study)) >= 2L,
      descriptive_eligible = nrow(retained) > 0L, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Choose fixed source/SCT sentinels without topology or outcome values.
mv07d_select_omitted_sentinels_v1 <- function(samples, seed = 20260805L) {
  omitted <- samples[samples$corpus_class ==
                       "single_study_tissue_descriptive_only", , drop = FALSE]
  if (nrow(omitted) != 34L || length(unique(omitted$tissue)) != 3L ||
      any(omitted$post_qc_cells < 384L)) {
    stop("The 34-sample omitted descriptive stratum is not feasible at 384 cells.",
         call. = FALSE)
  }
  rows <- lapply(sort(unique(omitted$tissue), method = "radix"), function(tissue) {
    x <- omitted[omitted$tissue == tissue, , drop = FALSE]
    x <- x[order(x$post_qc_cells, x$sample_id, method = "radix"), , drop = FALSE]
    selected <- x[c(1L, nrow(x)), , drop = FALSE]
    selected$selection_boundary <- c("minimum_post_qc_cells",
                                     "maximum_post_qc_cells")
    selected
  })
  result <- do.call(rbind, rows)
  result <- result[c("sample_id", "study", "tissue", "approach",
                     "post_qc_cells", "selection_boundary")]
  result$seed <- as.integer(seed)
  result$selected_cells <- 384L
  result$selection_basis <- "tissue_stratified_depth_extremes_no_topology_results"
  result$outcome_label_use <- "design_stratification_only"
  result$ph_authorized <- FALSE
  rownames(result) <- NULL
  result
}

#' Return the estimand-specific population registry.
mv07d_estimand_populations_v1 <- function() {
  data.frame(
    contract_id = "mv07d_estimand_populations_v1", population_order = 1:5,
    population_id = c("historical_candidate", "historical_reproduction",
      "corrected_primary_cross_study", "corrected_full_corpus_descriptive",
      "below_threshold_sensitivity"),
    samples = c(127L, 124L, 90L, 124L, 3L),
    scientific_role = c("sample_flow_denominator", "legacy_reproduction_only",
      "primary_blocked_tissue_retrieval", "descriptive_topology_and_clustering",
      "threshold_sensitivity_only"),
    current_computation = c("metadata_only", "legacy_outputs_orientation_ineligible",
      "corrected_complete", "corrected_90_complete_34_missing",
      "not_authorized"),
    primary_claim_eligible = c(FALSE, FALSE, TRUE, FALSE, FALSE),
    note = c("includes three samples removed before PH",
      "legacy feature-as-point outputs cannot substitute for corrected results",
      "five tissues represented in at least two studies",
      "adds three single-study tissues without cross-study inference",
      "cannot use the fixed 384-cell depth and cannot affect the primary endpoint"),
    stringsAsFactors = FALSE)
}

#' Preserve the accepted landscape definition for any future expansion.
mv07d_landscape_contract_v1 <- function() {
  data.frame(
    contract_id = "mv07d_landscape_contract_v1",
    item = c("finite_intervals", "essential_h0", "level_policy", "integration",
      "dimension_policy", "grid_policy", "level_cap_policy", "streaming"),
    required_state = c("all_finite_positive_persistence_intervals",
      "exclude_infinite_interval", "all_consecutive_active_levels",
      "exact_or_error_controlled_squared_l2_on_dimension_support",
      "h0_h1_separate", "no_universal_fixed_grid", "no_universal_level_cap",
      "stream_or_chunk_without_dense_landscape_materialization"),
    applies_to_full_corpus_expansion = TRUE, changed_by_mv07d = FALSE,
    stringsAsFactors = FALSE)
}

#' Return the prospective full-corpus gate disposition.
mv07d_expansion_gate_v1 <- function(samples, sentinels, source_coverage) {
  if (nrow(samples) != 127L || nrow(sentinels) != 6L ||
      nrow(source_coverage) != 3L ||
      !all(source_coverage$source_kind %in%
             c("candidate_sparse_rdata", "individual_seurat_rds",
               "accepted_corrected_artifacts")) ||
      any(source_coverage$observed_samples != source_coverage$expected_samples)) {
    stop("MV7-D structural expansion gate failed.", call. = FALSE)
  }
  data.frame(
    contract_id = "mv07d_expansion_gate_v1",
    decision = "authorize_six_sample_source_sct_feasibility_only",
    next_stage = "MV7-D1",
    primary_90_recalculation_authorized = FALSE,
    omitted_34_ph_authorized = FALSE,
    omitted_34_outcome_authorized = FALSE,
    below_250_sensitivity_authorized = FALSE,
    new_data_authorized = FALSE,
    manuscript_claim_promotion_authorized = FALSE,
    owner_input_required = FALSE,
    stringsAsFactors = FALSE)
}
