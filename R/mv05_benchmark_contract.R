# Internal MV-05 statistical-plan and benchmark-contract helpers.

.mv05_required_metadata <- c(
  "cohort", "sample_id", "study", "tissue", "approach", "filtered_cells"
)

mv05_endpoint_registry_v1 <- function() {
  data.frame(
    endpoint_id = c(
      "cross_study_tissue_mrr_v1",
      "cross_study_tissue_1nn_balanced_accuracy_v1",
      "cross_study_tissue_distance_contrast_v1",
      "within_tissue_study_distance_contrast_v1",
      "within_study_approach_distance_contrast_v1",
      "integration_biological_change_v1",
      "integration_technical_change_v1",
      "label_free_cluster_stability_v1",
      "oracle_tissue_ari_v0",
      "bone_approach_retrieval_v1"
    ),
    family = c(
      "biological_conservation", "biological_conservation",
      "biological_conservation", "technical_signal", "technical_signal",
      "integration_change", "integration_change", "clustering",
      "historical_sensitivity", "technical_signal"
    ),
    role = c(
      "primary", "supportive", "supportive", "co_primary_diagnostic",
      "secondary", "primary", "co_primary_diagnostic", "secondary",
      "historical_sensitivity", "bone_technical_only"
    ),
    analysis_unit = c(
      "sample_macro_averaged_within_tissue", "sample_macro_averaged_within_tissue",
      "study_pair_summary", "study_pair_summary", "sample_pair_summary",
      "paired_outer_fold_summary", "paired_outer_fold_summary",
      "sample_partition_across_subsample_seeds", "sample", "sample"
    ),
    split_contract = c(
      rep("large_leave_one_study_out_v1", 7L),
      "all_samples_label_free_repeated_subsamples_v1",
      "oracle_k_historical_only_v0", "bone_within_study_technical_v1"
    ),
    direction = c(
      "higher_is_better", "higher_is_better", "higher_is_better",
      "lower_is_better", "lower_is_better", "higher_is_better",
      "lower_is_better", "higher_is_better", "higher_is_better",
      "report_without_biological_direction"
    ),
    primary_summary = c(
      "macro mean reciprocal rank", "macro balanced accuracy",
      "between-tissue minus within-tissue cross-study distance",
      "cross-study minus same-study distance within tissue",
      "different-approach minus same-approach distance within study and tissue",
      "integrated minus SCT change in biological endpoint",
      "integrated minus SCT change in residual technical endpoint",
      "mean pairwise adjusted Rand stability with jackknife Monte Carlo SE",
      "sample-level adjusted Rand index", "approach-stratified retrieval summary"
    ),
    multiplicity_family = c(
      "F1_primary_retrieval", "none_supportive", "none_supportive",
      "F2_primary_technical", "none_secondary", "F3_integration_change",
      "F3_integration_change", "none_label_free", "none_historical",
      "none_bone_technical"
    ),
    stringsAsFactors = FALSE
  )
}

mv05_method_registry_v1 <- function() {
  data.frame(
    method_id = c(
      "cell_landscape_h0_v1", "cell_landscape_h1_v1",
      "gene_landscape_h0_v1", "gene_landscape_h1_v1",
      "cell_distribution_energy_shared_pca_v1",
      "pseudobulk_shared_panel_euclidean_v1",
      "gene_correlation_frobenius_v1",
      "composition_hellinger_conditional_v1",
      "pam_stability_k_v1", "hclust_average_v1",
      "spectral_gaussian_laplacian_v1", "oracle_k_v0",
      "cell_level_seurat_louvain_direct_comparison_v0",
      "kmeans_on_distance_matrix_v0"
    ),
    layer = c(
      rep("sample_distance", 8L), rep("sample_clustering", 6L)
    ),
    view = c(
      "cell", "cell", "gene", "gene", "cell", "expression", "gene",
      "composition", "both", "both", "both", "both", "cell", "both"
    ),
    role = c(
      "confirmatory", "confirmatory", "secondary", "secondary",
      "matched_cell_baseline", "context_baseline", "matched_gene_baseline",
      "conditional_sensitivity", "primary", "sensitivity", "sensitivity",
      "historical_sensitivity", "excluded", "excluded"
    ),
    frozen_parameters = c(
      "all levels; exact/error-controlled; H0 separate",
      "all levels; exact/error-controlled; H1 separate",
      "all levels; exact/error-controlled; H0 separate",
      "all levels; exact/error-controlled; H1 separate",
      "energy distance on the same 384 cells in training-fit shared 30-PC space",
      "Euclidean distance between same-panel sample means after training-fit scaling",
      "RMS Frobenius distance between same-gene correlation matrices",
      "Hellinger distance; execute only with harmonized externally frozen cell labels",
      "k=2:min(10,n-1); five seeds; one-SE stability rule",
      "average linkage; same selected k; general dissimilarity",
      "median-bandwidth Gaussian graph; normalized Laplacian; deterministic seed",
      "known class count; never selects method or primary k",
      "cell and sample units differ", "distance matrices are not feature tables"
    ),
    execution_status = c(
      rep("requires_fold_fitted_full_existing_data", 8L),
      "requires_complete_five_seed_matrices", "planned_sensitivity",
      "requires_new_tested_implementation", "historical_only", "prohibited",
      "prohibited"
    ),
    stringsAsFactors = FALSE
  )
}

mv05_label_use_registry_v1 <- function() {
  data.frame(
    stage = c(
      "sample_eligibility_audit", "outer_split_assignment", "fit_transform",
      "distance_normalization", "select_method", "select_k",
      "write_predictions", "evaluate_biological", "evaluate_technical",
      "oracle_k_sensitivity"
    ),
    allowed_labels = c(
      "study;tissue;approach", "study", "", "", "", "",
      "", "tissue", "study;approach", "tissue;study;approach"
    ),
    required_boundary = c(
      "design_only_no_outcomes", "study_blocks_only", "training_samples_only",
      "training_distances_only", "prespecified_registry_only",
      "repeated_distance_stability_only", "predictions_hashed_before_labels_open",
      "predictions_and_matrices_immutable", "predictions_and_matrices_immutable",
      "historical_sensitivity_never_primary"
    ),
    stringsAsFactors = FALSE
  )
}

mv05_multiplicity_registry_v1 <- function() {
  data.frame(
    family_id = c(
      "F1_primary_retrieval", "F2_primary_technical", "F3_integration_change"
    ),
    included_contrasts = c(
      "cell H0/H1 versus matched cell-energy baseline within SCT and integrated representations",
      "cell H0/H1 residual study contrast within SCT and integrated representations",
      "cell H0/H1 integrated-minus-SCT biological and technical changes"
    ),
    adjustment = "Holm",
    interval_policy = "paired study-block bootstrap; 2000 replicates; 95% percentile interval",
    monte_carlo_policy = "optional paired study-level sign flip; B=9999; p=(b+1)/(B+1)",
    stringsAsFactors = FALSE
  )
}

new_mv05_benchmark_contract_v1 <- function() {
  identity <- list(
    contract_id = "mv05_statistical_benchmark_contract_v1",
    endpoint_registry = mv05_endpoint_registry_v1(),
    method_registry = mv05_method_registry_v1(),
    label_use_registry = mv05_label_use_registry_v1(),
    multiplicity_registry = mv05_multiplicity_registry_v1(),
    outer_split = "leave one study out; large cohort only",
    confirmatory_tissue_rule = "at least two studies per tissue",
    subsample_seeds = 20260805:20260809,
    bootstrap_replicates = 2000L,
    monte_carlo_replicates = 9999L,
    integration_policy = paste(
      "inductive held-out mapping required for confirmatory integration claims;",
      "current all-sample integrated artifacts are descriptive only"
    ),
    pilot_policy = paste(
      "MV-04 descriptive pilot distances may be used for technical smoke tests;",
      "they are prohibited from confirmatory endpoint estimation"
    )
  )
  structure(
    c(identity, list(cache_key = paste0(
      "mv05_contract_v1:",
      digest::digest(identity, algo = "sha256", serialize = TRUE)
    ))),
    class = "scph_mv05_benchmark_contract_v1"
  )
}

.mv05_validate_metadata <- function(metadata) {
  if (!is.data.frame(metadata) ||
      !all(.mv05_required_metadata %in% names(metadata))) {
    stop("metadata lacks required MV-05 sample fields.", call. = FALSE)
  }
  result <- metadata[, .mv05_required_metadata, drop = FALSE]
  for (field in setdiff(.mv05_required_metadata, "filtered_cells")) {
    result[[field]] <- trimws(as.character(result[[field]]))
    if (anyNA(result[[field]]) || any(!nzchar(result[[field]]))) {
      stop("metadata contains missing or empty ", field, ".", call. = FALSE)
    }
  }
  result$filtered_cells <- as.integer(result$filtered_cells)
  if (anyNA(result$filtered_cells) || any(result$filtered_cells <= 0L) ||
      anyDuplicated(result$sample_id)) {
    stop("metadata has invalid cell counts or duplicated sample IDs.",
         call. = FALSE)
  }
  result
}

mv05_design_feasibility_v1 <- function(metadata) {
  metadata <- .mv05_validate_metadata(metadata)
  groups <- split(metadata, interaction(metadata$cohort, metadata$tissue,
                                        drop = TRUE, lex.order = TRUE))
  result <- do.call(rbind, lapply(groups, function(values) {
    studies <- sort(unique(values$study), method = "radix")
    approaches <- sort(unique(values$approach), method = "radix")
    eligible <- identical(values$cohort[[1L]], "large") &&
      length(studies) >= 2L
    data.frame(
      cohort = values$cohort[[1L]], tissue = values$tissue[[1L]],
      sample_count = nrow(values), study_count = length(studies),
      approach_count = length(approaches), study_ids = paste(studies, collapse = ";"),
      approaches = paste(approaches, collapse = ";"),
      cross_study_tissue_eligible = eligible,
      disposition = if (eligible) "confirmatory_existing_data_candidate" else if (
        identical(values$cohort[[1L]], "bone")
      ) "technical_approach_only" else "single_study_tissue_descriptive_only",
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result[order(result$cohort, result$tissue), , drop = FALSE]
}

mv05_loso_fold_summary_v1 <- function(metadata) {
  metadata <- .mv05_validate_metadata(metadata)
  metadata <- metadata[metadata$cohort == "large", , drop = FALSE]
  if (!nrow(metadata)) stop("No large-cohort samples are available.", call. = FALSE)
  feasibility <- mv05_design_feasibility_v1(metadata)
  eligible_tissues <- feasibility$tissue[feasibility$cross_study_tissue_eligible]
  studies <- sort(unique(metadata$study), method = "radix")
  result <- do.call(rbind, lapply(studies, function(study) {
    test <- metadata[metadata$study == study, , drop = FALSE]
    train <- metadata[metadata$study != study, , drop = FALSE]
    eligible <- test$tissue %in% eligible_tissues & test$tissue %in% train$tissue
    data.frame(
      split_id = paste0("large_loso_v1:", study), held_out_study = study,
      training_samples = nrow(train), test_samples = nrow(test),
      eligible_test_samples = sum(eligible),
      ineligible_test_samples = sum(!eligible),
      test_tissues = paste(sort(unique(test$tissue), method = "radix"),
                           collapse = ";"),
      confirmatory_tissues = paste(sort(unique(test$tissue[eligible]),
                                        method = "radix"), collapse = ";"),
      no_study_overlap = !any(train$study %in% test$study),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}

mv05_validate_label_use_v1 <- function(stage, labels = character()) {
  registry <- mv05_label_use_registry_v1()
  stage <- as.character(stage)
  if (length(stage) != 1L || !stage %in% registry$stage) {
    stop("Unknown MV-05 label-use stage.", call. = FALSE)
  }
  labels <- sort(unique(tolower(as.character(labels))), method = "radix")
  labels <- labels[nzchar(labels)]
  allowed <- registry$allowed_labels[registry$stage == stage]
  allowed <- if (!nzchar(allowed)) character() else
    strsplit(allowed, ";", fixed = TRUE)[[1L]]
  prohibited <- setdiff(labels, allowed)
  if (length(prohibited)) {
    stop("Prohibited label use at stage ", stage, ": ",
         paste(prohibited, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

mv05_monte_carlo_p_v1 <- function(exceedances, replicates) {
  exceedances <- as.numeric(exceedances)
  replicates <- as.numeric(replicates)
  if (length(exceedances) != 1L || length(replicates) != 1L ||
      !is.finite(exceedances) || !is.finite(replicates) ||
      exceedances < 0 || replicates < 1 || exceedances > replicates ||
      exceedances != floor(exceedances) || replicates != floor(replicates)) {
    stop("exceedances and replicates must be valid integer counts.",
         call. = FALSE)
  }
  (exceedances + 1) / (replicates + 1)
}

mv05_select_stable_k_v1 <- function(assignments) {
  required <- c("seed", "k", "sample_id", "cluster")
  if (!is.data.frame(assignments) || !all(required %in% names(assignments))) {
    stop("assignments lacks seed, k, sample_id, or cluster.", call. = FALSE)
  }
  assignments$seed <- as.character(assignments$seed)
  assignments$k <- as.integer(assignments$k)
  assignments$sample_id <- as.character(assignments$sample_id)
  assignments$cluster <- as.character(assignments$cluster)
  expected_seeds <- as.character(20260805:20260809)
  seeds <- sort(unique(assignments$seed), method = "radix")
  if (!identical(seeds, expected_seeds)) {
    return(list(status = "no_stable_k", selected_k = NA_integer_,
                summary = data.frame()))
  }
  rows <- lapply(sort(unique(assignments$k)), function(k) {
    part <- assignments[assignments$k == k, , drop = FALSE]
    by_seed <- split(part, part$seed)
    if (!identical(sort(names(by_seed), method = "radix"), seeds)) return(NULL)
    sample_axes <- lapply(by_seed, function(x) sort(x$sample_id, method = "radix"))
    if (!all(vapply(sample_axes, identical, logical(1L), sample_axes[[1L]])) ||
        any(vapply(by_seed, function(x) anyDuplicated(x$sample_id), integer(1L)))) {
      return(NULL)
    }
    clusterings <- lapply(by_seed, function(x) {
      x$cluster[match(sample_axes[[1L]], x$sample_id)]
    })
    if (any(vapply(clusterings, function(x) length(unique(x)) < 2L,
                   logical(1L)))) return(NULL)
    pair_values <- utils::combn(seq_along(clusterings), 2L, function(index) {
      mclust::adjustedRandIndex(clusterings[[index[[1L]]]],
                               clusterings[[index[[2L]]]])
    })
    mean_stability <- mean(pair_values)
    leave_one_out <- vapply(seq_along(clusterings), function(index) {
      kept <- clusterings[-index]
      mean(utils::combn(seq_along(kept), 2L, function(pair) {
        mclust::adjustedRandIndex(kept[[pair[[1L]]]], kept[[pair[[2L]]]])
      }))
    }, numeric(1L))
    jackknife_se <- sqrt((length(seeds) - 1) / length(seeds) *
                           sum((leave_one_out - mean(leave_one_out)) ^ 2))
    data.frame(k = k, mean_stability = mean_stability,
               monte_carlo_se = jackknife_se, pair_count = length(pair_values),
               stringsAsFactors = FALSE)
  })
  summary <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(summary) || !nrow(summary)) {
    return(list(status = "no_stable_k", selected_k = NA_integer_,
                summary = data.frame()))
  }
  best <- summary[which.max(summary$mean_stability), , drop = FALSE]
  threshold <- best$mean_stability - best$monte_carlo_se
  eligible <- summary$k[summary$mean_stability >= threshold]
  list(status = "selected", selected_k = min(eligible), summary = summary,
       threshold = threshold)
}

mv05_validate_matrix_record_v1 <- function(distance_matrix, analysis_role,
                                            fit_scope_id) {
  distance_matrix <- .validate_distance_matrix_v1(distance_matrix)
  analysis_role <- as.character(analysis_role)
  fit_scope_id <- as.character(fit_scope_id)
  if (length(analysis_role) != 1L || length(fit_scope_id) != 1L ||
      !nzchar(analysis_role) || !nzchar(fit_scope_id)) {
    stop("analysis_role and fit_scope_id must be non-empty strings.",
         call. = FALSE)
  }
  if (analysis_role == "confirmatory_blocked" &&
      grepl("descriptive_all|pilot", fit_scope_id, ignore.case = TRUE)) {
    stop("Descriptive/pilot fit scopes are prohibited for blocked confirmation.",
         call. = FALSE)
  }
  invisible(TRUE)
}
