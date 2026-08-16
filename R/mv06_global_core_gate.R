.mv06c_rank_variances <- function(variances) {
  apply(variances, 2L, function(value) {
    rank(-value, ties.method = "min", na.last = "keep")
  })
}

.mv06c_panel_sha256 <- function(panel) {
  payload <- data.frame(
    panel_order = as.integer(panel$panel_order),
    feature_id = as.character(panel$feature_id),
    gene = as.character(panel$gene),
    stringsAsFactors = FALSE
  )
  rownames(payload) <- NULL
  digest::digest(payload, algo = "sha256", serialize = TRUE)
}

#' Build the frozen MV6-C global-core panel from aggregate variances
#'
#' @keywords internal
mv06c_build_global_core_panel_v1 <- function(feature_ids, variances, seeds,
                                              panel_size = 500L,
                                              samples_per_seed = 90L) {
  feature_ids <- as.character(feature_ids)
  seeds <- as.integer(seeds)
  panel_size <- as.integer(panel_size)
  samples_per_seed <- as.integer(samples_per_seed)
  if (!is.matrix(variances) || !is.numeric(variances) ||
      nrow(variances) != length(feature_ids) ||
      ncol(variances) != length(seeds) ||
      !length(feature_ids) || anyNA(feature_ids) || any(!nzchar(feature_ids)) ||
      anyDuplicated(feature_ids) || anyNA(seeds) ||
      !identical(sort(unique(seeds)), 20260805:20260809) ||
      length(samples_per_seed) != 1L || is.na(samples_per_seed) ||
      samples_per_seed < 1L ||
      any(tabulate(match(seeds, 20260805:20260809), nbins = 5L) !=
          samples_per_seed) ||
      length(panel_size) != 1L || is.na(panel_size) || panel_size < 2L) {
    stop("MV6-C feature, variance, seed, or panel axes are invalid.",
         call. = FALSE)
  }
  canonical <- canonical_mv03_gene_ids(feature_ids)
  category <- mv03_feature_category(feature_ids)
  finite <- is.finite(variances)
  positive <- finite & variances > .Machine$double.eps
  finite_all <- rowSums(finite) == ncol(variances)
  positive_all <- rowSums(positive) == ncol(variances)
  ranks <- .mv06c_rank_variances(variances)
  median_rank <- apply(ranks, 1L, stats::median, na.rm = TRUE)
  minimum_variance <- apply(variances, 1L, function(value) {
    if (all(is.finite(value))) min(value) else NA_real_
  })
  candidate <- data.frame(
    feature_id = feature_ids, gene = canonical, category = category,
    finite_cache_count = rowSums(finite),
    positive_cache_count = rowSums(positive),
    finite_all_caches = finite_all,
    positive_all_caches = positive_all,
    median_variance_rank = median_rank,
    minimum_variance = minimum_variance,
    stringsAsFactors = FALSE
  )
  retained <- candidate[
    candidate$category == "retained_candidate" &
      candidate$finite_all_caches & candidate$positive_all_caches,
    , drop = FALSE
  ]
  retained <- retained[order(
    retained$median_variance_rank, retained$gene, retained$feature_id,
    method = "radix"
  ), , drop = FALSE]
  before_deduplication <- nrow(retained)
  retained <- retained[!duplicated(retained$gene), , drop = FALSE]
  retained$global_candidate_order <- seq_len(nrow(retained))
  selected_n <- min(panel_size, nrow(retained))
  panel <- retained[seq_len(selected_n), , drop = FALSE]
  panel$panel_order <- seq_len(nrow(panel))
  panel$selected_all_cache_core <- TRUE

  seed_rows <- lapply(20260805:20260809, function(seed) {
    columns <- seeds == seed
    seed_positive <- rowSums(positive[, columns, drop = FALSE]) == sum(columns)
    seed_rank <- apply(ranks[, columns, drop = FALSE], 1L,
                       stats::median, na.rm = TRUE)
    seed_candidate <- data.frame(
      feature_id = feature_ids, gene = canonical, category = category,
      positive = seed_positive, median_rank = seed_rank,
      stringsAsFactors = FALSE
    )
    seed_candidate <- seed_candidate[
      seed_candidate$category == "retained_candidate" &
        seed_candidate$positive,
      , drop = FALSE
    ]
    seed_candidate <- seed_candidate[order(
      seed_candidate$median_rank, seed_candidate$gene,
      seed_candidate$feature_id, method = "radix"
    ), , drop = FALSE]
    seed_candidate <- seed_candidate[!duplicated(seed_candidate$gene),,
                                     drop = FALSE]
    seed_top <- head(seed_candidate$feature_id, panel_size)
    overlap <- length(intersect(panel$feature_id, seed_top))
    union <- length(union(panel$feature_id, seed_top))
    shared_rank_axis <- match(retained$feature_id, seed_candidate$feature_id)
    correlation <- suppressWarnings(stats::cor(
      retained$median_variance_rank,
      seed_candidate$median_rank[shared_rank_axis],
      method = "spearman", use = "complete.obs"
    ))
    data.frame(
      contract_id = "mv06c_seed_panel_stability_v1",
      seed = seed, cache_count = sum(columns),
      eligible_unique_genes = nrow(seed_candidate),
      seed_top_panel_size = length(seed_top),
      overlap_with_global_panel = overlap,
      jaccard_with_global_panel = if (union) overlap / union else NA_real_,
      spearman_global_vs_seed_candidate_ranks = correlation,
      stringsAsFactors = FALSE
    )
  })
  seed_stability <- do.call(rbind, seed_rows)

  eligibility <- data.frame(
    contract_id = "mv06c_global_core_eligibility_summary_v1",
    common_features = length(feature_ids),
    retained_category_features = sum(category == "retained_candidate"),
    technical_category_features = sum(category != "retained_candidate"),
    nonfinite_any_cache_features = sum(!finite_all),
    nonpositive_any_cache_features = sum(finite_all & !positive_all),
    eligible_features_before_canonical_deduplication = before_deduplication,
    duplicate_canonical_features_removed = before_deduplication - nrow(retained),
    eligible_unique_canonical_genes = nrow(retained),
    requested_panel_size = panel_size,
    selected_panel_size = nrow(panel),
    eligibility_margin = nrow(retained) - panel_size,
    stringsAsFactors = FALSE
  )
  decision <- if (nrow(retained) >= panel_size) {
    "go_bounded_matched_sct_profile"
  } else {
    "stop_global_core_insufficient"
  }
  list(
    panel = panel, candidates = retained, eligibility = eligibility,
    seed_stability = seed_stability, decision = decision,
    panel_sha256 = .mv06c_panel_sha256(panel)
  )
}

#' Build the frozen MV6-C future-workload inventory
#'
#' @keywords internal
mv06c_future_workload_v1 <- function() {
  data.frame(
    contract_id = "mv06c_future_matched_sct_workload_v1",
    matched_cell_views = 6750L,
    matched_gene_views = 6750L,
    h0_h1_diagram_components = 27000L,
    directed_query_training_pairs = 35350L,
    four_component_landscape_distances = 141400L,
    training_fitted_component_scales = 300L,
    five_weight_fusion_pair_rows = 176750L,
    execution_authorized = FALSE,
    next_allowed_scope = "bounded_prespecified_fold_seed_profile_only",
    stringsAsFactors = FALSE
  )
}
