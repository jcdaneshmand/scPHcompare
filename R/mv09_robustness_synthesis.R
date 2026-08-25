.mv09_metrics_v1 <- function() c(
  "pearson", "spearman", "nonnegative_scale", "relative_stress",
  "left_median_distance", "right_median_distance",
  "median_abs_median_scaled_change", "p95_abs_median_scaled_change",
  "mean_neighbor_jaccard", "median_neighbor_jaccard",
  "p10_neighbor_jaccard"
)

.mv09_quantile <- function(x, p) as.numeric(stats::quantile(
  x, p, names = FALSE, type = 7
))

mv09_build_robustness_synthesis_v1 <- function(value) {
  required <- c(
    "comparison_order", "dataset_scope", "contrast_id", "seed",
    "homology_dimension", "left_stack", "right_stack", "units",
    "unordered_pairs", "distance_transform", .mv09_metrics_v1(),
    "interpretation", "outcome_label_state", "biological_outcomes_computed",
    "clustering_jobs", "fusion_jobs", "label_jobs", "outcome_jobs"
  )
  if (!is.data.frame(value) || nrow(value) != 40L ||
      !all(required %in% names(value)) ||
      !identical(as.integer(value$comparison_order), 1:40) ||
      anyDuplicated(value$comparison_order) ||
      sum(value$dataset_scope == "internal124") != 30L ||
      sum(value$dataset_scope == "external8") != 10L ||
      sum(value$homology_dimension == "H0") != 20L ||
      sum(value$homology_dimension == "H1") != 20L ||
      any(value$outcome_label_state != "closed") ||
      any(as.logical(value$biological_outcomes_computed)) ||
      any(value$clustering_jobs != 0L | value$fusion_jobs != 0L |
            value$label_jobs != 0L | value$outcome_jobs != 0L)) {
    stop("MV9 aggregate comparison contract is invalid", call. = FALSE)
  }
  metrics <- .mv09_metrics_v1()
  for (metric in metrics) {
    x <- as.numeric(value[[metric]])
    if (anyNA(x) || any(!is.finite(x))) stop("MV9 metric is nonfinite")
  }
  long <- do.call(rbind, lapply(metrics, function(metric) data.frame(
    comparison_order = value$comparison_order,
    dataset_scope = value$dataset_scope, contrast_id = value$contrast_id,
    seed = as.integer(value$seed),
    homology_dimension = value$homology_dimension,
    metric = metric, value = as.numeric(value[[metric]]),
    stringsAsFactors = FALSE
  )))
  long <- long[order(long$dataset_scope, long$contrast_id,
                     long$homology_dimension, long$seed, long$metric,
                     method = "radix"), , drop = FALSE]
  row.names(long) <- NULL
  internal <- long[long$dataset_scope == "internal124", , drop = FALSE]
  keys <- unique(internal[c("contrast_id", "homology_dimension", "metric")])
  keys <- keys[order(keys$contrast_id, keys$homology_dimension, keys$metric,
                     method = "radix"), , drop = FALSE]
  internal_summary <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    hit <- internal$contrast_id == key$contrast_id &
      internal$homology_dimension == key$homology_dimension &
      internal$metric == key$metric
    x <- internal$value[hit]
    if (length(x) != 5L) stop("MV9 internal seed cardinality drift")
    data.frame(
      contract_id = "mv09_internal_seed_summary_v1",
      dataset_scope = "internal124", contrast_id = key$contrast_id,
      homology_dimension = key$homology_dimension, metric = key$metric,
      seeds = 5L, minimum = min(x), q25 = .mv09_quantile(x, .25),
      median = stats::median(x), q75 = .mv09_quantile(x, .75),
      maximum = max(x), inference = "none_descriptive_only",
      stringsAsFactors = FALSE
    )
  }))
  external <- long[long$dataset_scope == "external8", , drop = FALSE]
  external$contract_id <- "mv09_external_singleton_v1"
  external$replication_units <- 1L
  external$interpretation <- "single_external_cohort_no_variability_or_generalization"
  external <- external[c("contract_id", "dataset_scope", "contrast_id",
                         "seed", "homology_dimension", "metric", "value",
                         "replication_units", "interpretation")]

  identifiers <- unique(value[c("dataset_scope", "contrast_id", "seed")])
  identifiers <- identifiers[order(identifiers$dataset_scope,
                                   identifiers$contrast_id, identifiers$seed,
                                   method = "radix"), , drop = FALSE]
  deltas <- do.call(rbind, lapply(seq_len(nrow(identifiers)), function(i) {
    key <- identifiers[i, , drop = FALSE]
    hit <- value$dataset_scope == key$dataset_scope &
      value$contrast_id == key$contrast_id & value$seed == key$seed
    block <- value[hit, , drop = FALSE]
    if (nrow(block) != 2L || !setequal(block$homology_dimension, c("H0", "H1"))) {
      stop("MV9 H0/H1 pairing drift", call. = FALSE)
    }
    h0 <- block[block$homology_dimension == "H0", , drop = FALSE]
    h1 <- block[block$homology_dimension == "H1", , drop = FALSE]
    data.frame(
      contract_id = "mv09_dimension_delta_v1",
      dataset_scope = key$dataset_scope, contrast_id = key$contrast_id,
      seed = as.integer(key$seed), metric = metrics,
      h0_value = as.numeric(h0[1, metrics]),
      h1_value = as.numeric(h1[1, metrics]),
      h1_minus_h0 = as.numeric(h1[1, metrics]) - as.numeric(h0[1, metrics]),
      interpretation = "paired_descriptive_dimension_difference",
      stringsAsFactors = FALSE
    )
  }))
  delta_keys <- unique(deltas[c("dataset_scope", "contrast_id", "metric")])
  delta_summary <- do.call(rbind, lapply(seq_len(nrow(delta_keys)), function(i) {
    key <- delta_keys[i, , drop = FALSE]
    hit <- deltas$dataset_scope == key$dataset_scope &
      deltas$contrast_id == key$contrast_id & deltas$metric == key$metric
    x <- deltas$h1_minus_h0[hit]
    data.frame(
      contract_id = "mv09_dimension_delta_summary_v1",
      dataset_scope = key$dataset_scope, contrast_id = key$contrast_id,
      metric = key$metric, replication_units = length(x),
      minimum = min(x), q25 = .mv09_quantile(x, .25),
      median = stats::median(x), q75 = .mv09_quantile(x, .75),
      maximum = max(x), inference = "none_descriptive_only",
      stringsAsFactors = FALSE
    )
  }))
  list(
    aggregate = value,
    plot_data = long,
    internal_seed_summary = internal_summary,
    external_singleton = external,
    dimension_delta = deltas,
    dimension_delta_summary = delta_summary
  )
}
