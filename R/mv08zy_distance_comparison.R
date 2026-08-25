# Internal helpers for label-closed MV8-ZY distance comparisons.

.mv08zy_validate_pairs <- function(value, expected_units = NULL) {
  required <- c("first_unit_id", "second_unit_id", "distance")
  if (!is.data.frame(value) || !all(required %in% names(value)) || !nrow(value)) {
    stop("distance pairs lack the required columns", call. = FALSE)
  }
  first <- as.character(value$first_unit_id)
  second <- as.character(value$second_unit_id)
  distance <- as.numeric(value$distance)
  if (anyNA(first) || anyNA(second) || any(!nzchar(first)) ||
      any(!nzchar(second)) || any(first == second) || anyNA(distance) ||
      any(!is.finite(distance)) || any(distance < 0)) {
    stop("distance pairs contain invalid identities or values", call. = FALSE)
  }
  left <- ifelse(first < second, first, second)
  right <- ifelse(first < second, second, first)
  key <- paste(left, right, sep = "\r")
  if (anyDuplicated(key)) stop("duplicated unordered pair", call. = FALSE)
  units <- sort(unique(c(left, right)), method = "radix")
  if (!is.null(expected_units) && !identical(units, sort(as.character(
    expected_units), method = "radix"))) {
    stop("distance-pair unit axis differs from the frozen axis", call. = FALSE)
  }
  expected <- length(units) * (length(units) - 1L) / 2L
  if (nrow(value) != expected) stop("distance pairs are incomplete", call. = FALSE)
  ordering <- order(left, right, method = "radix")
  data.frame(
    first_unit_id = left[ordering], second_unit_id = right[ordering],
    pair_key = key[ordering], distance = distance[ordering],
    stringsAsFactors = FALSE
  )
}

.mv08zy_neighbor_sets <- function(pairs, k) {
  units <- sort(unique(c(pairs$first_unit_id, pairs$second_unit_id)),
                method = "radix")
  result <- setNames(vector("list", length(units)), units)
  for (unit in units) {
    selected <- pairs$first_unit_id == unit | pairs$second_unit_id == unit
    neighbor <- ifelse(pairs$first_unit_id[selected] == unit,
                       pairs$second_unit_id[selected],
                       pairs$first_unit_id[selected])
    distance <- pairs$distance[selected]
    ordering <- order(distance, neighbor, method = "radix")
    result[[unit]] <- neighbor[ordering][seq_len(min(k, length(ordering)))]
  }
  result
}

mv08zy_compare_distance_pairs_v1 <- function(left, right, comparison_id) {
  left <- .mv08zy_validate_pairs(left)
  right <- .mv08zy_validate_pairs(right)
  if (!identical(left[c("first_unit_id", "second_unit_id", "pair_key")],
                 right[c("first_unit_id", "second_unit_id", "pair_key")])) {
    stop("comparison inputs do not have identical unordered-pair axes",
         call. = FALSE)
  }
  comparison_id <- as.character(comparison_id)
  if (length(comparison_id) != 1L || is.na(comparison_id) ||
      !nzchar(comparison_id)) stop("comparison_id must be non-empty")
  x <- left$distance
  y <- right$distance
  left_median <- stats::median(x)
  right_median <- stats::median(y)
  if (!is.finite(left_median) || !is.finite(right_median) ||
      left_median <= sqrt(.Machine$double.eps) ||
      right_median <= sqrt(.Machine$double.eps) ||
      sum(x ^ 2) <= .Machine$double.eps ||
      sum(y ^ 2) <= .Machine$double.eps ||
      stats::sd(x) <= sqrt(.Machine$double.eps) ||
      stats::sd(y) <= sqrt(.Machine$double.eps)) {
    stop("comparison distance stack is degenerate", call. = FALSE)
  }
  scale <- max(0, sum(x * y) / sum(x ^ 2))
  stress <- sqrt(sum((y - scale * x) ^ 2) / sum(y ^ 2))
  scaled_change <- abs(y / right_median - x / left_median)
  units <- sort(unique(c(left$first_unit_id, left$second_unit_id)),
                method = "radix")
  k <- min(10L, length(units) - 1L)
  left_sets <- .mv08zy_neighbor_sets(left, k)
  right_sets <- .mv08zy_neighbor_sets(right, k)
  overlap <- vapply(units, function(unit) {
    shared <- intersect(left_sets[[unit]], right_sets[[unit]])
    union <- union(left_sets[[unit]], right_sets[[unit]])
    length(shared) / length(union)
  }, numeric(1L))
  neighbor <- data.frame(
    comparison_id = comparison_id, unit_id = units, k = k,
    neighbor_jaccard = overlap, stringsAsFactors = FALSE
  )
  summary <- data.frame(
    contract_id = "mv08zy_distance_comparison_summary_v1",
    comparison_id = comparison_id, units = length(units),
    unordered_pairs = nrow(left), distance_transform = "sqrt_exact_squared_L2",
    pearson = stats::cor(x, y, method = "pearson"),
    spearman = stats::cor(x, y, method = "spearman"),
    nonnegative_scale = scale, relative_stress = stress,
    left_median_distance = left_median,
    right_median_distance = right_median,
    median_abs_median_scaled_change = stats::median(scaled_change),
    p95_abs_median_scaled_change = as.numeric(stats::quantile(
      scaled_change, 0.95, names = FALSE, type = 7
    )),
    neighbor_k = k, mean_neighbor_jaccard = mean(overlap),
    median_neighbor_jaccard = stats::median(overlap),
    p10_neighbor_jaccard = as.numeric(stats::quantile(
      overlap, 0.10, names = FALSE, type = 7
    )),
    interpretation = "descriptive_no_equivalence_or_biological_claim",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
    outcome_jobs = 0L, stringsAsFactors = FALSE
  )
  list(summary = summary, neighbor = neighbor,
       pair_axis = left[c("first_unit_id", "second_unit_id", "pair_key")])
}
