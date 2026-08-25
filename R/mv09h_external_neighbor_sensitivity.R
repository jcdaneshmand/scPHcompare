# Label-closed external-neighborhood sensitivity for the eight-unit HCA cohort.

mv09h_neighbor_degeneracy_v1 <- function(n_units = 8L, k = 7L) {
  n_units <- as.integer(n_units)
  k <- as.integer(k)
  if (length(n_units) != 1L || is.na(n_units) || n_units < 3L ||
      length(k) != 1L || is.na(k) || k < 1L || k > n_units - 1L) {
    stop("invalid neighborhood-degeneracy arguments", call. = FALSE)
  }
  data.frame(
    contract_id = "mv09h_neighbor_degeneracy_v1",
    units = n_units,
    k = k,
    k_equals_all_other_units = k == n_units - 1L,
    possible_neighbor_sets_per_unit = if (k == n_units - 1L) 1 else
      choose(n_units - 1L, k),
    jaccard_for_any_two_complete_rankings = if (k == n_units - 1L) 1 else NA_real_,
    informative_for_neighborhood_preservation = k < n_units - 1L,
    disposition = if (k == n_units - 1L)
      "structurally_noninformative_exclude_from_interpretation" else
        "potentially_informative",
    stringsAsFactors = FALSE
  )
}

mv09h_external_neighbor_sensitivity_v1 <- function(
    left, right, comparison_id, ks = c(2L, 3L)) {
  if (!exists(".mv08zy_validate_pairs", mode = "function") ||
      !exists(".mv08zy_neighbor_sets", mode = "function")) {
    stop("MV8-ZY distance-comparison helpers must be loaded", call. = FALSE)
  }
  left <- .mv08zy_validate_pairs(left)
  right <- .mv08zy_validate_pairs(right)
  if (!identical(left[c("first_unit_id", "second_unit_id", "pair_key")],
                 right[c("first_unit_id", "second_unit_id", "pair_key")])) {
    stop("external sensitivity inputs have different pair axes", call. = FALSE)
  }
  units <- sort(unique(c(left$first_unit_id, left$second_unit_id)),
                method = "radix")
  ks <- as.integer(ks)
  comparison_id <- as.character(comparison_id)
  if (length(units) != 8L || !identical(ks, c(2L, 3L)) ||
      length(comparison_id) != 1L || is.na(comparison_id) ||
      !nzchar(comparison_id)) {
    stop("external neighborhood contract drift", call. = FALSE)
  }
  rows <- lapply(ks, function(k) {
    left_sets <- .mv08zy_neighbor_sets(left, k)
    right_sets <- .mv08zy_neighbor_sets(right, k)
    overlap <- vapply(units, function(unit) {
      shared <- intersect(left_sets[[unit]], right_sets[[unit]])
      combined <- union(left_sets[[unit]], right_sets[[unit]])
      length(shared) / length(combined)
    }, numeric(1L))
    data.frame(
      contract_id = "mv09h_external_neighbor_unit_v1",
      comparison_id = comparison_id,
      unit_id = units,
      k = k,
      neighbor_jaccard = overlap,
      stringsAsFactors = FALSE
    )
  })
  unit <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(unit, unit$k), function(value) {
    data.frame(
      contract_id = "mv09h_external_neighbor_summary_v1",
      comparison_id = comparison_id,
      units = length(units),
      unordered_pairs = nrow(left),
      k = unique(value$k),
      mean_neighbor_jaccard = mean(value$neighbor_jaccard),
      median_neighbor_jaccard = stats::median(value$neighbor_jaccard),
      p10_neighbor_jaccard = as.numeric(stats::quantile(
        value$neighbor_jaccard, 0.10, names = FALSE, type = 7
      )),
      replication_units = 1L,
      interpretation = "descriptive_sensitivity_no_equivalence_or_biology",
      outcome_label_state = "closed",
      biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  rownames(unit) <- NULL
  list(summary = summary, unit = unit,
       pair_axis = left[c("first_unit_id", "second_unit_id", "pair_key")],
       degeneracy = mv09h_neighbor_degeneracy_v1(length(units),
                                                  length(units) - 1L))
}
