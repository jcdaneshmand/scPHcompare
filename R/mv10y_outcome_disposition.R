# Transparent, value-aware descriptive disposition for closed MV10 outcomes.

.mv10y_group_summary <- function(data, columns, prefix, expected_groups) {
  key <- interaction(data[columns], drop = TRUE, lex.order = TRUE)
  result <- do.call(rbind, lapply(split(data, key), function(x) {
    out <- x[1L, columns, drop = FALSE]
    out$values <- nrow(x)
    out$minimum <- min(x$seed_mean)
    out$median <- stats::median(x$seed_mean)
    out$maximum <- max(x$seed_mean)
    out$range <- out$maximum - out$minimum
    out
  }))
  rownames(result) <- NULL
  if (nrow(result) != expected_groups) {
    stop(prefix, " disposition cardinality drift", call. = FALSE)
  }
  names(result)[names(result) %in% c("minimum", "median", "maximum", "range")] <-
    paste0(prefix, "_", c("minimum", "median", "maximum", "range"))
  result
}

mv10y_build_outcome_disposition_v1 <- function(complete_summary,
                                                endpoint_coverage) {
  if (!is.data.frame(complete_summary) || nrow(complete_summary) != 300L ||
      !is.data.frame(endpoint_coverage) || nrow(endpoint_coverage) != 6L ||
      any(!is.finite(complete_summary$seed_mean)) ||
      any(complete_summary$inference_performed) ||
      any(complete_summary$ranking_performed) ||
      any(complete_summary$biological_claim)) {
    stop("MV10-Y source contract failed", call. = FALSE)
  }
  primary <- complete_summary[
    complete_summary$method_id == "pam_dissimilarity_v1", , drop = FALSE
  ]
  primary_envelope <- .mv10y_group_summary(
    primary, c("endpoint_id", "population_id", "label_axis",
               "homology_dimension", "metric_id"),
    "representation", 20L
  )
  method_envelope <- .mv10y_group_summary(
    complete_summary,
    c("endpoint_id", "population_id", "label_axis", "stack_id",
      "homology_dimension", "metric_id"), "method", 60L
  )
  identity <- c("stack_id", "homology_dimension", "method_id", "metric_id")
  contrast <- function(left, right, left_name, right_name, axes) {
    first <- complete_summary[
      complete_summary$endpoint_id == left, c(identity, "seed_mean"), drop = FALSE
    ]
    second <- complete_summary[
      complete_summary$endpoint_id == right, c(identity, "seed_mean"), drop = FALSE
    ]
    names(first)[names(first) == "seed_mean"] <- "left_mean"
    names(second)[names(second) == "seed_mean"] <- "right_mean"
    merged <- merge(first, second, by = identity, sort = TRUE)
    merged$left_endpoint <- left_name
    merged$right_endpoint <- right_name
    merged$contrast_axis <- axes
    merged$descriptive_difference <- merged$left_mean - merged$right_mean
    merged$inference_performed <- FALSE
    merged$causal_interpretation_allowed <- FALSE
    merged
  }
  context_contrast <- rbind(
    contrast(
      "primary90_context_restriction__tissue", "full124_descriptive__tissue",
      "primary90_tissue", "full124_tissue", "population_restriction"
    ),
    contrast(
      "primary90_context_restriction__study", "full124_descriptive__study",
      "primary90_study", "full124_study", "population_restriction"
    )
  )
  approach_contrast <- rbind(
    contrast(
      "full124_descriptive__canonical_approach", "full124_descriptive__tissue",
      "full124_approach", "full124_tissue", "confounded_proxy_minus_tissue"
    ),
    contrast(
      "full124_descriptive__canonical_approach", "full124_descriptive__study",
      "full124_approach", "full124_study", "confounded_proxy_minus_study"
    )
  )
  if (nrow(context_contrast) != 120L || nrow(approach_contrast) != 120L) {
    stop("MV10-Y contrast cardinality drift", call. = FALSE)
  }
  disposition <- data.frame(
    contract_id = "mv10z_outcome_disposition_v1",
    decision = "retain_frozen_PAM_without_outcome_tuning",
    outcome_role = "descriptive_external_alignment_only",
    primary_PAM_rows = nrow(primary), complete_method_rows = nrow(complete_summary),
    magnitude_threshold_applied = FALSE, p_values_computed = FALSE,
    method_selection_executed = FALSE, representation_ranking_executed = FALSE,
    H0_H1_combined = FALSE, approach_causal_interpretation = FALSE,
    biological_claims = FALSE, manuscript_claims = FALSE,
    next_stage = "cross_view_descriptive_synthesis_prefreeze_or_stop",
    stringsAsFactors = FALSE
  )
  list(
    primary_envelope = primary_envelope,
    method_envelope = method_envelope,
    context_contrast = context_contrast,
    approach_contrast = approach_contrast,
    disposition = disposition
  )
}
