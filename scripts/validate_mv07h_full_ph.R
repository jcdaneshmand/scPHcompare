#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "ripserr", "TDA")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: validate_mv07h_full_ph.R PREFREEZE PRIVATE_ROOT PUBLIC_DIR CHECKS CROSS_ENGINE DECISION")
}
prefreeze <- args[[1L]]
private_root <- args[[2L]]
public_dir <- args[[3L]]
checks_output <- args[[4L]]
cross_output <- args[[5L]]
decision_output <- args[[6L]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")

axis <- read.csv(file.path(prefreeze, "mv07h-sample-seed-axis.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
source_queue <- read.csv(file.path(prefreeze, "mv07h-source-queue.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
ph_queue <- read.csv(file.path(prefreeze, "mv07h-ph-queue.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
sentinel_axis <- read.csv(file.path(prefreeze, "mv07h-sentinel-axis.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
source_metrics <- read.csv(file.path(public_dir, "mv07h-source-metrics.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
ph_metrics <- read.csv(file.path(public_dir, "mv07h-ph-metrics.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
engine_attempts <- read.csv(file.path(public_dir,
  "mv07h-ph-engine-attempts.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
equivalence <- read.csv(
  file.path(public_dir, "mv07h-sentinel-equivalence.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
repeat_validation <- read.csv(
  file.path(public_dir, "mv07h-repeat-validation.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
production_decision <- read.csv(
  file.path(public_dir, "mv07h-full-ph-decision.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)

source_ok <- logical(nrow(source_queue))
for (index in seq_len(nrow(source_queue))) {
  row <- source_queue[index, , drop = FALSE]
  path <- file.path(private_root, row$output_file)
  record <- readRDS(path)
  mv07h_validate_source_record_v1(record)
  metric <- source_metrics[source_metrics$seed == row$seed, , drop = FALSE]
  source_ok[[index]] <- nrow(metric) == 1L && length(record$views) == 124L &&
    record$identity$seed == row$seed &&
    metric$output_sha256 == .mv07h_sha256(path) &&
    metric$output_bytes == as.numeric(file.info(path)$size)
  rm(record); invisible(gc(FALSE))
}

ph_ok <- logical(nrow(ph_queue))
h0_ok <- logical(nrow(ph_queue))
h1_ok <- logical(nrow(ph_queue))
orientation_ok <- logical(nrow(ph_queue))
for (index in seq_len(nrow(ph_queue))) {
  row <- ph_queue[index, , drop = FALSE]
  path <- file.path(private_root, row$output_file)
  record <- readRDS(path)
  mv07h_validate_ph_record_v1(record)
  metric <- ph_metrics[ph_metrics$job_id == row$job_id, , drop = FALSE]
  expected_points <- if (row$view_id == "cell_topology_v1") 384L else 500L
  ph_ok[[index]] <- nrow(metric) == 1L &&
    metric$output_sha256 == .mv07h_sha256(path) &&
    metric$output_bytes == as.numeric(file.info(path)$size) &&
      record$identity$seed == row$seed &&
      record$identity$sample_id == row$sample_id &&
      record$identity$view_id == row$view_id &&
      metric$ph_engine == record$topology_result$provenance$ph_engine &&
      metric$ph_engine_version ==
        record$topology_result$provenance$ph_engine_version
  orientation_ok[[index]] <- record$identity$point_count == expected_points &&
    record$h0_mst_oracle$finite_h0_intervals == expected_points - 1L
  h0_ok[[index]] <- isTRUE(record$h0_mst_oracle$passed)
  h1 <- record$topology_result$diagram[
    record$topology_result$diagram[, "dimension"] == 1, , drop = FALSE]
  h1_ok[[index]] <- nrow(h1) > 0L && all(is.finite(h1[, "death"])) &&
    all(h1[, "death"] > h1[, "birth"])
  rm(record)
  if (index %% 100L == 0L) invisible(gc(FALSE))
}

normalize_diagram <- function(value) {
  result <- as.matrix(value); storage.mode(result) <- "double"
  colnames(result) <- c("dimension", "birth", "death"); result
}
remove_capped_essential <- function(diagram, point_count) {
  h0 <- which(diagram[, "dimension"] == 0)
  removed <- FALSE
  if (length(h0) == point_count) {
    diagram <- diagram[-h0[[which.max(diagram[h0, "death"])]],, drop = FALSE]
    removed <- TRUE
  }
  list(diagram = diagram, removed = removed)
}
compare_dimension <- function(first, second, dimension, tolerance = 1e-6) {
  first <- first[first[, "dimension"] == dimension,
                 c("birth", "death"), drop = FALSE]
  second <- second[second[, "dimension"] == dimension,
                   c("birth", "death"), drop = FALSE]
  order_rows <- function(value) value[order(
    value[, "birth"], value[, "death"], method = "radix"),, drop = FALSE]
  first <- order_rows(first); second <- order_rows(second)
  error <- if (nrow(first) == nrow(second) && nrow(first)) {
    max(abs(first - second))
  } else if (!nrow(first) && !nrow(second)) 0 else Inf
  list(first_count = nrow(first), second_count = nrow(second),
       maximum_absolute_error = error,
       passed = is.finite(error) && error <= tolerance)
}
cross_rows <- list()
for (seed in source_queue$seed) {
  source_record <- readRDS(file.path(
    private_root, "source", paste0("mv07h__", seed, "__source.rds")))
  excluded <- sentinel_axis$sample_id[sentinel_axis$seed == seed]
  sample_id <- sort(setdiff(names(source_record$views), excluded),
                    method = "radix")[[1L]]
  for (view_id in .mv07h_views) {
    view <- source_record$views[[sample_id]][[view_id]]
    if (view_id == "cell_topology_v1") {
      input <- view$payload[seq_len(32L),, drop = FALSE]
      maximum_scale <- max(stats::dist(input)); ripser_input <- input
      tda_input <- input; tda_dist <- "euclidean"
      subset_id <- "first_32_ordered_cells"
    } else {
      distances <- as.matrix(view$payload)[seq_len(32L), seq_len(32L),
                                           drop = FALSE]
      maximum_scale <- max(distances); ripser_input <- stats::as.dist(distances)
      tda_input <- distances; tda_dist <- "arbitrary"
      subset_id <- "first_32_ordered_genes"
    }
    ripser <- remove_capped_essential(normalize_diagram(
      ripserr::vietoris_rips(ripser_input, max_dim = 1L, threshold = -1,
                             p = 2L, return_format = "mat")), 32L)
    gudhi <- remove_capped_essential(normalize_diagram(TDA::ripsDiag(
      X = tda_input, maxdimension = 1L, maxscale = maximum_scale,
      dist = tda_dist, library = "GUDHI", location = FALSE,
      printProgress = FALSE)$diagram), 32L)
    for (dimension in 0:1) {
      comparison <- compare_dimension(ripser$diagram, gudhi$diagram, dimension)
      cross_rows[[length(cross_rows) + 1L]] <- data.frame(
        contract_id = "mv07h_new_sample_cross_engine_validation_v1",
        seed = seed, sample_id = sample_id, view_id = view_id,
        deterministic_subset = subset_id, subset_points = 32L,
        homology_dimension = paste0("H", dimension),
        ripserr_intervals = comparison$first_count,
        gudhi_intervals = comparison$second_count,
        ripserr_capped_essential_h0_removed = ripser$removed,
        gudhi_capped_essential_h0_removed = gudhi$removed,
        maximum_absolute_error = comparison$maximum_absolute_error,
        tolerance = 1e-6, passed = comparison$passed,
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  rm(source_record); invisible(gc(FALSE))
}
cross_engine <- do.call(rbind, cross_rows)
fallback_metrics <- ph_metrics[
  ph_metrics$ph_engine == "TDA_ripsDiag_GUDHI",, drop = FALSE]
default_repeat_jobs <- c(
  paste0("source__", min(source_queue$seed)),
  ph_queue$job_id[
    ph_queue$seed == min(source_queue$seed) &
      ph_queue$sample_id %in% sort(unique(sentinel_axis$sample_id[
        sentinel_axis$seed == min(source_queue$seed)]), method = "radix")]
)
expected_repeat_jobs <- unique(c(default_repeat_jobs,
                                 fallback_metrics$job_id))
fallback_ok <- if (!nrow(fallback_metrics)) TRUE else all(vapply(
  fallback_metrics$job_id, function(job_id) {
    primary <- engine_attempts[
      engine_attempts$job_id == job_id &
        engine_attempts$attempt_scope == "primary",, drop = FALSE]
    fallback <- engine_attempts[
      engine_attempts$job_id == job_id &
        engine_attempts$attempt_scope == "fallback",, drop = FALSE]
    repeated <- repeat_validation[repeat_validation$job_id == job_id,,
                                  drop = FALSE]
    nrow(primary) == 1L && primary$ph_engine == "ripserr" &&
      primary$disposition == "rss_cap_exceeded" &&
      primary$peak_process_tree_rss_bytes > primary$rss_cap_bytes &&
      nrow(fallback) == 1L &&
      fallback$ph_engine == "TDA_ripsDiag_GUDHI" &&
      fallback$disposition == "completed" &&
      fallback$peak_process_tree_rss_bytes <= fallback$rss_cap_bytes &&
      fallback$rss_cap_bytes == 12 * 1024^3 &&
      nrow(repeated) == 1L && repeated$sha256_equal
  }, logical(1L)))
checks <- data.frame(
  contract_id = "mv07h_full_ph_independent_validation_v1",
  category = c("axis", "source_bundles", "typed_views", "ph_records",
               "point_orientation", "h0_mst", "h1_intervals",
               "sentinel_equivalence", "deterministic_repeat",
               "new_sample_cross_engine", "exact_resource_fallback",
               "resource_firewall"),
  passed = c(
    nrow(axis) == 620L && nrow(source_queue) == 5L && nrow(ph_queue) == 1240L,
    length(source_ok) == 5L && all(source_ok),
    sum(source_metrics$typed_views) == 1240L,
    length(ph_ok) == 1240L && all(ph_ok),
    all(orientation_ok), all(h0_ok), all(h1_ok),
    nrow(equivalence) == 60L && all(as.logical(equivalence$passed)),
    nrow(repeat_validation) == length(expected_repeat_jobs) &&
      setequal(repeat_validation$job_id, expected_repeat_jobs) &&
      all(as.logical(repeat_validation$bytes_equal)) &&
      all(as.logical(repeat_validation$sha256_equal)),
    nrow(cross_engine) == 20L && all(cross_engine$passed),
    all(ph_metrics$ph_engine %in% c("ripserr", "TDA_ripsDiag_GUDHI")) &&
      all(ph_metrics$ph_engine[ph_metrics$view_id == "cell_topology_v1"] ==
            "ripserr") && fallback_ok &&
      production_decision$gudhi_fallback_records == nrow(fallback_metrics) &&
      production_decision$rss_triggered_fallbacks == nrow(fallback_metrics),
    production_decision$decision ==
      "full_PH_complete_await_independent_validation" &&
      production_decision$aggregate_elapsed_seconds <=
        production_decision$aggregate_elapsed_cap_seconds &&
      production_decision$private_bytes <=
        production_decision$private_storage_cap_bytes &&
      production_decision$landscape_jobs == 0L &&
      production_decision$clustering_jobs == 0L &&
      production_decision$label_jobs == 0L &&
      production_decision$outcome_jobs == 0L
  ),
  detail = c("124 by five, two views", "five exact reused transforms",
             "1,240 typed views", "1,240 corrected diagrams",
             "384 cells or 500 genes", "all finite H0 deaths match MST",
             "all H1 finite positive and nonempty", "60 exact MV7-G views",
             "source sentinels plus every fallback", "20 Ripserr/GUDHI checks",
             "Ripserr primary; exact GUDHI only after RSS failure",
             "within caps; landscapes labels outcomes still closed"),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) stop("MV7-H full-PH validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "))
decision <- data.frame(
  contract_id = "mv07h_full_ph_validation_decision_v1",
  decision = "authorize_one_MV7H_landscape_stress_group",
  source_jobs = 5L, typed_views = 1240L, ph_jobs = 1240L,
  cross_engine_checks = 20L, repeat_artifacts = nrow(repeat_validation),
  ripserr_selected_records = sum(ph_metrics$ph_engine == "ripserr"),
  gudhi_fallback_records = nrow(fallback_metrics),
  landscape_groups_authorized = 1L, landscape_groups_closed = 19L,
  clustering_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write.csv(checks, checks_output, row.names = FALSE, na = "")
write.csv(cross_engine, cross_output, row.names = FALSE, na = "")
write.csv(decision, decision_output, row.names = FALSE, na = "")
message("MV7-H full-PH independent validation: 12/12 pass")
