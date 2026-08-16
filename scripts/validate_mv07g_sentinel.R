#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "ripserr", "TDA")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: validate_mv07g_sentinel.R PREFREEZE PRIVATE_ROOT PUBLIC_DIR CHECKS CROSS_ENGINE DECISION")
}
prefreeze <- args[[1]]
private_root <- args[[2]]
public_dir <- args[[3]]
checks_output <- args[[4]]
cross_output <- args[[5]]
decision_output <- args[[6]]
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
queue <- read.csv(file.path(prefreeze, "mv07g-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
axis <- read.csv(file.path(prefreeze, "mv07g-sentinel-axis.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
source_metrics <- read.csv(file.path(public_dir, "mv07g-source-metrics.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
ph_metrics <- read.csv(file.path(public_dir, "mv07g-ph-metrics.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
repeat_validation <- read.csv(
  file.path(public_dir, "mv07g-repeat-validation.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
projection <- read.csv(file.path(public_dir, "mv07g-full-ph-projection.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
production_decision <- read.csv(file.path(public_dir, "mv07g-decision.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE)
source_records <- lapply(source_metrics$seed, function(seed) {
  path <- file.path(private_root, "source",
                    paste0("mv07g__", seed, "__source.rds"))
  if (!file.exists(path)) stop("Missing MV7-G source bundle: ", seed)
  record <- readRDS(path)
  mv07g_validate_source_record_v1(record)
  record
})
names(source_records) <- as.character(source_metrics$seed)
source_hash_ok <- vapply(seq_len(nrow(source_metrics)), function(index) {
  path <- file.path(private_root, "source",
    paste0("mv07g__", source_metrics$seed[[index]], "__source.rds"))
  sha(path) == source_metrics$output_sha256[[index]] &&
    as.numeric(file.info(path)$size) == source_metrics$output_bytes[[index]]
}, logical(1L))
ph_records <- lapply(seq_len(nrow(ph_metrics)), function(index) {
  row <- ph_metrics[index,]
  path <- file.path(private_root, "ph", paste0(
    "mv07g__", row$seed, "__", row$sample_id, "__", row$view_id,
    "__ph.rds"))
  if (!file.exists(path)) stop("Missing MV7-G PH record: ", row$job_id)
  record <- readRDS(path)
  source_record <- source_records[[as.character(row$seed)]]
  view <- source_record$views[[row$sample_id]][[row$view_id]]
  mv07g_validate_ph_record_v1(record, view)
  if (sha(path) != row$output_sha256 ||
      as.numeric(file.info(path)$size) != row$output_bytes) {
    stop("MV7-G PH public/private artifact mismatch: ", row$job_id)
  }
  record
})
normalize_diagram <- function(value) {
  result <- as.matrix(value)
  storage.mode(result) <- "double"
  colnames(result) <- c("dimension", "birth", "death")
  result
}
remove_capped_essential <- function(diagram, point_count) {
  h0 <- which(diagram[, "dimension"] == 0)
  removed <- FALSE
  if (length(h0) == point_count) {
    diagram <- diagram[-h0[[which.max(diagram[h0, "death"])]], , drop = FALSE]
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
    value[, "birth"], value[, "death"], method = "radix"), , drop = FALSE]
  first <- order_rows(first)
  second <- order_rows(second)
  error <- if (nrow(first) == nrow(second) && nrow(first)) {
    max(abs(first - second))
  } else if (!nrow(first) && !nrow(second)) 0 else Inf
  list(first_count = nrow(first), second_count = nrow(second),
       maximum_absolute_error = error,
       passed = is.finite(error) && error <= tolerance)
}
cross_rows <- list()
cross_index <- 0L
cross_seed <- min(axis$seed)
source_record <- source_records[[as.character(cross_seed)]]
for (sample_id in source_record$identity$sentinel_ids) {
  for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
    view <- source_record$views[[sample_id]][[view_id]]
    if (view_id == "cell_topology_v1") {
      input <- view$payload[seq_len(32L), , drop = FALSE]
      maximum_scale <- max(stats::dist(input))
      ripser_input <- input
      tda_input <- input
      tda_dist <- "euclidean"
      subset_id <- "first_32_ordered_cells"
    } else {
      distances <- as.matrix(view$payload)[seq_len(32L), seq_len(32L),
                                           drop = FALSE]
      maximum_scale <- max(distances)
      ripser_input <- stats::as.dist(distances)
      tda_input <- distances
      tda_dist <- "arbitrary"
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
      comparison <- compare_dimension(ripser$diagram, gudhi$diagram,
                                      dimension)
      cross_index <- cross_index + 1L
      cross_rows[[cross_index]] <- data.frame(
        contract_id = "mv07g_reduced_cross_engine_validation_v1",
        seed = cross_seed, sample_id = sample_id, view_id = view_id,
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
}
cross_engine <- do.call(rbind, cross_rows)
full_hours <- sum(projection$projected_worker_hours)
full_bytes <- sum(projection$projected_private_bytes)
typed_ok <- all(vapply(source_records, function(record) {
  length(record$views) == 6L &&
    all(vapply(record$views, function(pair) {
      identical(sort(names(pair)),
                sort(c("cell_topology_v1", "gene_topology_v1")))
    }, logical(1L)))
}, logical(1L)))
point_ok <- all(vapply(ph_records, function(record) {
  expected <- if (record$identity$view_id == "cell_topology_v1") 384L else 500L
  record$identity$point_count == expected &&
    record$h0_mst_oracle$finite_h0_intervals == expected - 1L
}, logical(1L)))
h1_ok <- all(vapply(ph_records, function(record) {
  h1 <- record$topology_result$diagram[
    record$topology_result$diagram[, "dimension"] == 1, , drop = FALSE]
  all(is.finite(h1[, "death"])) && all(h1[, "death"] > h1[, "birth"])
}, logical(1L)))
checks <- data.frame(
  contract_id = "mv07g_independent_validation_v1",
  category = c("queue_axis", "source_bundles", "typed_views", "ph_records",
               "point_orientation", "h0_mst", "h1_intervals",
               "cross_engine", "deterministic_repeat", "resource_projection",
               "label_landscape_firewall"),
  passed = c(
    nrow(queue) == 65L && nrow(axis) == 30L && all(table(axis$seed) == 6L),
    length(source_records) == 5L && all(source_hash_ok),
    typed_ok, length(ph_records) == 60L,
    point_ok,
    all(vapply(ph_records, function(record) isTRUE(record$h0_mst_oracle$passed),
               logical(1L))),
    h1_ok, nrow(cross_engine) == 24L && all(cross_engine$passed),
    nrow(repeat_validation) == 13L &&
      all(repeat_validation$bytes_equal) && all(repeat_validation$sha256_equal),
    full_hours <= 24 && full_bytes <= 8 * 1024^3 &&
      production_decision$aggregate_elapsed_seconds <=
        production_decision$aggregate_elapsed_cap_seconds &&
      production_decision$private_bytes <=
        production_decision$private_storage_cap_bytes,
    production_decision$landscape_jobs == 0L &&
      production_decision$distance_jobs == 0L &&
      production_decision$clustering_jobs == 0L &&
      production_decision$label_jobs == 0L &&
      production_decision$outcome_jobs == 0L &&
      production_decision$outcome_label_state == "closed"
  ),
  detail = c("5 fits plus 60 PH jobs", "five exact global transforms",
             "60 typed cell/gene views", "60 corrected diagrams",
             "384 cells and 500 genes", "all finite H0 deaths match MST",
             "finite positive H1", "24 Ripserr/GUDHI checks",
             "one seed, 13 exact artifacts", "under 24h and 8GiB projection",
             "landscapes labels outcomes remain closed"),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) {
  stop("MV7-G validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "))
}
decision <- data.frame(
  contract_id = "mv07g_validation_decision_v1",
  decision = "authorize_MV7H_full_PH_landscape_prefreeze",
  panel_sha256 = unique(vapply(source_records, function(record)
    record$identity$panel_sha256, character(1L))),
  source_jobs = 5L, typed_views = 60L, ph_jobs = 60L,
  cross_engine_checks = 24L, repeat_artifacts = 13L,
  projected_full_worker_hours = full_hours,
  projected_full_private_bytes = full_bytes,
  landscape_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write.csv(checks, checks_output, row.names = FALSE, na = "")
write.csv(cross_engine, cross_output, row.names = FALSE, na = "")
write.csv(decision, decision_output, row.names = FALSE, na = "")
message("MV7-G independent validation: 11/11 pass")
