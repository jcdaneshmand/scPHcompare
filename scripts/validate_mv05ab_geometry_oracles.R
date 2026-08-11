#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("usage: validate_mv05ab_geometry_oracles.R QUEUE D1_ROOT G_ROOT",
             "RESULT_ROOT OUTPUT_DIR TOLERANCE"), call. = FALSE)
}
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
d1_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
g_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
result_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[5L]]
tolerance <- as.numeric(args[[6L]])
if (nrow(queue) != 150L || !is.finite(tolerance) || tolerance <= 0) {
  stop("MV5-AB independent geometry inputs are invalid.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
write_once <- function(x, name) {
  path <- file.path(output_dir, name)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
close <- function(x, y) length(x) == length(y) &&
  all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y)))
within_mean <- function(x) 2 * sum(stats::dist(x)) / nrow(x)^2
energy_direct <- function(x, y) {
  squared <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  cross <- mean(sqrt(pmax(squared, 0)))
  sqrt(max(0, 2 * cross - within_mean(x) - within_mean(y)))
}
training <- unique(queue[c("fold_id", "training_samples")])
training <- training[order(training$training_samples, training$fold_id), ]
positions <- c(1L, floor((nrow(training) + 1L) / 2L), nrow(training))
selected_folds <- training$fold_id[positions]

group_rows <- vector("list", nrow(queue))
energy_rows <- list()
energy_cursor <- 0L
total_views <- 0L
total_h0_edges <- 0L
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  fold <- sub("^large_loso_v1:", "", unit$fold_id)
  seed <- as.integer(unit$seed)
  d1_path <- file.path(
    d1_root, paste0(fold, "__", seed, "__sct_cell_fold.rds"))
  g_stem <- paste0("mv05g_group__", fold, "__", seed)
  g_path <- file.path(g_root, "results", g_stem, paste0(g_stem, ".rds"))
  d1 <- readRDS(d1_path)
  g <- readRDS(g_path)
  raw <- if (unit$representation == "sct_whole") {
    lapply(d1$payload$cell_views, function(view) view$payload)
  } else {
    g$payload$coordinates
  }
  path <- file.path(result_root, safe(unit$robustness_group_id))
  metrics <- read_csv(file.path(path, "view_metrics.csv"))
  intervals <- read_csv(file.path(path, "finite_intervals.csv"))
  pairs <- read_csv(file.path(path, "pair_scope.csv"))
  energy <- read_csv(file.path(path, "energy_pairs.csv"))
  ids <- sort(names(raw), method = "radix")
  if (length(ids) != 90L || !setequal(ids, metrics$sample_id) ||
      !setequal(ids, unique(intervals$sample_id))) {
    stop("MV5-AB independent sample axis failed at group ", index,
         call. = FALSE)
  }
  normalized <- vector("list", length(ids)); names(normalized) <- ids
  maximum_mst_error <- 0
  minimum_source_norm <- Inf
  for (sample_id in ids) {
    x <- as.matrix(raw[[sample_id]])
    norms <- sqrt(rowSums(x^2))
    if (!identical(dim(x), c(384L, 30L)) || any(!is.finite(x)) ||
        any(!is.finite(norms)) || any(norms <= 0)) {
      stop("MV5-AB independent normalization failed: ", sample_id,
           call. = FALSE)
    }
    minimum_source_norm <- min(minimum_source_norm, norms)
    x <- x / norms
    normalized[[sample_id]] <- x
    if (max(abs(sqrt(rowSums(x^2)) - 1)) > tolerance) {
      stop("MV5-AB independent unit norm failed: ", sample_id,
           call. = FALSE)
    }
    expected <- sort(stats::hclust(stats::dist(x), method = "single")$height)
    observed <- sort(intervals$death[
      intervals$sample_id == sample_id &
        intervals$homology_dimension == "H0"])
    if (!close(expected, observed)) {
      stop("MV5-AB independent H0 MST failed: ", sample_id,
           call. = FALSE)
    }
    maximum_mst_error <- max(maximum_mst_error,
      max(abs(expected - observed)))
    total_views <- total_views + 1L
    total_h0_edges <- total_h0_edges + length(observed)
  }
  metric_min <- metrics$minimum_row_norm[match(ids, metrics$sample_id)]
  metric_max <- metrics$maximum_row_norm[match(ids, metrics$sample_id)]
  if (any(abs(metric_min - 1) > tolerance) ||
      any(abs(metric_max - 1) > tolerance)) {
    stop("MV5-AB stored unit-norm metric failed at group ", index,
         call. = FALSE)
  }
  if (unit$fold_id %in% selected_folds) {
    pair <- pairs[order(pairs$pair_ordinal), , drop = FALSE][1L, ]
    observed <- energy$distance[
      energy$pair_request_id == pair$pair_request_id]
    expected <- energy_direct(normalized[[pair$first_sample_id]],
                              normalized[[pair$second_sample_id]])
    if (length(observed) != 1L || !close(expected, observed)) {
      stop("MV5-AB independent energy oracle failed at group ", index,
           call. = FALSE)
    }
    energy_cursor <- energy_cursor + 1L
    energy_rows[[energy_cursor]] <- data.frame(
      contract_id = "mv05ab_independent_energy_oracle_v1",
      robustness_group_id = unit$robustness_group_id,
      fold_id = unit$fold_id, seed = seed,
      representation = unit$representation,
      pair_request_id = pair$pair_request_id,
      expected_distance = expected, observed_distance = observed,
      absolute_error = abs(expected - observed), passed = TRUE,
      production_scientific_helpers_called = FALSE,
      labels_opened = FALSE, rankings_computed = FALSE,
      outcomes_computed = FALSE, stringsAsFactors = FALSE)
  }
  group_rows[[index]] <- data.frame(
    contract_id = "mv05ab_independent_geometry_group_v1",
    robustness_group_id = unit$robustness_group_id,
    fold_id = unit$fold_id, seed = seed,
    representation = unit$representation, views = length(ids),
    h0_mst_edges = sum(intervals$homology_dimension == "H0"),
    minimum_source_row_norm = minimum_source_norm,
    maximum_mst_absolute_error = maximum_mst_error,
    normalization_passed = TRUE, all_view_h0_mst_passed = TRUE,
    production_scientific_helpers_called = FALSE,
    labels_opened = FALSE, rankings_computed = FALSE,
    outcomes_computed = FALSE, stringsAsFactors = FALSE)
  if (index %% 10L == 0L) message("MV5-AB independent geometry groups: ", index)
}
groups <- do.call(rbind, group_rows)
energies <- do.call(rbind, energy_rows)
if (nrow(groups) != 150L || total_views != 13500L ||
    total_h0_edges != 13500L * 383L || nrow(energies) != 30L ||
    length(unique(energies$fold_id)) != 3L ||
    length(unique(energies$seed)) != 5L ||
    length(unique(energies$representation)) != 2L) {
  stop("MV5-AB independent oracle totals failed.", call. = FALSE)
}
summary <- data.frame(
  contract_id = "mv05ab_independent_geometry_summary_v1",
  groups = nrow(groups), normalized_views = total_views,
  h0_mst_edges = total_h0_edges,
  maximum_mst_absolute_error = max(groups$maximum_mst_absolute_error),
  energy_oracles = nrow(energies), energy_fold_strata = 3L,
  energy_representations = 2L, energy_seeds = 5L,
  maximum_energy_absolute_error = max(energies$absolute_error),
  production_scientific_helpers_called = FALSE,
  labels_opened = FALSE, rankings_computed = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE)
write_once(groups, "mv05ab-independent-geometry-groups-2026-08-11.csv")
write_once(energies, "mv05ab-independent-energy-oracles-2026-08-11.csv")
write_once(summary, "mv05ab-independent-geometry-summary-2026-08-11.csv")
message("MV5-AB independent geometry passed 13,500 views and 30 energy oracles.")
