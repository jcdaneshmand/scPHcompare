#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "ripserr")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-AJ independent validation.",
         call. = FALSE)
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("usage: validate_mv05aj_geometry_oracles.R QUEUE D1_ROOT G_ROOT",
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
  stop("MV5-AJ independent geometry inputs are invalid.", call. = FALSE)
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
point_order <- function(sample_id, seed, point_ids) {
  hashes <- vapply(point_ids, function(point_id) {
    digest::digest(list(
      contract_id = "mv05t_nested_point_order_v1",
      sample_id = sample_id, seed = as.integer(seed), point_id = point_id
    ), algo = "sha256", serialize = TRUE)
  }, character(1L))
  point_ids[order(hashes, point_ids, method = "radix")]
}
scientific_sha <- function(value) digest::digest(
  value, algo = "sha256", serialize = TRUE
)
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
h1_rows <- list()
h1_cursor <- 0L
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
    stop("MV5-AJ independent sample axis failed at group ", index,
         call. = FALSE)
  }
  selected <- vector("list", length(ids)); names(selected) <- ids
  maximum_mst_error <- 0
  maximum_norm_error <- 0
  selected_pair_ids <- if (unit$fold_id %in% selected_folds) {
    pair <- pairs[order(pairs$pair_ordinal), , drop = FALSE][1L, ]
    c(pair$first_sample_id, pair$second_sample_id)
  } else character()
  for (sample_id in ids) {
    x <- as.matrix(raw[[sample_id]])
    if (!identical(dim(x), c(384L, 30L)) || any(!is.finite(x)) ||
        !identical(rownames(x), rownames(raw[[sample_id]]))) {
      stop("MV5-AJ independent source view failed: ", sample_id,
           call. = FALSE)
    }
    ordering <- point_order(sample_id, seed, rownames(x))
    nested_192_ids <- ordering[seq_len(192L)]
    selected_ids <- ordering[seq_len(256L)]
    if (!identical(selected_ids[seq_len(192L)], nested_192_ids)) {
      stop("MV5-AJ independent nested inclusion failed: ", sample_id,
           call. = FALSE)
    }
    x <- x[selected_ids, seq_len(30L), drop = FALSE]
    selected[[sample_id]] <- x
    metric <- metrics[metrics$sample_id == sample_id, , drop = FALSE]
    norms <- sqrt(rowSums(x^2))
    if (nrow(metric) != 1L ||
        metric$selected_point_ids_sha256 != scientific_sha(selected_ids) ||
        metric$nested192_point_ids_sha256 != scientific_sha(nested_192_ids) ||
        !as.logical(metric$nested192_is_exact_prefix) ||
        metric$transformed_payload_sha256 != scientific_sha(x) ||
        metric$point_selection != "sha256_sample_seed_cell_nested_v1") {
      stop("MV5-AJ independent selected identity failed: ", sample_id,
           call. = FALSE)
    }
    maximum_norm_error <- max(maximum_norm_error,
      abs(metric$minimum_row_norm - min(norms)),
      abs(metric$maximum_row_norm - max(norms)))
    expected <- sort(stats::hclust(stats::dist(x), method = "single")$height)
    observed <- sort(intervals$death[
      intervals$sample_id == sample_id &
        intervals$homology_dimension == "H0"])
    if (!close(expected, observed)) {
      stop("MV5-AJ independent H0 MST failed: ", sample_id,
           call. = FALSE)
    }
    maximum_mst_error <- max(maximum_mst_error,
      max(abs(expected - observed)))
    total_views <- total_views + 1L
    total_h0_edges <- total_h0_edges + length(observed)
    if (sample_id %in% selected_pair_ids) {
      direct <- as.matrix(ripserr::vietoris_rips(
        dataset = x, max_dim = 1L, threshold = -1, p = 2L,
        return_format = "mat"
      ))
      if (ncol(direct) != 3L) {
        stop("MV5-AJ independent H1 diagram shape failed.", call. = FALSE)
      }
      colnames(direct) <- c("dimension", "birth", "death")
      expected_h1 <- direct[
        direct[, "dimension"] == 1 & is.finite(direct[, "death"]),
        c("birth", "death"), drop = FALSE
      ]
      observed_h1 <- as.matrix(intervals[
        intervals$sample_id == sample_id &
          intervals$homology_dimension == "H1",
        c("birth", "death"), drop = FALSE
      ])
      expected_h1 <- expected_h1[order(expected_h1[, 1L], expected_h1[, 2L]),
                                 , drop = FALSE]
      observed_h1 <- observed_h1[order(observed_h1[, 1L], observed_h1[, 2L]),
                                 , drop = FALSE]
      if (!identical(dim(expected_h1), dim(observed_h1)) ||
          !close(as.numeric(expected_h1), as.numeric(observed_h1))) {
        stop("MV5-AJ independent H1 diagram failed: ", sample_id,
             call. = FALSE)
      }
      h1_cursor <- h1_cursor + 1L
      h1_rows[[h1_cursor]] <- data.frame(
        contract_id = "mv05aj_independent_h1_diagram_oracle_v1",
        robustness_group_id = unit$robustness_group_id,
        fold_id = unit$fold_id, seed = seed,
        representation = unit$representation, sample_id = sample_id,
        h1_intervals = nrow(expected_h1), exact_match = TRUE,
        production_scientific_helpers_called = FALSE,
        labels_opened = FALSE, rankings_computed = FALSE,
        outcomes_computed = FALSE, stringsAsFactors = FALSE)
    }
  }
  metric_min <- metrics$minimum_row_norm[match(ids, metrics$sample_id)]
  metric_max <- metrics$maximum_row_norm[match(ids, metrics$sample_id)]
  if (any(!is.finite(metric_min)) || any(!is.finite(metric_max)) ||
      maximum_norm_error > tolerance) {
    stop("MV5-AJ stored Euclidean norm metric failed at group ", index,
         call. = FALSE)
  }
  if (unit$fold_id %in% selected_folds) {
    pair <- pairs[order(pairs$pair_ordinal), , drop = FALSE][1L, ]
    observed <- energy$distance[
      energy$pair_request_id == pair$pair_request_id]
    expected <- energy_direct(selected[[pair$first_sample_id]],
                              selected[[pair$second_sample_id]])
    if (length(observed) != 1L || !close(expected, observed)) {
      stop("MV5-AJ independent energy oracle failed at group ", index,
           call. = FALSE)
    }
    energy_cursor <- energy_cursor + 1L
    energy_rows[[energy_cursor]] <- data.frame(
      contract_id = "mv05aj_independent_energy_oracle_v1",
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
    contract_id = "mv05aj_independent_geometry_group_v1",
    robustness_group_id = unit$robustness_group_id,
    fold_id = unit$fold_id, seed = seed,
    representation = unit$representation, views = length(ids),
    h0_mst_edges = sum(intervals$homology_dimension == "H0"),
    selected_points_per_view = 256L,
    nested_192_exact_prefix_verified = TRUE,
    nested_256_subset_of_384_verified = TRUE,
    selected_identity_verified = TRUE,
    maximum_mst_absolute_error = maximum_mst_error,
    maximum_norm_metric_absolute_error = maximum_norm_error,
    all_view_h0_mst_passed = TRUE,
    production_scientific_helpers_called = FALSE,
    labels_opened = FALSE, rankings_computed = FALSE,
    outcomes_computed = FALSE, stringsAsFactors = FALSE)
  if (index %% 10L == 0L) message("MV5-AJ independent geometry groups: ", index)
}
groups <- do.call(rbind, group_rows)
energies <- do.call(rbind, energy_rows)
h1_oracles <- do.call(rbind, h1_rows)
if (nrow(groups) != 150L || total_views != 13500L ||
    total_h0_edges != 13500L * 255L || nrow(energies) != 30L ||
    nrow(h1_oracles) != 60L ||
    length(unique(energies$fold_id)) != 3L ||
    length(unique(energies$seed)) != 5L ||
    length(unique(energies$representation)) != 2L) {
  stop("MV5-AJ independent oracle totals failed.", call. = FALSE)
}
summary <- data.frame(
  contract_id = "mv05aj_independent_geometry_summary_v1",
  groups = nrow(groups), nested_views = total_views,
  h0_mst_edges = total_h0_edges,
  maximum_mst_absolute_error = max(groups$maximum_mst_absolute_error),
  energy_oracles = nrow(energies), energy_fold_strata = 3L,
  energy_representations = 2L, energy_seeds = 5L,
  maximum_energy_absolute_error = max(energies$absolute_error),
  h1_diagram_oracles = nrow(h1_oracles),
  h1_fold_strata = length(unique(h1_oracles$fold_id)),
  nested_192_exact_prefix_verified = TRUE,
  nested_256_subset_of_384_verified = TRUE,
  production_scientific_helpers_called = FALSE,
  labels_opened = FALSE, rankings_computed = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE)
write_once(groups, "mv05aj-independent-geometry-groups-2026-08-11.csv")
write_once(energies, "mv05aj-independent-energy-oracles-2026-08-11.csv")
write_once(h1_oracles, "mv05aj-independent-h1-diagram-oracles-2026-08-11.csv")
write_once(summary, "mv05aj-independent-geometry-summary-2026-08-11.csv")
message("MV5-AJ independent geometry passed 13,500 views and 30 energy oracles.")
