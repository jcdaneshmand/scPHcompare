#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "cluster", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv07i_label_closed_entry.R MV7I_PREFREEZE ",
       "MV7H_PREFREEZE MV7H_COMPLETE_VALIDATION LANDSCAPE_ROOT OUTPUT_DIR",
       call. = FALSE)
}
source("R/mv07h_full_topology.R")
source("R/mv05_benchmark_contract.R")
source("R/mv05n_clustering_gate.R")
source("R/mv05q_clustering_artifacts.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
mv07i <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
mv07h <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
mv07h_validation <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
landscape_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
output_dir <- args[[5L]]
decision <- read_csv(file.path(mv07i, "mv07i-decision.csv"))
registry <- read_csv(file.path(mv07i, "mv07i-representation-registry.csv"))
clustering_contract <- read_csv(file.path(mv07i, "mv07i-clustering-contract.csv"))
queue <- read_csv(file.path(mv07h, "mv07h-landscape-queue.csv"))
axis <- read_csv(file.path(mv07h, "mv07h-sample-seed-axis.csv"))
inventory <- read_csv(file.path(
  mv07h_validation, "mv07h-landscape-complete-group-inventory.csv"))
sample_ids <- sort(unique(axis$sample_id), method = "radix")
seeds <- sort(unique(as.integer(axis$seed)))
if (decision$decision !=
      "authorize_label_closed_matrix_and_clustering_production_only" ||
    as.logical(decision$labels_authorized) ||
    as.logical(decision$outcomes_authorized) || nrow(registry) != 6L ||
    any(registry$outcome_label_state != "closed") ||
    nrow(clustering_contract) != 6L ||
    any(as.logical(clustering_contract$label_values_used)) ||
    length(sample_ids) != 124L || !identical(seeds, 20260805:20260809) ||
    nrow(queue) != 20L || nrow(inventory) != 20L ||
    !setequal(queue$group_id, inventory$group_id) ||
    any(!as.logical(inventory$row_contract_passed))) {
  stop("MV7-I label-closed production admission is stale.", call. = FALSE)
}
safe <- function(value) gsub(":", "_", value, fixed = TRUE)
if (dir.exists(output_dir)) {
  status_path <- file.path(output_dir, "status.csv")
  if (!file.exists(status_path)) stop("Existing MV7-I output is incomplete.")
  status <- read_csv(status_path)
  named_files <- c(
    matrix_bundle = "matrix-bundle.rds",
    pair_summary = "pair-seed-summary.csv",
    h1_summary = "h1-contribution-summary.csv",
    candidate_partitions = "candidate-pam-partitions.csv",
    stability = "stability-summary.csv",
    selected_partitions = "selected-partitions.csv",
    provenance = "provenance.csv")
  paths <- file.path(output_dir, named_files)
  expected_hashes <- unname(unlist(status[paste0(
    names(named_files), "_sha256")], use.names = FALSE))
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      !all(file.exists(paths)) ||
      !identical(tolower(unname(vapply(paths, .mv07h_sha256, character(1L)))),
                 tolower(expected_hashes))) {
    stop("Existing MV7-I output is stale.", call. = FALSE)
  }
  message("Reused complete MV7-I label-closed artifacts.")
  quit(save = "no", status = 0L)
}
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = "mv07i__partial__", tmpdir = dirname(output_dir))
dir.create(partial)
started <- proc.time()[["elapsed"]]

matrix_from_rows <- function(rows) {
  expected <- utils::combn(sample_ids, 2L)
  if (nrow(rows) != 7626L ||
      !identical(as.character(rows$first_sample_id),
                 as.character(expected[1L, ])) ||
      !identical(as.character(rows$second_sample_id),
                 as.character(expected[2L, ])) ||
      any(!is.finite(rows$distance)) || any(rows$distance < 0)) {
    stop("MV7-I distance rows do not form the canonical matrix.")
  }
  result <- matrix(0, 124L, 124L, dimnames = list(sample_ids, sample_ids))
  index <- cbind(match(rows$first_sample_id, sample_ids),
                 match(rows$second_sample_id, sample_ids))
  result[index] <- rows$distance
  result[cbind(index[, 2L], index[, 1L])] <- rows$distance
  .mv05n_validate_distance_matrix(result)
}
components <- list()
source_rows <- list()
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  path <- file.path(landscape_root, safe(unit$group_id), "distances.csv")
  rows <- read_csv(path)
  expected <- inventory[inventory$group_id == unit$group_id, , drop = FALSE]
  observed_hash <- .mv07h_sha256(path)
  if (nrow(expected) != 1L ||
      !identical(tolower(observed_hash), tolower(expected$distances_sha256))) {
    stop("MV7-I source landscape hash drifted for ", unit$group_id, ".",
         call. = FALSE)
  }
  key <- paste(unit$view_id, unit$homology_dimension, unit$seed, sep = "__")
  components[[key]] <- matrix_from_rows(rows)
  source_rows[[index]] <- data.frame(
    group_id = unit$group_id, seed = unit$seed, view_id = unit$view_id,
    homology_dimension = unit$homology_dimension,
    distances_sha256 = observed_hash, stringsAsFactors = FALSE)
}
source_inventory <- do.call(rbind, source_rows)

matrix_sets <- list()
for (view in c("cell", "gene")) {
  view_id <- paste0(view, "_topology_v1")
  h0 <- stats::setNames(lapply(seeds, function(seed) components[[paste(
    view_id, "H0", seed, sep = "__")]]), as.character(seeds))
  h1 <- stats::setNames(lapply(seeds, function(seed) components[[paste(
    view_id, "H1", seed, sep = "__")]]), as.character(seeds))
  combined <- stats::setNames(Map(mv05q_combine_h0_h1_v1, h0, h1),
                              as.character(seeds))
  matrix_sets[[paste0(view, "_H0")]] <- h0
  matrix_sets[[paste0(view, "_H1")]] <- h1
  matrix_sets[[paste0(view, "_H0_H1_secondary")]] <- combined
}
if (!identical(names(matrix_sets), registry$representation_id)) {
  stop("MV7-I representation order drifted.", call. = FALSE)
}
upper <- upper.tri(matrix_sets[[1L]][[1L]])
pairs <- utils::combn(sample_ids, 2L)
pair_summaries <- list()
for (index in seq_along(matrix_sets)) {
  representation_id <- names(matrix_sets)[[index]]
  values <- vapply(matrix_sets[[index]], function(value) value[upper],
                   numeric(7626L))
  pair_summaries[[index]] <- data.frame(
    contract_id = "mv07i_pair_seed_summary_v1",
    representation_id = representation_id,
    first_sample_id = pairs[1L, ], second_sample_id = pairs[2L, ],
    seed_20260805 = values[, 1L], seed_20260806 = values[, 2L],
    seed_20260807 = values[, 3L], seed_20260808 = values[, 4L],
    seed_20260809 = values[, 5L],
    median = apply(values, 1L, stats::median),
    minimum = apply(values, 1L, min), maximum = apply(values, 1L, max),
    iqr = apply(values, 1L, stats::IQR),
    raw_mad = apply(values, 1L, stats::mad, constant = 1),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
pair_summary <- do.call(rbind, pair_summaries)
rownames(pair_summary) <- NULL

h1_summaries <- list()
for (view in c("cell", "gene")) {
  h0 <- matrix_sets[[paste0(view, "_H0")]]
  h1 <- matrix_sets[[paste0(view, "_H1")]]
  fractions <- vapply(seq_along(seeds), function(index) {
    h0_sq <- h0[[index]][upper] ^ 2
    h1_sq <- h1[[index]][upper] ^ 2
    denominator <- h0_sq + h1_sq
    ifelse(denominator == 0, 0, h1_sq / denominator)
  }, numeric(7626L))
  h1_summaries[[view]] <- data.frame(
    contract_id = "mv07i_h1_contribution_summary_v1", view_id = view,
    first_sample_id = pairs[1L, ], second_sample_id = pairs[2L, ],
    seed_20260805 = fractions[, 1L], seed_20260806 = fractions[, 2L],
    seed_20260807 = fractions[, 3L], seed_20260808 = fractions[, 4L],
    seed_20260809 = fractions[, 5L],
    median = apply(fractions, 1L, stats::median),
    minimum = apply(fractions, 1L, min), maximum = apply(fractions, 1L, max),
    iqr = apply(fractions, 1L, stats::IQR),
    raw_mad = apply(fractions, 1L, stats::mad, constant = 1),
    zero_rule = "zero_only_when_H0_and_H1_are_both_zero",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
h1_summary <- do.call(rbind, h1_summaries)
rownames(h1_summary) <- NULL

candidate_rows <- list(); stability_rows <- list(); selected_rows <- list()
for (index in seq_along(matrix_sets)) {
  representation_id <- names(matrix_sets)[[index]]
  candidates <- mv05n_fit_five_seed_partitions_v1(
    matrix_sets[[index]], method = "pam")
  selection <- mv05_select_stable_k_v1(candidates)
  if (selection$status != "selected" || is.na(selection$selected_k)) {
    stop("MV7-I returned no_stable_k for ", representation_id,
         "; preserve artifacts and prefreeze a reporting rebind.",
         call. = FALSE)
  }
  candidates$contract_id <- "mv07i_candidate_pam_partition_v1"
  candidates$representation_id <- representation_id
  candidates$outcome_label_state <- "closed"
  candidates$biological_outcomes_computed <- FALSE
  candidate_rows[[index]] <- candidates
  stability <- selection$summary
  names(stability)[names(stability) == "monte_carlo_se"] <- "jackknife_se"
  stability$contract_id <- "mv07i_stability_summary_v1"
  stability$representation_id <- representation_id
  stability$selected_k <- selection$selected_k
  stability$one_se_threshold <- selection$threshold
  stability$outcome_label_state <- "closed"
  stability$biological_outcomes_computed <- FALSE
  stability_rows[[index]] <- stability
  partitions <- list(); cursor <- 0L
  for (seed in seeds) for (algorithm in c("pam_stability_k_v1",
                                          "hclust_average_v1")) {
    fit <- if (algorithm == "pam_stability_k_v1") {
      mv05n_pam_partition_v1(
        matrix_sets[[index]][[as.character(seed)]], selection$selected_k)
    } else {
      mv05n_average_partition_v1(
        matrix_sets[[index]][[as.character(seed)]], selection$selected_k)
    }
    cursor <- cursor + 1L
    partitions[[cursor]] <- data.frame(
      contract_id = "mv07i_selected_partition_v1",
      representation_id = representation_id, seed = seed,
      algorithm_id = algorithm, selected_k = selection$selected_k,
      sample_id = fit$sample_id, cluster = fit$cluster,
      is_medoid = if ("is_medoid" %in% names(fit)) fit$is_medoid else FALSE,
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE)
  }
  selected_rows[[index]] <- do.call(rbind, partitions)
}
candidate_partitions <- do.call(rbind, candidate_rows)
stability <- do.call(rbind, stability_rows)
selected_partitions <- do.call(rbind, selected_rows)
for (value in list(candidate_partitions, stability, selected_partitions)) {
  if (any(value$outcome_label_state != "closed") ||
      any(as.logical(value$biological_outcomes_computed))) {
    stop("MV7-I label firewall opened.", call. = FALSE)
  }
}
median_matrices <- lapply(matrix_sets, function(values) {
  result <- apply(simplify2array(values), c(1L, 2L), stats::median)
  dimnames(result) <- dimnames(values[[1L]])
  .mv05n_validate_distance_matrix(result)
})
matrix_bundle <- list(
  contract_id = "mv07i_matrix_bundle_v1", sample_ids = sample_ids,
  seeds = seeds, seed_specific = matrix_sets, median = median_matrices,
  source_inventory = source_inventory,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE)
provenance <- data.frame(
  contract_id = "mv07i_label_closed_provenance_v1",
  source_groups = nrow(source_inventory), source_component_rows = 152520L,
  representations = length(matrix_sets), seed_specific_matrices = 30L,
  candidate_pam_fits = 270L, selected_average_linkage_fits = 30L,
  pair_summary_rows = nrow(pair_summary),
  h1_contribution_rows = nrow(h1_summary),
  candidate_partition_rows = nrow(candidate_partitions),
  stability_rows = nrow(stability),
  selected_partition_rows = nrow(selected_partitions),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  labels_loaded = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE)

paths <- c(
  matrix_bundle = file.path(partial, "matrix-bundle.rds"),
  pair_summary = file.path(partial, "pair-seed-summary.csv"),
  h1_summary = file.path(partial, "h1-contribution-summary.csv"),
  candidate_partitions = file.path(partial, "candidate-pam-partitions.csv"),
  stability = file.path(partial, "stability-summary.csv"),
  selected_partitions = file.path(partial, "selected-partitions.csv"),
  provenance = file.path(partial, "provenance.csv"))
saveRDS(matrix_bundle, paths[["matrix_bundle"]], version = 3L)
write.csv(pair_summary, paths[["pair_summary"]], row.names = FALSE, na = "")
write.csv(h1_summary, paths[["h1_summary"]], row.names = FALSE, na = "")
write.csv(candidate_partitions, paths[["candidate_partitions"]],
          row.names = FALSE, na = "")
write.csv(stability, paths[["stability"]], row.names = FALSE, na = "")
write.csv(selected_partitions, paths[["selected_partitions"]],
          row.names = FALSE, na = "")
write.csv(provenance, paths[["provenance"]], row.names = FALSE, na = "")
elapsed <- proc.time()[["elapsed"]] - started
status <- data.frame(
  contract_id = "mv07i_label_closed_status_v1", completion_state = "complete",
  elapsed_seconds = elapsed,
  matrix_bundle_sha256 = .mv07h_sha256(paths[["matrix_bundle"]]),
  pair_summary_sha256 = .mv07h_sha256(paths[["pair_summary"]]),
  h1_summary_sha256 = .mv07h_sha256(paths[["h1_summary"]]),
  candidate_partitions_sha256 = .mv07h_sha256(paths[["candidate_partitions"]]),
  stability_sha256 = .mv07h_sha256(paths[["stability"]]),
  selected_partitions_sha256 = .mv07h_sha256(paths[["selected_partitions"]]),
  provenance_sha256 = .mv07h_sha256(paths[["provenance"]]),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  labels_loaded = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE)
write.csv(status, file.path(partial, "status.csv"), row.names = FALSE, na = "")
if (!file.rename(partial, output_dir)) {
  unlink(partial, recursive = TRUE)
  stop("Failed to atomically publish MV7-I label-closed artifacts.")
}
message("Completed MV7-I label-closed matrix and clustering artifacts.")
