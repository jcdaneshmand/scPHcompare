#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "cluster", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: run_mv08g_comparison.R PREFREEZE MV07H_LANDSCAPE_ROOT MV08G_LANDSCAPE_ROOT MATCHED_SHIFT_ROOT OUTPUT")
}
prefreeze <- args[[1L]]; root500 <- args[[2L]]; root475 <- args[[3L]]
shift_root <- args[[4L]]; output <- args[[5L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G comparison output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv05n_clustering_gate.R")
source("R/mv05_benchmark_contract.R")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
queue475 <- read.csv(file.path(prefreeze, "mv08g-landscape-queue.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
shift_queue <- read.csv(file.path(prefreeze, "mv08g-matched-shift-queue.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
inventory500 <- read.csv(file.path(prefreeze,
  "mv08g-accepted500-landscape-inventory.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
if (nrow(queue475) != 20L || nrow(shift_queue) != 20L ||
    nrow(inventory500) != 20L) stop("MV8-G comparison queues are incomplete.")

component_id <- function(view, dimension) paste0(
  ifelse(view == "cell_topology_v1", "cell", "gene"), "_", dimension)
pair_key <- function(value) paste(value$first_sample_id,
                                  value$second_sample_id, sep = "\r")
to_matrix <- function(value) {
  ids <- sort(unique(c(value$first_sample_id, value$second_sample_id)),
              method = "radix")
  if (length(ids) != 124L || nrow(value) != 7626L ||
      anyDuplicated(pair_key(value)) || any(!is.finite(value$distance)) ||
      any(value$distance < 0)) stop("MV8-G distance group is malformed.")
  matrix <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  first <- match(value$first_sample_id, ids); second <- match(value$second_sample_id, ids)
  matrix[cbind(first, second)] <- value$distance
  matrix[cbind(second, first)] <- value$distance
  matrix
}
load500 <- function(seed, view, dimension) {
  row <- inventory500[inventory500$seed == seed &
    inventory500$view_id == view &
    inventory500$homology_dimension == dimension, , drop = FALSE]
  if (nrow(row) != 1L) stop("Accepted 500-gene group is not unique.")
  path <- file.path(root500, gsub(":", "_", row$group_id, fixed = TRUE),
                    "distances.csv")
  if (!file.exists(path) || sha(path) != row$distances_sha256) {
    stop("Accepted 500-gene distance identity drift.")
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
load475 <- function(seed, view, dimension) {
  row <- queue475[queue475$seed == seed & queue475$view_id == view &
    queue475$homology_dimension == dimension, , drop = FALSE]
  path <- file.path(root475, row$group_id, "distances.csv")
  status <- read.csv(file.path(root475, row$group_id, "status.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(row) != 1L || nrow(status) != 1L ||
      status$completion_state != "complete" ||
      sha(path) != status$distances_sha256) stop("475-gene distance identity drift.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
load_shift <- function(seed, view, dimension) {
  row <- shift_queue[shift_queue$seed == seed & shift_queue$view_id == view &
    shift_queue$homology_dimension == dimension, , drop = FALSE]
  path <- file.path(shift_root, row$group_id, "distances.csv")
  status <- read.csv(file.path(shift_root, row$group_id, "status.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(row) != 1L || nrow(status) != 1L ||
      status$completion_state != "complete" ||
      sha(path) != status$distances_sha256) stop("Matched-shift identity drift.")
  value <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(value) != 124L || anyDuplicated(value$sample_id) ||
      any(!is.finite(value$distance)) || any(value$distance < 0)) {
    stop("Matched-shift rows are malformed.")
  }
  value
}

matrices <- list(panel500 = setNames(vector("list", 4L), .mv08g_components),
                 panel475 = setNames(vector("list", 4L), .mv08g_components))
for (panel in names(matrices)) for (component in .mv08g_components) {
  matrices[[panel]][[component]] <- setNames(vector("list", 5L),
                                              as.character(.mv08g_seeds))
}
seed_metrics <- list(); neighbor_rows <- list(); normalized_shift_rows <- list()
for (seed in .mv08g_seeds) for (view in .mv08g_views) for (dimension in .mv08g_dimensions) {
  component <- component_id(view, dimension)
  value500 <- load500(seed, view, dimension)
  value475 <- load475(seed, view, dimension)
  key500 <- pair_key(value500); key475 <- pair_key(value475)
  matched <- match(key500, key475)
  if (anyNA(matched) || !identical(key500, key475[matched])) {
    stop("500- and 475-gene pair axes differ.")
  }
  value475 <- value475[matched, , drop = FALSE]
  matrix500 <- to_matrix(value500); matrix475 <- to_matrix(value475)
  matrices$panel500[[component]][[as.character(seed)]] <- matrix500
  matrices$panel475[[component]][[as.character(seed)]] <- matrix475
  stress <- mv08g_nonnegative_scale_stress_v1(value500$distance,
                                               value475$distance)
  neighbors <- mv08g_top_k_neighbor_overlap_v1(matrix500, matrix475, k = 10L)
  neighbors$contract_id <- "mv08g_top10_neighbor_overlap_v1"
  neighbors$seed <- seed; neighbors$component_id <- component
  neighbors$outcome_label_state <- "closed"
  neighbors$biological_outcomes_computed <- FALSE
  neighbor_rows[[length(neighbor_rows) + 1L]] <- neighbors
  shift <- load_shift(seed, view, dimension)
  denominator <- stats::median(value500$distance[value500$distance > 0])
  if (!is.finite(denominator) || denominator <= 0) {
    stop("MV8-G 500-gene distance scale is invalid.")
  }
  shift$contract_id <- "mv08g_normalized_matched_shift_v1"
  shift$component_id <- component
  shift$median_nonzero_panel500_between_sample_distance <- denominator
  shift$normalized_shift <- shift$distance / denominator
  normalized_shift_rows[[length(normalized_shift_rows) + 1L]] <- shift
  seed_metrics[[length(seed_metrics) + 1L]] <- data.frame(
    contract_id = "mv08g_seed_distance_comparison_v1", seed = seed,
    component_id = component,
    spearman = stats::cor(value500$distance, value475$distance,
                          method = "spearman"),
    kendall = stats::cor(value500$distance, value475$distance,
                         method = "kendall"),
    nonnegative_scale = stress$scale,
    normalized_stress = stress$normalized_stress,
    mean_top10_overlap = mean(neighbors$overlap),
    median_normalized_matched_shift = stats::median(shift$normalized_shift),
    p90_normalized_matched_shift = as.numeric(stats::quantile(
      shift$normalized_shift, 0.90, names = FALSE, type = 8)),
    maximum_normalized_matched_shift = max(shift$normalized_shift),
    pair_count = nrow(value500), sample_count = nrow(shift),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
}
seed_metrics <- do.call(rbind, seed_metrics)
neighbors <- do.call(rbind, neighbor_rows)
normalized_shifts <- do.call(rbind, normalized_shift_rows)

candidate_rows <- list()
for (panel in names(matrices)) for (component in .mv08g_components) {
  for (seed in .mv08g_seeds) for (k in 2:10) {
    fit <- mv05n_pam_partition_v1(
      matrices[[panel]][[component]][[as.character(seed)]], k)
    fit$contract_id <- "mv08g_candidate_pam_partition_v1"
    fit$panel_id <- panel; fit$component_id <- component; fit$seed <- seed
    candidate_rows[[length(candidate_rows) + 1L]] <- fit
  }
}
candidates <- do.call(rbind, candidate_rows)
selection_rows <- list(); stability_rows <- list()
for (panel in names(matrices)) for (component in .mv08g_components) {
  part <- candidates[candidates$panel_id == panel &
    candidates$component_id == component, , drop = FALSE]
  selected <- mv05_select_stable_k_v1(part[c("seed", "k", "sample_id", "cluster")])
  if (selected$status != "selected" || is.na(selected$selected_k)) {
    stop("MV8-G label-free k selection failed.")
  }
  selection_rows[[length(selection_rows) + 1L]] <- data.frame(
    contract_id = "mv08g_panel_selected_k_v1", panel_id = panel,
    component_id = component, selected_k = selected$selected_k,
    selection_rule = "five_seed_smallest_k_within_one_SE_of_maximum_mean_stability",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  summary <- selected$summary
  summary$contract_id <- "mv08g_pam_stability_summary_v1"
  summary$panel_id <- panel; summary$component_id <- component
  summary$selected_k <- selected$selected_k; summary$threshold <- selected$threshold
  stability_rows[[length(stability_rows) + 1L]] <- summary
}
selections <- do.call(rbind, selection_rows)
stability <- do.call(rbind, stability_rows)

pam_k_rows <- list()
for (component in .mv08g_components) for (seed in .mv08g_seeds) for (k in 2:10) {
  left <- candidates[candidates$panel_id == "panel500" &
    candidates$component_id == component & candidates$seed == seed &
    candidates$k == k, ]
  right <- candidates[candidates$panel_id == "panel475" &
    candidates$component_id == component & candidates$seed == seed &
    candidates$k == k, ]
  right <- right[match(left$sample_id, right$sample_id), ]
  pam_k_rows[[length(pam_k_rows) + 1L]] <- data.frame(
    contract_id = "mv08g_pam_k_panel_agreement_v1", component_id = component,
    seed = seed, k = k,
    adjusted_rand_index = mclust::adjustedRandIndex(left$cluster, right$cluster),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE)
}
pam_k <- do.call(rbind, pam_k_rows)

fixed_rows <- list(); fixed_partition_rows <- list(); k_agreement_rows <- list()
for (component in .mv08g_components) {
  k500 <- selections$selected_k[selections$panel_id == "panel500" &
    selections$component_id == component]
  k475 <- selections$selected_k[selections$panel_id == "panel475" &
    selections$component_id == component]
  k_agreement_rows[[length(k_agreement_rows) + 1L]] <- data.frame(
    contract_id = "mv08g_panel_selected_k_agreement_v1",
    component_id = component, panel500_selected_k = k500,
    panel475_selected_k = k475, exact_k_agreement = k500 == k475,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE)
  for (seed in .mv08g_seeds) for (algorithm in c("pam", "average")) {
    fits <- lapply(names(matrices), function(panel) {
      matrix <- matrices[[panel]][[component]][[as.character(seed)]]
      if (algorithm == "pam") mv05n_pam_partition_v1(matrix, k500) else
        mv05n_average_partition_v1(matrix, k500)
    })
    names(fits) <- names(matrices)
    right <- fits$panel475[match(fits$panel500$sample_id,
                                 fits$panel475$sample_id), ]
    fixed_rows[[length(fixed_rows) + 1L]] <- data.frame(
      contract_id = "mv08g_fixed500k_panel_agreement_v1",
      component_id = component, seed = seed, algorithm = algorithm,
      fixed_panel500_selected_k = k500,
      adjusted_rand_index = mclust::adjustedRandIndex(
        fits$panel500$cluster, right$cluster),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE)
    for (panel in names(fits)) {
      value <- fits[[panel]]
      value$contract_id <- "mv08g_fixed500k_partition_v1"
      value$panel_id <- panel; value$component_id <- component
      value$seed <- seed; value$algorithm <- algorithm
      fixed_partition_rows[[length(fixed_partition_rows) + 1L]] <- value
    }
  }
}
fixed <- do.call(rbind, fixed_rows)
fixed_partitions <- do.call(rbind, fixed_partition_rows)
k_agreement <- do.call(rbind, k_agreement_rows)

component_summary <- do.call(rbind, lapply(.mv08g_components, function(component) {
  metric <- seed_metrics[seed_metrics$component_id == component, ]
  pam <- fixed[fixed$component_id == component & fixed$algorithm == "pam", ]
  data.frame(
    contract_id = "mv08g_component_summary_v1", component_id = component,
    median_spearman = stats::median(metric$spearman),
    median_kendall = stats::median(metric$kendall),
    median_normalized_stress = stats::median(metric$normalized_stress),
    median_top10_overlap = stats::median(metric$mean_top10_overlap),
    median_normalized_matched_shift = stats::median(
      metric$median_normalized_matched_shift),
    median_fixed_k_pam_ari = stats::median(pam$adjusted_rand_index),
    panel_selected_k_agreement = k_agreement$exact_k_agreement[
      k_agreement$component_id == component],
    seeds = nrow(metric), outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE)
}))
classification <- mv08g_harmonization_class_v1(component_summary)
decision <- data.frame(
  contract_id = "mv08g_harmonization_decision_v1", decision = classification,
  complete_components = nrow(component_summary), complete_seeds = 5L,
  raw_read_reprocessing_recommendation = if (classification ==
    "material_panel_sensitivity") "recommend_before_external_topology_claim" else
    if (classification == "high_harmonization_stability")
      "optional_strengthening_analysis_owner_preference_recorded" else
        "owner_decision_after_component_specific_review",
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE, biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)

outputs <- list(
  "mv08g-seed-distance-comparison.csv" = seed_metrics,
  "mv08g-top10-neighbor-overlap.csv" = neighbors,
  "mv08g-normalized-matched-shifts.csv" = normalized_shifts,
  "mv08g-candidate-pam-partitions.csv" = candidates,
  "mv08g-pam-stability-summary.csv" = stability,
  "mv08g-panel-selected-k.csv" = selections,
  "mv08g-panel-selected-k-agreement.csv" = k_agreement,
  "mv08g-pam-k-panel-agreement.csv" = pam_k,
  "mv08g-fixed500k-partitions.csv" = fixed_partitions,
  "mv08g-fixed500k-panel-agreement.csv" = fixed,
  "mv08g-component-summary.csv" = component_summary,
  "mv08g-decision.csv" = decision)
paths <- vapply(names(outputs), function(name) {
  path <- file.path(output, name); write_provenance_csv(outputs[[name]], path); path
}, character(1L))
manifest <- data.frame(
  contract_id = "mv08g_comparison_artifact_manifest_v1", file = basename(paths),
  bytes = as.numeric(file.info(paths)$size), sha256 = vapply(paths, sha, character(1L)),
  contains_expression = FALSE, contains_cell_barcode = FALSE,
  contains_absolute_private_path = FALSE, contains_biological_label = FALSE,
  stringsAsFactors = FALSE)
write_provenance_csv(manifest, file.path(output, "mv08g-artifact-manifest.csv"))
message("MV8-G comparison complete: ", classification,
        "; raw reads remain unauthorized")
