#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "cluster", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop("usage: validate_mv08g_comparison.R PREFREEZE VALIDATION_PREFREEZE PRIMARY_PREFREEZE MV07H_LANDSCAPE_ROOT MV08G_LANDSCAPE_ROOT MATCHED_SHIFT_ROOT RESULT OUTPUT")
}
prefreeze <- args[[1L]]; validation_prefreeze <- args[[2L]]
primary <- args[[3L]]; root500 <- args[[4L]]; root475 <- args[[5L]]
shift_root <- args[[6L]]; result <- args[[7L]]; output <- args[[8L]]
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G comparison validation output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv05n_clustering_gate.R")
source("R/mv05_benchmark_contract.R")
source("R/mv08g_panel_sensitivity.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
validation_contract <- read.csv(file.path(validation_prefreeze,
  "mv08g-comparison-validation-contract.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
validation_decision <- read.csv(file.path(validation_prefreeze,
  "mv08g-comparison-validation-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
validation_freeze <- read.csv(file.path(validation_prefreeze,
  "mv08g-comparison-validation-source-freeze.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
execution_contract <- read.csv(file.path(prefreeze, "mv08g-comparison-contract.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (nrow(validation_contract) != 1L || nrow(validation_decision) != 1L ||
    validation_contract$validator_head != head ||
    validation_contract$execution_head != execution_contract$accepted_head ||
    validation_decision$decision !=
      "authorize_one_independent_comparison_validation_after_representation_closure" ||
    validation_decision$validation_jobs_authorized != 1L ||
    validation_decision$comparison_jobs_authorized != 0L ||
    any(!file.exists(validation_freeze$artifact_locator)) ||
    any(vapply(validation_freeze$artifact_locator, sha, character(1L)) !=
          validation_freeze$sha256)) {
  stop("MV8-G comparison validation-only prefreeze is stale.")
}
same_numeric <- function(first, second, tolerance = 1e-12) {
  length(first) == length(second) && all(is.finite(first)) && all(is.finite(second)) &&
    all(abs(as.numeric(first) - as.numeric(second)) <= tolerance *
      pmax(1, abs(as.numeric(first)), abs(as.numeric(second))))
}
read_result <- function(name) read.csv(file.path(result, name),
  stringsAsFactors = FALSE, check.names = FALSE)
manifest <- read_result("mv08g-artifact-manifest.csv")
if (nrow(manifest) != 12L ||
    sha(file.path(result, "mv08g-artifact-manifest.csv")) !=
      validation_contract$result_manifest_sha256 ||
    any(vapply(file.path(result, manifest$file),
  sha, character(1L)) != manifest$sha256) || any(truth(manifest$contains_expression)) ||
  any(truth(manifest$contains_cell_barcode)) ||
  any(truth(manifest$contains_absolute_private_path)) ||
  any(truth(manifest$contains_biological_label))) {
  stop("MV8-G comparison artifact manifest failed.")
}
seed_metrics <- read_result("mv08g-seed-distance-comparison.csv")
neighbors <- read_result("mv08g-top10-neighbor-overlap.csv")
shifts <- read_result("mv08g-normalized-matched-shifts.csv")
candidates <- read_result("mv08g-candidate-pam-partitions.csv")
stability <- read_result("mv08g-pam-stability-summary.csv")
selections <- read_result("mv08g-panel-selected-k.csv")
k_agreement <- read_result("mv08g-panel-selected-k-agreement.csv")
pam_agreement <- read_result("mv08g-pam-k-panel-agreement.csv")
fixed_partitions <- read_result("mv08g-fixed500k-partitions.csv")
fixed_agreement <- read_result("mv08g-fixed500k-panel-agreement.csv")
component_summary <- read_result("mv08g-component-summary.csv")
result_decision <- read_result("mv08g-decision.csv")
expected_rows <- c(seed_metrics = 20L, neighbors = 2480L, shifts = 2480L,
  candidates = 44640L, stability = 72L, selections = 8L,
  k_agreement = 4L, pam_agreement = 180L, fixed_partitions = 9920L,
  fixed_agreement = 40L, component_summary = 4L, result_decision = 1L)
actual_rows <- c(nrow(seed_metrics), nrow(neighbors), nrow(shifts),
  nrow(candidates), nrow(stability), nrow(selections), nrow(k_agreement),
  nrow(pam_agreement), nrow(fixed_partitions), nrow(fixed_agreement),
  nrow(component_summary), nrow(result_decision))
if (!identical(unname(actual_rows), unname(expected_rows)) ||
    anyDuplicated(candidates[c("panel_id", "component_id", "seed", "k", "sample_id")]) ||
    anyDuplicated(fixed_partitions[c("panel_id", "component_id", "seed",
                                     "algorithm", "sample_id")]) ||
    truth(result_decision$hca_fastq_download_authorized) ||
    truth(result_decision$raw_reprocessing_authorized) ||
    truth(result_decision$label_access_authorized)) {
  stop("MV8-G comparison result axes or stop boundary failed.")
}
queue475 <- read.csv(file.path(primary, "mv08g-landscape-queue.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
shift_queue <- read.csv(file.path(primary, "mv08g-matched-shift-queue.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
inventory500 <- read.csv(file.path(primary,
  "mv08g-accepted500-landscape-inventory.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
component_id <- function(view, dimension) paste0(
  ifelse(view == "cell_topology_v1", "cell", "gene"), "_", dimension)
pair_key <- function(value) paste(value$first_sample_id,
                                  value$second_sample_id, sep = "\r")
to_matrix <- function(value) {
  ids <- sort(unique(c(value$first_sample_id, value$second_sample_id)),
              method = "radix")
  matrix <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  first <- match(value$first_sample_id, ids); second <- match(value$second_sample_id, ids)
  matrix[cbind(first, second)] <- value$distance
  matrix[cbind(second, first)] <- value$distance
  matrix
}
load500 <- function(seed, view, dimension) {
  row <- inventory500[inventory500$seed == seed & inventory500$view_id == view &
    inventory500$homology_dimension == dimension, , drop = FALSE]
  path <- file.path(root500, gsub(":", "_", row$group_id, fixed = TRUE),
                    "distances.csv")
  if (nrow(row) != 1L || sha(path) != row$distances_sha256) stop("500 distance drift.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
load475 <- function(seed, view, dimension) {
  row <- queue475[queue475$seed == seed & queue475$view_id == view &
    queue475$homology_dimension == dimension, , drop = FALSE]
  path <- file.path(root475, row$group_id, "distances.csv")
  status <- read.csv(file.path(root475, row$group_id, "status.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(row) != 1L || sha(path) != status$distances_sha256) stop("475 distance drift.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
load_shift <- function(seed, view, dimension) {
  row <- shift_queue[shift_queue$seed == seed & shift_queue$view_id == view &
    shift_queue$homology_dimension == dimension, , drop = FALSE]
  path <- file.path(shift_root, row$group_id, "distances.csv")
  status <- read.csv(file.path(shift_root, row$group_id, "status.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(row) != 1L || sha(path) != status$distances_sha256) stop("Shift drift.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
matrices <- list(panel500 = setNames(vector("list", 4L), .mv08g_components),
                 panel475 = setNames(vector("list", 4L), .mv08g_components))
for (panel in names(matrices)) for (component in .mv08g_components) {
  matrices[[panel]][[component]] <- setNames(vector("list", 5L),
                                              as.character(.mv08g_seeds))
}
reconstructed_seed <- list(); reconstructed_neighbors <- list()
reconstructed_shifts <- list()
for (seed in .mv08g_seeds) for (view in .mv08g_views) for (dimension in .mv08g_dimensions) {
  component <- component_id(view, dimension)
  d500 <- load500(seed, view, dimension); d475 <- load475(seed, view, dimension)
  matched <- match(pair_key(d500), pair_key(d475)); d475 <- d475[matched, ]
  m500 <- to_matrix(d500); m475 <- to_matrix(d475)
  matrices$panel500[[component]][[as.character(seed)]] <- m500
  matrices$panel475[[component]][[as.character(seed)]] <- m475
  stress <- mv08g_nonnegative_scale_stress_v1(d500$distance, d475$distance)
  overlap <- mv08g_top_k_neighbor_overlap_v1(m500, m475, 10L)
  overlap$seed <- seed; overlap$component_id <- component
  reconstructed_neighbors[[length(reconstructed_neighbors) + 1L]] <- overlap
  shift <- load_shift(seed, view, dimension)
  denominator <- stats::median(d500$distance[d500$distance > 0])
  shift$component_id <- component; shift$normalized_shift <- shift$distance / denominator
  reconstructed_shifts[[length(reconstructed_shifts) + 1L]] <- shift
  reconstructed_seed[[length(reconstructed_seed) + 1L]] <- data.frame(
    seed = seed, component_id = component,
    spearman = stats::cor(d500$distance, d475$distance, method = "spearman"),
    kendall = stats::cor(d500$distance, d475$distance, method = "kendall"),
    nonnegative_scale = stress$scale, normalized_stress = stress$normalized_stress,
    mean_top10_overlap = mean(overlap$overlap),
    median_normalized_matched_shift = stats::median(shift$normalized_shift),
    p90_normalized_matched_shift = as.numeric(stats::quantile(
      shift$normalized_shift, 0.90, names = FALSE, type = 8)),
    maximum_normalized_matched_shift = max(shift$normalized_shift))
}
reconstructed_seed <- do.call(rbind, reconstructed_seed)
order_seed <- order(seed_metrics$component_id, seed_metrics$seed)
matched_seed <- match(paste(seed_metrics$component_id, seed_metrics$seed),
  paste(reconstructed_seed$component_id, reconstructed_seed$seed))
metric_fields <- c("spearman", "kendall", "nonnegative_scale", "normalized_stress",
  "mean_top10_overlap", "median_normalized_matched_shift",
  "p90_normalized_matched_shift", "maximum_normalized_matched_shift")
if (anyNA(matched_seed) || !all(vapply(metric_fields, function(field)
  same_numeric(seed_metrics[[field]], reconstructed_seed[[field]][matched_seed]),
  logical(1L)))) stop("MV8-G independent seed-metric reconstruction failed.")
reconstructed_neighbors <- do.call(rbind, reconstructed_neighbors)
neighbor_match <- match(paste(neighbors$component_id, neighbors$seed,
  neighbors$sample_id), paste(reconstructed_neighbors$component_id,
  reconstructed_neighbors$seed, reconstructed_neighbors$sample_id))
if (anyNA(neighbor_match) || !same_numeric(neighbors$overlap,
  reconstructed_neighbors$overlap[neighbor_match])) stop("Neighbor reconstruction failed.")
reconstructed_shifts <- do.call(rbind, reconstructed_shifts)
shift_match <- match(paste(shifts$component_id, shifts$seed, shifts$sample_id),
  paste(reconstructed_shifts$component_id, reconstructed_shifts$seed,
        reconstructed_shifts$sample_id))
if (anyNA(shift_match) || !same_numeric(shifts$normalized_shift,
  reconstructed_shifts$normalized_shift[shift_match])) stop("Shift reconstruction failed.")

reconstructed_candidates <- list()
for (panel in names(matrices)) for (component in .mv08g_components) {
  for (seed in .mv08g_seeds) for (k in 2:10) {
    fit <- mv05n_pam_partition_v1(matrices[[panel]][[component]][[as.character(seed)]], k)
    fit$panel_id <- panel; fit$component_id <- component; fit$seed <- seed
    reconstructed_candidates[[length(reconstructed_candidates) + 1L]] <- fit
  }
}
reconstructed_candidates <- do.call(rbind, reconstructed_candidates)
candidate_match <- match(paste(candidates$panel_id, candidates$component_id,
  candidates$seed, candidates$k, candidates$sample_id),
  paste(reconstructed_candidates$panel_id, reconstructed_candidates$component_id,
        reconstructed_candidates$seed, reconstructed_candidates$k,
        reconstructed_candidates$sample_id))
if (anyNA(candidate_match) || !identical(as.integer(candidates$cluster),
  as.integer(reconstructed_candidates$cluster[candidate_match]))) {
  stop("MV8-G candidate PAM partitions did not independently reproduce.")
}
reconstructed_selection <- list(); reconstructed_stability <- list()
for (panel in names(matrices)) for (component in .mv08g_components) {
  part <- reconstructed_candidates[reconstructed_candidates$panel_id == panel &
    reconstructed_candidates$component_id == component, ]
  selected <- mv05_select_stable_k_v1(part[c("seed", "k", "sample_id", "cluster")])
  reconstructed_selection[[length(reconstructed_selection) + 1L]] <- data.frame(
    panel_id = panel, component_id = component, selected_k = selected$selected_k)
  value <- selected$summary
  value$panel_id <- panel; value$component_id <- component
  value$selected_k <- selected$selected_k; value$threshold <- selected$threshold
  reconstructed_stability[[length(reconstructed_stability) + 1L]] <- value
}
reconstructed_selection <- do.call(rbind, reconstructed_selection)
reconstructed_stability <- do.call(rbind, reconstructed_stability)
selection_match <- match(paste(selections$panel_id, selections$component_id),
  paste(reconstructed_selection$panel_id, reconstructed_selection$component_id))
if (anyNA(selection_match) || !identical(as.integer(selections$selected_k),
  as.integer(reconstructed_selection$selected_k[selection_match]))) {
  stop("MV8-G selected k did not independently reproduce.")
}
stability_match <- match(paste(stability$panel_id, stability$component_id,
  stability$k), paste(reconstructed_stability$panel_id,
  reconstructed_stability$component_id, reconstructed_stability$k))
for (field in c("mean_stability", "monte_carlo_se", "threshold")) {
  if (anyNA(stability_match) || !same_numeric(stability[[field]],
      reconstructed_stability[[field]][stability_match])) {
    stop("MV8-G stability summary did not reproduce: ", field)
  }
}
recomputed_k_agreement <- vapply(k_agreement$component_id, function(component) {
  first <- reconstructed_selection$selected_k[
    reconstructed_selection$panel_id == "panel500" &
      reconstructed_selection$component_id == component]
  second <- reconstructed_selection$selected_k[
    reconstructed_selection$panel_id == "panel475" &
      reconstructed_selection$component_id == component]
  first == second
}, logical(1L))
if (!identical(unname(truth(k_agreement$exact_k_agreement)),
               unname(recomputed_k_agreement))) {
  stop("MV8-G selected-k agreement did not reproduce.")
}
recomputed_pam_ari <- vapply(seq_len(nrow(pam_agreement)), function(index) {
  row <- pam_agreement[index, ]
  first <- reconstructed_candidates[
    reconstructed_candidates$panel_id == "panel500" &
      reconstructed_candidates$component_id == row$component_id &
      reconstructed_candidates$seed == row$seed &
      reconstructed_candidates$k == row$k, ]
  second <- reconstructed_candidates[
    reconstructed_candidates$panel_id == "panel475" &
      reconstructed_candidates$component_id == row$component_id &
      reconstructed_candidates$seed == row$seed &
      reconstructed_candidates$k == row$k, ]
  second <- second[match(first$sample_id, second$sample_id), ]
  mclust::adjustedRandIndex(first$cluster, second$cluster)
}, numeric(1L))
if (!same_numeric(pam_agreement$adjusted_rand_index, recomputed_pam_ari)) {
  stop("MV8-G complete k=2:10 PAM agreement did not reproduce.")
}
reconstructed_fixed <- list()
for (panel in names(matrices)) for (component in .mv08g_components) {
  k500 <- reconstructed_selection$selected_k[
    reconstructed_selection$panel_id == "panel500" &
      reconstructed_selection$component_id == component]
  for (seed in .mv08g_seeds) for (algorithm in c("pam", "average")) {
    matrix <- matrices[[panel]][[component]][[as.character(seed)]]
    fit <- if (algorithm == "pam") mv05n_pam_partition_v1(matrix, k500) else
      mv05n_average_partition_v1(matrix, k500)
    fit <- mv08g_fixed_partition_schema_v1(fit, algorithm)
    fit$panel_id <- panel; fit$component_id <- component; fit$seed <- seed
    fit$algorithm <- algorithm
    reconstructed_fixed[[length(reconstructed_fixed) + 1L]] <- fit
  }
}
reconstructed_fixed <- do.call(rbind, reconstructed_fixed)
fixed_match <- match(paste(fixed_partitions$panel_id, fixed_partitions$component_id,
  fixed_partitions$seed, fixed_partitions$algorithm, fixed_partitions$sample_id),
  paste(reconstructed_fixed$panel_id, reconstructed_fixed$component_id,
        reconstructed_fixed$seed, reconstructed_fixed$algorithm,
        reconstructed_fixed$sample_id))
if (anyNA(fixed_match) || !identical(as.integer(fixed_partitions$cluster),
  as.integer(reconstructed_fixed$cluster[fixed_match]))) {
  stop("MV8-G fixed-k partitions did not independently reproduce.")
}
recomputed_fixed_ari <- vapply(seq_len(nrow(fixed_agreement)), function(index) {
  row <- fixed_agreement[index, ]
  first <- reconstructed_fixed[reconstructed_fixed$panel_id == "panel500" &
    reconstructed_fixed$component_id == row$component_id &
    reconstructed_fixed$seed == row$seed &
    reconstructed_fixed$algorithm == row$algorithm, ]
  second <- reconstructed_fixed[reconstructed_fixed$panel_id == "panel475" &
    reconstructed_fixed$component_id == row$component_id &
    reconstructed_fixed$seed == row$seed &
    reconstructed_fixed$algorithm == row$algorithm, ]
  second <- second[match(first$sample_id, second$sample_id), ]
  mclust::adjustedRandIndex(first$cluster, second$cluster)
}, numeric(1L))
if (!same_numeric(fixed_agreement$adjusted_rand_index, recomputed_fixed_ari)) {
  stop("MV8-G fixed-k ARI reconstruction failed.")
}
recomputed_summary <- do.call(rbind, lapply(.mv08g_components, function(component) {
  metric <- reconstructed_seed[reconstructed_seed$component_id == component, ]
  pam <- fixed_agreement[fixed_agreement$component_id == component &
    fixed_agreement$algorithm == "pam", ]
  data.frame(component_id = component,
    median_spearman = stats::median(metric$spearman),
    median_top10_overlap = stats::median(metric$mean_top10_overlap),
    median_fixed_k_pam_ari = stats::median(pam$adjusted_rand_index))
}))
summary_match <- match(component_summary$component_id, recomputed_summary$component_id)
if (anyNA(summary_match) || !same_numeric(component_summary$median_spearman,
    recomputed_summary$median_spearman[summary_match]) ||
    !same_numeric(component_summary$median_top10_overlap,
      recomputed_summary$median_top10_overlap[summary_match]) ||
    !same_numeric(component_summary$median_fixed_k_pam_ari,
      recomputed_summary$median_fixed_k_pam_ari[summary_match])) {
  stop("MV8-G component summary reconstruction failed.")
}
classification <- mv08g_harmonization_class_v1(recomputed_summary)
if (result_decision$decision != classification) {
  stop("MV8-G harmonization classification did not independently reproduce.")
}
checks <- data.frame(
  contract_id = "mv08g_comparison_independent_validation_v1",
  check_id = c("artifact_identity", "row_axes", "distance_metrics",
    "neighbors", "matched_shifts", "candidate_PAM", "stability_and_selected_k",
    "pam_k_agreement", "fixed_partitions", "fixed_ARI", "component_summary", "classification",
    "label_and_raw_read_firewall"), passed = TRUE,
  detail = c("12 result files match manifest", "all frozen row counts exact",
    "20 seed-component metric rows reconstructed", "2,480 rows reconstructed",
    "2,480 normalized shifts reconstructed", "44,640 rows refit",
    "72 stability rows and eight selections reproduced", "180 ARIs reproduced",
    "9,920 rows refit",
    "40 ARIs reproduced", "four summaries reproduced", classification,
    "no HCA download, raw processing, labels, or outcomes authorized"),
  stringsAsFactors = FALSE)
decision <- data.frame(
  contract_id = "mv08g_comparison_validation_decision_v1",
  decision = "comparison_exact_ready_for_raw_read_owner_decision",
  harmonization_classification = classification,
  independent_checks_passed = nrow(checks),
  hca_fastq_download_authorized = FALSE, raw_reprocessing_authorized = FALSE,
  label_access_authorized = FALSE,
  next_gate = "owner_review_raw_read_reprocessing_decision",
  stringsAsFactors = FALSE)
write_provenance_csv(checks,
  file.path(output, "mv08g-comparison-independent-validation.csv"))
write_provenance_csv(decision,
  file.path(output, "mv08g-comparison-validation-decision.csv"))
for (name in manifest$file) {
  value <- read_result(name)
  write_provenance_csv(value, file.path(output, name))
}
published <- c("mv08g-comparison-independent-validation.csv",
               "mv08g-comparison-validation-decision.csv", manifest$file)
published_paths <- file.path(output, published)
out_manifest <- data.frame(
  contract_id = "mv08g_comparison_validation_artifact_manifest_v1",
  file = published, bytes = as.numeric(file.info(published_paths)$size),
  sha256 = vapply(published_paths, sha, character(1L)), contains_expression = FALSE,
  contains_cell_barcode = FALSE, contains_absolute_private_path = FALSE,
  contains_biological_label = FALSE, stringsAsFactors = FALSE)
write_provenance_csv(out_manifest,
  file.path(output, "mv08g-comparison-validation-artifact-manifest.csv"))
message("MV8-G comparison independently reproduced: ", classification)
