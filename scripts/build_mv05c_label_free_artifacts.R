#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv05c_label_free_artifacts.R AUDIT_DIR PRIVATE_DIR",
       call. = FALSE)
}
audit_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
private_dir <- args[[2L]]
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)

source("R/provenance_utils.R")
source("R/topological_distance_contract.R")

landscape <- utils::read.csv(
  file.path(audit_dir, "mv05c-landscape-pairs-2026-08-06.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
baseline <- utils::read.csv(
  file.path(audit_dir, "mv05c-baseline-pairs-2026-08-06.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
samples <- utils::read.csv(
  file.path(audit_dir, "mv05c-label-closed-sample-manifest-2026-08-06.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
resources <- utils::read.csv(
  file.path(audit_dir, "mv05c-job-resources-2026-08-06.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
landscape_environment <- utils::read.csv(
  file.path(audit_dir, "mv05c-landscape-environment-2026-08-06.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

if (any(c("tissue", "approach") %in% names(samples)) ||
    any(samples$outcome_label_state != "closed") ||
    any(as.logical(samples$biological_outcomes_computed)) ||
    nrow(landscape) != 750L || any(landscape$status != "completed") ||
    any(!is.finite(landscape$distance)) ||
    any(tolower(landscape$exact) != "true") ||
    any(tolower(landscape$all_active_levels) != "true")) {
  stop("Label or exact-distance preflight failed.")
}

matrix_from_pairs <- function(rows, ids) {
  result <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (index in seq_len(nrow(rows))) {
    left <- match(rows$first_sample_id[[index]], ids)
    right <- match(rows$second_sample_id[[index]], ids)
    if (is.na(left) || is.na(right) || left == right) {
      stop("Pair rows contain invalid sample identities.")
    }
    result[left, right] <- result[right, left] <- rows$distance[[index]]
  }
  if (nrow(rows) != choose(length(ids), 2L) ||
      any(!is.finite(result)) || !identical(result, t(result)) ||
      any(diag(result) != 0)) {
    stop("Distance matrix failed completeness or metric-shape checks.")
  }
  result
}

landscape_method <- function(view_id, dimension) {
  paste0(if (view_id == "cell_topology_v1") "cell" else "gene",
         "_landscape_", tolower(dimension), "_v1")
}

matrices <- list()
index_rows <- list()
add_matrix <- function(fold_id, seed, representation, method_id, matrix,
                       source_contract_id) {
  key <- paste(fold_id, seed, representation, method_id, sep = "__")
  if (!is.null(matrices[[key]])) stop("Duplicate matrix identity: ", key)
  sha <- digest::digest(matrix, algo = "sha256", serialize = TRUE)
  matrices[[key]] <<- matrix
  index_rows[[length(index_rows) + 1L]] <<- data.frame(
    matrix_id = key, fold_id = fold_id, seed = as.integer(seed),
    representation = representation, method_id = method_id,
    sample_count = nrow(matrix), pair_count = choose(nrow(matrix), 2L),
    minimum_distance = min(matrix[upper.tri(matrix)]),
    maximum_distance = max(matrix[upper.tri(matrix)]),
    distance_matrix_sha256 = sha, source_contract_id = source_contract_id,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

sample_ids <- sort(samples$sample_id, method = "radix")
landscape_groups <- split(
  seq_len(nrow(landscape)),
  interaction(
    landscape$fold_id, landscape$seed, landscape$representation,
    landscape$view_id, landscape$homology_dimension,
    drop = TRUE, lex.order = TRUE
  )
)
for (indices in landscape_groups) {
  rows <- landscape[indices, , drop = FALSE]
  add_matrix(
    rows$fold_id[[1L]], rows$seed[[1L]], rows$representation[[1L]],
    landscape_method(rows$view_id[[1L]], rows$homology_dimension[[1L]]),
    matrix_from_pairs(rows, sample_ids), rows$contract_id[[1L]]
  )
}

baseline_groups <- split(
  seq_len(nrow(baseline)),
  interaction(
    baseline$fold_id, baseline$seed, baseline$representation,
    baseline$method_id, drop = TRUE, lex.order = TRUE
  )
)
for (indices in baseline_groups) {
  rows <- baseline[indices, , drop = FALSE]
  add_matrix(
    rows$fold_id[[1L]], rows$seed[[1L]], rows$representation[[1L]],
    rows$method_id[[1L]], matrix_from_pairs(rows, sample_ids),
    rows$contract_id[[1L]]
  )
}
matrix_index <- do.call(rbind, index_rows)
if (nrow(matrix_index) != 85L) stop("Expected exactly 85 completed matrices.")

neighbor_rows <- list()
for (index in seq_len(nrow(matrix_index))) {
  info <- matrix_index[index, , drop = FALSE]
  matrix <- matrices[[info$matrix_id]]
  held_out_study <- sub("^mv05c_loso_v1:", "", info$fold_id)
  query_ids <- samples$sample_id[samples$study == held_out_study]
  training_ids <- samples$sample_id[samples$study != held_out_study]
  for (query_id in sort(query_ids, method = "radix")) {
    values <- as.numeric(matrix[query_id, training_ids, drop = FALSE])
    names(values) <- training_ids
    ordered <- order(values, names(values), method = "radix")
    for (rank in seq_along(ordered)) {
      neighbor_rows[[length(neighbor_rows) + 1L]] <- data.frame(
        matrix_id = info$matrix_id, fold_id = info$fold_id,
        seed = info$seed, representation = info$representation,
        method_id = info$method_id, query_sample_id = query_id,
        neighbor_rank = rank,
        training_sample_id = names(values)[ordered[[rank]]],
        distance = unname(values[ordered[[rank]]]),
        prediction_payload = "ranked_training_sample_ids_only_no_labels",
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
}
neighbors <- do.call(rbind, neighbor_rows)
if (nrow(neighbors) != 425L) stop("Held-out neighbor ranking is incomplete.")

adjusted_rand <- function(first, second) {
  table <- table(first, second)
  choose2 <- function(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(table))
  sum_rows <- sum(choose2(rowSums(table)))
  sum_cols <- sum(choose2(colSums(table)))
  total <- choose2(sum(table))
  expected <- if (total == 0) 0 else sum_rows * sum_cols / total
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) return(1)
  (sum_cells - expected) / (maximum - expected)
}

method_groups <- split(
  seq_len(nrow(matrix_index)),
  interaction(
    matrix_index$fold_id, matrix_index$representation,
    matrix_index$method_id, drop = TRUE, lex.order = TRUE
  )
)
partition_rows <- list()
stability_rows <- list()
summary_rows <- list()
selected_rows <- list()
for (indices in method_groups) {
  info <- matrix_index[indices, , drop = FALSE]
  if (!identical(sort(info$seed), 20260805:20260809)) {
    stop("Label-free stability requires all five frozen seeds.")
  }
  assignments <- list()
  for (k in 2:5) {
    assignments[[as.character(k)]] <- list()
    for (row in seq_len(nrow(info))) {
      matrix <- matrices[[info$matrix_id[[row]]]]
      cluster <- cluster::pam(
        stats::as.dist(matrix), k = k, diss = TRUE,
        cluster.only = TRUE, do.swap = TRUE
      )
      cluster <- cluster[order(names(cluster), method = "radix")]
      assignments[[as.character(k)]][[as.character(info$seed[[row]])]] <- cluster
      partition_rows[[length(partition_rows) + 1L]] <- data.frame(
        fold_id = info$fold_id[[1L]],
        representation = info$representation[[1L]],
        method_id = info$method_id[[1L]], seed = info$seed[[row]], k = k,
        sample_id = names(cluster), pam_cluster = as.integer(cluster),
        outcome_label_state = "closed", biological_outcomes_computed = FALSE,
        stringsAsFactors = FALSE
      )
    }
    seeds <- sort(as.integer(names(assignments[[as.character(k)]])))
    values <- numeric()
    for (left in seq_len(length(seeds) - 1L)) {
      for (right in seq.int(left + 1L, length(seeds))) {
        ari <- adjusted_rand(
          assignments[[as.character(k)]][[as.character(seeds[[left]])]],
          assignments[[as.character(k)]][[as.character(seeds[[right]])]]
        )
        values <- c(values, ari)
        stability_rows[[length(stability_rows) + 1L]] <- data.frame(
          fold_id = info$fold_id[[1L]],
          representation = info$representation[[1L]],
          method_id = info$method_id[[1L]], k = k,
          first_seed = seeds[[left]], second_seed = seeds[[right]],
          adjusted_rand = ari, outcome_label_state = "closed",
          biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
        )
      }
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      fold_id = info$fold_id[[1L]],
      representation = info$representation[[1L]],
      method_id = info$method_id[[1L]], k = k,
      mean_pairwise_ari = mean(values),
      pairwise_ari_se = stats::sd(values) / sqrt(length(values)),
      seed_pairs = length(values), outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }
}
partitions <- do.call(rbind, partition_rows)
stability <- do.call(rbind, stability_rows)
stability_summary <- do.call(rbind, summary_rows)

selection_groups <- split(
  seq_len(nrow(stability_summary)),
  interaction(
    stability_summary$fold_id, stability_summary$representation,
    stability_summary$method_id, drop = TRUE, lex.order = TRUE
  )
)
for (indices in selection_groups) {
  rows <- stability_summary[indices, , drop = FALSE]
  best <- which.max(rows$mean_pairwise_ari)
  threshold <- rows$mean_pairwise_ari[[best]] - rows$pairwise_ari_se[[best]]
  selected_k <- min(rows$k[rows$mean_pairwise_ari >= threshold])
  selected_rows[[length(selected_rows) + 1L]] <- data.frame(
    fold_id = rows$fold_id[[1L]], representation = rows$representation[[1L]],
    method_id = rows$method_id[[1L]], selected_k = selected_k,
    one_se_threshold = threshold, best_mean_k = rows$k[[best]],
    best_mean_pairwise_ari = rows$mean_pairwise_ari[[best]],
    selection_rule = "smallest_k_within_one_se_of_max_mean_pairwise_ari",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
selections <- do.call(rbind, selected_rows)

hclust_rows <- list()
for (index in seq_len(nrow(selections))) {
  choice <- selections[index, , drop = FALSE]
  matches <- matrix_index$fold_id == choice$fold_id &
    matrix_index$representation == choice$representation &
    matrix_index$method_id == choice$method_id
  info <- matrix_index[matches, , drop = FALSE]
  for (row in seq_len(nrow(info))) {
    matrix <- matrices[[info$matrix_id[[row]]]]
    cluster <- stats::cutree(
      stats::hclust(stats::as.dist(matrix), method = "average"),
      k = choice$selected_k
    )
    hclust_rows[[length(hclust_rows) + 1L]] <- data.frame(
      fold_id = choice$fold_id, representation = choice$representation,
      method_id = choice$method_id, seed = info$seed[[row]],
      selected_k = choice$selected_k, sample_id = names(cluster),
      average_linkage_cluster = as.integer(cluster),
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
hclust_partitions <- do.call(rbind, hclust_rows)

private_bundle <- list(
  contract_id = "mv05c_label_free_distance_bundle_v1",
  matrices = matrices, matrix_index = matrix_index,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE
)
private_path <- file.path(private_dir, "mv05c-label-free-distance-bundle.rds")
if (file.exists(private_path)) stop("Refusing to overwrite private distance bundle.")
saveRDS(private_bundle, private_path, compress = FALSE, version = 3)
private_sha <- digest::digest(
  file = private_path, algo = "sha256", serialize = FALSE
)
matrix_index$private_bundle_file <- basename(private_path)
matrix_index$private_bundle_sha256 <- private_sha

observed_worker_seconds <- sum(resources$elapsed_seconds)
observed_job_median <- stats::median(resources$elapsed_seconds)
observed_job_max <- max(resources$elapsed_seconds)
observed_peak_bytes <- max(resources$maximum_resident_kbytes) * 1024
full_samples <- 90L
pilot_samples <- 6L
full_folds <- 15L
seeds <- 5L
full_jobs <- full_folds * seeds
sample_linear_factor <- full_samples / pilot_samples
full_job_seconds_linear <- observed_job_median * sample_linear_factor
full_worker_hours_linear <- full_job_seconds_linear * full_jobs / 3600
full_pairs_per_matrix <- choose(full_samples, 2L)
pilot_pairs_per_matrix <- choose(pilot_samples, 2L)
full_landscape_rows <- full_folds * seeds * 6L * full_pairs_per_matrix
landscape_seconds_per_pair <- landscape_environment$elapsed_seconds[[1L]] /
  nrow(landscape)
full_landscape_hours_linear <- full_landscape_rows *
  landscape_seconds_per_pair / 3600
projection <- data.frame(
  contract_id = "mv05c_to_mv05d_resource_projection_v1",
  observed_jobs = nrow(resources),
  observed_worker_seconds = observed_worker_seconds,
  observed_job_median_seconds = observed_job_median,
  observed_job_max_seconds = observed_job_max,
  observed_peak_process_tree_rss_bytes = observed_peak_bytes,
  observed_samples = pilot_samples, projected_samples = full_samples,
  projected_folds = full_folds, projected_seeds = seeds,
  projected_fold_seed_jobs = full_jobs,
  sample_linear_factor = sample_linear_factor,
  projected_job_seconds_sample_linear = full_job_seconds_linear,
  projected_worker_hours_sample_linear = full_worker_hours_linear,
  projected_landscape_pair_rows = full_landscape_rows,
  projected_landscape_hours_pair_linear = full_landscape_hours_linear,
  current_per_job_30_minute_cap_passes_projection =
    full_job_seconds_linear <= 1800,
  current_24_worker_hour_cap_passes_projection =
    full_worker_hours_linear + full_landscape_hours_linear <= 24,
  disposition = "full_mv05d_requires_cached_normalization_chunking_and_pair_scope_reduction",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

write_public <- function(value, filename) utils::write.csv(
  value, file.path(audit_dir, filename), row.names = FALSE, na = ""
)
write_public(matrix_index, "mv05c-distance-matrix-index-2026-08-06.csv")
write_public(neighbors, "mv05c-label-free-neighbor-rankings-2026-08-06.csv")
write_public(partitions, "mv05c-pam-partitions-2026-08-06.csv")
write_public(stability, "mv05c-pam-stability-detail-2026-08-06.csv")
write_public(stability_summary, "mv05c-pam-stability-summary-2026-08-06.csv")
write_public(selections, "mv05c-pam-one-se-selection-2026-08-06.csv")
write_public(hclust_partitions, "mv05c-average-linkage-partitions-2026-08-06.csv")
write_public(projection, "mv05c-to-mv05d-resource-projection-2026-08-06.csv")

stopifnot(
  nrow(matrix_index) == 85L, nrow(neighbors) == 425L,
  nrow(stability_summary) == 68L, nrow(selections) == 17L,
  all(!matrix_index$biological_outcomes_computed),
  all(!neighbors$biological_outcomes_computed),
  all(!partitions$biological_outcomes_computed),
  all(!hclust_partitions$biological_outcomes_computed),
  !projection$current_per_job_30_minute_cap_passes_projection,
  !projection$current_24_worker_hour_cap_passes_projection
)
message("Built MV5-C label-free matrices, rankings, clustering, and MV5-D projection.")
