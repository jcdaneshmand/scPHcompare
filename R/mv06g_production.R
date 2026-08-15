# Internal MV6-G complete label-closed scale/ranking production helpers.

mv06g_fit_group_scales_v1 <- function(training_distances,
                                       expected_training_pairs) {
  expected_training_pairs <- as.integer(expected_training_pairs)
  required <- c("group_id", "component_id", "distance",
                "outcome_label_state", "biological_outcomes_computed")
  if (!is.data.frame(training_distances) ||
      !all(required %in% names(training_distances)) ||
      nrow(training_distances) != 4L * expected_training_pairs ||
      any(!is.finite(training_distances$distance)) ||
      any(training_distances$distance < 0) ||
      any(training_distances$outcome_label_state != "closed") ||
      any(as.logical(training_distances$biological_outcomes_computed))) {
    stop("MV6-G production training distances are invalid.", call. = FALSE)
  }
  components <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
  rows <- lapply(seq_along(components), function(index) {
    component <- components[[index]]
    values <- training_distances$distance[
      training_distances$component_id == component
    ]
    scale <- stats::median(values)
    if (length(values) != expected_training_pairs || !is.finite(scale) ||
        scale <= sqrt(.Machine$double.eps)) {
      stop("MV6-G production scale is degenerate: ", component,
           call. = FALSE)
    }
    data.frame(
      contract_id = "mv06g_training_component_scale_v1",
      group_id = unique(training_distances$group_id),
      component_order = index, component_id = component,
      training_pairs = length(values), scale_statistic = "median",
      scale_value = scale, minimum_distance = min(values),
      maximum_distance = max(values), query_rows_used = 0L,
      labels_used = 0L, outcome_label_state = "closed",
      biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (nrow(result) != 4L || anyDuplicated(result$component_id)) {
    stop("MV6-G production did not fit four scales.", call. = FALSE)
  }
  result
}

mv06g_build_group_rankings_v1 <- function(query_distances, scales,
                                           expected_training_samples,
                                           expected_query_pairs) {
  required <- c("group_id", "fold_id", "seed", "query_sample_id",
                "training_sample_id", "view_id", "homology_dimension",
                "distance", "outcome_label_state",
                "biological_outcomes_computed")
  if (!is.data.frame(query_distances) ||
      nrow(query_distances) != 4L * expected_query_pairs ||
      !all(required %in% names(query_distances)) ||
      !is.data.frame(scales) || nrow(scales) != 4L ||
      any(!is.finite(query_distances$distance)) ||
      any(query_distances$distance < 0) ||
      any(query_distances$outcome_label_state != "closed") ||
      any(as.logical(query_distances$biological_outcomes_computed))) {
    stop("MV6-G production query distances or scales are invalid.",
         call. = FALSE)
  }
  query_distances$component_id <- mapply(
    mv06g_component_id_v1, query_distances$view_id,
    query_distances$homology_dimension, USE.NAMES = FALSE
  )
  scale_map <- stats::setNames(scales$scale_value, scales$component_id)
  query_distances$normalized_distance <- query_distances$distance /
    unname(scale_map[query_distances$component_id])
  pair_key <- paste(query_distances$query_sample_id,
                    query_distances$training_sample_id, sep = "\r")
  methods <- mv06g_method_panel_v1()
  rows <- lapply(split(query_distances, pair_key), function(component_rows) {
    component_rows <- component_rows[match(
      c("cell_H0", "cell_H1", "gene_H0", "gene_H1"),
      component_rows$component_id
    ), , drop = FALSE]
    if (nrow(component_rows) != 4L || anyNA(component_rows$component_id) ||
        anyDuplicated(component_rows$component_id)) {
      stop("MV6-G production query pair lacks four components.", call. = FALSE)
    }
    z <- stats::setNames(component_rows$normalized_distance,
                         component_rows$component_id)
    cell <- mean(z[c("cell_H0", "cell_H1")])
    gene <- mean(z[c("gene_H0", "gene_H1")])
    values <- c(z, cell_composite = cell,
                fusion_gene_weight_025 = 0.75 * cell + 0.25 * gene,
                fusion_gene_weight_050 = 0.50 * cell + 0.50 * gene,
                fusion_gene_weight_075 = 0.25 * cell + 0.75 * gene,
                gene_composite = gene)
    data.frame(
      contract_id = "mv06g_label_closed_ranking_v1",
      group_id = component_rows$group_id[[1L]],
      fold_id = component_rows$fold_id[[1L]],
      seed = as.integer(component_rows$seed[[1L]]),
      query_sample_id = component_rows$query_sample_id[[1L]],
      training_sample_id = component_rows$training_sample_id[[1L]],
      method_order = methods$method_order, method_id = methods$method_id,
      gene_weight = methods$gene_weight,
      normalized_distance = unname(values[methods$method_id]),
      scale_fit_partition = "fold_training_pairs_only",
      tie_break = "canonical_training_sample_id",
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      fusion_evaluations = 0L, outcome_jobs = 0L, stringsAsFactors = FALSE
    )
  })
  rankings <- do.call(rbind, rows)
  blocks <- split(seq_len(nrow(rankings)), paste(rankings$query_sample_id,
                                                 rankings$method_id, sep = "\r"))
  rankings$rank <- integer(nrow(rankings))
  for (indices in blocks) {
    ordered <- order(rankings$normalized_distance[indices],
                     rankings$training_sample_id[indices], method = "radix")
    rankings$rank[indices[ordered]] <- seq_along(indices)
  }
  rankings <- rankings[order(rankings$query_sample_id, rankings$method_order,
                             rankings$rank, rankings$training_sample_id,
                             method = "radix"), , drop = FALSE]
  rownames(rankings) <- NULL
  blocks <- split(seq_len(nrow(rankings)), paste(rankings$query_sample_id,
                                                 rankings$method_id, sep = "\r"))
  expected_blocks <- expected_query_pairs / expected_training_samples * 9L
  if (nrow(rankings) != 9L * expected_query_pairs ||
      length(blocks) != expected_blocks ||
      any(vapply(blocks, function(indices) {
        !identical(sort(as.integer(rankings$rank[indices])),
                   seq_len(expected_training_samples))
      }, logical(1L)))) {
    stop("MV6-G production ranking axis drifted.", call. = FALSE)
  }
  rankings
}

mv06g_production_source_paths_v1 <- function() {
  c(
    "R/mv06g_production.R",
    "scripts/build_mv06g_production_prefreeze.R",
    "scripts/validate_mv06g_production_prefreeze.R",
    "scripts/validate_mv06g_production_prefreeze_repeat.R",
    "scripts/run_mv06g_group.R",
    "scripts/run_mv06g_production_rebind_monitor.R",
    "scripts/validate_mv06g_production_rebind.R",
    "tests/testthat/test-mv06g-production.R"
  )
}

mv06g_validate_production_group_v1 <- function(path, queue_row, parent,
                                                policy, source_group) {
  files <- c("training-distances.csv", "scales.csv", "rankings.csv",
             "metrics.csv", "status.csv")
  paths <- file.path(path, files)
  expected_pairs <- queue_row$training_samples *
    (queue_row$training_samples - 1L) / 2L
  expected_training_rows <- 4L * expected_pairs
  expected_rankings <- 9L * queue_row$biological_pairs
  if (!dir.exists(path) || !all(file.exists(paths))) {
    stop("MV6-G production group is incomplete.", call. = FALSE)
  }
  status <- utils::read.csv(paths[[5L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
  hashes <- vapply(paths[1:4], .mv06f_sha256, character(1L))
  bytes <- unname(file.info(paths[1:4])$size)
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      status$group_id != queue_row$group_id ||
      status$production_implementation_root_sha256 !=
        policy$production_implementation_root_sha256 ||
      status$parent_contract_sha256 != policy$parent_contract_sha256 ||
      status$rust_library_sha256 != parent$rust_library_sha256 ||
      status$source_diagrams_sha256 != source_group$diagrams_sha256 ||
      status$source_distances_sha256 != source_group$distances_sha256 ||
      !identical(unname(hashes), as.character(unlist(status[c(
        "training_distances_sha256", "scales_sha256", "rankings_sha256",
        "metrics_sha256"
      )], use.names = FALSE))) ||
      !identical(as.numeric(bytes), as.numeric(unlist(status[c(
        "training_distances_bytes", "scales_bytes", "rankings_bytes",
        "metrics_bytes"
      )], use.names = FALSE)))) {
    stop("MV6-G production group status is stale.", call. = FALSE)
  }
  training <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
  scales <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
  rankings <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
  if (nrow(training) != expected_training_rows ||
      any(table(training$component_id) != expected_pairs) ||
      anyDuplicated(training$pair_id) || nrow(scales) != 4L ||
      any(scales$training_pairs != expected_pairs) ||
      any(scales$scale_value <= sqrt(.Machine$double.eps)) ||
      nrow(rankings) != expected_rankings ||
      any(table(rankings$query_sample_id, rankings$method_id) !=
            queue_row$training_samples) ||
      any(training$outcome_label_state != "closed") ||
      any(scales$outcome_label_state != "closed") ||
      any(rankings$outcome_label_state != "closed") ||
      any(as.logical(training$biological_outcomes_computed)) ||
      any(as.logical(rankings$biological_outcomes_computed)) ||
      any(rankings$fusion_evaluations != 0L) || any(rankings$outcome_jobs != 0L)) {
    stop("MV6-G production group scientific artifacts are invalid.",
         call. = FALSE)
  }
  invisible(status)
}
