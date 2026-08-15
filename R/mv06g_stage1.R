# Internal MV6-G stage-one label-closed training-scale/ranking helpers.

mv06g_component_id_v1 <- function(view_id, homology_dimension) {
  view_id <- match.arg(view_id, c("cell_topology_v1", "gene_topology_v1"))
  homology_dimension <- match.arg(homology_dimension, c("H0", "H1"))
  paste(if (view_id == "cell_topology_v1") "cell" else "gene",
        homology_dimension, sep = "_")
}

mv06g_training_pair_id_v1 <- function(group_id, first_sample_id,
                                       second_sample_id, view_id,
                                       homology_dimension) {
  samples <- sort(c(as.character(first_sample_id),
                    as.character(second_sample_id)), method = "radix")
  if (length(unique(samples)) != 2L) {
    stop("MV6-G training pairs require two distinct samples.", call. = FALSE)
  }
  identity <- list(
    contract_id = "mv06g_training_landscape_pair_identity_v1",
    group_id = as.character(group_id), first_sample_id = samples[[1L]],
    second_sample_id = samples[[2L]],
    view_id = match.arg(view_id,
                        c("cell_topology_v1", "gene_topology_v1")),
    homology_dimension = match.arg(homology_dimension, c("H0", "H1"))
  )
  paste0("mv06g_training_pair_v1:", digest::digest(
    identity, algo = "sha256", serialize = TRUE
  ))
}

mv06g_fit_component_scales_v1 <- function(training_distances) {
  required <- c("group_id", "view_id", "homology_dimension", "component_id",
                "distance", "outcome_label_state",
                "biological_outcomes_computed")
  if (!is.data.frame(training_distances) ||
      !all(required %in% names(training_distances)) ||
      nrow(training_distances) != 8320L ||
      any(!is.finite(training_distances$distance)) ||
      any(training_distances$distance < 0) ||
      any(training_distances$outcome_label_state != "closed") ||
      any(as.logical(training_distances$biological_outcomes_computed))) {
    stop("MV6-G training distances are incomplete or invalid.", call. = FALSE)
  }
  components <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
  rows <- lapply(seq_along(components), function(index) {
    component <- components[[index]]
    values <- training_distances$distance[
      training_distances$component_id == component
    ]
    scale <- stats::median(values)
    if (length(values) != 2080L || !is.finite(scale) ||
        scale <= sqrt(.Machine$double.eps)) {
      stop("MV6-G component scale is absent or degenerate: ", component,
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
    stop("MV6-G did not fit exactly four component scales.", call. = FALSE)
  }
  result
}

mv06g_build_rankings_v1 <- function(query_distances, scales) {
  required <- c("group_id", "fold_id", "seed", "query_sample_id",
                "training_sample_id", "view_id", "homology_dimension",
                "distance", "outcome_label_state",
                "biological_outcomes_computed")
  if (!is.data.frame(query_distances) || nrow(query_distances) != 6500L ||
      !all(required %in% names(query_distances)) ||
      !is.data.frame(scales) || nrow(scales) != 4L ||
      any(!is.finite(query_distances$distance)) ||
      any(query_distances$distance < 0) ||
      any(query_distances$outcome_label_state != "closed") ||
      any(as.logical(query_distances$biological_outcomes_computed))) {
    stop("MV6-G query distances or training scales are invalid.", call. = FALSE)
  }
  query_distances$component_id <- mapply(
    mv06g_component_id_v1, query_distances$view_id,
    query_distances$homology_dimension, USE.NAMES = FALSE
  )
  scale_map <- stats::setNames(scales$scale_value, scales$component_id)
  if (!setequal(names(scale_map), c("cell_H0", "cell_H1", "gene_H0",
                                    "gene_H1"))) {
    stop("MV6-G scale component axis is invalid.", call. = FALSE)
  }
  query_distances$normalized_distance <- query_distances$distance /
    unname(scale_map[query_distances$component_id])
  pair_key <- paste(query_distances$query_sample_id,
                    query_distances$training_sample_id, sep = "\r")
  split_rows <- split(query_distances, pair_key)
  methods <- mv06g_method_panel_v1()
  ranking_rows <- lapply(split_rows, function(rows) {
    rows <- rows[match(c("cell_H0", "cell_H1", "gene_H0", "gene_H1"),
                       rows$component_id), , drop = FALSE]
    if (nrow(rows) != 4L || anyNA(rows$component_id) ||
        anyDuplicated(rows$component_id)) {
      stop("MV6-G query pair lacks four unique components.", call. = FALSE)
    }
    z <- stats::setNames(rows$normalized_distance, rows$component_id)
    cell <- mean(z[c("cell_H0", "cell_H1")])
    gene <- mean(z[c("gene_H0", "gene_H1")])
    values <- c(z, cell_composite = cell,
                fusion_gene_weight_025 = 0.75 * cell + 0.25 * gene,
                fusion_gene_weight_050 = 0.50 * cell + 0.50 * gene,
                fusion_gene_weight_075 = 0.25 * cell + 0.75 * gene,
                gene_composite = gene)
    values <- unname(values[methods$method_id])
    data.frame(
      contract_id = "mv06g_label_closed_ranking_v1",
      group_id = rows$group_id[[1L]], fold_id = rows$fold_id[[1L]],
      seed = as.integer(rows$seed[[1L]]),
      query_sample_id = rows$query_sample_id[[1L]],
      training_sample_id = rows$training_sample_id[[1L]],
      method_order = methods$method_order, method_id = methods$method_id,
      gene_weight = methods$gene_weight, normalized_distance = values,
      scale_fit_partition = "fold_training_pairs_only",
      tie_break = "canonical_training_sample_id",
      outcome_label_state = "closed", biological_outcomes_computed = FALSE,
      fusion_evaluations = 0L, outcome_jobs = 0L,
      stringsAsFactors = FALSE
    )
  })
  rankings <- do.call(rbind, ranking_rows)
  blocks <- split(seq_len(nrow(rankings)), paste(rankings$query_sample_id,
                                                 rankings$method_id, sep = "\r"))
  rankings$rank <- integer(nrow(rankings))
  for (indices in blocks) {
    order_index <- order(rankings$normalized_distance[indices],
                         rankings$training_sample_id[indices],
                         method = "radix")
    rankings$rank[indices[order_index]] <- seq_along(indices)
  }
  rankings <- rankings[order(rankings$query_sample_id, rankings$method_order,
                             rankings$rank, rankings$training_sample_id,
                             method = "radix"), , drop = FALSE]
  rownames(rankings) <- NULL
  if (nrow(rankings) != 14625L || length(blocks) != 225L ||
      any(!is.finite(rankings$normalized_distance)) ||
      any(vapply(blocks, function(indices) {
        !identical(sort(as.integer(rankings$rank[indices])), 1:65)
      }, logical(1L)))) {
    stop("MV6-G label-closed ranking cardinality or ranks drifted.",
         call. = FALSE)
  }
  rankings
}

mv06g_stage1_source_paths_v1 <- function() {
  c(
    "R/mv06g_stage1.R",
    "scripts/build_mv06g_stage1_launch_prefreeze.R",
    "scripts/validate_mv06g_stage1_launch_prefreeze.R",
    "scripts/validate_mv06g_stage1_launch_repeat.R",
    "scripts/run_mv06g_stage1_group.R",
    "scripts/run_mv06g_stage1_monitor.R",
    "scripts/validate_mv06g_stage1.R",
    "scripts/validate_mv06g_stage1_persim.py",
    "scripts/check_mv06g_stage1_resume.R",
    "tests/testthat/test-mv06g-stage1.R"
  )
}

mv06g_validate_stage1_group_v1 <- function(path, queue_row, parent_contract,
                                            launch, source_group) {
  required_files <- c("training-distances.csv", "scales.csv", "rankings.csv",
                      "metrics.csv", "status.csv")
  paths <- file.path(path, required_files)
  if (!dir.exists(path) || !all(file.exists(paths)) ||
      !is.data.frame(queue_row) || nrow(queue_row) != 1L ||
      !is.data.frame(parent_contract) || nrow(parent_contract) != 1L ||
      !is.data.frame(launch) || nrow(launch) != 1L ||
      !is.data.frame(source_group) || nrow(source_group) != 1L) {
    stop("MV6-G stage-one directory or identities are incomplete.", call. = FALSE)
  }
  status <- utils::read.csv(paths[[5L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
  artifact_hashes <- vapply(paths[1:4], .mv06f_sha256, character(1L))
  artifact_bytes <- unname(file.info(paths[1:4])$size)
  if (nrow(status) != 1L || status$completion_state != "complete" ||
      status$group_id != queue_row$group_id ||
      status$parent_contract_sha256 != launch$parent_contract_sha256 ||
      status$stage1_implementation_root_sha256 !=
        launch$stage1_implementation_root_sha256 ||
      status$rust_library_sha256 != parent_contract$rust_library_sha256 ||
      status$source_diagrams_sha256 != source_group$diagrams_sha256 ||
      status$source_distances_sha256 != source_group$distances_sha256 ||
      status$outcome_label_state != "closed" ||
      as.logical(status$biological_outcomes_computed) ||
      any(unlist(status[c("fusion_evaluations", "outcome_jobs")]) != 0) ||
      !identical(unname(artifact_hashes), as.character(unlist(status[c(
        "training_distances_sha256", "scales_sha256", "rankings_sha256",
        "metrics_sha256"
      )], use.names = FALSE))) ||
      !identical(as.numeric(artifact_bytes), as.numeric(unlist(status[c(
        "training_distances_bytes", "scales_bytes", "rankings_bytes",
        "metrics_bytes"
      )], use.names = FALSE)))) {
    stop("MV6-G stage-one status or artifact identity is stale.", call. = FALSE)
  }
  training <- utils::read.csv(paths[[1L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
  scales <- utils::read.csv(paths[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
  rankings <- utils::read.csv(paths[[3L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
  component_counts <- table(training$component_id)
  rank_counts <- table(rankings$query_sample_id, rankings$method_id)
  if (nrow(training) != 8320L || anyDuplicated(training$pair_id) ||
      length(component_counts) != 4L || any(component_counts != 2080L) ||
      any(!as.logical(training$exact)) ||
      any(!as.logical(training$all_active_levels)) ||
      any(as.logical(training$level_cap_applied)) ||
      any(training$rust_status != 0L) || nrow(scales) != 4L ||
      any(scales$scale_value <= sqrt(.Machine$double.eps)) ||
      any(scales$query_rows_used != 0L) || any(scales$labels_used != 0L) ||
      nrow(rankings) != 14625L || any(rank_counts != 65L) ||
      !setequal(rankings$method_id, mv06g_method_panel_v1()$method_id) ||
      any(training$outcome_label_state != "closed") ||
      any(scales$outcome_label_state != "closed") ||
      any(rankings$outcome_label_state != "closed") ||
      any(as.logical(training$biological_outcomes_computed)) ||
      any(as.logical(scales$biological_outcomes_computed)) ||
      any(as.logical(rankings$biological_outcomes_computed))) {
    stop("MV6-G stage-one scientific artifacts are invalid.", call. = FALSE)
  }
  invisible(status)
}
