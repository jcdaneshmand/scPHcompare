# MV5-AL prediction-locked nested-256 cell-depth retrieval robustness helpers.

.mv05al_seeds <- 20260805:20260809
.mv05al_endpoints <- c(
  "cross_study_tissue_mrr_v1",
  "cross_study_tissue_1nn_balanced_accuracy_v1")

.mv05al_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05al_with_rng <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  force(expression)
}

.mv05al_boundary_exceedances <- function(null, observed) {
  tolerance <- 64 * .Machine$double.eps *
    pmax(1, abs(null), abs(observed))
  sum(abs(null) + tolerance >= abs(observed))
}

#' Construct canonical label-closed nested256 rankings for one frozen group.
#'
#' @param method_rows The accepted four-method MV5-AJ rows for one group.
#' @param group_row One MV5-AK queue row.
#' @return Canonically ordered prediction-side rankings.
mv05al_rank_group_v1 <- function(method_rows, group_row) {
  required <- c("robustness_group_id", "query_sample_id",
                "training_sample_id", "method_id", "distance",
                "outcome_label_state", "outcomes_computed")
  if (!is.data.frame(method_rows) || !all(required %in% names(method_rows)) ||
      !is.data.frame(group_row) || nrow(group_row) != 1L ||
      group_row$configuration_id != "nested_cells_256_pc30_euclidean_v1" ||
      !isTRUE(as.logical(group_row$coordinate_source_identity_exact)) ||
      !isTRUE(as.logical(group_row$nested192_prefix_identity_exact)) ||
      !isTRUE(as.logical(group_row$nested256_subset_384_identity_exact)) ||
      nrow(method_rows) != group_row$method_rows ||
      anyDuplicated(method_rows[c("query_sample_id", "training_sample_id",
                                  "method_id")]) ||
      any(!is.finite(method_rows$distance)) || any(method_rows$distance < 0) ||
      any(method_rows$outcome_label_state != "closed") ||
      any(tolower(as.character(method_rows$outcomes_computed)) == "true")) {
    stop("MV5-AL group method rows are malformed or not outcome-closed.",
         call. = FALSE)
  }
  expected_methods <- mv05ak_method_map_v1()
  expected_methods <- expected_methods[
    expected_methods$representation == group_row$representation,
    "nested256_method_id"]
  if (!setequal(method_rows$method_id, expected_methods)) {
    stop("MV5-AL group method axis drifted.", call. = FALSE)
  }
  keys <- interaction(method_rows$method_id, method_rows$query_sample_id,
                      drop = TRUE, lex.order = TRUE)
  result <- do.call(rbind, lapply(split(method_rows, keys), function(part) {
    part <- part[order(part$distance, part$training_sample_id,
                       method = "radix"), , drop = FALSE]
    runs <- rle(part$distance)
    part$neighbor_rank <- seq_len(nrow(part))
    part$distance_tie_size <- rep(runs$lengths, runs$lengths)
    part$distance_tied <- part$distance_tie_size > 1L
    part$tie_break_policy <-
      "ascending_distance_then_canonical_training_sample_id_radix_v1"
    part
  }))
  result$contract_id <- "mv05al_nested256_ranking_row_v1"
  result$fold_id <- group_row$fold_id
  result$seed <- as.integer(group_row$seed)
  result$representation <- group_row$representation
  result$configuration_id <- group_row$configuration_id
  result$prediction_locked <- TRUE
  result$labels_opened <- FALSE
  result$outcomes_computed <- FALSE
  result <- result[order(
    result$representation, result$fold_id, result$seed, result$method_id,
    result$query_sample_id, result$neighbor_rank, method = "radix"), , drop = FALSE]
  result$ranking_id <- vapply(seq_len(nrow(result)), function(index) {
    paste0("mv05al_ranking_v1:", .mv05al_digest(list(
      group = result$robustness_group_id[[index]],
      method = result$method_id[[index]],
      query = result$query_sample_id[[index]],
      training = result$training_sample_id[[index]],
      rank = result$neighbor_rank[[index]])))
  }, character(1L))
  rownames(result) <- NULL
  if (nrow(result) != group_row$method_rows ||
      anyDuplicated(result$ranking_id) ||
      any(table(result$method_id, result$query_sample_id) == 0L)) {
    stop("MV5-AL canonical group ranking completion failed.", call. = FALSE)
  }
  result[c(
    "contract_id", "ranking_id", "robustness_group_id", "fold_id", "seed",
    "representation", "configuration_id", "method_id", "query_sample_id",
    "training_sample_id", "distance", "neighbor_rank", "distance_tie_size",
    "distance_tied", "tie_break_policy", "prediction_locked",
    "labels_opened", "outcomes_computed")]
}

#' Evaluate one immutable NESTED256 ranking group after its prediction lock.
mv05al_evaluate_group_v1 <- function(rankings, labels,
                                     method_map = mv05ak_method_map_v1()) {
  required <- c("fold_id", "seed", "representation", "method_id",
                "query_sample_id", "training_sample_id", "distance",
                "neighbor_rank", "distance_tie_size", "distance_tied",
                "tie_break_policy", "prediction_locked", "labels_opened",
                "outcomes_computed")
  if (!is.data.frame(rankings) || !all(required %in% names(rankings)) ||
      !is.data.frame(labels) ||
      !all(c("sample_id", "study", "tissue") %in% names(labels)) ||
      anyDuplicated(labels$sample_id) || anyNA(labels) ||
      any(!rankings$prediction_locked) || any(rankings$labels_opened) ||
      any(rankings$outcomes_computed) ||
      any(rankings$tie_break_policy !=
            "ascending_distance_then_canonical_training_sample_id_radix_v1")) {
    stop("MV5-AL ranking or label input is malformed.", call. = FALSE)
  }
  representation <- unique(rankings$representation)
  if (length(representation) != 1L) stop("Mixed representation group.", call. = FALSE)
  map <- method_map[method_map$representation == representation, , drop = FALSE]
  label_index <- match(rankings$training_sample_id, labels$sample_id)
  query_index <- match(rankings$query_sample_id, labels$sample_id)
  if (anyNA(label_index) || anyNA(query_index)) {
    stop("MV5-AL ranking-to-label key join failed.", call. = FALSE)
  }
  rankings$.training_tissue <- labels$tissue[label_index]
  rankings$.training_study <- labels$study[label_index]
  heldout <- sub("^large_loso_v1:", "", rankings$fold_id)
  if (any(labels$study[query_index] != heldout) ||
      any(rankings$.training_study == heldout)) {
    stop("MV5-AL held-out study isolation failed.", call. = FALSE)
  }
  groups <- split(seq_len(nrow(rankings)), interaction(
    rankings$method_id, rankings$query_sample_id,
    drop = TRUE, lex.order = TRUE))
  observations <- do.call(rbind, lapply(groups, function(index) {
    part <- rankings[index, , drop = FALSE]
    part <- part[order(part$neighbor_rank), , drop = FALSE]
    if (!identical(as.integer(part$neighbor_rank), seq_len(nrow(part))) ||
        !identical(part$training_sample_id,
                   part$training_sample_id[order(
                     part$distance, part$training_sample_id,
                     method = "radix")])) {
      stop("MV5-AL immutable ranking sequence drifted.", call. = FALSE)
    }
    q <- match(part$query_sample_id[[1L]], labels$sample_id)
    q_tissue <- labels$tissue[[q]]
    q_study <- labels$study[[q]]
    same <- which(part$.training_tissue == q_tissue)
    status <- if (!q_tissue %in% .mv05e_eligible_tissues) {
      "single_study_tissue_not_estimable"
    } else if (!length(same)) {
      "training_tissue_absent_not_estimable"
    } else "estimable"
    first <- if (length(same)) same[[1L]] else NA_integer_
    method <- part$method_id[[1L]]
    family <- map$family_id[match(method, map$nested256_method_id)]
    data.frame(
      contract_id = "mv05al_nested256_query_method_outcome_v1",
      fold_id = part$fold_id[[1L]], held_out_study = q_study,
      seed = as.integer(part$seed[[1L]]),
      representation = representation, family_id = family,
      method_id = method, query_sample_id = part$query_sample_id[[1L]],
      query_tissue = q_tissue, training_samples = nrow(part),
      training_studies = length(unique(part$.training_study)),
      first_same_tissue_sample_id = if (is.na(first)) NA_character_ else
        part$training_sample_id[[first]],
      first_same_tissue_rank = if (is.na(first)) NA_integer_ else
        as.integer(part$neighbor_rank[[first]]),
      reciprocal_rank = if (is.na(first)) NA_real_ else
        1 / as.numeric(part$neighbor_rank[[first]]),
      nearest_sample_id = part$training_sample_id[[1L]],
      nearest_tissue = part$.training_tissue[[1L]],
      one_nn_correct = if (status == "estimable")
        part$.training_tissue[[1L]] == q_tissue else NA,
      nearest_distance_tied = as.logical(part$distance_tied[[1L]]),
      endpoint_status = status,
      labels_opened_after_prediction_lock = TRUE,
      upstream_refit = FALSE, reranked_after_label_open = FALSE,
      outcomes_computed = TRUE, evaluation_executed = TRUE,
      stringsAsFactors = FALSE)
  }))
  observations <- observations[order(
    observations$representation, observations$family_id,
    observations$fold_id, observations$seed, observations$query_sample_id,
    method = "radix"), , drop = FALSE]
  rownames(observations) <- NULL
  expected <- length(unique(rankings$query_sample_id)) * 4L
  if (nrow(observations) != expected ||
      any(observations$endpoint_status != "estimable")) {
    stop("MV5-AL query-method outcome completion failed.", call. = FALSE)
  }
  observations
}

#' Pair accepted 384-cell and new nested-256 query-method endpoints.
mv05al_pair_baseline_v1 <- function(nested256, baseline_sct, baseline_integrated,
                                    method_map = mv05ak_method_map_v1()) {
  required <- c("fold_id", "held_out_study", "seed", "method_id",
                "query_sample_id", "query_tissue", "training_samples",
                "training_studies", "reciprocal_rank", "one_nn_correct",
                "endpoint_status", "upstream_refit", "reranked_after_label_open")
  if (!all(required %in% names(nested256)) ||
      !all(required %in% names(baseline_sct)) ||
      !all(required %in% names(baseline_integrated))) {
    stop("MV5-AL paired endpoint schema is incomplete.", call. = FALSE)
  }
  baseline <- rbind(
    transform(baseline_sct, representation = "sct_whole"),
    transform(baseline_integrated, representation = "inductive_integrated"))
  baseline$family_id <- NA_character_
  for (index in seq_len(nrow(method_map))) {
    hit <- baseline$representation == method_map$representation[[index]] &
      baseline$method_id == method_map$baseline_method_id[[index]]
    baseline$family_id[hit] <- method_map$family_id[[index]]
  }
  baseline <- baseline[!is.na(baseline$family_id), , drop = FALSE]
  keys <- c("fold_id", "held_out_study", "seed", "representation",
            "family_id", "query_sample_id", "query_tissue")
  key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
  nested256_key <- key(nested256); baseline_key <- key(baseline)
  if (nrow(nested256) != 3600L || nrow(baseline) != 3600L ||
      anyDuplicated(nested256_key) || anyDuplicated(baseline_key) ||
      !setequal(nested256_key, baseline_key)) {
    stop("MV5-AL Euclidean/nested256 query endpoint pairing failed.", call. = FALSE)
  }
  baseline <- baseline[match(nested256_key, baseline_key), , drop = FALSE]
  common <- c("held_out_study", "query_tissue", "training_samples",
              "training_studies", "endpoint_status")
  if (!all(vapply(common, function(name) identical(nested256[[name]], baseline[[name]]),
                  logical(1L))) ||
      any(nested256$upstream_refit) || any(nested256$reranked_after_label_open) ||
      any(baseline$upstream_refit) || any(baseline$reranked_after_label_open)) {
    stop("MV5-AL paired identity, denominator, or no-refit contract failed.",
         call. = FALSE)
  }
  data.frame(
    contract_id = "mv05al_paired_query_method_endpoint_v1",
    nested256[keys],
    baseline_method_id = baseline$method_id,
    nested256_method_id = nested256$method_id,
    baseline_reciprocal_rank = as.numeric(baseline$reciprocal_rank),
    nested256_reciprocal_rank = as.numeric(nested256$reciprocal_rank),
    direct_reciprocal_rank_change =
      as.numeric(nested256$reciprocal_rank - baseline$reciprocal_rank),
    baseline_one_nn_correct = as.logical(baseline$one_nn_correct),
    nested256_one_nn_correct = as.logical(nested256$one_nn_correct),
    direct_one_nn_change = as.numeric(nested256$one_nn_correct) -
      as.numeric(baseline$one_nn_correct),
    endpoint_status = nested256$endpoint_status,
    outcomes_computed = TRUE, evaluation_executed = TRUE,
    method_selection_executed = FALSE,
    stringsAsFactors = FALSE)
}

#' Construct every frozen query-level nested256 robustness estimand.
mv05al_query_estimands_v1 <- function(
    paired, registry = mv05ak_estimand_registry_v1()) {
  keys <- c("fold_id", "held_out_study", "seed", "representation",
            "query_sample_id", "query_tissue")
  direct <- do.call(rbind, lapply(seq_len(nrow(paired)), function(index) {
    part <- paired[index, , drop = FALSE]
    rbind(
      data.frame(part[keys], family_id = part$family_id,
        endpoint_id = .mv05al_endpoints[[1L]],
        estimand_type = "direct_nested256_cells_minus_384_cells",
        value = part$direct_reciprocal_rank_change, stringsAsFactors = FALSE),
      data.frame(part[keys], family_id = part$family_id,
        endpoint_id = .mv05al_endpoints[[2L]],
        estimand_type = "direct_nested256_cells_minus_384_cells",
        value = part$direct_one_nn_change, stringsAsFactors = FALSE))
  }))
  did <- do.call(rbind, lapply(c("h0", "h1"), function(family) {
    topology <- paired[paired$family_id == family, , drop = FALSE]
    energy <- paired[paired$family_id == "energy", , drop = FALSE]
    match_keys <- c("fold_id", "held_out_study", "seed", "representation",
                    "query_sample_id", "query_tissue")
    k <- function(x) do.call(paste, c(x[match_keys], sep = "\r"))
    energy <- energy[match(k(topology), k(energy)), , drop = FALSE]
    if (anyNA(energy$query_sample_id)) stop("DID energy pairing failed.", call. = FALSE)
    rbind(
      data.frame(topology[match_keys], family_id = family,
        endpoint_id = .mv05al_endpoints[[1L]],
        estimand_type =
          "topology_increment_nested256_minus_384_cells_difference_in_differences",
        value = topology$direct_reciprocal_rank_change -
          energy$direct_reciprocal_rank_change, stringsAsFactors = FALSE),
      data.frame(topology[match_keys], family_id = family,
        endpoint_id = .mv05al_endpoints[[2L]],
        estimand_type =
          "topology_increment_nested256_minus_384_cells_difference_in_differences",
        value = topology$direct_one_nn_change - energy$direct_one_nn_change,
        stringsAsFactors = FALSE))
  }))
  result <- rbind(direct, did)
  lookup <- paste(registry$estimand_type, registry$representation,
                  registry$family_id, registry$endpoint_id, sep = "\r")
  observed <- paste(result$estimand_type, result$representation,
                    result$family_id, result$endpoint_id, sep = "\r")
  result$estimand_id <- registry$estimand_id[match(observed, lookup)]
  result$estimand_order <- registry$estimand_order[match(observed, lookup)]
  result$estimand_role <- registry$estimand_role[match(observed, lookup)]
  if (nrow(result) != 10800L || anyNA(result$estimand_id) ||
      any(!is.finite(result$value))) {
    stop("MV5-AL query estimand construction failed.", call. = FALSE)
  }
  result$contract_id <- "mv05al_query_estimand_v1"
  result <- result[order(result$estimand_order, result$query_tissue,
                         result$held_out_study, result$query_sample_id,
                         result$seed, method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  result[c("contract_id", "estimand_id", "estimand_order", "estimand_type",
           "estimand_role", "endpoint_id", "representation", "family_id",
           "fold_id", "held_out_study", "seed", "query_sample_id",
           "query_tissue", "value")]
}

#' Aggregate query estimands in the frozen seed/sample/tissue order.
mv05al_summarize_estimands_v1 <- function(query) {
  grouping <- c("estimand_id", "estimand_order", "estimand_type",
                "estimand_role", "endpoint_id", "representation", "family_id")
  sample <- stats::aggregate(
    value ~ estimand_id + estimand_order + estimand_type + estimand_role +
      endpoint_id + representation + family_id + query_sample_id +
      query_tissue + held_out_study,
    data = query, FUN = mean)
  counts <- stats::aggregate(
    seed ~ estimand_id + query_sample_id, data = query, FUN = length)
  names(counts)[names(counts) == "seed"] <- "completed_seeds"
  sample <- merge(sample, counts, by = c("estimand_id", "query_sample_id"),
                  sort = FALSE)
  names(sample)[names(sample) == "value"] <- "estimate"
  sample$contract_id <- "mv05al_sample_estimand_v1"
  tissue <- stats::aggregate(
    estimate ~ estimand_id + estimand_order + estimand_type + estimand_role +
      endpoint_id + representation + family_id + query_tissue,
    data = sample, FUN = mean)
  tissue_counts <- stats::aggregate(
    query_sample_id ~ estimand_id + query_tissue, data = sample, FUN = length)
  names(tissue_counts)[names(tissue_counts) == "query_sample_id"] <- "samples"
  tissue <- merge(tissue, tissue_counts,
                  by = c("estimand_id", "query_tissue"), sort = FALSE)
  tissue$study_blocks <- vapply(seq_len(nrow(tissue)), function(index) {
    length(unique(sample$held_out_study[
      sample$estimand_id == tissue$estimand_id[[index]] &
        sample$query_tissue == tissue$query_tissue[[index]]]))
  }, integer(1L))
  tissue$contract_id <- "mv05al_tissue_estimand_v1"
  macro <- stats::aggregate(
    estimate ~ estimand_id + estimand_order + estimand_type + estimand_role +
      endpoint_id + representation + family_id,
    data = tissue, FUN = mean)
  macro$contract_id <- "mv05al_macro_estimand_v1"
  macro$tissues <- 5L; macro$samples <- 90L; macro$study_blocks <- 15L
  macro$completed_seeds <- 5L
  order_rows <- function(x) {
    extra <- intersect(c("query_tissue", "held_out_study", "query_sample_id"),
                       names(x))
    args <- c(list(x$estimand_order), x[extra], list(method = "radix"))
    x <- x[do.call(order, args), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  sample <- order_rows(sample); tissue <- order_rows(tissue); macro <- order_rows(macro)
  if (nrow(sample) != 2160L || nrow(tissue) != 120L || nrow(macro) != 24L ||
      any(sample$completed_seeds != 5L)) {
    stop("MV5-AL frozen aggregation completion failed.", call. = FALSE)
  }
  list(query = query, sample = sample, tissue = tissue, macro = macro)
}

.mv05al_macro <- function(base, values, weights = NULL) {
  tissues <- sort(unique(base$query_tissue), method = "radix")
  tissue_values <- vapply(tissues, function(tissue) {
    hit <- base$query_tissue == tissue
    if (is.null(weights)) mean(values[hit]) else
      stats::weighted.mean(values[hit], weights[hit])
  }, numeric(1L))
  mean(tissue_values)
}

#' Apply the frozen paired blocked uncertainty and four-test family.
mv05al_block_inference_v1 <- function(
    sample_summary, registry = mv05ak_estimand_registry_v1(),
    bootstrap_replicates = 2000L, bootstrap_seed = 20260814L,
    randomization_replicates = 9999L, randomization_seed = 20260815L) {
  sample_ids <- sort(unique(sample_summary$query_sample_id), method = "radix")
  first_id <- registry$estimand_id[order(registry$estimand_order)][[1L]]
  base <- sample_summary[sample_summary$estimand_id == first_id,
                         c("query_sample_id", "query_tissue", "held_out_study")]
  base <- base[match(sample_ids, base$query_sample_id), , drop = FALSE]
  estimand_ids <- registry$estimand_id[order(registry$estimand_order)]
  values <- vapply(estimand_ids, function(id) {
    part <- sample_summary[sample_summary$estimand_id == id, , drop = FALSE]
    part$estimate[match(sample_ids, part$query_sample_id)]
  }, numeric(length(sample_ids)))
  colnames(values) <- estimand_ids
  if (nrow(base) != 90L || anyNA(base) || any(!is.finite(values))) {
    stop("MV5-AL inference sample matrix is incomplete.", call. = FALSE)
  }
  studies <- sort(unique(base$held_out_study), method = "radix")
  study_tissue <- unique(base[c("held_out_study", "query_tissue")])
  if (anyDuplicated(study_tissue$held_out_study)) {
    stop("MV5-AL bootstrap requires tissue-homogeneous study blocks.",
         call. = FALSE)
  }
  count_matrix <- .mv05al_with_rng(bootstrap_seed, {
    counts <- matrix(0L, bootstrap_replicates, length(studies),
                     dimnames = list(NULL, studies))
    for (tissue in sort(unique(base$query_tissue), method = "radix")) {
      ids <- sort(study_tissue$held_out_study[
        study_tissue$query_tissue == tissue], method = "radix")
      for (replicate_id in seq_len(bootstrap_replicates)) {
        draw <- sample(ids, length(ids), replace = TRUE)
        counts[replicate_id, match(names(table(draw)), studies)] <-
          as.integer(table(draw))
      }
    }
    counts
  })
  bootstrap <- matrix(NA_real_, bootstrap_replicates, ncol(values),
                      dimnames = list(NULL, estimand_ids))
  study_index <- match(base$held_out_study, studies)
  for (replicate_id in seq_len(bootstrap_replicates)) {
    weights <- count_matrix[replicate_id, study_index]
    for (column in seq_len(ncol(values))) {
      bootstrap[replicate_id, column] <-
        .mv05al_macro(base, values[, column], weights)
    }
  }
  intervals <- do.call(rbind, lapply(seq_along(estimand_ids), function(index) {
    row <- registry[match(estimand_ids[[index]], registry$estimand_id), ]
    ci <- stats::quantile(bootstrap[, index], c(0.025, 0.975),
                          names = FALSE, type = 7)
    data.frame(
      contract_id = "mv05al_estimand_interval_v1",
      estimand_id = estimand_ids[[index]], estimand_order = row$estimand_order,
      estimand_type = row$estimand_type, estimand_role = row$estimand_role,
      endpoint_id = row$endpoint_id, representation = row$representation,
      family_id = row$family_id,
      estimate = .mv05al_macro(base, values[, index]),
      ci_lower = ci[[1L]], ci_upper = ci[[2L]],
      bootstrap_replicates = bootstrap_replicates,
      bootstrap_seed = bootstrap_seed, tissues = 5L, study_blocks = 15L,
      samples = 90L, completed_seeds = 5L,
      inference_status = "estimable", stringsAsFactors = FALSE)
  }))
  primary <- intervals[intervals$estimand_role ==
                         "confirmatory_nested256_sensitivity", , drop = FALSE]
  primary$contract_id <- "mv05al_primary_contrast_v1"
  primary$multiplicity_family <- "F1_nested256_topology_increment_mrr_four_tests"
  primary$randomization_replicates <- randomization_replicates
  primary$randomization_seed <- randomization_seed
  signs <- .mv05al_with_rng(randomization_seed, matrix(
    sample(c(-1, 1), randomization_replicates * length(studies), replace = TRUE),
    nrow = randomization_replicates, ncol = length(studies),
    dimnames = list(NULL, studies)))
  primary$raw_p_value <- NA_real_; primary$holm_p_value <- NA_real_
  null_hash <- character(nrow(primary)); exceedances <- integer(nrow(primary))
  for (index in seq_len(nrow(primary))) {
    column <- match(primary$estimand_id[[index]], estimand_ids)
    null <- vapply(seq_len(randomization_replicates), function(replicate_id) {
      signed <- values[, column] * signs[replicate_id, study_index]
      .mv05al_macro(base, signed)
    }, numeric(1L))
    exceedances[[index]] <- .mv05al_boundary_exceedances(
      null, primary$estimate[[index]])
    primary$raw_p_value[[index]] <-
      (exceedances[[index]] + 1) / (randomization_replicates + 1)
    null_hash[[index]] <- .mv05al_digest(null)
  }
  primary$holm_p_value <- stats::p.adjust(primary$raw_p_value, method = "holm")
  primary <- primary[order(primary$estimand_order), , drop = FALSE]
  bootstrap_audit <- data.frame(
    contract_id = "mv05al_bootstrap_audit_v1",
    bootstrap_replicates = bootstrap_replicates,
    bootstrap_seed = bootstrap_seed, tissues = 5L, study_blocks = 15L,
    estimands = 24L, block_count_matrix_sha256 = .mv05al_digest(count_matrix),
    replicate_matrix_sha256 = .mv05al_digest(bootstrap),
    resampling_unit = "paired_tissue_stratified_heldout_study_block",
    same_draw_all_configurations_methods_estimands = TRUE,
    seeds_treated_as_independent = FALSE, stringsAsFactors = FALSE)
  randomization_audit <- data.frame(
    contract_id = "mv05al_randomization_audit_v1",
    multiplicity_family = "F1_nested256_topology_increment_mrr_four_tests",
    estimand_id = primary$estimand_id,
    randomization_replicates = randomization_replicates,
    randomization_seed = randomization_seed, study_blocks = length(studies),
    sign_matrix_sha256 = .mv05al_digest(signs),
    null_distribution_sha256 = null_hash,
    exceedance_count = exceedances,
    boundary_policy = "absolute_two_sided_64eps_ties_count_as_exceedance_v1",
    epsilon_multiplier = 64L, two_sided = TRUE,
    exact_zero_p_value_possible = FALSE, stringsAsFactors = FALSE)
  list(intervals = intervals[order(intervals$estimand_order), , drop = FALSE],
       primary = primary, bootstrap_audit = bootstrap_audit,
       randomization_audit = randomization_audit,
       bootstrap_counts = count_matrix, bootstrap_estimates = bootstrap,
       sign_matrix = signs)
}
