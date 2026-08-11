# MV5-S prediction-locked clustering-outcome execution helpers.

.mv05s_digest <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE)
}

.mv05s_comb2 <- function(value) {
  value * (value - 1) / 2
}

.mv05s_assert_partition_pair <- function(first, second) {
  if (length(first) != length(second) || length(first) < 2L ||
      anyNA(first) || anyNA(second)) {
    stop("MV5-S partition vectors must be complete, paired, and nontrivial.",
         call. = FALSE)
  }
  invisible(TRUE)
}

mv05s_adjusted_rand_index_v1 <- function(first, second) {
  .mv05s_assert_partition_pair(first, second)
  contingency <- table(first, second)
  pairs <- sum(.mv05s_comb2(contingency))
  row_pairs <- sum(.mv05s_comb2(rowSums(contingency)))
  column_pairs <- sum(.mv05s_comb2(colSums(contingency)))
  total_pairs <- .mv05s_comb2(sum(contingency))
  expected <- row_pairs * column_pairs / total_pairs
  denominator <- 0.5 * (row_pairs + column_pairs) - expected
  if (!is.finite(denominator) || abs(denominator) <= .Machine$double.eps) {
    return(NA_real_)
  }
  (pairs - expected) / denominator
}

mv05s_nmi_max_v1 <- function(first, second) {
  .mv05s_assert_partition_pair(first, second)
  contingency <- table(first, second)
  probabilities <- contingency / sum(contingency)
  row_probabilities <- rowSums(probabilities)
  column_probabilities <- colSums(probabilities)
  entropy <- function(values) {
    values <- values[values > 0]
    -sum(values * log(values))
  }
  first_entropy <- entropy(row_probabilities)
  second_entropy <- entropy(column_probabilities)
  joint_entropy <- entropy(as.numeric(probabilities))
  denominator <- max(first_entropy, second_entropy)
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    return(NA_real_)
  }
  (first_entropy + second_entropy - joint_entropy) / denominator
}

mv05s_metric_v1 <- function(metric_id, partition, labels) {
  value <- switch(
    metric_id,
    adjusted_rand_index = mv05s_adjusted_rand_index_v1(partition, labels),
    normalized_mutual_information = mv05s_nmi_max_v1(partition, labels),
    stop("Unsupported MV5-S metric: ", metric_id, call. = FALSE)
  )
  list(
    estimate = value,
    status = if (is.finite(value)) "completed" else
      "degenerate_label_or_partition_metric_not_identifiable")
}

mv05s_jackknife_se_v1 <- function(values) {
  if (length(values) < 2L || any(!is.finite(values))) return(NA_real_)
  leave_one_out <- vapply(seq_along(values), function(index) {
    mean(values[-index])
  }, numeric(1L))
  sqrt((length(values) - 1) / length(values) *
         sum((leave_one_out - mean(leave_one_out))^2))
}

.mv05s_queue_value <- function(queue_row, column) {
  if (!is.data.frame(queue_row) || nrow(queue_row) != 1L ||
      !column %in% names(queue_row)) {
    stop("MV5-S queue row is malformed.", call. = FALSE)
  }
  queue_row[[column]][[1L]]
}

.mv05s_join_label <- function(sample_ids, labels, label_axis) {
  if (!all(c("sample_id", "study", "tissue", "approach") %in% names(labels)) ||
      anyNA(labels) || anyDuplicated(labels$sample_id) ||
      !label_axis %in% c("study", "tissue", "approach")) {
    stop("MV5-S labels are malformed.", call. = FALSE)
  }
  index <- match(sample_ids, labels$sample_id)
  if (anyNA(index)) stop("MV5-S sample-label join is incomplete.", call. = FALSE)
  labels[[label_axis]][index]
}

mv05s_evaluate_training_unit_v1 <- function(queue_row, selected, labels) {
  if (.mv05s_queue_value(queue_row, "evaluation_scope") !=
      "overlapping_training_partition_alignment") {
    stop("MV5-S training evaluator received a held-out unit.", call. = FALSE)
  }
  group_id <- .mv05s_queue_value(queue_row, "analysis_group_id")
  algorithm_id <- .mv05s_queue_value(queue_row, "algorithm_id")
  expected_n <- as.integer(.mv05s_queue_value(queue_row, "training_samples"))
  part <- selected[selected$analysis_group_id == group_id &
                     selected$algorithm_id == algorithm_id, , drop = FALSE]
  seeds <- sort(unique(as.integer(part$seed)))
  if (!identical(seeds, 20260805:20260809)) {
    stop("MV5-S training unit has an incomplete seed axis.", call. = FALSE)
  }
  rows <- lapply(seeds, function(seed_id) {
    seed_part <- part[as.integer(part$seed) == seed_id, , drop = FALSE]
    if (nrow(seed_part) != expected_n || anyDuplicated(seed_part$sample_id) ||
        length(unique(seed_part$cluster)) != unique(seed_part$selected_k)) {
      stop("MV5-S training partition axis is malformed.", call. = FALSE)
    }
    label_values <- .mv05s_join_label(
      seed_part$sample_id, labels,
      .mv05s_queue_value(queue_row, "label_axis"))
    metric <- mv05s_metric_v1(
      .mv05s_queue_value(queue_row, "metric_id"),
      seed_part$cluster, label_values)
    data.frame(
      evaluation_unit_id = .mv05s_queue_value(queue_row, "evaluation_unit_id"),
      analysis_group_id = group_id,
      fold_id = .mv05s_queue_value(queue_row, "fold_id"),
      representation = .mv05s_queue_value(queue_row, "representation"),
      distance_id = .mv05s_queue_value(queue_row, "distance_id"),
      algorithm_id = algorithm_id,
      algorithm_role = .mv05s_queue_value(queue_row, "algorithm_role"),
      endpoint_id = .mv05s_queue_value(queue_row, "endpoint_id"),
      label_axis = .mv05s_queue_value(queue_row, "label_axis"),
      metric_id = .mv05s_queue_value(queue_row, "metric_id"),
      seed = seed_id,
      estimate = metric$estimate,
      samples = nrow(seed_part),
      label_classes = length(unique(label_values)),
      clusters = length(unique(seed_part$cluster)),
      status = metric$status,
      p_value_computed = FALSE,
      method_selection_executed = FALSE,
      stringsAsFactors = FALSE)
  })
  seed_rows <- do.call(rbind, rows)
  values <- seed_rows$estimate[seed_rows$status == "completed"]
  summary <- data.frame(
    evaluation_unit_id = .mv05s_queue_value(queue_row, "evaluation_unit_id"),
    analysis_group_id = group_id,
    fold_id = .mv05s_queue_value(queue_row, "fold_id"),
    representation = .mv05s_queue_value(queue_row, "representation"),
    distance_id = .mv05s_queue_value(queue_row, "distance_id"),
    algorithm_id = algorithm_id,
    algorithm_role = .mv05s_queue_value(queue_row, "algorithm_role"),
    endpoint_id = .mv05s_queue_value(queue_row, "endpoint_id"),
    label_axis = .mv05s_queue_value(queue_row, "label_axis"),
    metric_id = .mv05s_queue_value(queue_row, "metric_id"),
    seed_mean = if (length(values) == 5L) mean(values) else NA_real_,
    seed_jackknife_se = if (length(values) == 5L)
      mv05s_jackknife_se_v1(values) else NA_real_,
    completed_seeds = length(values),
    expected_seeds = 5L,
    samples_per_seed = expected_n,
    status = if (length(values) == 5L) "completed" else
      "degenerate_label_or_partition_metric_not_identifiable",
    p_value_computed = FALSE,
    method_selection_executed = FALSE,
    stringsAsFactors = FALSE)
  list(seed = seed_rows, summary = summary)
}

mv05s_balanced_accuracy_v1 <- function(correct, labels) {
  if (length(correct) != length(labels) || !length(correct) ||
      anyNA(correct) || anyNA(labels)) {
    stop("MV5-S balanced-accuracy inputs are malformed.", call. = FALSE)
  }
  class_values <- tapply(as.numeric(correct), labels, mean)
  list(estimate = mean(class_values), class_values = class_values)
}

mv05s_evaluate_heldout_unit_v1 <- function(queue_row, selected, heldout, labels) {
  if (.mv05s_queue_value(queue_row, "evaluation_scope") !=
      "heldout_label_prediction_from_frozen_training_cluster") {
    stop("MV5-S held-out evaluator received a training unit.", call. = FALSE)
  }
  group_id <- .mv05s_queue_value(queue_row, "analysis_group_id")
  algorithm_id <- .mv05s_queue_value(queue_row, "algorithm_id")
  label_axis <- .mv05s_queue_value(queue_row, "label_axis")
  expected_n <- as.integer(.mv05s_queue_value(queue_row, "training_samples"))
  training <- selected[selected$analysis_group_id == group_id &
                         selected$algorithm_id == algorithm_id, , drop = FALSE]
  query <- heldout[heldout$analysis_group_id == group_id &
                     heldout$algorithm_id == algorithm_id, , drop = FALSE]
  seeds <- sort(unique(as.integer(training$seed)))
  if (!identical(seeds, 20260805:20260809) ||
      !identical(seeds, sort(unique(as.integer(query$seed))))) {
    stop("MV5-S held-out unit has an incomplete seed axis.", call. = FALSE)
  }
  private_rows <- lapply(seeds, function(seed_id) {
    training_seed <- training[as.integer(training$seed) == seed_id, , drop = FALSE]
    query_seed <- query[as.integer(query$seed) == seed_id, , drop = FALSE]
    if (nrow(training_seed) != expected_n || anyDuplicated(training_seed$sample_id) ||
        !nrow(query_seed) || anyDuplicated(query_seed$query_sample_id)) {
      stop("MV5-S held-out training/query axis is malformed.", call. = FALSE)
    }
    training_labels <- data.frame(
      sample_id = training_seed$sample_id,
      label_value = .mv05s_join_label(training_seed$sample_id, labels, label_axis),
      stringsAsFactors = FALSE)
    mapping <- mv05r_plurality_map_v1(
      training_seed[c("sample_id", "cluster")], training_labels, "label_value")
    map_index <- match(as.character(query_seed$cluster), as.character(mapping$cluster))
    if (anyNA(map_index)) stop("MV5-S held-out cluster mapping is incomplete.",
                               call. = FALSE)
    label_index <- match(query_seed$query_sample_id, labels$sample_id)
    if (anyNA(label_index)) stop("MV5-S held-out label join is incomplete.",
                                 call. = FALSE)
    true_label <- labels[[label_axis]][label_index]
    predicted_label <- mapping$predicted_label[map_index]
    data.frame(
      evaluation_unit_id = .mv05s_queue_value(queue_row, "evaluation_unit_id"),
      analysis_group_id = group_id,
      fold_id = .mv05s_queue_value(queue_row, "fold_id"),
      representation = .mv05s_queue_value(queue_row, "representation"),
      distance_id = .mv05s_queue_value(queue_row, "distance_id"),
      algorithm_id = algorithm_id,
      algorithm_role = .mv05s_queue_value(queue_row, "algorithm_role"),
      endpoint_id = .mv05s_queue_value(queue_row, "endpoint_id"),
      label_axis = label_axis,
      seed = seed_id,
      query_sample_id = query_seed$query_sample_id,
      heldout_study = labels$study[label_index],
      cluster = query_seed$cluster,
      true_label = true_label,
      predicted_label = predicted_label,
      correct = true_label == predicted_label,
      plurality_count = mapping$plurality_count[map_index],
      plurality_tie_size = mapping$plurality_tie_size[map_index],
      stringsAsFactors = FALSE)
  })
  private <- do.call(rbind, private_rows)
  public_seed <- do.call(rbind, lapply(split(private, private$seed), function(part) {
    score <- mv05s_balanced_accuracy_v1(part$correct, part$true_label)
    data.frame(
      evaluation_unit_id = part$evaluation_unit_id[[1L]],
      analysis_group_id = part$analysis_group_id[[1L]],
      fold_id = part$fold_id[[1L]],
      representation = part$representation[[1L]],
      distance_id = part$distance_id[[1L]],
      algorithm_id = part$algorithm_id[[1L]],
      algorithm_role = part$algorithm_role[[1L]],
      endpoint_id = part$endpoint_id[[1L]],
      label_axis = part$label_axis[[1L]],
      seed = part$seed[[1L]],
      estimate = score$estimate,
      heldout_samples = nrow(part),
      heldout_label_classes = length(score$class_values),
      status = "completed",
      p_value_computed = FALSE,
      method_selection_executed = FALSE,
      stringsAsFactors = FALSE)
  }))
  values <- public_seed$estimate
  summary <- data.frame(
    evaluation_unit_id = .mv05s_queue_value(queue_row, "evaluation_unit_id"),
    analysis_group_id = group_id,
    fold_id = .mv05s_queue_value(queue_row, "fold_id"),
    representation = .mv05s_queue_value(queue_row, "representation"),
    distance_id = .mv05s_queue_value(queue_row, "distance_id"),
    algorithm_id = algorithm_id,
    algorithm_role = .mv05s_queue_value(queue_row, "algorithm_role"),
    endpoint_id = .mv05s_queue_value(queue_row, "endpoint_id"),
    label_axis = label_axis,
    seed_mean = mean(values),
    seed_jackknife_se = mv05s_jackknife_se_v1(values),
    completed_seeds = length(values),
    expected_seeds = 5L,
    heldout_samples_per_seed = nrow(private) / length(values),
    status = "completed",
    p_value_computed = FALSE,
    method_selection_executed = FALSE,
    stringsAsFactors = FALSE)
  list(seed = public_seed, summary = summary, private = private)
}

.mv05s_with_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(expression)
}

mv05s_bootstrap_counts_v1 <- function(labels, label_axis,
                                       replicates = 2000L,
                                       seed = 20260810L) {
  labels <- unique(labels[c("sample_id", "study", label_axis)])
  studies <- sort(unique(labels$study), method = "radix")
  if (length(studies) != 15L || replicates < 1L || anyNA(labels)) {
    stop("MV5-S bootstrap label design is invalid.", call. = FALSE)
  }
  counts <- matrix(0L, nrow = replicates, ncol = length(studies),
                   dimnames = list(NULL, studies))
  .mv05s_with_seed(seed, {
    if (label_axis == "tissue") {
      study_label <- unique(labels[c("study", "tissue")])
      if (anyDuplicated(study_label$study)) {
        stop("MV5-S tissue bootstrap requires tissue-homogeneous studies.",
             call. = FALSE)
      }
      strata <- split(study_label$study, study_label$tissue)
      if (length(strata) != 5L || any(lengths(strata) < 2L)) {
        stop("MV5-S tissue bootstrap strata lack study support.", call. = FALSE)
      }
      for (replicate_id in seq_len(replicates)) {
        sampled <- unlist(lapply(strata, function(ids) {
          sample(ids, length(ids), replace = TRUE)
        }), use.names = FALSE)
        counts[replicate_id, ] <- tabulate(match(sampled, studies),
                                           nbins = length(studies))
      }
    } else if (label_axis == "approach") {
      approaches_per_study <- tapply(labels$approach, labels$study,
                                     function(x) length(unique(x)))
      if (sum(approaches_per_study > 1L) != 3L) {
        stop("MV5-S approach bootstrap mixed-study design drifted.",
             call. = FALSE)
      }
      for (replicate_id in seq_len(replicates)) {
        sampled <- sample(studies, length(studies), replace = TRUE)
        counts[replicate_id, ] <- tabulate(match(sampled, studies),
                                           nbins = length(studies))
      }
    } else {
      stop("MV5-S bootstrap supports only tissue or approach.", call. = FALSE)
    }
  })
  counts
}

mv05s_apply_bootstrap_v1 <- function(sample_summary, label_axis, counts,
                                      minimum_estimable_fraction = 0.95) {
  required <- c("sample_id", "study", label_axis, "seed_mean_correct")
  if (!all(required %in% names(sample_summary)) || anyNA(sample_summary) ||
      anyDuplicated(sample_summary$sample_id) ||
      !all(sample_summary$study %in% colnames(counts))) {
    stop("MV5-S bootstrap inputs are malformed.", call. = FALSE)
  }
  classes <- sort(unique(sample_summary[[label_axis]]), method = "radix")
  estimates <- vapply(seq_len(nrow(counts)), function(index) {
    weights <- counts[index, match(sample_summary$study, colnames(counts))]
    class_estimates <- vapply(classes, function(class_id) {
      keep <- sample_summary[[label_axis]] == class_id
      if (sum(weights[keep]) == 0) return(NA_real_)
      stats::weighted.mean(sample_summary$seed_mean_correct[keep], weights[keep])
    }, numeric(1L))
    if (anyNA(class_estimates)) NA_real_ else mean(class_estimates)
  }, numeric(1L))
  estimable <- estimates[is.finite(estimates)]
  enough <- length(estimable) >= ceiling(nrow(counts) * minimum_estimable_fraction)
  interval <- if (enough) stats::quantile(
    estimable, c(0.025, 0.975), names = FALSE, type = 7) else c(NA_real_, NA_real_)
  list(
    estimates = estimates,
    ci_lower = interval[[1L]],
    ci_upper = interval[[2L]],
    estimable_replicates = length(estimable),
    requested_replicates = nrow(counts),
    status = if (enough) "completed" else "bootstrap_support_insufficient")
}
