mv10p_method_registry_v1 <- function() {
  data.frame(
    method_order = 1:5,
    method_id = c("pam_dissimilarity_v1", "hclust_average_v1",
                  "hclust_complete_v1", "diana_dissimilarity_v1",
                  "hclust_single_diagnostic_v1"),
    method_role = c("primary", "sensitivity", "sensitivity", "sensitivity",
                    "diagnostic"),
    stringsAsFactors = FALSE
  )
}

mv10p_select_outcome_partitions_v1 <- function(partitions) {
  required <- c("stack_id", "seed", "homology_dimension", "sample_id",
                "cluster", "k", "method_id", "outcome_label_state",
                "biological_outcomes_computed")
  if (!is.data.frame(partitions) || nrow(partitions) != 167400L ||
      !all(required %in% names(partitions)) || anyDuplicated(partitions[c(
        "stack_id", "seed", "homology_dimension", "method_id", "k",
        "sample_id")]) || any(partitions$outcome_label_state != "closed") ||
      any(as.logical(partitions$biological_outcomes_computed))) {
    stop("MV10-P partition source is malformed or not label closed.",
         call. = FALSE)
  }
  selected_k <- c(H0 = 2L, H1 = 3L)
  selected <- partitions[
    partitions$k == selected_k[partitions$homology_dimension], , drop = FALSE
  ]
  registry <- mv10p_method_registry_v1()
  expected_stacks <- c("allqc_data_exact500", "allqc_residual_exact500",
                       "existing_selectedfit_data_exact500")
  expected_seeds <- 20260805:20260809
  if (nrow(selected) != 18600L ||
      !setequal(unique(selected$stack_id), expected_stacks) ||
      !setequal(unique(selected$method_id), registry$method_id) ||
      !identical(sort(unique(selected$seed)), expected_seeds) ||
      !identical(sort(unique(selected$homology_dimension)), c("H0", "H1")) ||
      any(table(selected$stack_id, selected$homology_dimension,
                selected$method_id, selected$seed) != 124L)) {
    stop("MV10-P selected partition axis drifted.", call. = FALSE)
  }
  selected[order(match(selected$stack_id, expected_stacks),
                 selected$homology_dimension,
                 match(selected$method_id, registry$method_id),
                 selected$seed, selected$sample_id), , drop = FALSE]
}

mv10p_build_queue_v1 <- function(partition_registry, endpoints, metrics) {
  scheduled <- endpoints[endpoints$execution_status == "scheduled", ,
                         drop = FALSE]
  queue <- merge(
    scheduled[c("endpoint_order", "endpoint_id", "population_id",
                "label_axis")],
    partition_registry[c("partition_order", "stack_id",
                         "homology_dimension", "method_id", "method_role",
                         "selected_k", "partition_group_sha256")], all = TRUE
  )
  queue <- merge(queue, metrics[c("metric_order", "metric_id")], all = TRUE)
  queue <- queue[order(queue$endpoint_order, queue$partition_order,
                       queue$metric_order), , drop = FALSE]
  queue$contract_id <- "mv10p_outcome_execution_queue_v1"
  queue$execution_order <- seq_len(nrow(queue))
  queue$evaluation_unit_id <- vapply(seq_len(nrow(queue)), function(index) {
    paste0("mv10p_outcome_v1:", digest::digest(list(
      endpoint_id = queue$endpoint_id[[index]],
      stack_id = queue$stack_id[[index]],
      homology_dimension = queue$homology_dimension[[index]],
      method_id = queue$method_id[[index]],
      metric_id = queue$metric_id[[index]],
      partition_group_sha256 = queue$partition_group_sha256[[index]]
    ), algo = "sha256", serialize = TRUE))
  }, character(1L))
  queue$expected_seeds <- 5L
  queue$expected_samples_per_seed <- ifelse(
    queue$population_id == "full124_descriptive", 124L, 90L
  )
  queue$cluster_metadata_join_executed <- FALSE
  queue$association_computed <- FALSE
  queue$method_selection_executed <- FALSE
  queue[c("contract_id", "evaluation_unit_id", "execution_order",
          "endpoint_id", "population_id", "label_axis", "stack_id",
          "homology_dimension", "method_id", "method_role", "selected_k",
          "metric_id", "partition_group_sha256", "expected_seeds",
          "expected_samples_per_seed", "cluster_metadata_join_executed",
          "association_computed", "method_selection_executed")]
}

mv10p_evaluate_outcomes_v1 <- function(selected, metadata, queue) {
  required_metadata <- c("sample_id", "tissue", "study",
                         "canonical_approach", "corrected_primary_90")
  if (!is.data.frame(metadata) || nrow(metadata) != 124L ||
      !all(required_metadata %in% names(metadata)) ||
      anyDuplicated(metadata$sample_id) || anyNA(metadata[required_metadata]) ||
      !setequal(selected$sample_id, metadata$sample_id) || nrow(queue) != 300L) {
    stop("MV10-Q metadata or queue axis is malformed.", call. = FALSE)
  }
  seed_rows <- vector("list", nrow(queue) * 5L)
  summary_rows <- vector("list", nrow(queue))
  contingency_rows <- vector("list", nrow(queue) * 5L / 2L)
  seed_cursor <- 0L; contingency_cursor <- 0L
  seeds <- 20260805:20260809
  for (unit_index in seq_len(nrow(queue))) {
    unit <- queue[unit_index, , drop = FALSE]
    population_ids <- if (unit$population_id == "full124_descriptive") {
      metadata$sample_id
    } else metadata$sample_id[as.logical(metadata$corrected_primary_90)]
    metric_id <- if (unit$metric_id == "normalized_mutual_information_max") {
      "normalized_mutual_information"
    } else unit$metric_id
    estimates <- numeric(5L)
    for (seed_index in seq_along(seeds)) {
      seed <- seeds[[seed_index]]
      partition <- selected[
        selected$stack_id == unit$stack_id &
          selected$homology_dimension == unit$homology_dimension &
          selected$method_id == unit$method_id & selected$seed == seed,
        c("sample_id", "cluster"), drop = FALSE
      ]
      partition <- partition[partition$sample_id %in% population_ids, ,
                             drop = FALSE]
      partition <- partition[order(partition$sample_id), , drop = FALSE]
      metadata_index <- match(partition$sample_id, metadata$sample_id)
      labels <- metadata[[unit$label_axis]][metadata_index]
      if (nrow(partition) != unit$expected_samples_per_seed ||
          anyNA(metadata_index) || anyNA(labels) ||
          anyDuplicated(partition$sample_id)) {
        stop("MV10-Q outcome unit sample/label axis is malformed.",
             call. = FALSE)
      }
      metric <- mv05s_metric_v1(metric_id, partition$cluster, labels)
      estimates[[seed_index]] <- metric$estimate
      seed_cursor <- seed_cursor + 1L
      seed_rows[[seed_cursor]] <- data.frame(
        contract_id = "mv10q_outcome_seed_metric_v1",
        evaluation_unit_id = unit$evaluation_unit_id,
        execution_order = unit$execution_order, endpoint_id = unit$endpoint_id,
        population_id = unit$population_id, label_axis = unit$label_axis,
        stack_id = unit$stack_id,
        homology_dimension = unit$homology_dimension,
        method_id = unit$method_id, method_role = unit$method_role,
        selected_k = unit$selected_k, metric_id = unit$metric_id, seed = seed,
        estimate = metric$estimate, samples = nrow(partition),
        label_classes = length(unique(labels)),
        partition_clusters = length(unique(partition$cluster)),
        status = metric$status, p_value_computed = FALSE,
        method_selection_executed = FALSE, stringsAsFactors = FALSE
      )
      if (unit$metric_id == "adjusted_rand_index") {
        cells <- as.data.frame(table(cluster = partition$cluster,
                                     label_value = labels),
                              stringsAsFactors = FALSE)
        contingency_cursor <- contingency_cursor + 1L
        contingency_rows[[contingency_cursor]] <- data.frame(
          contract_id = "mv10q_outcome_private_contingency_v1",
          endpoint_id = unit$endpoint_id, population_id = unit$population_id,
          label_axis = unit$label_axis, stack_id = unit$stack_id,
          homology_dimension = unit$homology_dimension,
          method_id = unit$method_id, selected_k = unit$selected_k, seed = seed,
          cluster = cells$cluster, label_value = cells$label_value,
          samples = cells$Freq, method_selection_executed = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
    summary_rows[[unit_index]] <- data.frame(
      contract_id = "mv10q_outcome_unit_summary_v1",
      evaluation_unit_id = unit$evaluation_unit_id,
      execution_order = unit$execution_order, endpoint_id = unit$endpoint_id,
      population_id = unit$population_id, label_axis = unit$label_axis,
      stack_id = unit$stack_id,
      homology_dimension = unit$homology_dimension,
      method_id = unit$method_id, method_role = unit$method_role,
      selected_k = unit$selected_k, metric_id = unit$metric_id,
      seed_mean = mean(estimates), seed_median = stats::median(estimates),
      seed_minimum = min(estimates), seed_maximum = max(estimates),
      seed_jackknife_se = mv05s_jackknife_se_v1(estimates),
      completed_seeds = sum(is.finite(estimates)), expected_seeds = 5L,
      status = if (all(is.finite(estimates))) "completed" else
        "degenerate_label_or_partition_metric_not_identifiable",
      p_value_computed = FALSE, method_selection_executed = FALSE,
      stringsAsFactors = FALSE
    )
  }
  seed_metrics <- do.call(rbind, seed_rows)
  unit_summaries <- do.call(rbind, summary_rows)
  contingency <- do.call(rbind, contingency_rows[seq_len(contingency_cursor)])
  if (nrow(seed_metrics) != 1500L || nrow(unit_summaries) != 300L ||
      any(!is.finite(seed_metrics$estimate)) ||
      any(seed_metrics$status != "completed") ||
      any(unit_summaries$status != "completed") ||
      any(seed_metrics$p_value_computed) ||
      any(seed_metrics$method_selection_executed)) {
    stop("MV10-Q outcome result contract failed.", call. = FALSE)
  }
  list(seed_metrics = seed_metrics, unit_summaries = unit_summaries,
       contingency = contingency)
}
