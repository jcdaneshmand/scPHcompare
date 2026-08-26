test_that("MV10-P outcome helpers preserve fixed axes and complete reporting", {
  methods <- mv10p_method_registry_v1()
  expect_equal(nrow(methods), 5L)
  expect_identical(methods$method_role,
                   c("primary", "sensitivity", "sensitivity", "sensitivity",
                     "diagnostic"))
  registry <- expand.grid(
    stack_id = c("a", "b", "c"), homology_dimension = c("H0", "H1"),
    method_id = methods$method_id, stringsAsFactors = FALSE
  )
  registry$partition_order <- seq_len(nrow(registry))
  registry$method_role <- methods$method_role[
    match(registry$method_id, methods$method_id)
  ]
  registry$selected_k <- ifelse(registry$homology_dimension == "H0", 2L, 3L)
  registry$partition_group_sha256 <- sprintf("%064d", registry$partition_order)
  endpoints <- expand.grid(
    population_id = c("full124_descriptive", "primary90_context_restriction"),
    label_axis = c("tissue", "study", "canonical_approach"),
    stringsAsFactors = FALSE
  )
  endpoints$endpoint_order <- seq_len(nrow(endpoints))
  endpoints$endpoint_id <- paste(endpoints$population_id, endpoints$label_axis,
                                 sep = "__")
  endpoints$execution_status <- "scheduled"
  endpoints$execution_status[
    endpoints$population_id == "primary90_context_restriction" &
      endpoints$label_axis == "canonical_approach"
  ] <- "structurally_not_estimable_single_class"
  metrics <- data.frame(metric_order = 1:2,
                        metric_id = c("adjusted_rand_index",
                                      "normalized_mutual_information_max"))
  queue <- mv10p_build_queue_v1(registry, endpoints, metrics)
  expect_equal(nrow(queue), 300L)
  expect_equal(anyDuplicated(queue$evaluation_unit_id), 0L)
  expect_equal(length(unique(queue$stack_id)), 3L)
  expect_equal(length(unique(queue$homology_dimension)), 2L)
  expect_equal(length(unique(queue$method_id)), 5L)
  expect_equal(length(unique(queue$endpoint_id)), 5L)
  expect_equal(length(unique(queue$metric_id)), 2L)
  expect_true(all(queue$expected_seeds == 5L))
  expect_false(any(queue$association_computed))
  expect_false(any(queue$method_selection_executed))

  samples <- sprintf("sample_%03d", 1:124)
  selected <- expand.grid(
    stack_id = c("a", "b", "c"), homology_dimension = c("H0", "H1"),
    method_id = methods$method_id, seed = 20260805:20260809,
    sample_id = samples, stringsAsFactors = FALSE
  )
  selected$k <- ifelse(selected$homology_dimension == "H0", 2L, 3L)
  selected$cluster <- 1L + (match(selected$sample_id, samples) +
                             match(selected$method_id, methods$method_id) +
                             selected$seed) %% selected$k
  metadata <- data.frame(
    sample_id = samples,
    tissue = rep(paste0("tissue_", 1:8), length.out = 124L),
    study = rep(paste0("study_", 1:18), length.out = 124L),
    canonical_approach = rep(c("scRNA-seq", "snRNA-seq"),
                             length.out = 124L),
    corrected_primary_90 = seq_len(124L) <= 90L,
    stringsAsFactors = FALSE
  )
  evaluated <- mv10p_evaluate_outcomes_v1(selected, metadata, queue)
  expect_equal(nrow(evaluated$seed_metrics), 1500L)
  expect_equal(nrow(evaluated$unit_summaries), 300L)
  expect_true(nrow(evaluated$contingency) > 0L)
  expect_true(all(is.finite(evaluated$seed_metrics$estimate)))
  expect_true(all(evaluated$seed_metrics$status == "completed"))
  expect_false(any(evaluated$seed_metrics$p_value_computed))
  expect_false(any(evaluated$seed_metrics$method_selection_executed))
})
