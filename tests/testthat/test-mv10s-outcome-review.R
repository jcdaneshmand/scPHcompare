test_that("MV10-S builds complete outcome-review tables without ranking", {
  stacks <- names(.mv10s_stack_labels)
  dimensions <- c("H0", "H1")
  methods <- names(.mv10s_method_labels)
  endpoints <- names(.mv10s_endpoint_labels)[1:5]
  metrics <- names(.mv10s_metric_labels)
  grid <- expand.grid(
    endpoint_id = endpoints, stack_id = stacks,
    homology_dimension = dimensions, method_id = methods,
    metric_id = metrics, stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$endpoint_id, grid$stack_id,
                     grid$homology_dimension, grid$method_id,
                     grid$metric_id), ]
  grid$execution_order <- seq_len(nrow(grid))
  grid$evaluation_unit_id <- paste0("unit-", grid$execution_order)
  grid$population_id <- ifelse(grepl("^full124", grid$endpoint_id),
                               "full124_descriptive",
                               "primary90_context_restriction")
  grid$label_axis <- sub(".*__", "", grid$endpoint_id)
  grid$method_role <- ifelse(grid$method_id == "pam_dissimilarity_v1",
                             "primary", "sensitivity")
  grid$selected_k <- ifelse(grid$homology_dimension == "H0", 2L, 3L)
  seed <- grid[rep(seq_len(nrow(grid)), each = 5L), ]
  seed$seed <- rep(20260805:20260809, nrow(grid))
  seed$estimate <- ((seed$execution_order %% 17L) + seed$seed %% 5L) / 25
  seed$status <- "completed"
  split_seed <- split(seed, seed$evaluation_unit_id)
  summary <- do.call(rbind, lapply(split_seed, function(x) data.frame(
    x[1, c("evaluation_unit_id", "execution_order", "endpoint_id",
           "population_id", "label_axis", "stack_id", "homology_dimension",
           "method_id", "method_role", "selected_k", "metric_id")],
    seed_mean = mean(x$estimate), seed_median = median(x$estimate),
    seed_minimum = min(x$estimate), seed_maximum = max(x$estimate),
    seed_jackknife_se = .mv10s_jackknife_se(x$estimate),
    completed_seeds = 5L, expected_seeds = 5L, status = "completed"
  )))
  structural <- data.frame(
    endpoint_id = names(.mv10s_endpoint_labels)[[6L]],
    population_id = "primary90_context_restriction",
    label_axis = "canonical_approach",
    status = "structurally_not_estimable_single_class"
  )
  result <- mv10s_build_outcome_review_v1(seed, summary, structural)
  expect_equal(nrow(result$complete_summary), 300L)
  expect_equal(nrow(result$primary_summary), 60L)
  expect_equal(nrow(result$endpoint_coverage), 6L)
  expect_equal(length(unique(result$complete_summary$stack_id)), 3L)
  expect_equal(length(unique(result$complete_summary$method_id)), 5L)
  expect_equal(length(unique(result$complete_summary$endpoint_id)), 5L)
  expect_equal(length(unique(result$complete_summary$metric_id)), 2L)
  expect_true(all(!result$complete_summary$inference_performed))
  expect_true(all(!result$complete_summary$ranking_performed))
  expect_true(all(!result$complete_summary$biological_claim))
  expect_equal(tail(result$endpoint_coverage$execution_status, 1L),
               "structurally_not_estimable_single_class")
})
