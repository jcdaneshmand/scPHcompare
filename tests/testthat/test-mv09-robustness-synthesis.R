mv09_fixture <- function() {
  internal <- expand.grid(
    dataset_scope = "internal124",
    contrast_id = c("fit_scope_effect", "layer_effect", "net_migration_effect"),
    seed = 20260805:20260809,
    homology_dimension = c("H0", "H1"), stringsAsFactors = FALSE
  )
  external <- expand.grid(
    dataset_scope = "external8",
    contrast_id = c("fit_scope_effect", "layer_effect",
                    "residual_panel_effect", "common475_net_migration_effect",
                    "exact500_total_migration_effect"),
    seed = 20260805L, homology_dimension = c("H0", "H1"),
    stringsAsFactors = FALSE
  )
  value <- rbind(internal, external)
  value$comparison_order <- seq_len(nrow(value))
  value$left_stack <- "left"
  value$right_stack <- "right"
  value$units <- ifelse(value$dataset_scope == "internal124", 124L, 8L)
  value$unordered_pairs <- ifelse(value$units == 124L, 7626L, 28L)
  value$distance_transform <- "sqrt_exact_squared_L2"
  for (i in seq_along(.mv09_metrics_v1())) {
    value[[.mv09_metrics_v1()[[i]]]] <-
      seq_len(nrow(value)) / (100 + i)
  }
  value$interpretation <- "descriptive_no_equivalence_or_biological_claim"
  value$outcome_label_state <- "closed"
  value$biological_outcomes_computed <- FALSE
  value$clustering_jobs <- 0L
  value$fusion_jobs <- 0L
  value$label_jobs <- 0L
  value$outcome_jobs <- 0L
  value
}

test_that("MV9 builds deterministic descriptive robustness strata", {
  result <- mv09_build_robustness_synthesis_v1(mv09_fixture())
  expect_equal(nrow(result$aggregate), 40L)
  expect_equal(nrow(result$plot_data), 440L)
  expect_equal(nrow(result$internal_seed_summary), 66L)
  expect_true(all(result$internal_seed_summary$seeds == 5L))
  expect_equal(nrow(result$external_singleton), 110L)
  expect_true(all(result$external_singleton$replication_units == 1L))
  expect_equal(nrow(result$dimension_delta), 220L)
  expect_equal(nrow(result$dimension_delta_summary), 88L)
  expect_true(all(result$dimension_delta_summary$inference ==
                    "none_descriptive_only"))
  expect_false(any(grepl("label|outcome", names(result$plot_data))))
})

test_that("MV9 fails closed on nonfinite metrics and open firewalls", {
  value <- mv09_fixture()
  value$pearson[[1L]] <- NA_real_
  expect_error(mv09_build_robustness_synthesis_v1(value), "nonfinite")
  value <- mv09_fixture()
  value$label_jobs[[1L]] <- 1L
  expect_error(mv09_build_robustness_synthesis_v1(value), "invalid")
})
