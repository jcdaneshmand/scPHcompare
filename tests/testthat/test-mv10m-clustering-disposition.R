test_that("MV10-M disposition applies only the frozen methodological rules", {
  dimensions <- c("H0", "H1")
  stacks <- c("allqc_data_exact500", "allqc_residual_exact500",
              "existing_selectedfit_data_exact500")
  methods <- c("average", "complete", "diana", "pam_dissimilarity_v1", "single")
  stability <- expand.grid(
    stack_id = stacks, homology_dimension = dimensions,
    method_id = methods, k = 2:10, stringsAsFactors = FALSE
  )
  stability$contract_id <- "fixture"
  stability$mean_adjusted_rand <- 0.8
  stability$minimum_adjusted_rand <- 0.7
  stability$maximum_adjusted_rand <- 0.9
  stability$representation_label <- stability$stack_id
  stability$method_label <- stability$method_id
  quality <- stability
  quality$median_mean_silhouette <- 0.4
  agreement <- expand.grid(
    stack_id = stacks, homology_dimension = dimensions, k = 2:10,
    method_pair_id = paste0("pam_pair_", 1:10), stringsAsFactors = FALSE
  )
  agreement$contract_id <- "fixture"
  agreement$method_pair_label <- rep(c(
    paste0("PAM (primary) vs method ", 1:4), paste0("other pair ", 1:6)
  ), length.out = nrow(agreement))
  agreement$median_adjusted_rand <- 0.5
  selection <- data.frame(
    homology_dimension = dimensions, method_id = "pam_dissimilarity_v1",
    selected_k = c(2L, 3L), stringsAsFactors = FALSE
  )
  primary_stability <- expand.grid(
    homology_dimension = dimensions, k = 2:10, stringsAsFactors = FALSE
  )
  primary_stability$method_id <- "pam_dissimilarity_v1"
  primary_stability$mean_adjusted_rand <- 0.8
  primary_stability$minimum_adjusted_rand <- 0.7
  primary_stability$maximum_adjusted_rand <- 0.9
  primary_stability$threshold <- 0.75
  primary_quality <- expand.grid(
    homology_dimension = dimensions, seed = 1:5, k = 2:10,
    stringsAsFactors = FALSE
  )
  primary_quality$method_id <- "pam_dissimilarity_v1"
  primary_quality$mean_silhouette <- 0.4
  primary_quality$negative_silhouette_fraction <- 0.1
  primary_quality$minimum_cluster_size <- 3L
  primary_quality$singleton_clusters <- 0L
  result <- mv10m_build_methodological_disposition_v1(
    stability, quality, agreement, selection, primary_stability, primary_quality
  )
  expect_equal(nrow(result$selected_primary_seed_metrics), 10L)
  expect_equal(nrow(result$primary_summary), 2L)
  expect_equal(nrow(result$representation_context), 6L)
  expect_equal(nrow(result$method_sensitivity), 24L)
  expect_true(all(result$primary_summary$structurally_nondegenerate))
  expect_equal(result$disposition$selected_H0_k, 2L)
  expect_equal(result$disposition$selected_H1_k, 3L)
  expect_true(result$disposition$H0_H1_remain_separate)
  expect_true(result$disposition$all_selected_partitions_structurally_nondegenerate)
  expect_false(result$disposition$labels_used)
  expect_false(result$disposition$outcomes_used)
  expect_false(result$disposition$biological_interpretation)
  expect_false(result$disposition$manuscript_claims)

  primary_quality$singleton_clusters[
    primary_quality$homology_dimension == "H1" &
      primary_quality$seed == 1L & primary_quality$k == 3L
  ] <- 1L
  failed <- mv10m_build_methodological_disposition_v1(
    stability, quality, agreement, selection, primary_stability, primary_quality
  )
  expect_false(failed$disposition$all_selected_partitions_structurally_nondegenerate)
  expect_equal(failed$disposition$decision,
               "do_not_open_labels_due_to_structural_degeneracy")
})
