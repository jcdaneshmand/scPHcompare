test_that("MV5-AC freezes complete method, endpoint, and estimand panels", {
  methods <- mv05ac_method_map_v1()
  endpoints <- mv05ac_endpoint_registry_v1()
  estimands <- mv05ac_estimand_registry_v1(methods, endpoints)

  expect_equal(nrow(methods), 8L)
  expect_equal(as.integer(table(methods$representation)), c(4L, 4L))
  expect_equal(names(table(methods$representation)),
               c("inductive_integrated", "sct_whole"))
  expect_true(all(methods$baseline_coordinate_count == 30L))
  expect_true(all(methods$cosine_coordinate_count == 30L))
  expect_true(all(methods$baseline_point_metric == "euclidean"))
  expect_true(all(methods$cosine_point_metric ==
                    "euclidean_chord_after_row_unit_normalization"))
  expect_equal(nrow(endpoints), 2L)
  expect_equal(nrow(estimands), 24L)
  expect_equal(sum(estimands$estimand_role ==
                     "confirmatory_cosine_sensitivity"), 4L)
  expect_false(any(estimands$equivalence_or_noninferiority_claim_authorized))
  expect_silent(mv05ac_validate_estimand_registry_v1(estimands))
})

test_that("MV5-AC rejects pre-outcome result leakage", {
  expect_error(.mv05ac_assert_preoutcome(
    data.frame(reciprocal_rank = 1)), "result columns")
  expect_error(.mv05ac_assert_preoutcome(
    data.frame(outcomes_computed = TRUE)), "prohibited completed")
  expect_silent(.mv05ac_assert_preoutcome(
    data.frame(outcomes_computed = FALSE, execution_authorized = FALSE)))
})

test_that("MV5-AC queue is exact and remains unauthorized", {
  folds <- paste0("large_loso_v1:S", sprintf("%02d", 1:15))
  scope <- expand.grid(
    fold_id = folds, seed = 20260805:20260809,
    representation = c("sct_whole", "inductive_integrated"),
    stringsAsFactors = FALSE)
  study_sizes <- c(5L, 2L, 2L, 10L, 6L, 4L, 1L, 5L, 9L, 4L, 4L,
                   1L, 25L, 4L, 8L)
  heldout <- rep(rep(study_sizes, times = 5L), 2L)
  scope$robustness_group_id <- paste0("group_", seq_len(nrow(scope)))
  scope$configuration_id <- "cells384_pc30_cosine_chord_v1"
  scope$query_samples <- heldout
  scope$training_samples <- 90L - heldout
  scope$biological_pairs <- scope$query_samples * scope$training_samples
  scope$method_rows <- scope$biological_pairs * 4L
  scope$baseline_group_id <- paste0("baseline_", seq_len(nrow(scope)))
  scope$baseline_group_sha256 <- paste(rep("a", 64L), collapse = "")
  scope$cosine_group_manifest_sha256 <- paste(rep("b", 64L), collapse = "")
  scope$private_result_locator <- paste0("ignored_tmp:", seq_len(nrow(scope)))

  # Adjust the synthetic held-out counts to the accepted 35,350 pairs per view.
  delta <- 35350L - sum(scope$biological_pairs[scope$representation == "sct_whole"])
  expect_equal(delta, 0L)
  queue <- mv05ac_build_evaluation_queue_v1(
    scope, paste(rep("c", 64L), collapse = ""))
  expect_equal(nrow(queue), 150L)
  expect_equal(sum(queue$expected_query_endpoint_rows), 7200L)
  expect_false(any(queue$execution_authorized))
  expect_silent(mv05ac_validate_evaluation_queue_v1(queue))
})

test_that("MV5-AC queue validator detects authorization and axis drift", {
  folds <- paste0("large_loso_v1:S", sprintf("%02d", 1:15))
  scope <- expand.grid(
    fold_id = folds, seed = 20260805:20260809,
    representation = c("sct_whole", "inductive_integrated"),
    stringsAsFactors = FALSE)
  study_sizes <- c(5L, 2L, 2L, 10L, 6L, 4L, 1L, 5L, 9L, 4L, 4L,
                   1L, 25L, 4L, 8L)
  heldout <- rep(study_sizes, times = 5L)
  scope$robustness_group_id <- paste0("group_", seq_len(nrow(scope)))
  scope$configuration_id <- "cells384_pc30_cosine_chord_v1"
  scope$query_samples <- rep(heldout, 2L)
  scope$training_samples <- 90L - scope$query_samples
  scope$biological_pairs <- scope$query_samples * scope$training_samples
  scope$method_rows <- scope$biological_pairs * 4L
  scope$baseline_group_id <- paste0("baseline_", seq_len(nrow(scope)))
  scope$baseline_group_sha256 <- paste(rep("a", 64L), collapse = "")
  scope$cosine_group_manifest_sha256 <- paste(rep("b", 64L), collapse = "")
  scope$private_result_locator <- paste0("ignored_tmp:", seq_len(nrow(scope)))
  queue <- mv05ac_build_evaluation_queue_v1(
    scope, paste(rep("c", 64L), collapse = ""))
  queue$execution_authorized[[1L]] <- TRUE
  expect_error(mv05ac_validate_evaluation_queue_v1(queue),
               "prospective 150-group")
})
