mv05aj_fixture_queue <- function() {
  queue <- expand.grid(
    fold_id = paste0("large_loso_v1:S", sprintf("%02d", 1:15)),
    seed = 20260805:20260809,
    representation = c("inductive_integrated", "sct_whole"),
    stringsAsFactors = FALSE)
  queue <- queue[order(queue$fold_id, queue$representation, queue$seed), ]
  n <- nrow(queue)
  queue$contract_id <- "mv05aj_nested256_execution_queue_v1"
  queue$robustness_group_id <- paste0("g", seq_len(n))
  queue$configuration_id <- "nested_cells_256_pc30_euclidean_v1"
  queue$configuration_order <- 4L
  queue$cells <- 256L; queue$coordinates <- 30L
  queue$point_metric <- "euclidean"
  queue$view_count <- 90L
  fold_training <- stats::setNames(65:79,
    paste0("large_loso_v1:S", sprintf("%02d", 1:15)))
  queue$training_samples <- unname(fold_training[queue$fold_id])
  pairs <- rep(471L, n)
  pairs[seq_len(50L)] <- pairs[seq_len(50L)] + 1L
  queue$biological_pairs <- pairs
  queue$landscape_request_rows <- 2L * pairs
  queue$landscape_subchunks <- rep(4L, n)
  queue$landscape_subchunks[seq_len(120L)] <- 5L
  queue$energy_request_rows <- pairs
  queue$assembled_method_rows <- 4L * pairs
  hash <- paste(rep("a", 64L), collapse = "")
  queue$coordinate_source_sha256 <- hash
  queue$base_pair_axis_sha256 <- hash
  queue$source_freeze_sha256 <- hash
  queue$implementation_sha256 <- hash
  queue$prospective_head <- paste(rep("b", 40L), collapse = "")
  queue$python_executable_sha256 <- hash
  queue$python_version <- "3.10"; queue$persim_version <- "0.3.8"
  queue$numpy_version <- "2"; queue$scipy_version <- "1"
  queue$configuration_execution_order <- seq_len(n)
  queue$deterministic_repeat_required <- FALSE
  queue$deterministic_repeat_required[c(1L, n)] <- TRUE
  queue$representation[[1L]] <- "inductive_integrated"
  queue$representation[[n]] <- "sct_whole"
  queue$execution_authorized <- TRUE; queue$execution_completed <- FALSE
  queue$labels_opened <- FALSE; queue$rankings_computed <- FALSE
  queue$outcomes_computed <- FALSE
  queue
}

test_that("MV5-AJ accepts only the complete nested-256 axis", {
  queue <- mv05aj_fixture_queue()
  expect_invisible(mv05aj_validate_nested_queue_v1(queue))
  scope <- mv05aj_nested_scope_v1(queue)
  expect_equal(scope$groups, 150L)
  expect_equal(scope$biological_pairs, 70700L)
  expect_equal(scope$landscape_rows, 141400L)
  expect_equal(scope$method_rows, 282800L)
})

test_that("MV5-AJ rejects configuration, order, and firewall drift", {
  queue <- mv05aj_fixture_queue()
  broken <- queue; broken$configuration_id[[1L]] <- "cells384_pc20_euclidean_v1"
  expect_error(mv05aj_validate_nested_queue_v1(broken), "150 unopened nested-256")
  broken <- queue; broken$configuration_execution_order[[2L]] <- 1L
  expect_error(mv05aj_validate_nested_queue_v1(broken), "150 unopened nested-256")
  broken <- queue; broken$outcomes_computed[[1L]] <- TRUE
  expect_error(mv05aj_validate_nested_queue_v1(broken), "150 unopened nested-256")
  broken <- queue; broken$tissue <- "x"
  expect_error(mv05aj_validate_nested_queue_v1(broken), "150 unopened nested-256")
})

test_that("MV5-AJ rejects incomplete counts and resource-axis drift", {
  queue <- mv05aj_fixture_queue()
  expect_error(mv05aj_validate_nested_queue_v1(queue[-1L, ]),
               "150 unopened nested-256")
  broken <- queue; broken$biological_pairs[[1L]] <-
    broken$biological_pairs[[1L]] + 1L
  broken$landscape_request_rows[[1L]] <-
    2L * broken$biological_pairs[[1L]]
  broken$energy_request_rows[[1L]] <- broken$biological_pairs[[1L]]
  broken$assembled_method_rows[[1L]] <- 4L * broken$biological_pairs[[1L]]
  expect_error(mv05aj_validate_nested_queue_v1(broken),
               "150 unopened nested-256")
})

test_that("MV5-AJ scope freezes the deterministic nested point policy", {
  scope <- mv05aj_nested_scope_v1(mv05aj_fixture_queue())
  expect_identical(scope$point_selection,
                   "sha256_sample_seed_cell_nested_v1")
  expect_false(scope$other_configuration_execution_authorized)
  expect_equal(scope$configuration_elapsed_cap_seconds, 6 * 60 * 60)
})

test_that("MV5-AJ 192 selection is the exact prefix of nested 256", {
  point_ids <- sprintf("cell_%03d", seq_len(384L))
  nested192 <- mv05t_nested_point_ids_v1(
    "sample_a", 20260805L, point_ids, 192L)
  nested256 <- mv05t_nested_point_ids_v1(
    "sample_a", 20260805L, point_ids, 256L)
  expect_identical(nested256[seq_len(192L)], nested192)
  expect_equal(length(setdiff(nested256, nested192)), 64L)
  expect_true(all(nested256 %in% point_ids))
})
