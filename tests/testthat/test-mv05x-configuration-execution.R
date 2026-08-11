mv05x_test_queue <- function() {
  configurations <- c(
    "cells384_pc20_euclidean_v1", "cells384_pc30_cosine_chord_v1",
    "nested_cells_192_pc30_euclidean_v1",
    "nested_cells_256_pc30_euclidean_v1"
  )
  base <- expand.grid(
    representation = c("inductive_integrated", "sct_whole"),
    fold = sprintf("SRA%06d", 1:15), seed = 20260805:20260809,
    stringsAsFactors = FALSE
  )[seq_len(150), ]
  base$biological_pairs <- c(521L, rep(471L, 149L))
  base$landscape_subchunks <- c(rep(5L, 120L), rep(4L, 30L))
  rows <- do.call(rbind, lapply(seq_along(configurations), function(index) {
    data.frame(
      contract_id = "mv05x_pc20_execution_queue_v1",
      robustness_group_id = paste0(
        "mv05v_group_v1:",
        sprintf("%064x", seq_len(150) + 150L * (index - 1L))
      ),
      fold_id = paste0("large_loso_v1:", base$fold), seed = base$seed,
      representation = base$representation,
      configuration_id = configurations[[index]], cells = 384L,
      coordinates = 30L, point_metric = "euclidean",
      heldout_samples = 5L, training_samples = 85L, view_count = 90L,
      biological_pairs = base$biological_pairs,
      landscape_request_rows = 2L * base$biological_pairs,
      landscape_subchunks = base$landscape_subchunks,
      energy_request_rows = base$biological_pairs,
      assembled_method_rows = 4L * base$biological_pairs,
      deterministic_repeat_required = FALSE,
      coordinate_source_sha256 = strrep("a", 64),
      base_pair_axis_sha256 = strrep("b", 64),
      source_freeze_sha256 = strrep("c", 64),
      implementation_sha256 = strrep("d", 64),
      prospective_head = strrep("e", 40), outcome_label_state = "closed",
      outcomes_computed = FALSE, execution_authorized = index == 1L,
      execution_completed = FALSE,
      python_executable_sha256 = strrep("f", 64),
      python_version = "3", persim_version = "1",
      numpy_version = "2", scipy_version = "1",
      configuration_execution_order = if (index == 1L) 1:150 else NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
  rows$execution_order <- seq_len(nrow(rows))
  rows$deterministic_repeat_required[
    unlist(lapply(split(seq_len(nrow(rows)), paste(
      rows$configuration_id, rows$representation
    )), `[[`, 1L), use.names = FALSE)
  ] <- TRUE
  rows
}

test_that("MV5-X authorizes exactly the complete PC20 configuration", {
  queue <- mv05x_test_queue()
  expect_invisible(mv05x_validate_configuration_queue_v1(queue))
  scope <- mv05x_configuration_scope_v1(queue)
  expect_equal(scope$authorized_groups, 150L)
  expect_equal(scope$views, 13500L)
  expect_equal(scope$biological_pairs, 70700L)
  expect_equal(scope$landscape_rows, 141400L)
  expect_equal(scope$energy_rows, 70700L)
  expect_equal(scope$method_rows, 282800L)
  expect_false(scope$labels_opened)
  expect_false(scope$outcomes_computed)
})

test_that("MV5-X fails closed on authorization and ordering leakage", {
  queue <- mv05x_test_queue()
  queue$execution_authorized[[151L]] <- TRUE
  expect_error(mv05x_validate_configuration_queue_v1(queue),
               "exactly its 150")
  queue <- mv05x_test_queue()
  queue$configuration_execution_order[[1L]] <- 2L
  expect_error(mv05x_validate_configuration_queue_v1(queue),
               "exactly its 150")
  queue <- mv05x_test_queue()
  queue$outcomes_computed[[1L]] <- TRUE
  expect_error(mv05x_validate_configuration_queue_v1(queue),
               "exactly its 150")
})
