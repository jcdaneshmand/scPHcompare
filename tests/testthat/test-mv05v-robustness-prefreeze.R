mv05v_test_queue <- function() {
  configurations <- c(
    "cells384_pc20_euclidean_v1", "cells384_pc30_cosine_chord_v1",
    "nested_cells_192_pc30_euclidean_v1",
    "nested_cells_256_pc30_euclidean_v1"
  )
  base <- expand.grid(
    representation = c("inductive_integrated", "sct_whole"),
    fold = sprintf("SRA%06d", 1:15), seed = 20260805:20260809,
    stringsAsFactors = FALSE
  )
  base <- base[seq_len(150), ]
  base$biological_pairs <- c(521L, rep(471L, 149L))
  base$landscape_subchunks <- c(rep(5L, 120L), rep(4L, 30L))
  rows <- do.call(rbind, lapply(seq_along(configurations), function(index) {
    data.frame(
      contract_id = "mv05v_full_group_queue_v1",
      robustness_group_id = paste0("mv05v_group_v1:", sprintf("%064x", seq_len(150) + 150L * (index - 1L))),
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
      outcomes_computed = FALSE, execution_authorized = FALSE,
      execution_completed = FALSE, stringsAsFactors = FALSE
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

test_that("MV5-V accepts only the four frozen OFAT configurations", {
  configurations <- data.frame(
    configuration_id = c(
      "nested_cells_192_pc30_euclidean_v1",
      "nested_cells_256_pc30_euclidean_v1",
      "cells384_pc20_euclidean_v1",
      "cells384_pc30_cosine_chord_v1"
    ),
    candidate_id = c("nested", "nested", "pc20", "cosine"),
    cells = c(192L, 256L, 384L, 384L),
    coordinates = c(30L, 30L, 20L, 30L),
    point_metric = c("euclidean", "euclidean", "euclidean", "chord"),
    comparison_design = "one_factor_at_a_time_against_384cell_30pc_euclidean",
    baseline_reused = TRUE, role = "postoutcome_secondary_sensitivity",
    outcomes_authorized = FALSE, stringsAsFactors = FALSE
  )
  result <- mv05v_validate_configurations_v1(configurations)
  expect_equal(nrow(result), 4L)
  expect_identical(result$configuration_order, 1:4)
  configurations$outcomes_authorized[[1L]] <- TRUE
  expect_error(mv05v_validate_configurations_v1(configurations),
               "exactly the four")
})

test_that("MV5-V full queue contract freezes exact counts and label closure", {
  queue <- mv05v_test_queue()
  expect_invisible(mv05v_validate_group_queue_v1(queue))
  expect_equal(sum(queue$view_count), 54000L)
  expect_equal(sum(queue$landscape_request_rows), 565600L)
  expect_equal(sum(queue$energy_request_rows), 282800L)
  expect_equal(sum(queue$landscape_subchunks), 2880L)
  expect_equal(sum(queue$deterministic_repeat_required), 8L)
  queue$outcomes_computed[[1L]] <- TRUE
  expect_error(mv05v_validate_group_queue_v1(queue), "violates")
})

test_that("MV5-V projection keeps execution separate from prefreeze", {
  historical <- c(
    sct_ph = 3767.75852203369,
    integrated_ph = 3952.67263936996,
    sct_landscape = 4194.15161585808,
    integrated_landscape = 2072.89317178726,
    sct_assembly = 1461.44837808609,
    integrated_assembly = 1896.36905193329
  )
  projection <- mv05v_resource_projection_v1(historical)
  expect_equal(nrow(projection), 3L)
  expect_equal(attr(projection, "projected_worker_hours"),
               19.27254819896, tolerance = 1e-8)
  decision <- mv05v_prefreeze_decision_v1(mv05v_test_queue(), projection)
  expect_true(decision$gate_complete)
  expect_false(decision$full_execution_authorized)
  expect_false(decision$labels_opened)
  expect_false(decision$outcomes_computed)
})
