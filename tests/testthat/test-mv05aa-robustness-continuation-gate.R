test_that("MV5-AA recovers only the frozen canonical next configuration", {
  ids <- c("cells384_pc20_euclidean_v1",
           "cells384_pc30_cosine_chord_v1",
           "nested_cells_192_pc30_euclidean_v1",
           "nested_cells_256_pc30_euclidean_v1")
  queue <- do.call(rbind, lapply(seq_along(ids), function(i) data.frame(
    configuration_id = ids[[i]], configuration_order = i,
    execution_order = ((i - 1L) * 150L + 1L):(i * 150L))))
  order <- mv05aa_configuration_order_v1(queue)
  expect_identical(order$configuration_id, ids)
  expect_identical(order$position_state,
                   c("complete", "next_eligible", "later_closed", "later_closed"))
  broken <- queue
  broken$configuration_order[broken$configuration_order == 2L] <- 3L
  expect_error(mv05aa_configuration_order_v1(broken), "ambiguous or drifted")
})

test_that("MV5-AA binds complete PC20 reporting without result slicing", {
  production <- data.frame(prediction_groups = 150, ranking_rows = 282800,
    outcome_groups = 150, query_method_rows = 3600,
    long_query_endpoint_rows = 7200, macro_estimands = 24, intervals = 24,
    primary_tests = 4, clustering_executed = FALSE,
    other_configurations_executed = FALSE, method_selection_executed = FALSE,
    equivalence_claim_authorized = FALSE)
  evidence <- mv05aa_bind_pc20_evidence_v1(
    production, data.frame(x = 1:4), data.frame(x = 1:24),
    data.frame(x = 1:24))
  expect_true(evidence$complete_evidence_bound)
  expect_false(evidence$cosine_geometry_answered_by_pc20)
  production$macro_estimands <- 23
  expect_error(mv05aa_bind_pc20_evidence_v1(
    production, data.frame(x = 1:4), data.frame(x = 1:24),
    data.frame(x = 1:24)), "missing or drifted")
})

test_that("MV5-AA requires every frozen criterion and cannot skip cosine", {
  ids <- c("cells384_pc20_euclidean_v1",
           "cells384_pc30_cosine_chord_v1",
           "nested_cells_192_pc30_euclidean_v1",
           "nested_cells_256_pc30_euclidean_v1")
  order <- data.frame(contract_id = "x", configuration_order = 1:4,
    configuration_id = ids,
    position_state = c("complete", "next_eligible", "later_closed", "later_closed"))
  evidence <- data.frame(complete_evidence_bound = TRUE,
    cosine_geometry_answered_by_pc20 = FALSE)
  pass <- mv05aa_continuation_criteria_v1()["criterion_id"]
  pass$passed <- TRUE
  decision <- mv05aa_decide_v1(order, evidence, pass)
  expect_identical(decision$authorized_configuration_id,
                   "cells384_pc30_cosine_chord_v1")
  pass$passed[[7L]] <- FALSE
  expect_error(mv05aa_decide_v1(order, evidence, pass), "cannot be authorized")
})

test_that("MV5-AA cosine queue remains label and outcome closed", {
  id <- "cells384_pc30_cosine_chord_v1"
  queue <- expand.grid(fold_id = paste0("f", 1:15), seed = 1:5,
    representation = c("sct_whole", "inductive_integrated"),
    stringsAsFactors = FALSE)
  queue$configuration_id <- id
  queue$robustness_group_id <- paste0("g", seq_len(nrow(queue)))
  queue$execution_order <- 151:300
  queue$cells <- 384L
  queue$coordinates <- 30L
  queue$point_metric <- "euclidean_chord_after_row_unit_normalization"
  queue$outcomes_computed <- FALSE
  decision <- data.frame(authorized_configuration_id = id,
    cosine_calculation_executed = FALSE)
  result <- mv05aa_cosine_queue_v1(queue, decision)
  expect_equal(nrow(result), 150L)
  expect_true(all(result$later_label_closed_calculation_authorized))
  expect_false(any(result$labels_opened | result$outcomes_computed |
                   result$rankings_computed | result$execution_completed))
  queue$outcomes_computed[[1L]] <- TRUE
  expect_error(mv05aa_cosine_queue_v1(queue, decision), "contaminated")
})
