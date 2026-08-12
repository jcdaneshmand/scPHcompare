test_that("MV5-AI recovers nested 256 as the sole canonical next setting", {
  ids <- c("cells384_pc20_euclidean_v1",
           "cells384_pc30_cosine_chord_v1",
           "nested_cells_192_pc30_euclidean_v1",
           "nested_cells_256_pc30_euclidean_v1")
  queue <- do.call(rbind, lapply(seq_along(ids), function(i) data.frame(
    configuration_id = ids[[i]], configuration_order = i,
    execution_order = ((i - 1L) * 150L + 1L):(i * 150L))))
  order <- mv05ai_configuration_order_v1(queue)
  expect_identical(order$configuration_id, ids)
  expect_identical(order$position_state,
                   c("complete", "complete", "complete", "next_eligible"))
  broken <- queue; broken$configuration_order[broken$configuration_order == 3L] <- 4L
  expect_error(mv05ai_configuration_order_v1(broken), "ambiguous or drifted")
})

test_that("MV5-AI binds complete unsliced robustness panels", {
  registry <- mv05ac_estimand_registry_v1()
  production <- data.frame(
    prediction_groups = 150, ranking_rows = 282800, outcome_groups = 150,
    query_method_rows = 3600, long_query_endpoint_rows = 7200,
    macro_estimands = 24, intervals = 24, primary_tests = 4,
    clustering_executed = FALSE, nested_configurations_executed = FALSE,
    method_selection_executed = FALSE, equivalence_claim_authorized = FALSE)
  primary <- registry[registry$estimand_role == "confirmatory_cosine_sensitivity", ]
  evidence <- mv05ai_bind_complete_evidence_v1(
    "cosine_chord_vs_euclidean", production, primary,
    data.frame(x = 1:24), registry)
  expect_true(evidence$complete_evidence_bound)
  expect_false(evidence$nested256_sensitivity_answered)
  production$macro_estimands <- 23
  expect_error(mv05ai_bind_complete_evidence_v1(
    "cosine_chord_vs_euclidean", production, primary,
    data.frame(x = 1:24), registry), "missing or drifted")
})

test_that("MV5-AI requires every criterion and cannot skip nested 256", {
  ids <- c("cells384_pc20_euclidean_v1",
           "cells384_pc30_cosine_chord_v1",
           "nested_cells_192_pc30_euclidean_v1",
           "nested_cells_256_pc30_euclidean_v1")
  order <- data.frame(configuration_order = 1:4, configuration_id = ids,
    position_state = c("complete", "complete", "complete", "next_eligible"))
  evidence <- data.frame(
    analysis_id = c("pc20_vs_pc30", "cosine_chord_vs_euclidean",
                    "nested192_vs_384_cells"),
    complete_evidence_bound = TRUE, nested256_sensitivity_answered = FALSE)
  pass <- mv05ai_continuation_criteria_v1()["criterion_id"]; pass$passed <- TRUE
  decision <- mv05ai_decide_v1(order, evidence, pass)
  expect_identical(decision$authorized_configuration_id,
                   "nested_cells_256_pc30_euclidean_v1")
  expect_true(decision$nested_256_calculation_authorized)
  pass$passed[[7L]] <- FALSE
  expect_error(mv05ai_decide_v1(order, evidence, pass), "cannot be authorized")
})

test_that("MV5-AI nested-256 queue remains label, rank, and outcome closed", {
  id <- "nested_cells_256_pc30_euclidean_v1"
  queue <- expand.grid(fold_id = paste0("f", 1:15), seed = 1:5,
    representation = c("sct_whole", "inductive_integrated"),
    stringsAsFactors = FALSE)
  queue$configuration_id <- id
  queue$robustness_group_id <- paste0("g", seq_len(nrow(queue)))
  queue$execution_order <- 451:600; queue$cells <- 256L
  queue$coordinates <- 30L; queue$point_metric <- "euclidean"
  queue$outcomes_computed <- FALSE
  decision <- data.frame(authorized_configuration_id = id,
    nested_256_calculation_executed = FALSE)
  result <- mv05ai_nested_256_queue_v1(queue, decision)
  expect_equal(nrow(result), 150L)
  expect_true(all(result$later_label_closed_calculation_authorized))
  expect_false(any(result$labels_opened | result$rankings_computed |
                   result$outcomes_computed | result$execution_completed))
})
