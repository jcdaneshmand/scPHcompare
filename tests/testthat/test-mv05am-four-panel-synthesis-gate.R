test_that("MV5-AM registry exhausts exactly four configurations", {
  x <- mv05am_panel_registry_v1()
  expect_equal(nrow(x), 4L)
  expect_identical(x$panel_order, 1:4)
  expect_false(any(x$configuration_sequence_exhausted_after_panel[1:3]))
  expect_true(x$configuration_sequence_exhausted_after_panel[[4L]])
})

test_that("MV5-AM requires complete unsliced panels", {
  macro <- data.frame(estimand_id = paste0("e", 1:24), estimate = 1:24)
  intervals <- data.frame(estimand_id = paste0("e", 1:24), estimate = 1:24,
    ci_lower = 0:23, ci_upper = 2:25)
  primary <- data.frame(estimand_id = paste0("e", 1:4),
    raw_p_value = rep(.5, 4), holm_p_value = rep(1, 4))
  production <- as.data.frame(as.list(c(prediction_groups = 150,
    ranking_rows = 282800, outcome_groups = 150, query_method_rows = 3600,
    long_query_endpoint_rows = 7200, macro_estimands = 24, intervals = 24,
    primary_tests = 4, representations = 2, methods = 4, endpoints = 2,
    samples = 90, studies = 15, tissues = 5, seeds = 5)))
  production$clustering_executed <- FALSE
  bound <- mv05am_bind_complete_panel_v1(
    "pc20_vs_pc30", production, macro, intervals, primary)
  expect_true(bound$complete_panel_bound)
  expect_equal(bound$result_rows_used_for_continuation_decision, 0L)
  expect_error(mv05am_bind_complete_panel_v1(
    "pc20_vs_pc30", production, macro[-1, ], intervals, primary),
    "complete 24-estimand")
})

test_that("MV5-AM preserves corrected landscape semantics", {
  x <- mv05am_landscape_contract_v1()
  expect_equal(nrow(x), 8L)
  expect_true(all(x$preserved_in_all_four_panels))
  expect_false(any(x$changed_by_synthesis))
  expect_true("all_consecutive_active_levels" %in% x$required_state)
  expect_true("no_fixed_uniform_grid" %in% x$required_state)
})

test_that("MV5-AM decision cannot consume results or authorize execution", {
  registry <- mv05am_panel_registry_v1()
  bindings <- data.frame(panel_id = registry$panel_id,
    complete_panel_bound = TRUE, result_rows_used_for_continuation_decision = 0L)
  criteria <- mv05am_continuation_criteria_v1()["criterion_id"]
  criteria$passed <- TRUE
  decision <- mv05am_decide_v1(registry, bindings,
    mv05am_comparability_v1(), mv05am_evidence_gaps_v1(), criteria)
  expect_identical(decision$authorized_next_sprint, "MV5-AN")
  expect_true(decision$prefreeze_only)
  expect_false(any(unlist(decision[c("new_calculation_authorized",
    "clustering_authorized", "public_default_change_authorized")])))
  bindings$estimate <- 1:4
  expect_error(mv05am_decide_v1(registry, bindings,
    mv05am_comparability_v1(), mv05am_evidence_gaps_v1(), criteria),
    "prohibited result input")
})
