test_that("MV7-A preserves the corrected landscape contract", {
  x <- mv07a_landscape_contract_v1()
  expect_equal(nrow(x), 8L)
  expect_true(all(x$preserved))
  expect_false(any(x$changed_by_mv07a))
  expect_true("all_consecutive_active_levels" %in% x$required_state)
  expect_true("no_universal_fixed_grid" %in% x$required_state)
  expect_true("no_universal_level_cap" %in% x$required_state)
})

test_that("MV7-A coverage distinguishes completed and missing axes", {
  r <- mv07a_robustness_coverage_v1()
  c <- mv07a_confounding_coverage_v1()
  expect_equal(nrow(r), 14L)
  expect_equal(nrow(c), 10L)
  expect_equal(r$coverage[r$axis_id == "homology_dimension"], "complete")
  expect_equal(r$coverage[r$axis_id == "gene_panel_size"], "gap")
  expect_equal(c$coverage[c$axis_id == "library_size"], "unavailable")
  expect_false(any(c$causal_adjustment_authorized))
})

test_that("MV7-A authorizes only existing-artifact diagnostics", {
  d <- mv07a_decide_v1(mv07a_robustness_coverage_v1(),
    mv07a_confounding_coverage_v1(), mv07a_gap_registry_v1(),
    mv07a_landscape_contract_v1())
  expect_equal(d$authorized_next_sprint, "MV7-B")
  expect_true(d$prefreeze_only)
  expect_false(d$new_ph_authorized)
  expect_false(d$new_data_authorized)
  expect_false(d$method_or_weight_selection_authorized)
})

test_that("MV7-A decision rejects numerical-result interfaces", {
  r <- mv07a_robustness_coverage_v1()
  r$mrr <- 0
  expect_error(mv07a_decide_v1(r, mv07a_confounding_coverage_v1(),
    mv07a_gap_registry_v1(), mv07a_landscape_contract_v1()),
    "prohibited numerical result")
})
