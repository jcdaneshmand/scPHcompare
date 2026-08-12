test_that("MV5-AR freezes additive default-off integration", {
  tables <- mv05ar_prefreeze_tables_v1()
  expect_equal(nrow(tables$controls), 12L)
  expect_equal(nrow(tables$boundaries), 6L)
  expect_false(any(tables$boundaries$default_behavior_changed))
  expect_false(any(tables$boundaries$corrected_downstream_consumption))
  expect_equal(nrow(tables$artifacts), 7L)
  expect_false(any(tables$artifacts$overwrite_allowed))
  expect_false(any(tables$artifacts$legacy_filename_reused))
  expect_equal(tables$resources$value[tables$resources$policy ==
    "adaptive_pair_seconds"], 240)
  expect_equal(sum(tables$stages$authorized_now), 1L)
})

test_that("MV5-AR decision authorizes implementation only", {
  decision <- mv05ar_decide_v1(TRUE, TRUE, TRUE, TRUE, 0, 0, 0, 0)
  expect_true(decision$additive_implementation_authorized)
  expect_false(decision$corrected_downstream_consumption_authorized)
  expect_false(decision$workflow_default_change_authorized)
  expect_false(decision$legacy_artifact_rewrite_authorized)
  expect_false(decision$optimization_authorized)
  stopped <- mv05ar_decide_v1(TRUE, TRUE, TRUE, TRUE, 0, 1, 0, 0)
  expect_false(stopped$additive_implementation_authorized)
})
