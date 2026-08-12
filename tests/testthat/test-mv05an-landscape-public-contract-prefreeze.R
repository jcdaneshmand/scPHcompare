test_that("MV5-AN classifies every known landscape family", {
  expect_identical(mv05an_classify_landscape_function_v1(
    "landscape_reference.R", "landscape_reference_distance")$pathway_class,
    "corrected_exact_or_error_controlled_engine")
  expect_identical(mv05an_classify_landscape_function_v1(
    "PH_PostProcessing_andAnalysis.R", "ComputePersistenceLandscapes")$pathway_class,
    "legacy_grid_artifact_workflow")
  expect_true(mv05an_classify_landscape_function_v1(
    "unknown.R", "landscape_x")$ambiguous)
})

test_that("MV5-AN freezes the complete scientific contract", {
  x <- mv05an_target_contract_v1()
  expect_equal(nrow(x), 12L)
  expect_true(all(x$immutable_scientific_semantics))
  expect_true("all_consecutive_active_levels_zero_pad_missing_depth" %in%
                x$required_state)
  expect_true("no_universal_uniform_grid_or_level_cap" %in% x$required_state)
})

test_that("MV5-AN authorizes new APIs without redirecting defaults", {
  x <- mv05an_public_api_decision_v1()
  expect_equal(nrow(x), 2L)
  expect_true(all(x$export_in_later_sprint))
  expect_false(any(x$legacy_function_redirected | x$legacy_artifact_overwritten |
                     x$workflow_default_changed))
  migration <- mv05an_migration_plan_v1()
  expect_identical(which(migration$authorized_in_mv05ao), 1:4)
  expect_false(any(migration$behavior_change[1:4]))
})

test_that("MV5-AN decision fails on ambiguity and preserves stop boundary", {
  inventory <- data.frame(ambiguous = rep(FALSE, 45))
  entry <- data.frame(x = 1:6); artifacts <- data.frame(x = 1:6)
  validation <- data.frame(passed = TRUE)
  decision <- mv05an_decide_v1(inventory, entry, artifacts,
    mv05an_target_contract_v1(), mv05an_migration_plan_v1(), validation)
  expect_identical(decision$authorized_next_sprint, "MV5-AO")
  expect_false(decision$current_workflow_default_change_authorized)
  inventory$ambiguous[[1]] <- TRUE
  expect_error(mv05an_decide_v1(inventory, entry, artifacts,
    mv05an_target_contract_v1(), mv05an_migration_plan_v1(), validation),
    "incomplete inventory")
})
