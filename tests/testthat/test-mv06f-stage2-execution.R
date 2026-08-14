source(testthat::test_path("..", "..", "R", "mv06f_production.R"))
source(testthat::test_path("..", "..", "R", "mv06f_stage2_execution.R"))

testthat::test_that("MV6-F stage-two guards stop each frozen cap", {
  plan <- data.frame(
    guard = c(
      "group_elapsed_seconds", "group_process_tree_rss_bytes",
      "concurrent_process_tree_rss_bytes", "production_worker_seconds",
      "private_root_bytes", "maximum_workers"
    ), value = c(1800, 8 * 1024^3, 12 * 1024^3, 14.4 * 3600,
                 10 * 1024^3, 1)
  )
  pass <- mv06f_stage2_guard_v1(0, 10, 10, 100, 100, 100, plan)
  testthat::expect_true(pass$launch_authorized)
  cases <- list(
    c(0, 1801, 1801, 100, 100, 100),
    c(0, 10, 10, 8 * 1024^3 + 1, 100, 100),
    c(0, 10, 10, 100, 12 * 1024^3 + 1, 100),
    c(14.4 * 3600, 1, 1, 100, 100, 100),
    c(0, 10, 10, 100, 100, 10 * 1024^3 + 1)
  )
  testthat::expect_true(all(vapply(cases, function(value) {
    !mv06f_stage2_guard_v1(
      value[[1L]], value[[2L]], value[[3L]], value[[4L]], value[[5L]],
      value[[6L]], plan
    )$launch_authorized
  }, logical(1L))))
})

testthat::test_that("MV6-F stage-two admission binds every stage-one gate", {
  contract <- data.frame(
    queue_root_sha256 = strrep("a", 64L),
    implementation_root_sha256 = strrep("b", 64L),
    rust_library_sha256 = strrep("c", 64L)
  )
  admission <- data.frame(
    contract_id = "mv06f_stage2_admission_v1",
    queue_root_sha256 = contract$queue_root_sha256,
    implementation_root_sha256 = contract$implementation_root_sha256,
    rust_library_sha256 = contract$rust_library_sha256,
    stage1_scientific_equivalence_passed = TRUE,
    stage1_clean_repeat_passed = TRUE, stage1_r_oracles_passed = TRUE,
    stage1_persim_oracles_passed = TRUE, stage1_resume_passed = TRUE,
    stage2_authorized = TRUE, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, fusion_jobs = 0L,
    clustering_jobs = 0L, outcome_jobs = 0L
  )
  testthat::expect_silent(mv06f_validate_stage2_admission_v1(
    admission, contract
  ))
  admission$stage1_persim_oracles_passed <- FALSE
  testthat::expect_error(mv06f_validate_stage2_admission_v1(
    admission, contract
  ), "incomplete")
})
