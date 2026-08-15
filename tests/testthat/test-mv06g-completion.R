source(testthat::test_path("..", "..", "R", "mv06f_production.R"))
source(testthat::test_path("..", "..", "R", "mv06g_completion.R"))

testthat::test_that("MV6-G completion source inventory is exact", {
  paths <- mv06g_completion_source_paths_v1()
  testthat::expect_length(paths, 9L)
  testthat::expect_identical(anyDuplicated(paths), 0L)
  testthat::expect_true(all(file.exists(testthat::test_path("..", "..", paths))))
})

testthat::test_that("MV6-G completion guard fails every frozen cap", {
  policy <- data.frame(
    elapsed_cap_seconds_per_group = 1800,
    rss_cap_bytes_per_group = 12 * 1024^3,
    private_storage_cap_bytes = 5 * 1024^3,
    total_worker_cap_seconds = 12 * 3600
  )
  testthat::expect_true(mv06g_completion_guard_v1(0, 1, 1, 1, policy)$pass)
  cases <- list(
    c(0, 1801, 1, 1), c(0, 1, 12 * 1024^3 + 1, 1),
    c(0, 1, 1, 5 * 1024^3 + 1), c(12 * 3600, 1, 1, 1)
  )
  testthat::expect_true(all(vapply(cases, function(x) {
    !mv06g_completion_guard_v1(x[[1L]], x[[2L]], x[[3L]], x[[4L]],
                               policy)$pass
  }, logical(1L))))
})

testthat::test_that("MV6-G completion metric remains label closed", {
  policy <- data.frame(
    rss_cap_bytes_per_group = 100, total_worker_cap_seconds = 100,
    private_storage_cap_bytes = 100,
    scientific_implementation_root_sha256 = strrep("a", 64L),
    execution_implementation_root_sha256 = strrep("b", 64L),
    rust_library_sha256 = strrep("c", 64L)
  )
  row <- data.frame(group_id = "g", execution_order = 2L)
  metric <- data.frame(
    contract_id = "mv06g_serial_completion_metric_v1", group_id = "g",
    execution_order = 2L, disposition = "completed", exit_status = 0L,
    elapsed_seconds = 1, charged_worker_seconds = 1,
    peak_process_tree_rss_bytes = 1, cumulative_worker_seconds = 1,
    cumulative_private_bytes = 1, scientific_group_complete = TRUE,
    scientific_implementation_root_sha256 = strrep("a", 64L),
    execution_implementation_root_sha256 = strrep("b", 64L),
    rust_library_sha256 = strrep("c", 64L),
    runner_stdout_sha256 = strrep("d", 64L),
    runner_stderr_sha256 = strrep("e", 64L), outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, fusion_evaluations = 0L,
    outcome_jobs = 0L
  )
  testthat::expect_silent(mv06g_validate_completion_metric_v1(
    metric, row, policy
  ))
  metric$outcome_jobs <- 1L
  testthat::expect_error(mv06g_validate_completion_metric_v1(
    metric, row, policy
  ), "stale or failed")
})
