test_that("MV5-P frozen group projections reconcile exactly", {
  repo_path <- function(...) testthat::test_path("..", "..", ...)
  mv05p_group_projection_v1 <- getFromNamespace(
    "mv05p_group_projection_v1", "scPHcompare")
  mv05p_validate_group_queue_v1 <- getFromNamespace(
    "mv05p_validate_group_queue_v1", "scPHcompare")
  representation <- rep(c("sct_whole", "inductive_integrated"), each = 75L)
  groups <- data.frame(
    production_group_id = sprintf("synthetic_group_%03d", seq_len(150L)),
    execution_order = seq_len(150L), representation = representation,
    landscape_request_rows = 100L, energy_pair_rows = 50L,
    shared_pseudobulk_pair_rows = ifelse(representation == "sct_whole",
                                         50L, 0L),
    max_parallel_workers = 2L, per_unit_elapsed_cap_seconds = 900,
    per_process_rss_cap_bytes = 4294967296,
    stage_worker_hour_cap = 21.6,
    stage_private_storage_cap_bytes = 10737418240,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    clustering_jobs_executed = 0L, production_executed = FALSE,
    source_freeze_sha256 = strrep("a", 64L), stringsAsFactors = FALSE
  )
  combined <- data.frame(
    component = c(
      "sct_whole_energy", "inductive_integrated_energy",
      "shared_pseudobulk_reused_across_representations",
      "all_required_distance_matrices_with_reserve"
    ),
    projected_worker_hours = c(
      1.96094184027778, 2.00654513888889, 0.022801649305555,
      16.1170472529584
    ), stringsAsFactors = FALSE
  )
  landscape <- data.frame(
    representation = c("sct_whole", "inductive_integrated"),
    projected_worker_hours = c(7.77175293860758, 2.88981957197332),
    projected_output_bytes = c(304029896, 314815988),
    stringsAsFactors = FALSE
  )
  projection <- mv05p_group_projection_v1(groups, combined, landscape)
  expect_equal(nrow(projection), 150L)
  expect_equal(sum(projection$projected_worker_seconds) / 3600,
               16.1170472529584, tolerance = 1e-10)
  expect_equal(sum(projection$projected_private_bytes), 1277893355,
               tolerance = 1)
  expect_true(all(projection$projected_worker_seconds > 0))
  expect_true(all(projection$projected_private_bytes > 0))
  real_queue <- repo_path(
    "docs", "audits", "mv05o-production-group-queue-2026-08-10.csv")
  if (file.exists(real_queue)) {
    expect_silent(mv05p_validate_group_queue_v1(read.csv(
      real_queue, stringsAsFactors = FALSE, check.names = FALSE
    )))
  }
})

test_that("MV5-P launch gate closes at either frozen cap", {
  mv05p_launch_budget_v1 <- getFromNamespace(
    "mv05p_launch_budget_v1", "scPHcompare")
  pass <- mv05p_launch_budget_v1(0, 16.1170472529584 * 3600,
                                 0, 1277893355)
  expect_true(pass$launch_authorized)
  expect_true(pass$worker_cap_passed)
  expect_true(pass$storage_cap_passed)
  expect_false(mv05p_launch_budget_v1(
    21.6 * 3600, 1, 0, 0)$launch_authorized)
  expect_false(mv05p_launch_budget_v1(
    0, 0, 10737418240, 1)$launch_authorized)
  expect_error(mv05p_launch_budget_v1(-1, 0, 0, 0), "invalid values")
})

test_that("MV5-P orchestrators preserve label-closed scope", {
  worker_path <- testthat::test_path(
    "..", "..", "scripts", "run_mv05p_production_group.R")
  monitor_path <- testthat::test_path(
    "..", "..", "scripts", "monitor_mv05p_production.R")
  if (file.exists(worker_path) && file.exists(monitor_path)) {
    worker <- paste(readLines(worker_path, warn = FALSE), collapse = "\n")
    monitor <- paste(readLines(monitor_path, warn = FALSE), collapse = "\n")
    expect_match(worker, 'outcome_label_state = "closed"', fixed = TRUE)
    expect_match(worker, "clustering_jobs_executed = 0L")
    expect_match(monitor, "persim.__version__ == '0.3.8'", fixed = TRUE)
    expect_match(monitor, "length(running) < 2L", fixed = TRUE)
    expect_false(grepl("ARI|NMI|tissue.*outcome|method_selection",
                       paste(worker, monitor), ignore.case = TRUE))
  } else {
    expect_true(TRUE)
  }
})
