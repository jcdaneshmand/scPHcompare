test_that("MV11-B through MV11-F scripts parse and preserve prospective gates", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv11b_cell_benchmark_prefreeze.R",
    "run_mv11_cell_matrix_worker.R",
    "run_mv11c_cell_benchmark_sentinel.R",
    "build_mv11d_cell_benchmark_sentinel_closure.R",
    "run_mv11e_full_cell_benchmark.R",
    "build_mv11f_full_cell_benchmark_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  for (path in scripts) expect_silent(parse(path))
  text <- setNames(lapply(scripts, function(path) paste(
    readLines(path, warn = FALSE), collapse = "\n"
  )), basename(scripts))

  prefreeze <- text[["build_mv11b_cell_benchmark_prefreeze.R"]]
  expect_match(prefreeze, "partition_fits = 450L", fixed = TRUE)
  expect_match(prefreeze, "private_assignment_rows = 55800L", fixed = TRUE)
  expect_match(prefreeze,
               "sentinel_execution_authorized_after_commit = TRUE",
               fixed = TRUE)
  expect_match(prefreeze, "full_execution_authorized = FALSE", fixed = TRUE)
  expect_match(prefreeze, "automatic_retries = 0L", fixed = TRUE)
  expect_match(prefreeze, "cross_view_comparison_allowed <- FALSE",
               fixed = TRUE)
  validation_position <- regexpr("if (!all(validation$passed))", prefreeze,
                                 fixed = TRUE)[[1L]]
  create_position <- regexpr("dir.create(output", prefreeze,
                             fixed = TRUE)[[1L]]
  expect_gt(validation_position, 0L)
  expect_gt(create_position, validation_position)

  worker <- text[["run_mv11_cell_matrix_worker.R"]]
  expect_match(worker, "mv11_cell_matrix_v1", fixed = TRUE)
  expect_match(worker, "mv10_partition_grid_v1(matrix)", fixed = TRUE)
  expect_match(worker, "nrow(partitions) != 5580L", fixed = TRUE)
  sentinel <- text[["run_mv11c_cell_benchmark_sentinel.R"]]
  expect_match(sentinel, "process$kill_tree()", fixed = TRUE)
  expect_match(sentinel, "tree_rss(process$get_pid())", fixed = TRUE)
  expect_false(grepl("--resume", sentinel, fixed = TRUE))
  sentinel_closure <-
    text[["build_mv11d_cell_benchmark_sentinel_closure.R"]]
  expect_match(sentinel_closure, "aggregate_assignment_exact_repeat",
               fixed = TRUE)
  expect_match(sentinel_closure,
               "full_execution_authorized_after_commit = TRUE", fixed = TRUE)

  full <- text[["run_mv11e_full_cell_benchmark.R"]]
  expect_match(full, "for (i in seq_len(nrow(queue)))", fixed = TRUE)
  expect_match(full, "mv10_seed_stability_v1(assignments)", fixed = TRUE)
  expect_match(full, "mv10_method_agreement_v1(assignments)", fixed = TRUE)
  expect_match(full, "mv11_select_primary_k_v1(assignments)", fixed = TRUE)
  expect_false(grepl("--resume", full, fixed = TRUE))
  closure <- text[["build_mv11f_full_cell_benchmark_closure.R"]]
  expect_match(closure, "assignments_exact_repeat", fixed = TRUE)
  expect_match(closure, "cross_view_comparison_authorized_now = FALSE",
               fixed = TRUE)
})
