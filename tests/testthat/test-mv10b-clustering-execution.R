test_that("MV10-B through MV10-F execution scripts parse and preserve gates", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv10b_clustering_execution_prefreeze.R",
    "run_mv10_clustering_matrix_worker.R",
    "run_mv10c_clustering_sentinel.R",
    "build_mv10d_clustering_sentinel_closure.R",
    "run_mv10e_full_clustering_benchmark.R",
    "build_mv10f_full_clustering_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  for (path in scripts) expect_silent(parse(path))

  text <- setNames(lapply(scripts, function(path) paste(
    readLines(path, warn = FALSE), collapse = "\n"
  )), basename(scripts))
  prefreeze <- text[["build_mv10b_clustering_execution_prefreeze.R"]]
  expect_match(prefreeze, "sentinel_execution_authorized_after_commit = TRUE",
               fixed = TRUE)
  expect_match(prefreeze, "full_execution_authorized = FALSE", fixed = TRUE)
  expect_match(prefreeze, "corrected_primary_representation_H1", fixed = TRUE)
  expect_match(prefreeze, "process_tree_rss_cap_bytes = 4 * 1024^3",
               fixed = TRUE)
  expect_match(prefreeze, "automatic_retries = 0L", fixed = TRUE)
  validation_position <- regexpr(
    "if (!all(validation$passed))", prefreeze, fixed = TRUE
  )[[1L]]
  create_position <- regexpr(
    "if (!dir.exists(output)) dir.create", prefreeze, fixed = TRUE
  )[[1L]]
  expect_gt(validation_position, 0L)
  expect_gt(create_position, validation_position)

  worker <- text[["run_mv10_clustering_matrix_worker.R"]]
  expect_match(worker, "mv10_partition_grid_v1(matrix)", fixed = TRUE)
  expect_match(worker, "nrow(partitions) != 5580L", fixed = TRUE)
  expect_match(worker, "nrow(quality) != 45L", fixed = TRUE)
  sentinel <- text[["run_mv10c_clustering_sentinel.R"]]
  expect_match(sentinel, "process$kill_tree()", fixed = TRUE)
  expect_match(sentinel, "tree_rss(process$get_pid())", fixed = TRUE)
  expect_match(sentinel, "truth(decision$full_execution_authorized)",
               fixed = TRUE)
  expect_false(grepl("--resume", sentinel, fixed = TRUE))

  sentinel_closure <- text[["build_mv10d_clustering_sentinel_closure.R"]]
  expect_match(sentinel_closure, "mv10_partition_grid_v1(matrix)", fixed = TRUE)
  expect_match(sentinel_closure, "partition_identity", fixed = TRUE)
  expect_match(sentinel_closure,
               "full_execution_authorized_after_commit = TRUE", fixed = TRUE)

  runner <- text[["run_mv10e_full_clustering_benchmark.R"]]
  expect_match(runner, "for (i in seq_len(nrow(queue)))", fixed = TRUE)
  expect_match(runner, "process$kill_tree()", fixed = TRUE)
  expect_match(runner, "mv10_seed_stability_v1(assignments)", fixed = TRUE)
  expect_match(runner, "mv10_method_agreement_v1(assignments)", fixed = TRUE)
  expect_match(runner, "mv10_select_primary_k_v1(assignments)", fixed = TRUE)
  expect_false(grepl("--resume", runner, fixed = TRUE))

  full_closure <- text[["build_mv10f_full_clustering_closure.R"]]
  expect_match(full_closure, "mv08zy_read_distance_stack_v1", fixed = TRUE)
  expect_match(full_closure, "mv10_partition_grid_v1(matrix)", fixed = TRUE)
  expect_match(full_closure, "aggregate_assignment_exact_repeat", fixed = TRUE)
  expect_match(full_closure, "result_interpretation_state", fixed = TRUE)
})
