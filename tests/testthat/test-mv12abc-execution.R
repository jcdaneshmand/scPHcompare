test_that("MV12-A through MV12-C scripts parse and preserve gates", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv12a_historical_fusion_prefreeze.R",
    "run_mv12b_historical_fusion.R",
    "build_mv12c_historical_fusion_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  for (path in scripts) expect_silent(parse(path))
  text <- setNames(lapply(scripts, function(path) paste(
    readLines(path, warn = FALSE), collapse = "\n"
  )), basename(scripts))
  prefreeze <- text[["build_mv12a_historical_fusion_prefreeze.R"]]
  expect_match(prefreeze, "partition_fits = 500L", fixed = TRUE)
  expect_match(prefreeze, "private_assignment_rows = 62000L", fixed = TRUE)
  expect_match(prefreeze, "primary_gene_weight = 0.5", fixed = TRUE)
  expect_match(prefreeze, "option2_required_unless", fixed = TRUE)
  expect_match(prefreeze, "automatic_retries = 0L", fixed = TRUE)
  runner <- text[["run_mv12b_historical_fusion.R"]]
  expect_match(runner, "mv12_fit_fusion_grid_v1(bundle)", fixed = TRUE)
  expect_match(runner, "mv12_fusion_decision_v1(stability, consensus)",
               fixed = TRUE)
  expect_match(runner, "method_or_weight_selected = FALSE", fixed = TRUE)
  closure <- text[["build_mv12c_historical_fusion_closure.R"]]
  expect_match(closure, "all_scientific_artifacts_exact", fixed = TRUE)
  expect_match(closure, "option2_new_allqc_cell_topology_required", fixed = TRUE)
})
