test_that("MV11-H through MV11-J scripts parse and preserve firewalls", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root, "scripts", c(
    "build_mv11h_cross_view_prefreeze.R",
    "run_mv11i_cross_view_agreement.R",
    "build_mv11j_cross_view_closure.R"
  ))
  expect_true(all(file.exists(scripts)))
  for (path in scripts) expect_silent(parse(path))
  text <- setNames(lapply(scripts, function(path) paste(
    readLines(path, warn = FALSE), collapse = "\n"
  )), basename(scripts))
  prefreeze <- text[["build_mv11h_cross_view_prefreeze.R"]]
  expect_match(prefreeze, "common_k = \"2;3\"", fixed = TRUE)
  expect_match(prefreeze, "comparison_units = 100L", fixed = TRUE)
  expect_match(prefreeze, "precommit_values_emitted_or_inspected = FALSE",
               fixed = TRUE)
  expect_match(prefreeze, "view_ranking_allowed = FALSE", fixed = TRUE)
  runner <- text[["run_mv11i_cross_view_agreement.R"]]
  expect_match(runner, "mv11g_cross_view_agreement_v1(gene, cell)",
               fixed = TRUE)
  expect_match(runner, "view_ranking_performed = FALSE", fixed = TRUE)
  closure <- text[["build_mv11j_cross_view_closure.R"]]
  expect_match(closure, "seed_exact_repeat", fixed = TRUE)
  expect_match(closure, "summary_exact_repeat", fixed = TRUE)
  expect_match(closure, "view_ranking_authorized = FALSE", fixed = TRUE)
})
