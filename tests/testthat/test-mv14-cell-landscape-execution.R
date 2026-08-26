test_that("MV14 pair identities are deterministic and dimension-specific", {
  binding <- data.frame(
    axis_order = 1:4,
    unit_id = paste0("unit", 1:4),
    diagram_sha256 = sprintf("%064d", 1:4),
    stringsAsFactors = FALSE
  )
  h0 <- .mv14_group_pairs(binding, "cohort__panel__seed__H0")
  h1 <- .mv14_group_pairs(binding, "cohort__panel__seed__H1")
  expect_equal(nrow(h0), 6L)
  expect_equal(h0$pair_ordinal, 1:6)
  expect_false(any(h0$pair_identity_sha256 == h1$pair_identity_sha256))
  expect_identical(h0, .mv14_group_pairs(binding, "cohort__panel__seed__H0"))
})

test_that("MV14 freezes the complete label-closed cell landscape surface", {
  expect_equal(10L * choose(124L, 2L) + 4L * choose(8L, 2L), 76372L)
  expect_equal(10L * ceiling(choose(124L, 2L) / 250L) +
                 4L * ceiling(choose(8L, 2L) / 250L), 314L)
  root <- testthat::test_path("..", "..")
  files <- c(
    "build_mv14_cell_landscape_prefreeze.R",
    "run_mv14_cell_landscape_chunk.R",
    "run_mv14_cell_landscape_production.R",
    "build_mv14_cell_landscape_closure.R"
  )
  text <- setNames(lapply(files, function(file) paste(readLines(
    file.path(root, "scripts", file), warn = FALSE), collapse = "\n"
  )), files)
  freeze <- text[[1L]]; worker <- text[[2L]]
  runner <- text[[3L]]; close <- text[[4L]]
  expect_match(freeze, "all_consecutive_active_levels", fixed = TRUE)
  expect_match(freeze, "no_universal_fixed_grid", fixed = TRUE)
  expect_match(freeze, "no_universal_level_cap", fixed = TRUE)
  expect_match(freeze, "required_exact_R_oracles = 14L", fixed = TRUE)
  expect_match(worker, "landscape_rust_prototype_dimension", fixed = TRUE)
  expect_match(worker, "exact = TRUE, all_active_levels = TRUE", fixed = TRUE)
  expect_match(runner, "exact completed prefix", fixed = TRUE)
  expect_match(runner, "retries = 0L", fixed = TRUE)
  expect_match(close, "landscape_reference_exact_dimension", fixed = TRUE)
  expect_match(close, "exact_pair_ids == 76372L", fixed = TRUE)
  for (value in text[1:3]) {
    expect_match(value, "outcome_label_state = \"closed\"", fixed = TRUE)
  }
  expect_match(close, "labels_authorized = FALSE", fixed = TRUE)
  expect_match(close, "manuscript_claims_authorized = FALSE", fixed = TRUE)
})

test_that("MV14 standalone scripts remain parseable", {
  root <- testthat::test_path("..", "..")
  files <- file.path(root, "scripts", c(
    "build_mv14_cell_landscape_prefreeze.R",
    "run_mv14_cell_landscape_chunk.R",
    "run_mv14_cell_landscape_production.R",
    "build_mv14_cell_landscape_closure.R"
  ))
  for (file in files) expect_silent(parse(file))
})
