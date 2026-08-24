test_that("MV8-ZO freezes the kernel failure before any repair", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zo-landscape-kernel-failure-prefreeze-v1"
  )
  expect_true(dir.exists(root))
  manifest <- read.csv(file.path(root, "mv08zo-artifact-manifest.csv"),
                       check.names = FALSE)
  checks <- read.csv(file.path(root, "mv08zo-validation.csv"),
                     check.names = FALSE)
  decision <- read.csv(file.path(root, "mv08zo-decision.csv"),
                       check.names = FALSE)
  exposure <- read.csv(file.path(root, "mv08zo-prefix-exposure.csv"),
                       check.names = FALSE)
  fixture <- read.csv(file.path(root, "mv08zo-public-synthetic-fixture.csv"),
                      check.names = FALSE)
  expect_equal(nrow(checks), 24L)
  expect_true(all(checks$passed))
  expect_equal(nrow(manifest), 8L)
  expect_equal(decision$stopped_prefix_chunks_preserved, 502L)
  expect_false(decision$stopped_prefix_scientifically_accepted)
  expect_false(decision$old_root_resume_authorized)
  expect_false(decision$old_root_delete_authorized)
  expect_true(decision$candidate_rust_source_change_authorized)
  expect_false(decision$prior_output_reuse_authorized)
  expect_false(decision$fresh_production_authorized)
  expect_equal(exposure$known_prior_pairs_using_failed_diagram, 13L)
  expect_false(exposure$active_level_guard_sufficient)
  expect_equal(nrow(fixture), 4L)
  expect_equal(fixture$expected_active_levels[[1L]], 3L)
})

test_that("MV8-ZO builder cannot resume or destroy production", {
  path <- testthat::test_path(
    "..", "..", "scripts",
    "build_mv08zo_landscape_kernel_failure_prefreeze.R"
  )
  expect_silent(parse(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl("processx::process\\$new|system2\\(", text,
                     perl = TRUE))
  expect_false(grepl("unlink\\(|file.remove\\(", text, perl = TRUE))
  expect_match(text, "old_root_resume_authorized = FALSE", fixed = TRUE)
  expect_match(text, "fresh_production_authorized = FALSE", fixed = TRUE)
})
