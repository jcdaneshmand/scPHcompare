test_that("MV11-H prospectively freezes symmetric common-K agreement", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits", "mv11h-cross-view-prefreeze-v1")
  read_audit <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_audit("mv11h-contract.csv")
  sources <- read_audit("mv11h-source-bindings.csv")
  implementation <- read_audit("mv11h-implementation-bindings.csv")
  validation <- read_audit("mv11h-validation.csv")
  decision <- read_audit("mv11h-decision.csv")
  manifest <- read_audit("mv11h-artifact-manifest.csv")
  expect_equal(contract$execution_head,
               "2f87770ec9477f082936bd8d2aa6e687923a855d")
  expect_equal(c(contract$selected_rows_per_view, contract$comparison_units,
                 contract$public_seed_rows, contract$public_summary_rows),
               c(12400, 100, 100, 20))
  expect_equal(contract$common_k, "2;3")
  expect_equal(contract$metric, "adjusted_rand_index")
  expect_true(contract$precommit_real_source_dry_run)
  expect_false(contract$precommit_values_emitted_or_inspected)
  expect_equal(nrow(sources), 2L)
  source_paths <- file.path(root, sources$private_path_not_published)
  expect_equal(as.numeric(file.info(source_paths)$size), sources$bytes)
  expect_equal(unname(vapply(source_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), sources$sha256)
  paths <- file.path(root, implementation$file)
  expect_equal(as.numeric(file.info(paths)$size), implementation$bytes)
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256)
  expect_equal(nrow(validation), 21L)
  expect_true(all(validation$passed))
  expect_true(decision$fixed_execution_authorized_after_commit)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$view_ranking_authorized)
  expect_false(decision$fusion_authorized)
  artifacts <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(artifacts)$size), manifest$bytes)
  expect_equal(unname(vapply(artifacts, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
