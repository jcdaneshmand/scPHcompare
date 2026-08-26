test_that("MV13-E full cell-PH closure remains exact and immutable", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv13e-allqc-cell-full-closure-v1")
  read_at <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  validation <- read_at("mv13e-validation.csv")
  decision <- read_at("mv13e-decision.csv")
  repeat_rows <- read_at("mv13e-group-repeat.csv")
  resource <- read_at("mv13e-resource-closure.csv")
  manifest <- read_at("mv13e-artifact-manifest.csv")

  expect_equal(nrow(validation), 27L)
  expect_true(all(validation$passed))
  expect_true(decision$full_cell_PH_independently_closed)
  expect_true(decision$landscapes_prefreeze_eligible_next)
  expect_false(decision$landscapes_authorized_by_this_closure)
  expect_false(decision$comparisons_authorized)
  expect_false(decision$clustering_authorized)
  expect_false(decision$fusion_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)

  expect_equal(nrow(repeat_rows), 7L)
  expect_equal(sum(repeat_rows$units), 636L)
  expect_true(all(repeat_rows$axis_exact))
  expect_true(all(repeat_rows$model_rotation_exact))
  expect_equal(sum(repeat_rows$exact_view_payloads), 636L)
  expect_equal(sum(repeat_rows$exact_PH_diagrams), 636L)
  expect_equal(sum(repeat_rows$H0_MST_oracles_passed), 636L)
  expect_true(all(repeat_rows$all_exact))
  expect_false(any(repeat_rows$labels_used))
  expect_false(any(repeat_rows$outcomes_used))

  expect_equal(resource$GNU_time_exit_status, 0)
  expect_equal(resource$stderr_bytes, 0)
  expect_equal(resource$workers, 1)
  expect_equal(resource$retries, 0)
  expect_lte(resource$elapsed_seconds, resource$elapsed_cap_seconds)
  expect_lte(resource$peak_process_rss_bytes, resource$rss_cap_bytes)
  expect_lte(resource$private_bytes, resource$private_cap_bytes)
  expect_lte(resource$public_bytes, resource$public_cap_bytes)

  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
