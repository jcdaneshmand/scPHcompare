test_that("MV14 public prefreeze exactly binds the closed cell landscape queue", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv14-cell-landscape-prefreeze-v1")
  read_at <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  contract <- read_at("mv14-contract.csv")
  groups <- read_at("mv14-group-queue.csv")
  queue <- read_at("mv14-production-queue.csv")
  inputs <- read_at("mv14-input-bindings.csv")
  implementation <- read_at("mv14-implementation-bindings.csv")
  closure <- read_at("mv14-prospective-closure.csv")
  validation <- read_at("mv14-validation.csv")
  decision <- read_at("mv14-decision.csv")
  manifest <- read_at("mv14-artifact-manifest.csv")

  expect_equal(contract$execution_head,
               "12635568691705434c20b4b9fbc0791e0536b351")
  expect_equal(contract$scientific_engine_version, 2L)
  expect_equal(contract$landscape_groups, 14L)
  expect_equal(contract$production_chunks, 314L)
  expect_equal(contract$production_pairs, 76372L)
  expect_equal(contract$level_policy, "all_consecutive_active_levels")
  expect_equal(contract$grid_policy, "no_universal_fixed_grid")
  expect_equal(contract$level_cap_policy, "no_universal_level_cap")
  expect_equal(contract$workers, 1L)
  expect_equal(contract$automatic_retries, 0L)

  expect_equal(nrow(groups), 14L)
  expect_equal(as.integer(table(groups$homology_dimension)), c(7L, 7L))
  expect_equal(nrow(queue), 314L)
  expect_equal(sum(queue$pair_count), 76372L)
  expect_identical(as.integer(queue$production_order), 1:314)
  expect_lte(max(queue$pair_count), 250L)
  expect_false(any(c("unit_id", "first_unit_id", "second_unit_id") %in%
                     names(groups)))
  expect_false(any(c("unit_id", "first_unit_id", "second_unit_id") %in%
                     names(queue)))
  expect_equal(inputs$role,
               c("mv13e_manifest", "engine_v2_admission_manifest",
                 "private_axis_bindings", "engine_v2_library"))
  expect_equal(closure$required_exact_R_oracles, 14L)
  expect_equal(nrow(validation), 31L)
  expect_true(all(validation$passed))
  expect_true(decision$production_authorized_after_commit)
  expect_false(decision$comparisons_authorized)
  expect_false(decision$clustering_authorized)
  expect_false(decision$fusion_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$manuscript_claims_authorized)

  implementation_paths <- file.path(root, implementation$file)
  expect_equal(as.numeric(file.info(implementation_paths)$size),
               implementation$bytes)
  expect_equal(unname(vapply(implementation_paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))), implementation$sha256)
  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
