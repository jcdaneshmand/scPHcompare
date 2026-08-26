test_that("MV14 full cell landscapes remain independently closed and immutable", {
  root <- testthat::test_path("..", "..")
  audit <- file.path(root, "docs", "audits",
                     "mv14-cell-landscape-closure-v1")
  read_at <- function(name) read.csv(
    file.path(audit, name), stringsAsFactors = FALSE, check.names = FALSE
  )
  groups <- read_at("mv14-group-closure.csv")
  oracles <- read_at("mv14-exact-R-oracles.csv")
  resource <- read_at("mv14-resource-closure.csv")
  validation <- read_at("mv14-validation.csv")
  decision <- read_at("mv14-decision.csv")
  manifest <- read_at("mv14-artifact-manifest.csv")

  expect_equal(nrow(validation), 30L)
  expect_true(all(validation$passed))
  expect_equal(nrow(groups), 14L)
  expect_equal(sum(groups$chunks), 314L)
  expect_equal(sum(groups$pairs), 76372L)
  expect_equal(as.integer(table(groups$homology_dimension)), c(7L, 7L))
  expect_true(all(groups$artifact_rehashed))
  expect_true(all(groups$all_chunk_hashes))
  expect_true(all(groups$all_pair_identities))
  expect_true(all(groups$all_active_depths))
  expect_true(all(groups$all_exact))
  expect_true(all(groups$exact_R_oracle))
  expect_false(any(groups$labels_used))
  expect_false(any(groups$outcomes_used))
  expect_equal(sum(groups$downstream_jobs), 0L)

  expect_equal(nrow(oracles), 14L)
  expect_true(all(oracles$passed))
  expect_lte(max(abs(oracles$R_squared_distance -
                       oracles$Rust_squared_distance)), 2e-15)
  expect_equal(max(abs(oracles$Rust_squared_distance -
                         oracles$stored_squared_distance)), 0)

  expect_equal(resource$GNU_time_exit_status, 0)
  expect_equal(resource$runner_stderr_bytes, 0)
  expect_equal(resource$workers, 1L)
  expect_equal(resource$retries, 0L)
  expect_lte(resource$aggregate_child_seconds,
             resource$aggregate_elapsed_cap_seconds)
  expect_lte(resource$peak_child_tree_rss_bytes,
             resource$child_rss_cap_bytes)
  expect_lte(resource$private_bytes, resource$private_cap_bytes)
  expect_lte(resource$public_bytes, resource$public_cap_bytes)

  expect_true(decision$full_cell_landscapes_independently_closed)
  expect_true(decision$comparison_prefreeze_eligible_next)
  expect_false(decision$comparisons_authorized_by_this_closure)
  expect_false(decision$clustering_authorized)
  expect_false(decision$fusion_authorized)
  expect_false(decision$labels_authorized)
  expect_false(decision$outcomes_authorized)
  expect_false(decision$manuscript_claims_authorized)

  files <- file.path(audit, manifest$artifact)
  expect_equal(as.numeric(file.info(files)$size), manifest$bytes)
  expect_equal(unname(vapply(files, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)
})
