test_that("MV8-ZY prospectively binds all comparison inputs and firewalls", {
  root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zy-distance-comparison-execution-prefreeze-v1"
  )
  manifest <- read.csv(file.path(root, "mv08zy-artifact-manifest.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  paths <- file.path(root, manifest$artifact)
  expect_true(all(file.exists(paths)))
  expect_equal(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes))
  expect_equal(unname(vapply(paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), manifest$sha256)

  contract <- read.csv(file.path(root, "mv08zy-contract.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  stacks <- read.csv(file.path(root, "mv08zy-stack-bindings.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  queue <- read.csv(file.path(root, "mv08zy-comparison-queue.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  validation <- read.csv(file.path(root, "mv08zy-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  decision <- read.csv(file.path(root, "mv08zy-decision.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  implementation <- read.csv(file.path(root,
                                       "mv08zy-implementation-bindings.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(nrow(stacks), 38L)
  expect_equal(nrow(queue), 40L)
  expect_identical(as.integer(queue$comparison_order), 1:40)
  expect_equal(sum(queue$dataset_scope == "internal124"), 30L)
  expect_equal(sum(queue$dataset_scope == "external8"), 10L)
  expect_equal(sum(queue$homology_dimension == "H0"), 20L)
  expect_equal(sum(queue$homology_dimension == "H1"), 20L)
  expect_true(all(stacks$source_rehashed))
  expect_true(all(stacks$available))
  expect_true(all(nzchar(stacks$payload_set_sha256)))
  expect_true(all(nzchar(stacks$pair_axis_sha256)))
  expect_true(all(nzchar(queue$pair_axis_sha256)))
  expect_true(all(queue$input_ready))
  expect_true(all(queue$outcome_label_state == "closed"))
  expect_false(any(queue$biological_outcomes_computed))
  expect_true(all(queue$clustering_jobs == 0L))
  expect_true(all(queue$fusion_jobs == 0L))
  expect_equal(contract$workers, 1L)
  expect_equal(contract$retries, 0L)
  expect_true(contract$strict_prefix_resume)
  expect_true(contract$independent_recomputation_required)
  expect_true(all(validation$passed))
  expect_equal(nrow(validation), 20L)
  expect_true(decision$execution_authorized_after_commit)
  unchanged <- implementation$file !=
    "scripts/build_mv08zz_comparison_closure.R"
  implementation_paths <- testthat::test_path(
    "..", "..", implementation$file[unchanged]
  )
  expect_true(all(file.exists(implementation_paths)))
  expect_equal(unname(vapply(implementation_paths, function(path) digest::digest(
    file = path, algo = "sha256", serialize = FALSE
  ), character(1L))), implementation$sha256[unchanged])
  recovery_root <- testthat::test_path(
    "..", "..", "docs", "audits",
    "mv08zza-closure-serialization-recovery-prefreeze-v1"
  )
  recovery <- read.csv(file.path(recovery_root,
                                 "mv08zza-implementation-bindings.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  closure <- recovery$file == "scripts/build_mv08zz_comparison_closure.R"
  expect_equal(sum(closure), 1L)
  closure_path <- testthat::test_path("..", "..", recovery$file[closure])
  expect_equal(digest::digest(file = closure_path, algo = "sha256",
                              serialize = FALSE), recovery$sha256[closure])
})
