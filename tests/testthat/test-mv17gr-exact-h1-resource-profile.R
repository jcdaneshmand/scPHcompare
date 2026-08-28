test_that("MV17-GR cone threshold makes the Rips complex a cone", {
  set.seed(17101L)
  x <- matrix(rnorm(90L), 30L, 3L)
  contract <- mv17gr_distance_contract_v1(x)
  expect_equal(nrow(contract$distances), 30L)
  expect_true(contract$cone_threshold <= contract$diameter + 1e-12)
  expect_true(all(
    contract$distances[contract$cone_anchor_index, ] <=
      contract$cone_threshold
  ))
  expect_gt(contract$admitted_edges, 0L)
  expect_lte(contract$edge_fraction, 1)
})

test_that("MV17-GR finite cone cutoff preserves exact Ripserr H1", {
  fixtures <- list(
    random = {
      set.seed(17102L)
      matrix(rnorm(60L), 20L, 3L)
    },
    circle = cbind(
      cos(seq(0, 2 * pi, length.out = 21L)[-21L]),
      sin(seq(0, 2 * pi, length.out = 21L)[-21L])
    )
  )
  for (x in fixtures) {
    full <- mv17gr_run_exact_h1_v1(x, "ripserr_infinite_v1")
    cone <- mv17gr_run_exact_h1_v1(x, "ripserr_cone_exact_v1")
    expect_lte(mv17gr_maximum_h1_difference_v1(full$h1, cone$h1), 1e-10)
    expect_equal(full$h1_metrics$value, cone$h1_metrics$value, tolerance = 1e-10)
    expect_true(cone$exact_cone_argument)
  }
})

test_that("MV17-GR exact GUDHI and Ripserr cone H1 agree", {
  skip_if_not_installed("TDA")
  theta <- seq(0, 2 * pi, length.out = 17L)[-17L]
  x <- cbind(cos(theta), sin(theta))
  ripser <- mv17gr_run_exact_h1_v1(x, "ripserr_cone_exact_v1")
  gudhi <- mv17gr_run_exact_h1_v1(x, "gudhi_cone_exact_v1")
  expect_lte(mv17gr_maximum_h1_difference_v1(ripser$h1, gudhi$h1), 1e-6)
  expect_equal(ripser$h1_metrics$value, gudhi$h1_metrics$value, tolerance = 1e-6)
})

test_that("MV17-GR null generation precedes correlation-chord transform", {
  set.seed(17103L)
  x <- matrix(rnorm(60L), 20L, 3L)
  a <- mv17gr_materialize_gene_case_v1(x, "coordinate_permutation", 901L)
  b <- mv17gr_materialize_gene_case_v1(x, "coordinate_permutation", 901L)
  expect_identical(a, b)
  expect_lte(max(abs(rowMeans(a))), 1e-12)
  expect_lte(max(abs(rowSums(a ^ 2) - 1)), 1e-12)
  expect_error(
    mv17gr_materialize_gene_case_v1(x, "observed", 1L), "seed zero"
  )
})

test_that("MV17-GR worker atomically promotes a typed exact result", {
  root <- normalizePath(testthat::test_path("..", ".."))
  directory <- tempfile("mv17gr-worker-")
  dir.create(directory)
  matrix_path <- file.path(directory, "matrix.rds")
  output <- file.path(directory, "result.rds")
  set.seed(17104L)
  saveRDS(matrix(rnorm(60L), 20L, 3L), matrix_path, version = 3)
  status <- withr::with_dir(root, system2(
    "Rscript",
    c(
      "--vanilla", "scripts/run_mv17gr_exact_h1_profile_worker.R",
      matrix_path, "observed", "0", "ripserr_cone_exact_v1", output
    ), stdout = TRUE, stderr = TRUE
  ))
  exit <- attr(status, "status")
  if (is.null(exit)) exit <- 0L
  expect_equal(exit, 0L, info = paste(status, collapse = "\n"))
  expect_true(file.exists(output))
  expect_false(file.exists(paste0(output, ".partial")))
  result <- readRDS(output)
  expect_identical(result$contract_id, "mv17gr_exact_h1_result_v1")
  expect_identical(result$geometry, "euclidean_correlation_chord_v1")
  expect_true(result$finite)
  expect_false(result$labels_opened)
  expect_false(result$outcomes_opened)
})

test_that("MV17-GR scripts freeze profiling without production retry", {
  root <- testthat::test_path("..", "..")
  paths <- file.path(root, "scripts", c(
    "build_mv17gr_resource_profile_prefreeze.R",
    "run_mv17gr_exact_h1_profile_worker.R",
    "run_mv17gr_exact_h1_resource_profile.R",
    "build_mv17gr_resource_profile_closure.R"
  ))
  expect_true(all(file.exists(paths)))
  expect_silent(lapply(paths, parse))
  prefreeze <- paste(readLines(paths[[1L]], warn = FALSE), collapse = "\n")
  runner <- paste(readLines(paths[[3L]], warn = FALSE), collapse = "\n")
  closure <- paste(readLines(paths[[4L]], warn = FALSE), collapse = "\n")
  expect_match(prefreeze, "mv17gr_resource_profile_authorized_v1", fixed = TRUE)
  expect_match(prefreeze, "failed_production_retry_authorized = FALSE", fixed = TRUE)
  expect_match(prefreeze, "attempt_address_space_cap_bytes = 80 * 1024 ^ 3", fixed = TRUE)
  expect_match(prefreeze, "profile_queue) != 10L", fixed = TRUE)
  expect_match(runner, "prlimit", fixed = TRUE)
  expect_match(runner, "mv17g_parallel_thread_environment_v1", fixed = TRUE)
  expect_match(runner, "suppressWarnings(system2", fixed = TRUE)
  expect_match(runner, "workers != 1L", fixed = TRUE)
  expect_match(runner, "retries != 0L", fixed = TRUE)
  expect_match(closure, "production_retry_authorized = FALSE", fixed = TRUE)
  expect_match(closure, "exact_H1_infeasible", fixed = TRUE)
})

test_that("MV17-GR public audits are aggregate-only when present", {
  root <- normalizePath(testthat::test_path("..", ".."))
  audits <- c(
    "mv17gr-exact-h1-resource-profile-prefreeze-v1",
    "mv17gr-exact-h1-resource-profile-closure-v1"
  )
  for (name in audits) {
    path <- file.path(root, "docs", "audits", name)
    if (!dir.exists(path)) next
    text <- paste(
      unlist(lapply(list.files(path, full.names = TRUE), readLines,
                    warn = FALSE)), collapse = "\n"
    )
    expect_false(grepl(
      "source_job_order|unit_order|seed|artifact =|source_path|unit_id|barcode",
      text, ignore.case = TRUE
    ))
  }
})
