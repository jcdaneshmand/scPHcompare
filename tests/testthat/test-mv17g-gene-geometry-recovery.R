test_that("MV17-G v2 preserves cell geometry and enforces gene chord geometry", {
  set.seed(17017L)
  x <- matrix(rnorm(60L), 20L, 3L)
  rownames(x) <- sprintf("p%02d", seq_len(nrow(x)))

  cell <- mv17g_topology_coordinates_v2(x, "cell")
  gene <- mv17g_topology_coordinates_v2(x, "gene")
  expect_identical(cell, x)
  expect_lte(max(abs(rowMeans(gene))), 1e-12)
  expect_lte(max(abs(rowSums(gene ^ 2) - 1)), 1e-12)
  expect_gt(
    max(abs(as.matrix(stats::dist(x)) - as.matrix(stats::dist(gene)))),
    1e-3
  )
  correlation_chord <- sqrt(pmax(0, 2 * (1 - stats::cor(t(x)))))
  expect_equal(
    as.vector(as.matrix(stats::dist(gene))),
    as.vector(correlation_chord),
    tolerance = 1e-12
  )
})

test_that("MV17-G v2 exactly preserves admitted cell metrics", {
  set.seed(17018L)
  x <- matrix(rnorm(60L), 20L, 3L)
  legacy <- mv17g_group_metrics_v1(
    x, "coordinate_permutation", c(51L, 52L)
  )
  corrected <- mv17g_group_metrics_v2(
    x, "cell", "coordinate_permutation", c(51L, 52L)
  )
  expect_identical(
    legacy,
    corrected[setdiff(names(corrected), "geometry")]
  )
  expect_identical(unique(corrected$geometry), "euclidean_shared_pca_v1")
})

test_that("MV17-G v2 gene metrics are deterministic and fail closed", {
  set.seed(17019L)
  x <- matrix(rnorm(60L), 20L, 3L)
  first <- mv17g_group_metrics_v2(
    x, "gene", "coordinate_permutation", c(61L, 62L)
  )
  repeated <- mv17g_group_metrics_v2(
    x, "gene", "coordinate_permutation", c(61L, 62L)
  )
  observed <- mv17g_group_metrics_v2(x, "gene", "observed", 0L)
  legacy_observed <- mv17g_group_metrics_v1(x, "observed", 0L)

  expect_identical(first, repeated)
  expect_equal(nrow(first), 16L)
  expect_equal(nrow(observed), 8L)
  expect_identical(
    unique(first$geometry), "euclidean_correlation_chord_v1"
  )
  expect_false(identical(observed$value, legacy_observed$value))
  expect_setequal(first$summary_id, mv17c_summary_registry_v1()$summary_id)
  expect_setequal(first$homology_dimension, c("H0", "H1"))
  expect_true(all(first$h0_mst_maximum_absolute_error <= 1e-8))
  expect_error(
    mv17g_group_metrics_v2(x, "cell", "within_row_axis_shuffle", 61L),
    "incompatible"
  )
  constant <- x
  constant[1L, ] <- 1
  expect_error(
    mv17g_topology_coordinates_v2(constant, "gene"), "constant row"
  )
})

test_that("MV17-G worker v2 promotes a typed atomic synthetic payload", {
  root <- normalizePath(testthat::test_path("..", ".."))
  directory <- tempfile("mv17g-worker-v2-")
  dir.create(directory)
  matrix_path <- file.path(directory, "matrix.rds")
  output <- file.path(directory, "result.rds")
  set.seed(17020L)
  saveRDS(matrix(rnorm(60L), 20L, 3L), matrix_path, version = 3)
  status <- withr::with_dir(root, system2(
    "Rscript",
    c(
      "--vanilla", file.path(
        root, "scripts", "run_mv17g_calibration_group_worker_v2.R"
      ), matrix_path, "gene", "coordinate_permutation", "301", "2", output
    ), stdout = TRUE, stderr = TRUE
  ))
  exit <- attr(status, "status")
  if (is.null(exit)) exit <- 0L
  expect_equal(exit, 0L, info = paste(status, collapse = "\n"))
  expect_true(file.exists(output))
  expect_false(file.exists(paste0(output, ".partial")))
  result <- readRDS(output)
  expect_identical(result$contract_id, "mv17g_group_result_v2")
  expect_identical(result$view, "gene")
  expect_identical(result$geometry, "euclidean_correlation_chord_v1")
  expect_equal(result$replicate_count, 2L)
  expect_equal(c(result$seed_first, result$seed_last), c(301L, 302L))
  expect_equal(nrow(result$metrics), 16L)
  expect_true(result$finite)
  expect_false(result$labels_opened)
  expect_false(result$outcomes_opened)
})

test_that("MV17-G recovery prefix is relative to the frozen subqueue", {
  queue <- data.frame(
    job_order = 529:532,
    view = "gene",
    unit_order = 1L,
    null_family = "observed"
  )
  scan <- data.frame(
    job_order = 529:532,
    state = c("complete", "complete", "absent", "absent")
  )
  expect_identical(mv17g_complete_queue_prefix_v2(scan, queue), 2L)
  bad <- scan
  bad$state <- c("complete", "absent", "complete", "absent")
  expect_error(
    mv17g_complete_queue_prefix_v2(bad, queue), "not a queue prefix"
  )
  partial <- scan
  partial$state[3L] <- "partial"
  expect_error(mv17g_complete_queue_prefix_v2(partial, queue), "invalid")
})

test_that("MV17-G v2 recovery wave passes the explicit view to its worker", {
  skip_if(.Platform$OS.type != "unix")
  root <- normalizePath(testthat::test_path("..", ".."))
  directory <- tempfile("mv17g-wave-v2-")
  private <- file.path(directory, "private")
  matrix_root <- file.path(directory, "matrices-root")
  dir.create(file.path(private, "jobs"), recursive = TRUE)
  dir.create(file.path(private, "logs"), recursive = TRUE)
  dir.create(file.path(matrix_root, "matrices"), recursive = TRUE)
  set.seed(17021L)
  saveRDS(
    matrix(rnorm(60L), 20L, 3L),
    file.path(matrix_root, "matrices", "gene__001.rds"), version = 3
  )
  queue <- data.frame(
    job_order = 529L, view = "gene", unit_order = 1L,
    null_family = "observed", seed_first = 0L, seed_last = 0L,
    replicate_count = 1L, scientific_runs = 1L
  )
  contract <- data.frame(
    workers = 1L, child_timeout_seconds = 120L,
    child_RSS_cap_bytes = 1024 ^ 3
  )
  result <- withr::with_dir(root, mv17g_run_parallel_wave_v2(
    queue, private, matrix_root,
    normalizePath("scripts/run_mv17g_calibration_group_worker_v2.R"),
    contract
  ))
  expect_identical(result$exit_status, 0L)
  expect_identical(result$artifacts, 4L)
  paths <- mv17g_job_artifacts_v1(queue, private)
  payload <- readRDS(paths[["result"]])
  expect_identical(payload$view, "gene")
  expect_identical(payload$geometry, "euclidean_correlation_chord_v1")
  expect_identical(payload$contract_id, "mv17g_group_result_v2")
  expect_equal(file.info(paths[["stdout"]])$size, 0)
  expect_equal(file.info(paths[["stderr"]])$size, 0)
})

test_that("MV17-G geometry recovery runner is prefreeze-gated", {
  root <- testthat::test_path("..", "..")
  runner <- file.path(root, "scripts", "run_mv17g_gene_geometry_recovery.R")
  expect_true(file.exists(runner))
  expect_silent(parse(runner))
  text <- paste(readLines(runner, warn = FALSE), collapse = "\n")
  expect_match(text, "execution_authorized_after_commit")
  expect_match(text, "gene_primary_children != 660L", fixed = TRUE)
  expect_match(text, "gene_primary_scientific_runs != 52404L", fixed = TRUE)
  expect_match(text, "mv17g_group_result_v2", fixed = TRUE)
  expect_match(text, "euclidean_correlation_chord_v1", fixed = TRUE)
  expect_match(text, "retries != 0L", fixed = TRUE)
  expect_match(text, "labels_opened = FALSE", fixed = TRUE)
  expect_match(text, "outcomes_opened = FALSE", fixed = TRUE)
})

test_that("MV17-G recovery prefreeze binds salvage and rejection partitions", {
  root <- testthat::test_path("..", "..")
  builder <- file.path(
    root, "scripts", "build_mv17g_gene_geometry_recovery_prefreeze.R"
  )
  expect_true(file.exists(builder))
  expect_silent(parse(builder))
  text <- paste(readLines(builder, warn = FALSE), collapse = "\n")
  expect_match(text, "mv17g_gene_geometry_recovery_authorized_v1", fixed = TRUE)
  expect_match(text, "artifact_prefix < 528L", fixed = TRUE)
  expect_match(text, "salvageable_cell_prefix_v1", fixed = TRUE)
  expect_match(text, "rejected_raw_euclidean_gene_evidence_v1", fixed = TRUE)
  expect_match(
    text, "identical(as.integer(gene_queue$job_order), 529:1188)", fixed = TRUE
  )
  expect_match(
    text, "sum(gene_queue$scientific_runs) != 52404L", fixed = TRUE
  )
  expect_match(text, "controller_processes != 0L", fixed = TRUE)
  expect_match(text, "worker_processes != 0L", fixed = TRUE)
  expect_match(text, "execution_only_after_prefreeze_commit = TRUE", fixed = TRUE)
  expect_match(text, "localization_authorized = FALSE", fixed = TRUE)
})
