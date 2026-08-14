profiled_diagram <- function(h0_shift = 0, h1_shift = 0) {
  h0 <- cbind(0, seq(0, 3.82, length.out = 383) + h0_shift,
              seq(0, 3.82, length.out = 383) + h0_shift + 0.5)
  h1 <- cbind(1, seq(0, 0.78, length.out = 79) + h1_shift,
              seq(0, 0.78, length.out = 79) + h1_shift + 0.2)
  result <- rbind(h0, h1)
  colnames(result) <- c("dimension", "birth", "death")
  result
}

corrected_control <- function(max_wall_seconds = 120, max_pairs = 3L) list(
  contract_id = "scph_corrected_landscape_workflow_control_v1",
  enabled = TRUE, max_wall_seconds = max_wall_seconds, max_pairs = max_pairs
)

test_that("corrected landscape control is strict and default-off", {
  expect_null(.validate_corrected_landscape_control_v1(NULL))
  expect_error(.validate_corrected_landscape_control_v1(list()),
               "must be NULL or a uniquely named list")
  expect_error(.validate_corrected_landscape_control_v1(c(
    corrected_control(), list(method = "adaptive")
  )), "cannot be overridden")
  expect_error(.validate_corrected_landscape_control_v1(c(
    corrected_control(), list(workers = 2L)
  )), "workers=1")
  expect_error(.validate_corrected_landscape_control_v1(c(
    corrected_control(), list(max_rss_bytes = 1024^3)
  )), "1.5 GiB")
  value <- .validate_corrected_landscape_control_v1(corrected_control())
  expect_identical(value$method, "auto")
  expect_identical(value$exact_max_intervals, 500L)
  expect_identical(value$abs_tol, 1e-8)
  expect_identical(value$downstream_use, "artifacts_only")
})

test_that("resource admission refuses budgets and unprofiled inputs", {
  diagrams <- list(a = profiled_diagram(), b = profiled_diagram(0.01))
  control <- .validate_corrected_landscape_control_v1(corrected_control())
  plan <- .corrected_landscape_resource_plan_v1(diagrams, control, "iter")
  expect_equal(nrow(plan$pairs), 1L)
  expect_equal(plan$planned_wall_seconds, 60)
  expect_error(.corrected_landscape_resource_plan_v1(
    diagrams, .validate_corrected_landscape_control_v1(
      corrected_control(max_wall_seconds = 59)
    ), "iter"), "wall time")
  expect_error(.corrected_landscape_resource_plan_v1(
    diagrams, .validate_corrected_landscape_control_v1(
      corrected_control(max_pairs = 0L)
    ), "iter"), "positive integer")
  expect_error(.corrected_landscape_resource_plan_v1(
    list(a = matrix(c(0, 0, 1), ncol = 3)), control, "iter"),
    "outside the profiled")
})

test_that("additive artifacts are exact, resumable, and completion-last", {
  root <- tempfile("corrected-landscape-")
  dir.create(root)
  pd_path <- file.path(root, "pd.rds")
  diagrams <- list(
    sample_b = profiled_diagram(0.02, 0.01),
    sample_a = profiled_diagram(),
    sample_c = profiled_diagram(0.04, 0.02)
  )
  saveRDS(diagrams, pd_path)
  control <- corrected_control(max_wall_seconds = 120, max_pairs = 3L)
  expect_error(produce_corrected_landscape_artifacts_v1(
    pd_path, "Toy Iter", root, control, log_message = function(...) NULL,
    stop_after_pairs = 1L
  ), "intentional interruption")
  artifact_dirs <- list.dirs(file.path(root, "corrected_landscape_v1"),
                             recursive = FALSE, full.names = TRUE)
  expect_length(artifact_dirs, 1L)
  expect_false(file.exists(file.path(artifact_dirs, "completion-v1.csv")))
  expect_length(list.files(file.path(artifact_dirs, "pairs"), pattern = "rds$"), 1L)

  result <- produce_corrected_landscape_artifacts_v1(
    pd_path, "Toy Iter", root, control, log_message = function(...) NULL
  )
  expect_false(result$resumed)
  expect_identical(result$downstream_use, "artifacts_only")
  expect_true(file.exists(file.path(result$artifact_dir, "completion-v1.csv")))
  matrix_value <- readRDS(file.path(result$artifact_dir, "distance-matrix-v1.rds"))
  direct <- persistence_landscape_distance_matrix(diagrams, method = "auto")
  expect_identical(matrix_value$matrices, direct$matrices)
  expect_identical(matrix_value$cache_key, direct$cache_key)
  expect_identical(matrix_value$sample_ids, sort(names(diagrams)))
  expect_true(all(vapply(matrix_value$matrices, function(x) {
    isTRUE(all.equal(x, t(x))) && all(diag(x) == 0)
  }, logical(1L))))
  expect_false(matrix_value$provenance$legacy_reproduction)
  expect_false(matrix_value$provenance$existing_workflow_default_changed)

  before <- file.info(list.files(result$artifact_dir, recursive = TRUE,
                                 full.names = TRUE))
  resumed <- produce_corrected_landscape_artifacts_v1(
    pd_path, "Toy Iter", root, control, log_message = function(...) NULL
  )
  after <- file.info(rownames(before))
  expect_true(resumed$resumed)
  expect_identical(before$size, after$size)
  expect_identical(before$mtime, after$mtime)

  pair_index <- read.csv(file.path(result$artifact_dir, "pair-index-v1.csv"),
                         stringsAsFactors = FALSE)
  pair_path <- file.path(result$artifact_dir, pair_index$pair_artifact[[1L]])
  bad <- readRDS(pair_path)
  bad$distances[["H0"]] <- bad$distances[["H0"]] + 1
  saveRDS(bad, pair_path)
  expect_error(produce_corrected_landscape_artifacts_v1(
    pd_path, "Toy Iter", root, control, log_message = function(...) NULL
  ), "hash validation")
})

test_that("workflow signatures add only null corrected control", {
  expect_null(formals(run_unified_pipeline)$corrected_landscape_control)
  expect_null(formals(run_postprocessing_pipeline)$corrected_landscape_control)
  expect_false("corrected_landscape_control" %in% names(formals(
    process_iteration_calculate_matrices
  )))
  expect_false("corrected_landscape_control" %in% names(formals(
    run_modular_analysis
  )))
  expect_false("corrected_landscape_control" %in% names(formals(
    run_cross_iteration
  )))
})

test_that("postprocessing corrected-only mode returns sidecar without legacy use", {
  root <- tempfile("corrected-postprocessing-")
  dir.create(root)
  pd_path <- file.path(root, "pd.rds")
  saveRDS(list(sample_a = profiled_diagram(),
               sample_b = profiled_diagram(0.01, 0.01)), pd_path)
  iteration <- list(
    name = "Corrected Only", pd_list = pd_path,
    expr_list = list(sample_a = matrix(1), sample_b = matrix(2))
  )
  result <- run_postprocessing_pipeline(
    list(data_iterations = list(iteration)), results_dir = root,
    run_standard_seurat_clustering = FALSE,
    run_kmeans_clustering = FALSE,
    run_hierarchical_ph_clustering = FALSE,
    run_spectral_clustering = FALSE,
    run_visualizations = FALSE,
    run_sample_level_heatmap = FALSE,
    run_cluster = FALSE, run_betti = FALSE, run_cross_iteration = FALSE,
    metadata_path = NULL,
    corrected_landscape_control = corrected_control(
      max_wall_seconds = 60, max_pairs = 1L
    )
  )
  sidecar <- result$data_iterations[[1L]]$corrected_landscape_v1
  expect_identical(sidecar$downstream_use, "artifacts_only")
  expect_true(file.exists(file.path(sidecar$artifact_dir, "completion-v1.csv")))
  expect_null(result$data_iterations[[1L]]$landscape_l2_distance_matrix)
  expect_null(result$data_iterations[[1L]]$landscape_list)
})
