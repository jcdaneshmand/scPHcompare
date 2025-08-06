library(testthat)
library(scPHcompare)

# Setup fake objects for mocking
fake_seurat <- list(meta.data = data.frame(Tissue="T1", SRA="S1", Approach="A1"))

tmp_pd <- tempfile(fileext = ".rds")
saveRDS(list(), tmp_pd)

fake_iter <- list(
  name = "iter1",
  seurat_obj = fake_seurat,
  assay = "RNA",
  pd_list = tmp_pd
)

fake_results <- list(
  data_iterations = list(fake_iter),
  SRA_col = "SRA",
  Tissue_col = "Tissue",
  Approach_col = "Approach"
)

# Test run_postprocessing_pipeline -------------------------------------------------

test_that("run_postprocessing_pipeline calls run_modular_analysis", {
  results_dir <- tempfile(pattern = "results_")
  modular_called <- FALSE
  captured_results <- NULL

  mock_modular <- function(ph_results, ...) {
    modular_called <<- TRUE
    captured_results <<- ph_results
  }

  mockr::with_mock(
    process_iteration_calculate_matrices = function(...) {},
    assignRandomGroup = function(obj, ...) obj,
    apply_all_clustering_methods = function(obj, ...) obj,
    generate_visualizations_for_iteration = function(obj, ...) obj,
    run_modular_analysis = mock_modular,
    {
      res <- run_postprocessing_pipeline(
        fake_results,
        results_dir = results_dir,
        num_cores = 1,
        run_standard_seurat_clustering = FALSE,
        run_kmeans_clustering = FALSE,
        run_hierarchical_ph_clustering = FALSE,
        run_spectral_clustering = FALSE,
        run_visualizations = FALSE,
        run_sample_level_heatmap = FALSE,
        run_cluster = TRUE,
        run_betti = FALSE,
        run_cross_iteration = FALSE,
        metadata_path = NULL
      )
    }
  )

  expect_true(dir.exists(results_dir))
  expect_true(modular_called)
  expect_identical(captured_results, fake_results)
  expect_identical(res, fake_results)
})

# Test run_modular_analysis --------------------------------------------------------

test_that("run_modular_analysis triggers selected modules", {
  results_dir <- tempfile(pattern = "results_")
  cluster_called <- FALSE
  betti_called <- 0
  cross_called <- FALSE

  mock_cluster <- function(data_iterations, ...) { cluster_called <<- TRUE }
  mock_betti <- function(...) { betti_called <<- betti_called + 1 }
  mock_cross <- function(data_iterations, ...) { cross_called <<- TRUE }

  mockr::with_mock(
    run_cluster_comparison = mock_cluster,
    compute_and_compare_betti_curves = mock_betti,
    run_cross_iteration = mock_cross,
    {
      run_modular_analysis(
        fake_results,
        results_dir = results_dir,
        run_cluster = TRUE,
        run_betti = TRUE,
        run_cross_iteration = TRUE
      )
    }
  )

  expect_true(dir.exists(results_dir))
  expect_true(cluster_called)
  expect_true(cross_called)
  expect_gt(betti_called, 0)
})

# Test run_cross_iteration ---------------------------------------------------------

test_that("run_cross_iteration forwards to cross_iteration_comparison_with_betti", {
  results_dir <- tempfile(pattern = "results_")
  captured_args <- NULL

  mock_cicwb <- function(data_iterations, group_by_col, output_folder, ...) {
    captured_args <<- list(data_iterations = data_iterations,
                           group_by_col = group_by_col,
                           output_folder = output_folder)
    "done"
  }

  res <- mockr::with_mock(
    cross_iteration_comparison_with_betti = mock_cicwb,
    {
      run_cross_iteration(
        fake_results$data_iterations,
        results_folder = results_dir,
        group_by_col = "Tissue"
      )
    }
  )

  expect_true(dir.exists(file.path(results_dir, "cross_iteration_comparisons")))
  expect_identical(captured_args$group_by_col, "Tissue")
  expect_identical(res, "done")
})
