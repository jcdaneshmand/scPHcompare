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
    generate_visualizations_for_iteration = function(seurat_obj, ...) seurat_obj,
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

test_that("assignRandomGroup is reproducible and balanced by sample identity", {
  counts <- Matrix::Matrix(
    matrix(
      1,
      nrow = 2,
      ncol = 8,
      dimnames = list(c("gene-a", "gene-b"), paste0("cell-", seq_len(8)))
    ),
    sparse = TRUE
  )
  obj <- Seurat::CreateSeuratObject(counts = counts, project = "fixture")
  obj@meta.data$orig.ident <- rep(c("sample-a", "sample-b", "sample-c", "sample-d"), each = 2)

  set.seed(987)
  rng_before <- .Random.seed
  grouped_one <- scPHcompare:::assignRandomGroup(obj, k = 2L, seed = 123L)
  grouped_two <- scPHcompare:::assignRandomGroup(obj, k = 2L, seed = 123L)

  expect_identical(.Random.seed, rng_before)
  expect_identical(grouped_one@meta.data$Random_Group,
                   grouped_two@meta.data$Random_Group)

  per_sample <- split(
    grouped_one@meta.data$Random_Group,
    grouped_one@meta.data$orig.ident
  )
  expect_true(all(vapply(per_sample, function(x) length(unique(x)) == 1L, logical(1))))
  group_sizes <- table(vapply(per_sample, `[[`, integer(1), 1L))
  expect_lte(max(group_sizes) - min(group_sizes), 1L)
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

test_that("cross-iteration comparison rejects an empty iteration set clearly", {
  expect_error(
    scPHcompare:::cross_iteration_comparison_with_betti(
      data_iterations = list(),
      output_folder = tempfile("cross-iteration-")
    ),
    "at least one iteration"
  )
})

test_that("spectral clustering constructs similarities from distances", {
  distance_matrix <- as.matrix(stats::dist(matrix(
    c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE
  )))
  rownames(distance_matrix) <- colnames(distance_matrix) <- c("a", "b", "c")

  clusters <- scPHcompare:::perform_spectral_clustering(distance_matrix, 2L)

  expect_length(clusters, 3L)
  expect_identical(names(clusters), c("a", "b", "c"))
  expect_error(
    scPHcompare:::perform_spectral_clustering(distance_matrix, 3L),
    "smaller than the matrix dimension"
  )
})
