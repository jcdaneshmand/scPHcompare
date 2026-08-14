library(testthat)
library(scPHcompare)

# Mocking process_datasets_PH to avoid heavy computation
fake_results <- list(
  data_iterations = list(),
  SRA_col = "SRA",
  Tissue_col = "Tissue",
  Approach_col = "Approach"
)

mock_pipeline <- function(metadata_path, ...) {
  fake_results
}


test_that("run_unified_pipeline returns list structure", {
  metadata <- data.frame(`File Path` = "dummy.RData")
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(metadata, tmp)

  mockr::with_mock(process_datasets_PH = mock_pipeline, {
    result <- run_unified_pipeline(tmp, results_dir = tempfile("results_"))
  })

  expect_type(result, "list")
  expect_named(result, c(names(fake_results), "provenance"))
  expect_true(file.exists(result$provenance$pipeline_metrics_file))
  expect_equal(
    result$provenance$pipeline_metrics$stage,
    c("metadata_load", "ph_processing", "pipeline_total")
  )
  expect_true(all(c(
    "rss_before_bytes", "rss_after_bytes", "rss_delta_bytes"
  ) %in% names(result$provenance$pipeline_metrics)))
})

test_that("run_unified_pipeline creates results_dir and triggers post processing", {
  metadata <- data.frame(`File Path` = "dummy.RData")
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(metadata, tmp)

  results_dir <- tempfile(pattern = "results_")
  postprocess_called <- FALSE
  captured_results <- NULL

  mock_postprocess <- function(ph_results, ...) {
    postprocess_called <<- TRUE
    captured_results <<- ph_results
  }

  mockr::with_mock(
    process_datasets_PH = mock_pipeline,
    run_postprocessing_pipeline = mock_postprocess,
    {
      run_unified_pipeline(tmp, results_dir = results_dir, run_cluster = TRUE)
    }
  )

  expect_true(dir.exists(results_dir))
  expect_true(postprocess_called)
  expect_identical(captured_results, fake_results)
})
