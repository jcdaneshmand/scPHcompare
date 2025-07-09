run_unified_pipeline <- function(metadata_path, results_dir,
                                  run_cluster = FALSE, ...) {
  source("PH_Calculation.R")
  source("cluster_comparison.R")

  metadata <- read.csv(metadata_path)
  data_iterations <- process_datasets_PH(metadata, ...)

  if (run_cluster) {
    run_cluster_comparison(data_iterations, results_dir, ...)
  }

  invisible(data_iterations)
}
