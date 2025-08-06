#' Generate toy datasets
#'
#' This function runs the helper script located at
#' `inst/scripts/generate_toy_data.R` to recreate 20 toy sparse
#' expression matrices and accompanying metadata under
#' `inst/extdata/toy`. The datasets span five tissues, two sequencing
#' approaches, and two SRA identifiers.
#'
#' @return A named list with `matrices` (a character vector of paths to
#'   the generated `.sparse.RData` files) and `metadata` (the path to the
#'   generated `metadata.csv`).
#' @export
generate_toy_data <- function() {
  script <- system.file("scripts", "generate_toy_data.R", package = "scPHcompare")
  if (script == "") {
    stop("generate_toy_data.R script not found")
  }

  env <- new.env(parent = emptyenv())
  sys.source(script, envir = env)

  files <- list(
    matrices = file.path(env$out_dir, paste0(env$samples$sample, ".sparse.RData")),
    metadata = env$metadata_path
  )

  files
}
