#' Generate toy datasets
#'
#' This function calls the script in `inst/scripts/generate_toy_data.R`
#' to regenerate the toy sparse expression matrices and accompanying
#' metadata located under `inst/extdata/toy`.
#'
#' @return Invisibly returns TRUE on success.
#' @export
generate_toy_data <- function() {
  script <- system.file("scripts", "generate_toy_data.R", package = "scPHcompare")
  if (script == "") {
    stop("generate_toy_data.R script not found")
  }
  source(script, local = FALSE)
  invisible(TRUE)
}
