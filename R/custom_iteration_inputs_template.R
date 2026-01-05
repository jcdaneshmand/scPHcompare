# Helpers for loading custom iteration inputs from a template file

#' Load the default custom iteration input template
#'
#' @return A named list of iteration overrides if the template exists, otherwise `NULL`.
#' @examples
#' # Load overrides defined in the package template (if any paths are populated)
#' load_custom_iteration_inputs_template()
#' @export
load_custom_iteration_inputs_template <- function() {
  template_path <- system.file("extdata", "custom_iteration_inputs_template.R", package = "scPHcompare")
  if (!nzchar(template_path) || !file.exists(template_path)) {
    return(NULL)
  }
  env <- new.env(parent = emptyenv())
  sys.source(template_path, envir = env)
  candidate_names <- c("custom_iteration_inputs_template", "custom_iteration_inputs")
  for (name in candidate_names) {
    if (exists(name, envir = env, inherits = FALSE)) {
      inputs <- get(name, envir = env, inherits = FALSE)
      has_existing_paths <- any(vapply(unlist(inputs, recursive = FALSE, use.names = FALSE), function(path) {
        is.character(path) && length(path) == 1 && nzchar(path) && file.exists(path)
      }, logical(1)))
      return(if (has_existing_paths) inputs else NULL)
    }
  }
  NULL
}

#' Import custom iteration inputs from a list or file path
#'
#' @param custom_iteration_inputs Either a named list of iteration overrides or a
#'   single file path pointing to an R script that defines a
#'   `custom_iteration_inputs` or `custom_iteration_inputs_template` object.
#'
#' @return A named list of custom iteration overrides, or `NULL` if none are
#'   available.
#' @examples
#' # Provide a named list directly
#' import_custom_iteration_inputs(list(Raw = list(ph_list_path = "./raw_pd_list.rds")))
#'
#' # Or point to a template file containing the list
#' # import_custom_iteration_inputs("./custom_iteration_inputs_template.R")
#' @export
import_custom_iteration_inputs <- function(custom_iteration_inputs) {
  if (is.null(custom_iteration_inputs)) {
    return(load_custom_iteration_inputs_template())
  }

  if (is.character(custom_iteration_inputs) && length(custom_iteration_inputs) == 1) {
    if (!file.exists(custom_iteration_inputs)) {
      stop(sprintf("Custom iteration input file does not exist: %s", custom_iteration_inputs))
    }
    env <- new.env(parent = emptyenv())
    sys.source(custom_iteration_inputs, envir = env)
    candidate_names <- c("custom_iteration_inputs", "custom_iteration_inputs_template")
    for (name in candidate_names) {
      if (exists(name, envir = env, inherits = FALSE)) {
        return(get(name, envir = env, inherits = FALSE))
      }
    }
    stop("No custom iteration inputs were found in the provided file. Define `custom_iteration_inputs` or `custom_iteration_inputs_template`.")
  }

  custom_iteration_inputs
}
