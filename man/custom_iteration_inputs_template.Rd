% Generated manually to expose custom iteration helpers
\name{load_custom_iteration_inputs_template}
\alias{load_custom_iteration_inputs_template}
\alias{import_custom_iteration_inputs}
\title{Load or import custom iteration inputs}
\description{Helper functions for loading the custom iteration input template shipped with the package or importing custom overrides from a list or file path.}
\usage{
load_custom_iteration_inputs_template()

import_custom_iteration_inputs(custom_iteration_inputs)
}
\arguments{
  \item{custom_iteration_inputs}{Either a named list of iteration overrides or a single file path pointing to an R script that defines a \code{custom_iteration_inputs} or \code{custom_iteration_inputs_template} object.}
}
\value{A named list of custom iteration overrides, or \code{NULL} if none are available.}
\examples{
# Load overrides defined in the package template (if any paths are populated)
load_custom_iteration_inputs_template()

# Provide a named list directly
import_custom_iteration_inputs(list(Raw = list(ph_list_path = "./raw_pd_list.rds")))

# Or point to a template file containing the list
# import_custom_iteration_inputs("./custom_iteration_inputs_template.R")
}
