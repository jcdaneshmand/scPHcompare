args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
setwd(normalizePath(root, mustWork = TRUE))
devtools::test(filter = "mv06-fusion", stop_on_failure = TRUE)
