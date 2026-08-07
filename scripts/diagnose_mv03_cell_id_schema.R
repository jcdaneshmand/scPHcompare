#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: diagnose_mv03_cell_id_schema.R <sct-rds> <integrated-rds> <sample-id>",
    call. = FALSE
  )
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required to inspect historical sparse objects.", call. = FALSE)
}

inspect_entry <- function(path, sample_id, label) {
  values <- readRDS(path)
  if (!is.list(values) || !(sample_id %in% names(values))) {
    stop(label, " does not contain requested sample: ", sample_id, call. = FALSE)
  }
  value <- values[[sample_id]]
  cat(label, " class: ", paste(class(value), collapse = ";"), "\n", sep = "")
  cat(label, " typeof: ", typeof(value), "\n", sep = "")
  cat(label, " dim: ", paste(dim(value), collapse = "x"), "\n", sep = "")
  if (methods::is(value, "Matrix")) {
    cat(label, " Dim slot: ", paste(methods::slot(value, "Dim"), collapse = "x"),
        "\n", sep = "")
  }
  cells <- colnames(value)
  genes <- rownames(value)
  cat(label, " gene head:\n", sep = "")
  dput(utils::head(genes, 3L))
  cat(label, " cell head:\n", sep = "")
  dput(utils::head(cells, 5L))
  cat(label, " cell tail:\n", sep = "")
  dput(utils::tail(cells, 5L))
  result <- list(cells = cells, genes = genes)
  rm(values, value)
  invisible(gc())
  result
}

sct <- inspect_entry(args[[1L]], args[[3L]], "sct")
integrated <- inspect_entry(args[[2L]], args[[3L]], "integrated")

candidate_normalizations <- list(
  exact = identity,
  remove_trailing_numeric_suffix = function(x) sub("_[0-9]+$", "", x),
  remove_leading_numeric_prefix = function(x) sub("^[0-9]+_", "", x),
  remove_first_field = function(x) sub("^[^_]+_", "", x),
  remove_last_field = function(x) sub("_[^_]+$", "", x),
  barcode_before_first_underscore = function(x) sub("_.*$", "", x),
  barcode_after_last_underscore = function(x) sub("^.*_", "", x)
)

for (name in names(candidate_normalizations)) {
  normalizer <- candidate_normalizations[[name]]
  first <- normalizer(sct$cells)
  second <- normalizer(integrated$cells)
  cat(
    "candidate ", name,
    ": shared=", length(intersect(first, second)),
    ", first_unique=", length(unique(first)),
    ", second_unique=", length(unique(second)),
    ", ordered_identical=", identical(first, second),
    "\n",
    sep = ""
  )
}
