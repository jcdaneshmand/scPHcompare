#!/usr/bin/env Rscript

options(warn = 2)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required to inspect historical sparse objects.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: audit_mv03_source_intermediates.R <historical-dir> <manifest-csv> ",
    "<detail-output-csv> <summary-output-csv>",
    call. = FALSE
  )
}

historical_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
detail_output <- args[[3L]]
summary_output <- args[[4L]]

manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
required_manifest <- c("cohort", "sample_id", "filtered_cells")
if (!all(required_manifest %in% names(manifest)) ||
    anyDuplicated(manifest$sample_id)) {
  stop("Pilot manifest has an incompatible schema or duplicate sample IDs.",
       call. = FALSE)
}

source_table <- data.frame(
  cohort = c("large", "large", "bone", "bone"),
  representation = c("sct_whole", "seurat_integration", "sct_whole", "integrated"),
  source_file = c(
    "expr_list_sctWhole.Rds",
    "expr_list_integrated.Rds",
    "expr_list_sctWhole_bonemarrow.Rds",
    "expr_list_integrated_bonemarrow.Rds"
  ),
  historical_extraction = c(
    "SCT_assay_data_slot",
    "integrated_assay_data_slot",
    "SCT_assay_data_slot",
    "integrated_assay_data_slot"
  ),
  stringsAsFactors = FALSE
)

matrix_values <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    return(methods::slot(x, "x"))
  }
  as.numeric(x)
}

detail_rows <- list()
source_summaries <- list()
identity <- list()

for (index in seq_len(nrow(source_table))) {
  source <- source_table[index, , drop = FALSE]
  path <- file.path(historical_dir, source$source_file)
  if (!file.exists(path)) {
    stop("Required historical intermediate is missing: ", path, call. = FALSE)
  }
  file_info <- file.info(path)
  started <- proc.time()[["elapsed"]]
  values <- readRDS(path)
  read_seconds <- proc.time()[["elapsed"]] - started
  if (!is.list(values) || is.null(names(values)) || any(!nzchar(names(values))) ||
      anyDuplicated(names(values))) {
    stop("Intermediate must be a uniquely named list: ", path, call. = FALSE)
  }
  requested <- manifest$sample_id[manifest$cohort == source$cohort]
  source_summaries[[length(source_summaries) + 1L]] <- data.frame(
    cohort = source$cohort,
    representation = source$representation,
    source_file = source$source_file,
    historical_extraction = source$historical_extraction,
    file_size_bytes = as.numeric(file_info$size),
    file_modified_utc = format(
      file_info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    ),
    read_seconds = read_seconds,
    list_samples = length(values),
    unique_list_sample_ids = length(unique(names(values))),
    requested_samples = length(requested),
    requested_present = sum(requested %in% names(values)),
    stringsAsFactors = FALSE
  )
  for (sample_id in requested) {
    present <- sample_id %in% names(values)
    if (!present) {
      detail_rows[[length(detail_rows) + 1L]] <- data.frame(
        cohort = source$cohort,
        representation = source$representation,
        source_file = source$source_file,
        historical_extraction = source$historical_extraction,
        sample_id = sample_id,
        present = FALSE,
        entry_status = "missing_sample",
        entry_type = "",
        matrix_class = "",
        genes = NA_integer_,
        cells = NA_integer_,
        unique_gene_ids = FALSE,
        unique_cell_ids = FALSE,
        finite_stored_values = FALSE,
        stored_value_min = NA_real_,
        stored_value_max = NA_real_,
        manifest_filtered_cells = manifest$filtered_cells[
          match(sample_id, manifest$sample_id)
        ],
        cell_count_matches_manifest = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    matrix_value <- values[[sample_id]]
    if (length(dim(matrix_value)) != 2L) {
      detail_rows[[length(detail_rows) + 1L]] <- data.frame(
        cohort = source$cohort,
        representation = source$representation,
        source_file = source$source_file,
        historical_extraction = source$historical_extraction,
        sample_id = sample_id,
        present = TRUE,
        entry_status = "not_matrix_like",
        entry_type = typeof(matrix_value),
        matrix_class = paste(class(matrix_value), collapse = ";"),
        genes = NA_integer_,
        cells = NA_integer_,
        unique_gene_ids = FALSE,
        unique_cell_ids = FALSE,
        finite_stored_values = FALSE,
        stored_value_min = NA_real_,
        stored_value_max = NA_real_,
        manifest_filtered_cells = manifest$filtered_cells[
          match(sample_id, manifest$sample_id)
        ],
        cell_count_matches_manifest = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    stored <- matrix_values(matrix_value)
    gene_ids <- rownames(matrix_value)
    cell_ids <- colnames(matrix_value)
    identity[[paste(source$cohort, source$representation, sample_id, sep = "::")]] <- list(
      genes = gene_ids,
      cells = cell_ids
    )
    expected_cells <- manifest$filtered_cells[match(sample_id, manifest$sample_id)]
    detail_rows[[length(detail_rows) + 1L]] <- data.frame(
      cohort = source$cohort,
      representation = source$representation,
      source_file = source$source_file,
      historical_extraction = source$historical_extraction,
      sample_id = sample_id,
      present = TRUE,
      entry_status = "matrix_like",
      entry_type = typeof(matrix_value),
      matrix_class = class(matrix_value)[[1L]],
      genes = nrow(matrix_value),
      cells = ncol(matrix_value),
      unique_gene_ids = !is.null(gene_ids) && !anyNA(gene_ids) &&
        all(nzchar(gene_ids)) && !anyDuplicated(gene_ids),
      unique_cell_ids = !is.null(cell_ids) && !anyNA(cell_ids) &&
        all(nzchar(cell_ids)) && !anyDuplicated(cell_ids),
      finite_stored_values = all(is.finite(stored)),
      stored_value_min = if (length(stored)) min(c(0, stored)) else 0,
      stored_value_max = if (length(stored)) max(c(0, stored)) else 0,
      manifest_filtered_cells = expected_cells,
      cell_count_matches_manifest = identical(
        as.integer(ncol(matrix_value)), as.integer(expected_cells)
      ),
      stringsAsFactors = FALSE
    )
  }
  rm(values)
  invisible(gc())
}

detail <- do.call(rbind, detail_rows)
summary <- do.call(rbind, source_summaries)

pair_rows <- list()
canonical_cell_id <- function(x) sub("^.*__", "", x)
for (cohort in unique(manifest$cohort)) {
  representations <- if (cohort == "large") {
    c("sct_whole", "seurat_integration")
  } else {
    c("sct_whole", "integrated")
  }
  sample_ids <- manifest$sample_id[manifest$cohort == cohort]
  common_genes <- NULL
  for (sample_id in sample_ids) {
    first <- identity[[paste(cohort, representations[[1L]], sample_id, sep = "::")]]
    second <- identity[[paste(cohort, representations[[2L]], sample_id, sep = "::")]]
    if (is.null(first) || is.null(second)) {
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        cohort = cohort,
        sample_id = sample_id,
        first_representation = representations[[1L]],
        second_representation = representations[[2L]],
        common_genes = NA_integer_,
        first_cells = if (is.null(first)) NA_integer_ else length(first$cells),
        second_cells = if (is.null(second)) NA_integer_ else length(second$cells),
        shared_cells = NA_integer_,
        identical_ordered_cells = FALSE,
        canonical_shared_cells = NA_integer_,
        canonical_first_unique = FALSE,
        canonical_second_unique = FALSE,
        canonical_ordered_identical = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    sample_common <- intersect(first$genes, second$genes)
    common_genes <- if (is.null(common_genes)) {
      sample_common
    } else {
      intersect(common_genes, sample_common)
    }
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      cohort = cohort,
      sample_id = sample_id,
      first_representation = representations[[1L]],
      second_representation = representations[[2L]],
      common_genes = length(sample_common),
      first_cells = length(first$cells),
      second_cells = length(second$cells),
      shared_cells = length(intersect(first$cells, second$cells)),
      identical_ordered_cells = identical(first$cells, second$cells),
      canonical_shared_cells = length(intersect(
        canonical_cell_id(first$cells), canonical_cell_id(second$cells)
      )),
      canonical_first_unique = !anyDuplicated(canonical_cell_id(first$cells)),
      canonical_second_unique = !anyDuplicated(canonical_cell_id(second$cells)),
      canonical_ordered_identical = identical(
        canonical_cell_id(first$cells), canonical_cell_id(second$cells)
      ),
      stringsAsFactors = FALSE
    )
  }
  summary$all_pilot_common_genes[
    summary$cohort == cohort
  ] <- length(common_genes)
}

pair_detail <- do.call(rbind, pair_rows)
pair_detail$record_type <- "paired_representation"
detail$record_type <- "source_entry"
for (column in setdiff(names(detail), names(pair_detail))) {
  pair_detail[[column]] <- NA
}
for (column in setdiff(names(pair_detail), names(detail))) {
  detail[[column]] <- NA
}
detail <- rbind(
  detail[, union(names(detail), names(pair_detail)), drop = FALSE],
  pair_detail[, union(names(detail), names(pair_detail)), drop = FALSE]
)

dir.create(dirname(detail_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(detail, detail_output, row.names = FALSE, na = "")
utils::write.csv(summary, summary_output, row.names = FALSE, na = "")
