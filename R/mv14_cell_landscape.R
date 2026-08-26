# Internal helpers for the MV14 all-QC cell-landscape contract.

.mv14_truth <- function(value) {
  if (is.logical(value)) return(!is.na(value) & value)
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}

.mv14_sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv14_sha256_text <- function(value) {
  digest::digest(paste(value, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

.mv14_read_csv <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

.mv14_atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) {
    stop("failed to atomically publish ", basename(path), call. = FALSE)
  }
  invisible(path)
}

.mv14_verify_manifest <- function(root, name) {
  path <- file.path(root, name)
  manifest <- .mv14_read_csv(path)
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(files, .mv14_sha256_file, character(1L)) ==
             manifest$sha256)) {
    stop("MV14 manifest drift: ", name, call. = FALSE)
  }
  list(rows = nrow(manifest), sha256 = .mv14_sha256_file(path))
}

.mv14_intervals <- function(record, homology_dimension) {
  dimension <- match(homology_dimension, c("H0", "H1")) - 1L
  if (is.na(dimension)) stop("MV14 dimension must be H0 or H1.", call. = FALSE)
  landscape_reference_intervals(record$result$diagram, dimension)
}

.mv14_active_depth <- function(intervals) {
  if (!nrow(intervals)) return(0L)
  points <- sort(unique(c(intervals[, "birth"], intervals[, "death"])))
  births <- tabulate(match(intervals[, "birth"], points), nbins = length(points))
  deaths <- tabulate(match(intervals[, "death"], points), nbins = length(points))
  as.integer(max(cumsum(births - deaths)))
}

.mv14_pair_identity <- function(group_id, pair_ordinal,
                                first_diagram_sha256,
                                second_diagram_sha256) {
  digest::digest(paste(
    "mv14_cell_landscape_pair_v1", group_id, as.integer(pair_ordinal),
    first_diagram_sha256, second_diagram_sha256, sep = "|"
  ), algo = "sha256", serialize = FALSE)
}

.mv14_group_pairs <- function(bindings, group_id) {
  bindings <- bindings[order(as.integer(bindings$axis_order)), , drop = FALSE]
  if (nrow(bindings) < 2L || anyDuplicated(bindings$unit_id) ||
      !identical(as.integer(bindings$axis_order), seq_len(nrow(bindings)))) {
    stop("MV14 private group axis is invalid.", call. = FALSE)
  }
  indices <- utils::combn(seq_len(nrow(bindings)), 2L)
  pairs <- data.frame(
    pair_ordinal = seq_len(ncol(indices)),
    first_axis_order = indices[1L, ],
    second_axis_order = indices[2L, ],
    first_unit_id = bindings$unit_id[indices[1L, ]],
    second_unit_id = bindings$unit_id[indices[2L, ]],
    first_diagram_sha256 = bindings$diagram_sha256[indices[1L, ]],
    second_diagram_sha256 = bindings$diagram_sha256[indices[2L, ]],
    stringsAsFactors = FALSE
  )
  pairs$pair_identity_sha256 <- vapply(seq_len(nrow(pairs)), function(index) {
    .mv14_pair_identity(
      group_id, pairs$pair_ordinal[[index]],
      pairs$first_diagram_sha256[[index]],
      pairs$second_diagram_sha256[[index]]
    )
  }, character(1L))
  pairs
}

.mv14_safe_group <- function(order) sprintf("group_%02d", as.integer(order))
.mv14_safe_chunk <- function(order) sprintf("chunk_%03d", as.integer(order))

.mv14_private_bytes <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(files)$isdir]
  sum(as.numeric(file.info(files)$size))
}
