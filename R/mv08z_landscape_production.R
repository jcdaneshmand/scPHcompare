# Internal helpers for the MV8-Z streamed production-landscape contract.

.mv08z_truth <- function(value) {
  if (is.logical(value)) return(!is.na(value) & value)
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}

.mv08z_sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.mv08z_sha256_text <- function(value) {
  digest::digest(paste(value, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

.mv08z_read_csv <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

.mv08z_atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) {
    stop("failed to atomically publish ", basename(path), call. = FALSE)
  }
  invisible(path)
}

.mv08z_verify_manifest <- function(root, name) {
  manifest_path <- file.path(root, name)
  manifest <- .mv08z_read_csv(manifest_path)
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !all(as.numeric(file.info(paths)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(paths, .mv08z_sha256_file, character(1L)) ==
             manifest$sha256)) {
    stop("MV8-Z bound manifest drift: ", name, call. = FALSE)
  }
  list(rows = nrow(manifest), sha256 = .mv08z_sha256_file(manifest_path))
}

.mv08z_pair_identity <- function(group_id, pair_ordinal,
                                 first_diagram_sha256,
                                 second_diagram_sha256) {
  payload <- paste(
    "mv08z_landscape_pair_v1", group_id, as.integer(pair_ordinal),
    first_diagram_sha256, second_diagram_sha256, sep = "|"
  )
  digest::digest(payload, algo = "sha256", serialize = FALSE)
}

.mv08z_group_pairs <- function(bindings) {
  bindings <- bindings[order(as.integer(bindings$axis_order)), , drop = FALSE]
  if (!identical(as.integer(bindings$axis_order), seq_len(nrow(bindings))) ||
      anyDuplicated(bindings$job_id) || nrow(bindings) < 2L) {
    stop("MV8-Z private group axis is invalid", call. = FALSE)
  }
  indices <- utils::combn(seq_len(nrow(bindings)), 2L)
  data.frame(
    pair_ordinal = seq_len(ncol(indices)),
    first_axis_order = indices[1L, ], second_axis_order = indices[2L, ],
    first_job_id = bindings$job_id[indices[1L, ]],
    second_job_id = bindings$job_id[indices[2L, ]],
    first_unit_id = bindings$unit_id[indices[1L, ]],
    second_unit_id = bindings$unit_id[indices[2L, ]],
    first_diagram_sha256 = bindings$diagram_sha256[indices[1L, ]],
    second_diagram_sha256 = bindings$diagram_sha256[indices[2L, ]],
    stringsAsFactors = FALSE
  )
}

.mv08z_add_pair_identities <- function(pairs, group_id) {
  pairs$pair_identity_sha256 <- vapply(seq_len(nrow(pairs)), function(index) {
    .mv08z_pair_identity(
      group_id, pairs$pair_ordinal[[index]],
      pairs$first_diagram_sha256[[index]],
      pairs$second_diagram_sha256[[index]]
    )
  }, character(1L))
  pairs
}

.mv08z_finite_intervals <- function(record, homology_dimension) {
  dimension <- as.integer(sub("H", "", homology_dimension, fixed = TRUE))
  if (is.na(dimension) || !dimension %in% 0:1) {
    stop("MV8-Z homology dimension must be H0 or H1", call. = FALSE)
  }
  landscape_reference_intervals(record$topology_result$diagram, dimension)
}

.mv08z_active_depth <- function(intervals) {
  if (!nrow(intervals)) return(0L)
  points <- sort(unique(c(intervals[, "birth"], intervals[, "death"])))
  births <- tabulate(match(intervals[, "birth"], points), nbins = length(points))
  deaths <- tabulate(match(intervals[, "death"], points), nbins = length(points))
  as.integer(max(cumsum(births - deaths)))
}

.mv08z_safe_group <- function(group_order) sprintf("group_%02d", as.integer(group_order))

.mv08z_safe_chunk <- function(chunk_order) sprintf("chunk_%03d", as.integer(chunk_order))
