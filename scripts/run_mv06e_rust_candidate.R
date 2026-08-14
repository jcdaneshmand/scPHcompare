#!/usr/bin/env Rscript

args <- getOption("mv06e.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 7L) {
  stop("usage: run_mv06e_rust_candidate.R INTERVALS PAIRS REFERENCES ",
       "RUST_LIBRARY OUTPUT_CSV METRICS_CSV EXPECTED_LIBRARY_SHA",
       call. = FALSE)
}
interval_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
pair_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
reference_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
library_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
output_path <- args[[5L]]
metrics_path <- args[[6L]]
expected_library_sha <- tolower(args[[7L]])
source("R/landscape_rust_prototype.R")
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- tempfile(pattern = paste0(basename(path), "."),
                      tmpdir = dirname(path))
  utils::write.csv(value, partial, row.names = FALSE, na = "")
  if (file.exists(path)) {
    unlink(partial)
    stop("Refusing to overwrite an MV6-E Rust artifact.", call. = FALSE)
  }
  if (!file.rename(partial, path)) {
    unlink(partial)
    stop("Failed to atomically publish an MV6-E Rust artifact.",
         call. = FALSE)
  }
}
if (!identical(file_sha(library_path), expected_library_sha)) {
  stop("MV6-E Rust library differs from the accepted binary hash.",
       call. = FALSE)
}
if (file.exists(output_path) || file.exists(metrics_path)) {
  if (!file.exists(output_path) || !file.exists(metrics_path)) {
    stop("MV6-E Rust candidate has an incomplete artifact pair.",
         call. = FALSE)
  }
  existing <- utils::read.csv(
    output_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  existing_metrics <- utils::read.csv(
    metrics_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (nrow(existing) == 240L &&
      all(existing$contract_id == "mv06e_rust_exact_result_v1") &&
      nrow(existing_metrics) == 1L &&
      identical(existing_metrics$contract_id, "mv06e_rust_metrics_v1")) {
    message("Reused validated MV6-E Rust output.")
    quit(status = 0L, save = "no")
  }
  stop("MV6-E Rust candidate has partial or stale output.", call. = FALSE)
}
started <- proc.time()[["elapsed"]]
interval_rows <- utils::read.csv(interval_path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
pairs <- utils::read.csv(pair_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
references <- utils::read.csv(reference_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
if (nrow(pairs) != 180L || nrow(references) != 20L) {
  stop("MV6-E Rust pair/reference cardinality differs from prefreeze.",
       call. = FALSE)
}
split_key <- paste(interval_rows$diagram_id,
                   interval_rows$homology_dimension, sep = "\r")
intervals <- split(interval_rows[c("birth", "death")], split_key)
intervals <- lapply(intervals, function(value) {
  result <- as.matrix(value)
  storage.mode(result) <- "double"
  result[order(result[, 1L], result[, 2L]), , drop = FALSE]
})
loaded_seconds <- proc.time()[["elapsed"]] - started

compute <- function(row_type, pair_id, view_id, dimension,
                    first_id, second_id) {
  key <- function(id) paste(id, dimension, sep = "\r")
  first <- intervals[[key(first_id)]]
  second <- intervals[[key(second_id)]]
  result <- landscape_rust_prototype_dimension(
    first, second, as.integer(sub("H", "", dimension)), library_path
  )
  if (!isTRUE(result$rust_used) || result$status != 0L ||
      result$engine_version != 1L) {
    stop("MV6-E Rust kernel returned a non-success status.", call. = FALSE)
  }
  data.frame(
    contract_id = "mv06e_rust_exact_result_v1",
    engine_id = "rust_scph_landscape_kernel_v1", row_type = row_type,
    pair_id = pair_id, view_id = view_id,
    homology_dimension = dimension, first_diagram_id = first_id,
    second_diagram_id = second_id,
    first_finite_intervals = as.integer(result$first_finite_intervals),
    second_finite_intervals = as.integer(result$second_finite_intervals),
    squared_distance = result$squared_distance,
    distance = sqrt(result$squared_distance),
    active_levels = as.integer(result$active_levels),
    event_segments = as.integer(result$event_segments),
    critical_points = NA_integer_, exact = TRUE, all_active_levels = TRUE,
    level_cap_applied = FALSE, status = result$status,
    engine_version = result$engine_version, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
}

primary_started <- proc.time()[["elapsed"]]
primary <- lapply(seq_len(nrow(pairs)), function(index) {
  row <- pairs[index, ]
  compute("primary", row$pair_id, row$view_id, row$homology_dimension,
          row$first_diagram_id, row$second_diagram_id)
})
primary_seconds <- proc.time()[["elapsed"]] - primary_started
primary_frame <- do.call(rbind, primary)
primary_by_id <- split(primary_frame, primary_frame$pair_id)

checks_started <- proc.time()[["elapsed"]]
reverse <- lapply(seq_len(nrow(references)), function(index) {
  original <- primary_by_id[[references$pair_id[[index]]]][1L, ]
  compute("reverse", original$pair_id, original$view_id,
          original$homology_dimension, original$second_diagram_id,
          original$first_diagram_id)
})
diagram_manifest <- unique(interval_rows[c(
  "diagram_id", "view_id", "homology_dimension"
)])
self <- lapply(seq_len(nrow(diagram_manifest)), function(index) {
  row <- diagram_manifest[index, ]
  compute("self", paste0("self:", row$diagram_id, ":",
                          row$homology_dimension), row$view_id,
          row$homology_dimension, row$diagram_id, row$diagram_id)
})
check_seconds <- proc.time()[["elapsed"]] - checks_started
output <- do.call(rbind, c(primary, reverse, self))
output <- output[order(output$row_type, output$view_id,
                       output$homology_dimension, output$pair_id,
                       output$first_diagram_id, output$second_diagram_id,
                       method = "radix"), , drop = FALSE]
rownames(output) <- NULL
if (nrow(output) != 240L) stop("MV6-E Rust output row count drifted.",
                               call. = FALSE)
identity <- list(
  contract_id = "mv06e_rust_execution_identity_v1",
  intervals_sha256 = file_sha(interval_path), pairs_sha256 = file_sha(pair_path),
  references_sha256 = file_sha(reference_path),
  library_sha256 = file_sha(library_path),
  shim_sha256 = file_sha("R/landscape_rust_prototype.R"), rows = nrow(output)
)
identity_sha <- digest::digest(identity, algo = "sha256", serialize = TRUE)
output$execution_identity_sha256 <- identity_sha
atomic_csv(output, output_path)
total_seconds <- proc.time()[["elapsed"]] - started
status_lines <- if (file.exists("/proc/self/status")) {
  readLines("/proc/self/status", warn = FALSE)
} else character()
peak_line <- grep("^VmHWM:", status_lines, value = TRUE)
peak_rss <- if (length(peak_line) == 1L) {
  as.numeric(gsub("[^0-9]", "", peak_line)) * 1024
} else NA_real_
metrics <- data.frame(
  contract_id = "mv06e_rust_metrics_v1",
  engine_id = "rust_scph_landscape_kernel_v1",
  diagram_dimensions_built = 0L, primary_pair_rows = nrow(pairs),
  reverse_rows = nrow(references), self_rows = nrow(diagram_manifest),
  load_seconds = loaded_seconds, landscape_build_seconds = 0,
  primary_pair_seconds = primary_seconds,
  reverse_self_seconds = check_seconds, total_seconds = total_seconds,
  peak_process_rss_bytes = peak_rss,
  rust_library_sha256 = file_sha(library_path), engine_version = 1L,
  execution_identity_sha256 = identity_sha, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic_csv(metrics, metrics_path)
message("Completed MV6-E Rust candidate: ", nrow(output),
        " scientific rows.")
