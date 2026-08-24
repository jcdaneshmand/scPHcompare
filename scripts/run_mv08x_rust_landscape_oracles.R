#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv08x_rust_landscape_oracles.R <private-selection.csv>",
  "<mv08s-private> <mv08v-private> <rust-library>",
  "<expected-library-sha256> <run-id> <output-dir>"
), call. = FALSE)

selection_path <- normalizePath(args[[1L]], mustWork = TRUE)
s_root <- normalizePath(args[[2L]], mustWork = TRUE)
v_root <- normalizePath(args[[3L]], mustWork = TRUE)
library_path <- normalizePath(args[[4L]], mustWork = TRUE)
expected_library_sha <- tolower(args[[5L]])
run_id <- args[[6L]]
output_dir <- normalizePath(args[[7L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-X oracle output", call. = FALSE)
if (!grepl("^[0-9a-f]{64}$", expected_library_sha) ||
    !run_id %in% c("a", "b")) {
  stop("invalid MV8-X oracle invocation identity", call. = FALSE)
}
for (package in "digest") if (!requireNamespace(package, quietly = TRUE)) {
  stop(package, " required", call. = FALSE)
}

source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")
source("R/mv08s_ph_sentinel.R")

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}

if (sha_file(library_path) != expected_library_sha) {
  stop("MV8-X Rust candidate SHA-256 mismatch", call. = FALSE)
}
selection <- read_csv(selection_path)
if (nrow(selection) != 28L ||
    !identical(as.integer(selection$oracle_order), 1:28) ||
    !all(selection$homology_dimension %in% c("H0", "H1")) ||
    any(selection$outcome_label_state != "closed") ||
    any(selection$biological_outcomes_computed)) {
  stop("MV8-X private oracle selection drift", call. = FALSE)
}

root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("unknown MV8-X PH source role", call. = FALSE)
}
record_cache <- new.env(parent = emptyenv())
load_record <- function(row, prefix) {
  role <- row[[paste0(prefix, "_source_role")]]
  relative <- row[[paste0(prefix, "_output_file")]]
  expected_sha <- row[[paste0(prefix, "_output_sha256")]]
  expected_bytes <- as.numeric(row[[paste0(prefix, "_output_bytes")]])
  key <- expected_sha
  if (exists(key, envir = record_cache, inherits = FALSE)) {
    return(get(key, envir = record_cache, inherits = FALSE))
  }
  path <- file.path(root_for(role), relative)
  if (!file.exists(path) || as.numeric(file.info(path)$size) != expected_bytes ||
      sha_file(path) != expected_sha) {
    stop("MV8-X selected PH artifact drift", call. = FALSE)
  }
  record <- readRDS(path)
  mv08s_validate_ph_record_v1(record)
  expected_diagram <- row[[paste0(prefix, "_diagram_sha256")]]
  if (record$topology_result$provenance$diagram_sha256 != expected_diagram) {
    stop("MV8-X selected diagram identity drift", call. = FALSE)
  }
  assign(key, record, envir = record_cache)
  record
}

reference_one <- function(first, second, dimension, row) {
  if (row$reference_route == "r_exact_breakpoint") {
    landscape_reference_exact_dimension(
      first, second, dimension, as.integer(row$exact_max_intervals)
    )
  } else if (row$reference_route == "r_adaptive_certified") {
    result <- landscape_reference_adaptive_dimension(
      first, second, dimension, as.numeric(row$absolute_tolerance),
      as.numeric(row$relative_tolerance), as.integer(row$subdivisions)
    )
    if (!isTRUE(result$within_requested_tolerance)) {
      stop("MV8-X canonical R adaptive oracle failed certification", call. = FALSE)
    }
    result
  } else stop("unknown MV8-X reference route", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
reference_seconds <- 0
rust_seconds <- 0
results <- lapply(seq_len(nrow(selection)), function(index) {
  row <- selection[index, , drop = FALSE]
  first_record <- load_record(row, "first")
  second_record <- load_record(row, "second")
  first_pd <- first_record$topology_result$diagram
  second_pd <- second_record$topology_result$diagram
  dimension <- as.integer(sub("H", "", row$homology_dimension))
  first_intervals <- landscape_reference_intervals(first_pd, dimension)
  second_intervals <- landscape_reference_intervals(second_pd, dimension)
  if (nrow(first_intervals) != as.integer(row$first_finite_intervals) ||
      nrow(second_intervals) != as.integer(row$second_finite_intervals)) {
    stop("MV8-X finite-interval count drift", call. = FALSE)
  }

  reference_started <- proc.time()[["elapsed"]]
  reference <- reference_one(first_pd, second_pd, dimension, row)
  reference_seconds <<- reference_seconds +
    proc.time()[["elapsed"]] - reference_started

  rust_started <- proc.time()[["elapsed"]]
  forward <- landscape_rust_prototype_dimension(
    first_intervals, second_intervals, dimension, library_path
  )
  reverse <- landscape_rust_prototype_dimension(
    second_intervals, first_intervals, dimension, library_path
  )
  first_self <- landscape_rust_prototype_dimension(
    first_intervals, first_intervals, dimension, library_path
  )
  second_self <- landscape_rust_prototype_dimension(
    second_intervals, second_intervals, dimension, library_path
  )
  rust_seconds <<- rust_seconds + proc.time()[["elapsed"]] - rust_started

  threshold <- if (isTRUE(reference$exact)) {
    1e-10 + 1e-10 * abs(reference$squared_distance)
  } else {
    reference$achieved_absolute_error_estimate +
      100 * .Machine$double.eps * max(1, abs(reference$squared_distance))
  }
  error <- abs(forward$squared_distance - reference$squared_distance)
  passed <- isTRUE(forward$rust_used) && forward$status == 0L &&
    forward$engine_version == 1L && is.finite(forward$squared_distance) &&
    forward$squared_distance >= 0 && error <= threshold &&
    identical(forward$squared_distance, reverse$squared_distance) &&
    forward$first_finite_intervals == reverse$second_finite_intervals &&
    forward$second_finite_intervals == reverse$first_finite_intervals &&
    forward$active_levels == reverse$active_levels &&
    forward$event_segments == reverse$event_segments &&
    identical(first_self$squared_distance, 0) &&
    identical(second_self$squared_distance, 0) &&
    forward$active_levels == max(nrow(first_intervals), nrow(second_intervals))
  data.frame(
    contract_id = "mv08x_oracle_result_v1",
    oracle_order = as.integer(row$oracle_order), oracle_id = row$oracle_id,
    pair_identity_sha256 = row$pair_identity_sha256,
    dataset_scope = row$dataset_scope,
    representation_id = row$representation_id, panel_id = row$panel_id,
    view_kind = row$view_kind,
    seed = as.integer(row$seed), homology_dimension = row$homology_dimension,
    group_units = as.integer(row$group_units),
    group_unordered_pairs = as.integer(row$group_unordered_pairs),
    reference_route = row$reference_route, reference_exact = reference$exact,
    reference_squared_distance = reference$squared_distance,
    reference_error_estimate = reference$achieved_absolute_error_estimate,
    candidate_squared_distance = forward$squared_distance,
    absolute_error = error, acceptance_threshold = threshold,
    first_finite_intervals = nrow(first_intervals),
    second_finite_intervals = nrow(second_intervals),
    active_levels = forward$active_levels, event_segments = forward$event_segments,
    status = forward$status, engine_version = forward$engine_version,
    reverse_bit_identical = identical(
      forward$squared_distance, reverse$squared_distance
    ),
    reverse_counts_swap =
      forward$first_finite_intervals == reverse$second_finite_intervals &&
      forward$second_finite_intervals == reverse$first_finite_intervals,
    reverse_diagnostics_match =
      forward$active_levels == reverse$active_levels &&
      forward$event_segments == reverse$event_segments,
    first_self_exact_zero = identical(first_self$squared_distance, 0),
    second_self_exact_zero = identical(second_self$squared_distance, 0),
    all_active_levels = forward$active_levels ==
      max(nrow(first_intervals), nrow(second_intervals)),
    passed = passed, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, results)
if (!all(results$passed)) stop("MV8-X canonical R oracle gate failed", call. = FALSE)

empty <- matrix(numeric(), nrow = 0L, ncol = 2L,
                dimnames = list(NULL, c("birth", "death")))
matrix2 <- function(values) {
  result <- matrix(values, ncol = 2L, byrow = TRUE)
  colnames(result) <- c("birth", "death")
  result
}
fixture_definitions <- list(
  single_tent = list(matrix2(c(0, 2)), empty, 2 / 3),
  sign_changing_tents = list(matrix2(c(0, 2)), matrix2(c(0.25, 2.25)), 7 / 64),
  narrow_feature = list(matrix2(c(0.499, 0.501)), empty, 0.002^3 / 12),
  duplicate_tents = list(matrix2(c(0, 2, 0, 2)), empty, 4 / 3),
  empty_pair = list(empty, empty, 0),
  unequal_depth = list(
    matrix2(c(0, 2, 0.5, 1.5)), matrix2(c(0, 2)), NA_real_
  ),
  multi_interval_self = list(
    matrix2(c(0, 2, 0.2, 2.7, 1.1, 1.4)),
    matrix2(c(0, 2, 0.2, 2.7, 1.1, 1.4)), 0
  )
)
fixtures <- lapply(names(fixture_definitions), function(case) {
  definition <- fixture_definitions[[case]]
  expected <- definition[[3L]]
  if (!is.finite(expected)) {
    first_pd <- cbind(0, definition[[1L]])
    second_pd <- cbind(0, definition[[2L]])
    expected <- landscape_reference_exact_dimension(
      first_pd, second_pd, 0L, exact_max_intervals = 500L
    )$squared_distance
  }
  candidate <- landscape_rust_prototype_dimension(
    definition[[1L]], definition[[2L]], 0L, library_path
  )
  error <- abs(candidate$squared_distance - expected)
  data.frame(
    contract_id = "mv08x_fixture_result_v1", case = case,
    expected_squared_distance = expected,
    observed_squared_distance = candidate$squared_distance,
    absolute_error = error, acceptance_threshold = 1e-12,
    status = candidate$status, engine_version = candidate$engine_version,
    rust_used = candidate$rust_used, fallback_used = FALSE,
    passed = isTRUE(candidate$rust_used) && candidate$status == 0L &&
      candidate$engine_version == 1L && error <= 1e-12,
    stringsAsFactors = FALSE
  )
})
fallback_first <- matrix2(c(0, 2))
missing <- landscape_rust_prototype_with_fallback(
  fallback_first, empty, 0L, reference_squared = function() 2 / 3,
  library = paste0(library_path, ".missing")
)
fixtures[[length(fixtures) + 1L]] <- data.frame(
  contract_id = "mv08x_fixture_result_v1", case = "missing_library_fallback",
  expected_squared_distance = 2 / 3, observed_squared_distance = missing$squared_distance,
  absolute_error = abs(missing$squared_distance - 2 / 3),
  acceptance_threshold = 0, status = missing$status,
  engine_version = missing$engine_version, rust_used = missing$rust_used,
  fallback_used = missing$fallback_used,
  passed = isTRUE(missing$fallback_used) && !missing$rust_used &&
    missing$status == 9001L && identical(missing$squared_distance, 2 / 3),
  stringsAsFactors = FALSE
)
corrupt <- tempfile("mv08x-corrupt-library-")
writeBin(charToRaw("not a dynamic library"), corrupt)
on.exit(unlink(corrupt), add = TRUE)
corrupt_result <- landscape_rust_prototype_with_fallback(
  fallback_first, empty, 0L, reference_squared = function() 2 / 3,
  library = corrupt
)
fixtures[[length(fixtures) + 1L]] <- data.frame(
  contract_id = "mv08x_fixture_result_v1", case = "corrupt_library_fallback",
  expected_squared_distance = 2 / 3,
  observed_squared_distance = corrupt_result$squared_distance,
  absolute_error = abs(corrupt_result$squared_distance - 2 / 3),
  acceptance_threshold = 0, status = corrupt_result$status,
  engine_version = corrupt_result$engine_version,
  rust_used = corrupt_result$rust_used,
  fallback_used = corrupt_result$fallback_used,
  passed = isTRUE(corrupt_result$fallback_used) && !corrupt_result$rust_used &&
    corrupt_result$status == 9001L &&
    identical(corrupt_result$squared_distance, 2 / 3),
  stringsAsFactors = FALSE
)
fixtures <- do.call(rbind, fixtures)
if (nrow(fixtures) != 9L || !all(fixtures$passed)) {
  stop("MV8-X analytical/fallback fixture gate failed", call. = FALSE)
}

status_lines <- if (file.exists("/proc/self/status")) {
  readLines("/proc/self/status", warn = FALSE)
} else character()
peak_line <- grep("^VmHWM:", status_lines, value = TRUE)
peak_rss <- if (length(peak_line) == 1L) {
  as.numeric(gsub("[^0-9]", "", peak_line)) * 1024
} else NA_real_
resource <- data.frame(
  contract_id = "mv08x_oracle_resource_v1", run_id = run_id,
  oracle_pairs = nrow(results), reference_seconds = reference_seconds,
  rust_seconds = rust_seconds,
  total_seconds = proc.time()[["elapsed"]] - started,
  peak_process_rss_bytes = peak_rss, elapsed_cap_seconds = 3600,
  rss_cap_bytes = 12 * 1024^3, workers = 1L, retries = 0L,
  rust_library_sha256 = sha_file(library_path),
  private_selection_sha256 = sha_file(selection_path),
  production_landscape_jobs = 0L, comparison_jobs = 0L,
  clustering_jobs = 0L, fusion_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
if (!is.finite(resource$peak_process_rss_bytes) ||
    resource$total_seconds > resource$elapsed_cap_seconds ||
    resource$peak_process_rss_bytes > resource$rss_cap_bytes) {
  stop("MV8-X oracle resource cap failed", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE)
atomic_csv(results, file.path(output_dir, "oracle-results.csv"))
atomic_csv(fixtures, file.path(output_dir, "fixture-results.csv"))
atomic_csv(resource, file.path(output_dir, "resource.csv"))
files <- list.files(output_dir, full.names = TRUE)
manifest <- data.frame(
  artifact = basename(files), bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "artifact-manifest.csv"))
cat("MV8-X oracle run ", run_id, ": ", sum(results$passed), "/",
    nrow(results), " pairs; fixtures ", sum(fixtures$passed), "/",
    nrow(fixtures), "\n", sep = "")
