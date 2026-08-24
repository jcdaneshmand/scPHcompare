#!/usr/bin/env Rscript

# Aggregate-only diagnosis for an MV8-X oracle-gate failure. This reruns only
# fast Rust invariants on the frozen pairs; it does not run the canonical R
# integrator, publish private locators, or authorize production landscapes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(paste(
  "usage: run_mv08xa_rust_invariant_diagnostic.R <private-selection.csv>",
  "<mv08s-private> <mv08v-private> <rust-library> <library-sha256>",
  "<execution-head> <output.csv>"
), call. = FALSE)

selection_path <- normalizePath(args[[1L]], mustWork = TRUE)
s_root <- normalizePath(args[[2L]], mustWork = TRUE)
v_root <- normalizePath(args[[3L]], mustWork = TRUE)
library_path <- normalizePath(args[[4L]], mustWork = TRUE)
expected_sha <- args[[5L]]
execution_head <- args[[6L]]
output_path <- normalizePath(args[[7L]], mustWork = FALSE)
if (file.exists(output_path)) stop("refusing to overwrite diagnostic", call. = FALSE)
if (!grepl("^[0-9a-f]{64}$", expected_sha) ||
    !grepl("^[0-9a-f]{40}$", execution_head)) {
  stop("invalid diagnostic identity", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_rust_prototype.R")

sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
if (sha_file(library_path) != expected_sha) stop("candidate hash drift", call. = FALSE)
selection <- utils::read.csv(
  selection_path, check.names = FALSE, stringsAsFactors = FALSE
)
if (nrow(selection) != 28L) stop("diagnostic selection drift", call. = FALSE)

root_for <- function(role) {
  if (role == "mv08s_private_v3") return(s_root)
  if (role == "mv08v_recovery_private_v2") return(v_root)
  stop("unknown source role", call. = FALSE)
}
load_diagram <- function(row, prefix) {
  path <- file.path(
    root_for(row[[paste0(prefix, "_source_role")]]),
    row[[paste0(prefix, "_output_file")]]
  )
  if (!file.exists(path) ||
      sha_file(path) != row[[paste0(prefix, "_output_sha256")]]) {
    stop("selected PH hash drift", call. = FALSE)
  }
  readRDS(path)$topology_result$diagram
}

checks <- lapply(seq_len(nrow(selection)), function(index) {
  row <- selection[index, , drop = FALSE]
  dimension <- as.integer(sub("H", "", row$homology_dimension))
  first <- landscape_reference_intervals(load_diagram(row, "first"), dimension)
  second <- landscape_reference_intervals(load_diagram(row, "second"), dimension)
  forward <- landscape_rust_prototype_dimension(first, second, dimension, library_path)
  reverse <- landscape_rust_prototype_dimension(second, first, dimension, library_path)
  first_self <- landscape_rust_prototype_dimension(first, first, dimension, library_path)
  second_self <- landscape_rust_prototype_dimension(second, second, dimension, library_path)
  c(
    engine = isTRUE(forward$rust_used) && forward$status == 0L &&
      forward$engine_version == 1L && is.finite(forward$squared_distance) &&
      forward$squared_distance >= 0,
    reverse_bit = identical(forward$squared_distance, reverse$squared_distance),
    reverse_counts = forward$first_finite_intervals == reverse$second_finite_intervals &&
      forward$second_finite_intervals == reverse$first_finite_intervals,
    reverse_diagnostics = forward$active_levels == reverse$active_levels &&
      forward$event_segments == reverse$event_segments,
    first_self_zero = identical(first_self$squared_distance, 0),
    second_self_zero = identical(second_self$squared_distance, 0),
    original_interval_count_assertion =
      forward$active_levels == max(nrow(first), nrow(second))
  )
})
checks <- do.call(rbind, checks)
receipt <- data.frame(
  contract_id = "mv08xa_rust_invariant_diagnostic_v1",
  execution_head = execution_head,
  diagnostic_runner_sha256 = sha_file(
    "scripts/run_mv08xa_rust_invariant_diagnostic.R"
  ),
  candidate_sha256 = expected_sha,
  private_selection_sha256 = sha_file(selection_path),
  pairs = nrow(checks),
  engine_passes = sum(checks[, "engine"]),
  reverse_bit_passes = sum(checks[, "reverse_bit"]),
  reverse_count_passes = sum(checks[, "reverse_counts"]),
  reverse_diagnostic_passes = sum(checks[, "reverse_diagnostics"]),
  first_self_zero_passes = sum(checks[, "first_self_zero"]),
  second_self_zero_passes = sum(checks[, "second_self_zero"]),
  all_active_level_passes = sum(checks[, "original_interval_count_assertion"]),
  production_landscape_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
partial <- paste0(output_path, ".partial")
utils::write.csv(receipt, partial, row.names = FALSE, quote = TRUE, na = "")
if (!file.rename(partial, output_path)) stop("failed to publish diagnostic", call. = FALSE)
cat("MV8-XA Rust invariants: engine/reverse/self=28/28; original active assertion=",
    receipt$all_active_level_passes, "/28\n", sep = "")
