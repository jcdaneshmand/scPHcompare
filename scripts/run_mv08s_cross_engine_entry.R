#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: run_mv08s_cross_engine_entry.R <mv08s-prefreeze> <cross-check-id>",
  "<source-rds> <common-panel> <subset32|full> <output-csv>"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
check_id <- args[[2L]]
source_path <- normalizePath(args[[3L]], mustWork = TRUE)
panel_path <- normalizePath(args[[4L]], mustWork = TRUE)
mode <- match.arg(args[[5L]], c("subset32", "full"))
output_path <- normalizePath(args[[6L]], mustWork = FALSE)
if (file.exists(output_path)) stop("refusing to overwrite MV8-S cross-engine output", call. = FALSE)
for (package in c("digest", "ripserr", "TDA")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")

contract <- utils::read.csv(
  file.path(prefreeze, "mv08s-cross-engine-contract.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
spec <- contract[contract$cross_check_id == check_id, , drop = FALSE]
if (nrow(spec) != 1L || spec$mode != mode || spec$workers != 1L ||
    spec$retries != 0L || spec$outcome_label_state != "closed" ||
    isTRUE(spec$biological_outcomes_computed)) {
  stop("MV8-S cross-engine check is not frozen", call. = FALSE)
}
queue <- utils::read.csv(
  file.path(prefreeze, "mv08s-ph-sentinel-queue.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
row <- queue[queue$job_id == spec$job_id, , drop = FALSE]
if (nrow(row) != 1L) stop("MV8-S cross-engine PH row missing", call. = FALSE)
panel <- utils::read.csv(panel_path, check.names = FALSE, stringsAsFactors = FALSE)
source_record <- readRDS(source_path)
if (row$execution_role == "source_produced_gene_ph") {
  view <- mv08s_residual_gene_view_v1(source_record, row, panel)
} else {
  mv08s_validate_baseline_record_v1(source_record)
  view <- source_record$views[[row$view_kind]]
}
validate_topology_view(view)

if (view$view_id == "cell_topology_v1") {
  value <- if (mode == "subset32") view$payload[seq_len(32L), , drop = FALSE] else
    view$payload
  point_count <- nrow(value)
  ripser_input <- value
  tda_input <- value
  tda_dist <- "euclidean"
  maximum_scale <- max(stats::dist(value))
} else {
  distance <- as.matrix(view$payload)
  if (mode == "subset32") {
    distance <- distance[seq_len(32L), seq_len(32L), drop = FALSE]
  }
  point_count <- nrow(distance)
  ripser_input <- stats::as.dist(distance)
  tda_input <- distance
  tda_dist <- "arbitrary"
  maximum_scale <- max(distance)
}
normalize_diagram <- function(value) {
  result <- as.matrix(value)
  storage.mode(result) <- "double"
  colnames(result) <- c("dimension", "birth", "death")
  result
}
remove_capped_essential <- function(diagram) {
  h0 <- which(diagram[, "dimension"] == 0)
  removed <- FALSE
  if (length(h0) == point_count) {
    diagram <- diagram[-h0[[which.max(diagram[h0, "death"])]], , drop = FALSE]
    removed <- TRUE
  }
  list(diagram = diagram, removed = removed)
}
compare_dimension <- function(first, second, dimension, tolerance = 1e-6) {
  first <- first[first[, "dimension"] == dimension,
                 c("birth", "death"), drop = FALSE]
  second <- second[second[, "dimension"] == dimension,
                   c("birth", "death"), drop = FALSE]
  order_rows <- function(value) value[order(
    value[, "birth"], value[, "death"], method = "radix"
  ), , drop = FALSE]
  first <- order_rows(first)
  second <- order_rows(second)
  error <- if (nrow(first) == nrow(second) && nrow(first)) {
    max(abs(first - second))
  } else if (!nrow(first) && !nrow(second)) 0 else Inf
  list(first_count = nrow(first), second_count = nrow(second),
       maximum_absolute_error = error,
       passed = is.finite(error) && error <= tolerance)
}
ripser <- remove_capped_essential(normalize_diagram(
  ripserr::vietoris_rips(
    ripser_input, max_dim = 1L, threshold = -1, p = 2L,
    return_format = "mat"
  )
))
gudhi <- remove_capped_essential(normalize_diagram(TDA::ripsDiag(
  X = tda_input, maxdimension = 1L, maxscale = maximum_scale,
  dist = tda_dist, library = "GUDHI", location = FALSE,
  printProgress = FALSE
)$diagram))
rows <- lapply(0:1, function(dimension) {
  comparison <- compare_dimension(ripser$diagram, gudhi$diagram, dimension)
  data.frame(
    contract_id = "mv08s_cross_engine_result_v1",
    cross_check_id = check_id,
    job_id = row$job_id,
    unit_id = row$unit_id,
    view_kind = row$view_kind,
    mode = mode,
    points = point_count,
    homology_dimension = paste0("H", dimension),
    ripserr_intervals = comparison$first_count,
    gudhi_intervals = comparison$second_count,
    ripserr_capped_essential_h0_removed = ripser$removed,
    gudhi_capped_essential_h0_removed = gudhi$removed,
    maximum_absolute_error = comparison$maximum_absolute_error,
    tolerance = 1e-6,
    passed = comparison$passed,
    outcome_label_state = "closed",
    biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
if (!all(result$passed)) stop("MV8-S cross-engine equivalence failed", call. = FALSE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
partial <- paste0(output_path, ".partial")
utils::write.csv(result, partial, row.names = FALSE, quote = TRUE, na = "")
if (!file.rename(partial, output_path)) {
  stop("failed to atomically publish MV8-S cross-engine result", call. = FALSE)
}
cat("MV8-S cross-engine check=", check_id, "; mode=", mode,
    "; points=", point_count, "; checks=2/2\n", sep = "")
