#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: stage_mv05ba_equivalence_inputs.R MV05AY_UNITS OUTPUT_DIR")
}
units <- normalizePath(args[[1L]], mustWork = TRUE)
out <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)

manifest <- utils::read.csv(
  "docs/audits/mv04-input-diagram-manifest-2026-08-05.csv",
  stringsAsFactors = FALSE
)
scope <- utils::read.csv(
  "docs/audits/mv05ax-complete-scope-2026-08-12.csv",
  stringsAsFactors = FALSE
)

verified <- vapply(seq_len(nrow(manifest)), function(index) {
  object <- readRDS(manifest$result_file[[index]])
  identical(digest::digest(manifest$result_file[[index]], algo = "sha256",
                           file = TRUE), manifest$result_file_sha256[[index]]) &&
    identical(digest::digest(object$diagram, algo = "sha256"),
              manifest$diagram_sha256[[index]]) &&
    identical(object$provenance$diagram_sha256,
              manifest$diagram_sha256[[index]])
}, logical(1))
if (!all(verified) || nrow(manifest) != 56L) stop("Source verification failed.")

interval_rows <- list()
reference_rows <- list()
diagram_index <- list()
for (row in seq_len(nrow(manifest))) {
  object <- readRDS(manifest$result_file[[row]])
  diagram <- .canonical_public_diagram_v1(object$diagram)
  finite <- is.finite(diagram[, "death"])
  values <- diagram[finite, , drop = FALSE]
  interval_rows[[row]] <- data.frame(
    contract_id = "mv05ba_private_interval_v1",
    stratum_id = manifest$stratum_id[[row]],
    diagram_id = manifest$diagram_id[[row]],
    diagram_sha256 = digest::digest(diagram, algo = "sha256", serialize = TRUE),
    interval_order = seq_len(nrow(values)),
    dimension = as.integer(values[, "dimension"]),
    birth = values[, "birth"], death = values[, "death"],
    stringsAsFactors = FALSE
  )
  diagram_index[[manifest$diagram_id[[row]]]] <- diagram
}

ref_index <- 0L
for (stratum_id in scope$stratum_id) {
  roots <- list.dirs(file.path(units, stratum_id, "corrected_landscape_v1"),
                     recursive = FALSE, full.names = TRUE)
  if (length(roots) != 1L) stop("Expected one MV5-AY sidecar: ", stratum_id)
  completion <- utils::read.csv(file.path(roots, "completion-v1.csv"),
                                stringsAsFactors = FALSE)
  .verify_completion_v1(roots, completion)
  pair_index <- utils::read.csv(file.path(roots, "pair-index-v1.csv"),
                                stringsAsFactors = FALSE)
  for (pair_row in seq_len(nrow(pair_index))) {
    path <- file.path(roots, pair_index$pair_artifact[[pair_row]])
    shard <- readRDS(path)
    if (!identical(digest::digest(path, algo = "sha256", file = TRUE),
                   pair_index$pair_sha256[[pair_row]])) {
      stop("Accepted shard hash drift: ", path)
    }
    for (dimension in c("H0", "H1")) {
      ref_index <- ref_index + 1L
      value <- shard$dimensions[[dimension]]
      reference_rows[[ref_index]] <- data.frame(
        contract_id = "mv05ba_private_reference_v1",
        stratum_id = stratum_id,
        pair_order = pair_index$pair_order[[pair_row]],
        pair_cache_key = shard$cache_key,
        pair_sha256 = pair_index$pair_sha256[[pair_row]],
        first_diagram_id = shard$provenance$first_source_id,
        second_diagram_id = shard$provenance$second_source_id,
        dimension = dimension,
        reference_distance = value$distance,
        reference_squared_distance = value$squared_distance,
        reference_method = value$method,
        reference_exact = value$exact,
        achieved_absolute_error_estimate =
          value$achieved_absolute_error_estimate,
        requested_absolute_tolerance = value$requested_absolute_tolerance,
        requested_relative_tolerance = value$requested_relative_tolerance,
        first_finite_intervals = value$first_finite_intervals,
        second_finite_intervals = value$second_finite_intervals,
        stringsAsFactors = FALSE
      )
    }
  }
}

intervals <- do.call(rbind, interval_rows)
references <- do.call(rbind, reference_rows)
references <- references[order(references$stratum_id, references$pair_order,
                               references$dimension, method = "radix"), ]
rownames(references) <- NULL
if (nrow(references) != 408L || sum(references$reference_exact) != 318L ||
    sum(!references$reference_exact) != 90L) {
  stop("Reference corpus counts changed.")
}

h1 <- references[references$dimension == "H1" & !references$reference_exact, ]
h1$interval_sum <- h1$first_finite_intervals + h1$second_finite_intervals
panel <- do.call(rbind, lapply(split(h1, h1$stratum_id), function(rows) {
  rows <- rows[order(-rows$interval_sum, rows$pair_order, method = "radix"), ]
  rows[seq_len(3L), c("stratum_id", "pair_order", "pair_cache_key",
                      "first_diagram_id", "second_diagram_id",
                      "first_finite_intervals", "second_finite_intervals",
                      "interval_sum")]
}))
panel <- panel[order(panel$stratum_id, panel$pair_order, method = "radix"), ]
rownames(panel) <- NULL
panel <- cbind(data.frame(contract_id = "mv05ba_speed_panel_v1",
                          panel_order = seq_len(nrow(panel))), panel)
if (nrow(panel) != 6L) stop("Expected six frozen high-depth benchmark pairs.")

utils::write.csv(intervals, file.path(out, "private-intervals.csv"),
                 row.names = FALSE, quote = TRUE)
saveRDS(diagram_index, file.path(out, "private-canonical-diagrams.rds"),
        version = 3)
utils::write.csv(references, file.path(out, "private-reference-distances.csv"),
                 row.names = FALSE, quote = TRUE)
saveRDS(references, file.path(out, "private-reference-distances.rds"),
        version = 3)
utils::write.csv(panel, file.path(out, "private-speed-panel.csv"),
                 row.names = FALSE, quote = TRUE)
speed_references <- do.call(rbind, lapply(seq_len(nrow(panel)), function(index) {
  references[references$stratum_id == panel$stratum_id[[index]] &
               references$pair_order == panel$pair_order[[index]], ]
}))
rownames(speed_references) <- NULL
if (nrow(speed_references) != 12L) stop("Expected twelve speed references.")
utils::write.csv(speed_references,
                 file.path(out, "private-speed-references.csv"),
                 row.names = FALSE, quote = TRUE)
summary <- data.frame(
  contract_id = "mv05ba_stage_summary_v1", diagrams = nrow(manifest),
  interval_rows = nrow(intervals), pair_dimension_rows = nrow(references),
  exact_reference_rows = sum(references$reference_exact),
  adaptive_reference_rows = sum(!references$reference_exact),
  speed_panel_pairs = nrow(panel), source_files_verified = all(verified),
  labels_opened = FALSE, outcomes_computed = FALSE, stringsAsFactors = FALSE
)
utils::write.csv(summary, file.path(out, "stage-summary.csv"),
                 row.names = FALSE, quote = TRUE)
cat("MV5-BA staged", nrow(references), "dimension references and",
    nrow(panel), "speed pairs\n")
