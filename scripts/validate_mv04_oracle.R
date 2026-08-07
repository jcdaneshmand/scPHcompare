#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: validate_mv04_oracle.R PAIRS INTERVALS OUTPUT")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
pairs <- read.csv(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
intervals <- read.csv(args[[2L]], stringsAsFactors = FALSE, check.names = FALSE)

make_diagram <- function(diagram_id) {
  part <- intervals[intervals$diagram_id == diagram_id, , drop = FALSE]
  result <- cbind(part$homology_dimension, part$birth, part$death)
  colnames(result) <- c("dimension", "Birth", "Death")
  result
}

# One eligible pair per stratum, choosing the smallest tractable diagram pair.
eligible <- pairs[pairs$homology_dimension == "H1" &
                    pmax(pairs$first_finite_intervals,
                         pairs$second_finite_intervals) <= 200L, , drop = FALSE]
eligible$size <- pmax(eligible$first_finite_intervals,
                      eligible$second_finite_intervals)
selected <- do.call(rbind, lapply(split(eligible, eligible$stratum_id), function(x) {
  x[order(x$size, x$pair_id), , drop = FALSE][1L, , drop = FALSE]
}))
# Add the smallest eligible H0 pair so both retained homology dimensions are
# exercised against the independent R implementation.
h0 <- pairs[pairs$homology_dimension == "H0", , drop = FALSE]
h0$size <- pmax(h0$first_finite_intervals, h0$second_finite_intervals)
selected <- rbind(selected, h0[order(h0$size, h0$pair_id), , drop = FALSE][1L, ])
if (!nrow(selected)) stop("No tractable eligible oracle cases were found.")

rows <- lapply(seq_len(nrow(selected)), function(index) {
  row <- selected[index, ]
  started <- proc.time()[["elapsed"]]
  dimension <- as.integer(sub("H", "", row$homology_dimension))
  observed <- landscape_reference_exact_dimension(
    make_diagram(row$first_diagram_id), make_diagram(row$second_diagram_id),
    dimension = dimension,
    exact_max_intervals = max(row$first_finite_intervals,
                              row$second_finite_intervals)
  )$distance
  elapsed <- proc.time()[["elapsed"]] - started
  tolerance <- max(1e-10, 1e-9 * abs(row$distance))
  data.frame(
    pair_id = row$pair_id, stratum_id = row$stratum_id,
    homology_dimension = row$homology_dimension,
    production_distance = row$distance,
    r_oracle_distance = observed, absolute_difference = abs(observed - row$distance),
    tolerance = tolerance, passed = abs(observed - row$distance) <= tolerance,
    oracle_seconds = elapsed, oracle_contract = "exact_breakpoint_stream_v1",
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
write.csv(result, args[[3L]], row.names = FALSE)
if (!all(result$passed)) stop("One or more eligible R oracle checks failed.")
message("Validated ", nrow(result), " eligible cases against the R exact oracle.")
