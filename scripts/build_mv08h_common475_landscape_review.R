#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv08h_common475_landscape_review.R <private-diagrams.rds> <output-dir>", call. = FALSE)
}

diagrams_path <- args[[1L]]
output_dir <- args[[2L]]
if (!file.exists(diagrams_path)) stop("private diagrams file not found: ", diagrams_path, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("R/landscape_contract.R")
source("R/landscape_reference.R")

diagrams <- readRDS(diagrams_path)
required_views <- c("cell_topology_v1", "gene_topology_v1")
if (!all(required_views %in% names(diagrams))) stop("private diagrams must contain both frozen views", call. = FALSE)

integrate_line_squared <- function(slope, intercept, left, right) {
  # Integral of (slope*x + intercept)^2 over one breakpoint partition.
  slope ^ 2 * (right ^ 3 - left ^ 3) / 3 +
    slope * intercept * (right ^ 2 - left ^ 2) +
    intercept ^ 2 * (right - left)
}

profile_dimension <- function(pd, dimension) {
  intervals <- landscape_reference_intervals(pd, dimension)
  if (!nrow(intervals)) {
    return(list(
      finite_intervals = 0L, active_levels = 0L, breakpoint_count = 0L,
      energy = numeric(), breakpoints = numeric()
    ))
  }
  births <- intervals[, "birth"]
  deaths <- intervals[, "death"]
  mids <- (births + deaths) / 2
  segments <- rbind(
    data.frame(interval = seq_len(nrow(intervals)), left = births, right = mids,
               slope = 1, intercept = -births),
    data.frame(interval = seq_len(nrow(intervals)), left = mids, right = deaths,
               slope = -1, intercept = deaths)
  )
  breaks <- landscape_reference_breakpoints(intervals)
  energy <- numeric()
  max_active <- 0L
  for (index in seq_len(length(breaks) - 1L)) {
    left <- breaks[[index]]
    right <- breaks[[index + 1L]]
    if (!is.finite(left) || !is.finite(right) || right <= left) next
    location <- (left + right) / 2
    active <- segments$left < location & location < segments$right
    if (!any(active)) next
    active_segments <- segments[active, , drop = FALSE]
    values <- active_segments$slope * location + active_segments$intercept
    order_index <- order(values, decreasing = TRUE, method = "radix")
    active_segments <- active_segments[order_index, , drop = FALSE]
    max_active <- max(max_active, nrow(active_segments))
    contributions <- vapply(seq_len(nrow(active_segments)), function(rank) {
      integrate_line_squared(
        active_segments$slope[[rank]], active_segments$intercept[[rank]],
        left, right
      )
    }, numeric(1))
    if (length(energy) < length(contributions)) {
      energy <- c(energy, numeric(length(contributions) - length(energy)))
    }
    energy[seq_along(contributions)] <- energy[seq_along(contributions)] + contributions
  }
  energy[!is.finite(energy)] <- 0
  list(
    finite_intervals = nrow(intervals), active_levels = max_active,
    breakpoint_count = length(breaks), energy = pmax(0, energy),
    breakpoints = breaks
  )
}

summary_rows <- list()
execution_rows <- list()
row_index <- 0L
for (view_id in required_views) {
  pd <- diagrams[[view_id]]
  diagram_hash <- digest::digest(pd, algo = "sha256")
  for (dimension in 0:1) {
    profile <- profile_dimension(pd, dimension)
    hdim <- paste0("H", dimension)
    if (length(profile$energy)) {
      for (level in seq_along(profile$energy)) {
        row_index <- row_index + 1L
        summary_rows[[row_index]] <- data.frame(
          view_id = view_id, homology_dimension = hdim, level = level,
          finite_positive_intervals = profile$finite_intervals,
          active_level_count = profile$active_levels,
          integrated_lambda_squared = profile$energy[[level]],
          integration_method = "exact_critical_breakpoint_v1",
          level_policy = "all_active_consecutive_levels",
          grid_policy = "none",
          infinite_interval_policy = "exclude_before_landscape_construction",
          diagram_sha256 = diagram_hash,
          labels_outcomes_opened = FALSE, fusion_computed = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
    execution_rows[[length(execution_rows) + 1L]] <- data.frame(
      view_id = view_id, homology_dimension = hdim,
      finite_positive_intervals = profile$finite_intervals,
      active_level_count = profile$active_levels,
      breakpoint_count = profile$breakpoint_count,
      total_integrated_lambda_squared = sum(profile$energy),
      integration_method = "exact_critical_breakpoint_v1",
      level_policy = "all_active_consecutive_levels", grid_policy = "none",
      infinite_interval_policy = "exclude_before_landscape_construction",
      diagram_sha256 = diagram_hash, labels_outcomes_opened = FALSE,
      fusion_computed = FALSE, stringsAsFactors = FALSE
    )
  }
}

summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
execution <- do.call(rbind, execution_rows)
write.csv(summary, file.path(output_dir, "mv08h-common475-landscape-level-summary.csv"), row.names = FALSE, quote = TRUE)
write.csv(execution, file.path(output_dir, "mv08h-common475-landscape-execution.csv"), row.names = FALSE, quote = TRUE)

manifest <- data.frame(
  artifact = c("mv08h-common475-landscape-level-summary.csv", "mv08h-common475-landscape-execution.csv"),
  purpose = c("one row per retained active landscape level", "one row per frozen view and homology dimension"),
  public = TRUE, contains_private_diagrams = FALSE, stringsAsFactors = FALSE
)
write.csv(manifest, file.path(output_dir, "mv08h-common475-landscape-artifact-manifest.csv"), row.names = FALSE, quote = TRUE)

report_lines <- c(
  paste0("# MV8-H common475 all-active-level landscape review (", format(Sys.Date()), ")"),
  "",
  "This is label-closed feasibility evidence for HCA_BM_002 only; it is not a biological result or validation claim.",
  "",
  "- Cell and gene topology views are reported separately.",
  "- H0 and H1 are integrated separately; no combined/fusion landscape was computed.",
  "- All active consecutive levels are retained; there is no fixed level cap.",
  "- Exact integration uses critical breakpoints (births, midpoints, deaths, and line crossings); no uniform grid is used.",
  "- Essential H0 intervals are excluded before landscape construction.",
  "- Labels, outcomes, other HCA units, and deletion remain closed.",
  "",
  paste0("Private input: ", basename(diagrams_path), " (not copied to public output).")
)
writeLines(report_lines, file.path(output_dir, paste0("MV08H_COMMON475_LANDSCAPE_", format(Sys.Date()), ".md")))

cat("completed common475 landscape review: ", nrow(execution), " view-dimension profiles; ",
    nrow(summary), " active levels\n", sep = "")
