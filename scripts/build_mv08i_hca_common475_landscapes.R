#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: build_mv08i_hca_common475_landscapes.R <private-diagram-dir> <output-dir>", call. = FALSE)
}
private_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source("R/landscape_contract.R")
source("R/landscape_reference.R")

files <- sort(list.files(private_dir, pattern = "-topology-diagrams\\.rds$", full.names = TRUE))
if (length(files) != 8L) stop("MV8-I requires eight private per-unit diagram files", call. = FALSE)

integrate_line_squared <- function(slope, intercept, left, right) {
  slope ^ 2 * (right ^ 3 - left ^ 3) / 3 +
    slope * intercept * (right ^ 2 - left ^ 2) +
    intercept ^ 2 * (right - left)
}

profile_dimension <- function(pd, dimension) {
  intervals <- landscape_reference_intervals(pd, dimension)
  if (!nrow(intervals)) return(list(finite_intervals = 0L, active_levels = 0L, breakpoint_count = 0L, energy = numeric()))
  births <- intervals[, "birth"]
  deaths <- intervals[, "death"]
  mids <- (births + deaths) / 2
  segments <- rbind(
    data.frame(interval = seq_len(nrow(intervals)), left = births, right = mids, slope = 1, intercept = -births),
    data.frame(interval = seq_len(nrow(intervals)), left = mids, right = deaths, slope = -1, intercept = deaths)
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
    order_index <- order(active_segments$slope * location + active_segments$intercept,
                         decreasing = TRUE, method = "radix")
    active_segments <- active_segments[order_index, , drop = FALSE]
    max_active <- max(max_active, nrow(active_segments))
    contributions <- vapply(seq_len(nrow(active_segments)), function(rank) {
      integrate_line_squared(active_segments$slope[[rank]], active_segments$intercept[[rank]], left, right)
    }, numeric(1))
    if (length(energy) < length(contributions)) energy <- c(energy, numeric(length(contributions) - length(energy)))
    energy[seq_along(contributions)] <- energy[seq_along(contributions)] + contributions
  }
  energy[!is.finite(energy)] <- 0
  list(finite_intervals = nrow(intervals), active_levels = max_active, breakpoint_count = length(breaks), energy = pmax(0, energy))
}

summary_rows <- list()
execution_rows <- list()
level_index <- 0L
execution_index <- 0L
for (path in files) {
  diagrams <- readRDS(path)
  unit_id <- diagrams$unit_id
  if (is.null(unit_id) || !all(c("cell_topology_v1", "gene_topology_v1") %in% names(diagrams))) {
    stop("private diagram file is missing the frozen unit/views: ", basename(path), call. = FALSE)
  }
  for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
    pd <- diagrams[[view_id]]
    diagram_hash <- digest::digest(pd, algo = "sha256")
    for (dimension in 0:1) {
      profile <- profile_dimension(pd, dimension)
      hdim <- paste0("H", dimension)
      execution_index <- execution_index + 1L
      execution_rows[[execution_index]] <- data.frame(
        unit_id = unit_id, view_id = view_id, homology_dimension = hdim,
        finite_positive_intervals = profile$finite_intervals, active_level_count = profile$active_levels,
        breakpoint_count = profile$breakpoint_count, total_integrated_lambda_squared = sum(profile$energy),
        integration_method = "exact_critical_breakpoint_v1", level_policy = "all_active_consecutive_levels",
        grid_policy = "none", infinite_interval_policy = "exclude_before_landscape_construction",
        diagram_sha256 = diagram_hash, labels_outcomes_opened = FALSE, fusion_computed = FALSE,
        stringsAsFactors = FALSE)
      if (length(profile$energy)) for (level in seq_along(profile$energy)) {
        level_index <- level_index + 1L
        summary_rows[[level_index]] <- data.frame(
          unit_id = unit_id, view_id = view_id, homology_dimension = hdim, level = level,
          finite_positive_intervals = profile$finite_intervals, active_level_count = profile$active_levels,
          integrated_lambda_squared = profile$energy[[level]], integration_method = "exact_critical_breakpoint_v1",
          level_policy = "all_active_consecutive_levels", grid_policy = "none",
          infinite_interval_policy = "exclude_before_landscape_construction", diagram_sha256 = diagram_hash,
          labels_outcomes_opened = FALSE, fusion_computed = FALSE, stringsAsFactors = FALSE)
      }
    }
  }
}
execution <- do.call(rbind, execution_rows)
levels <- do.call(rbind, summary_rows)
write.csv(execution, file.path(output_dir, "mv08i-landscape-execution.csv"), row.names = FALSE, quote = TRUE)
write.csv(levels, file.path(output_dir, "mv08i-landscape-level-summary.csv"), row.names = FALSE, quote = TRUE)
manifest <- data.frame(
  artifact = c("mv08i-landscape-execution.csv", "mv08i-landscape-level-summary.csv"),
  purpose = c("one row per unit, view, and homology dimension", "one row per retained active landscape level"),
  public = TRUE, contains_private_diagrams = FALSE, stringsAsFactors = FALSE)
write.csv(manifest, file.path(output_dir, "mv08i-landscape-artifact-manifest.csv"), row.names = FALSE, quote = TRUE)
writeLines(c(
  paste0("# MV8-I common475 HCA all-active-level landscapes (", format(Sys.Date()), ")"), "",
  "This is label-closed technical external-validation evidence for eight adult bone-marrow HCA units.", "",
  "- Cell and gene views are reported separately; H0 and H1 are integrated separately.",
  "- All active consecutive levels are retained; no universal level cap is used.",
  "- Exact critical-breakpoint integration uses births, midpoints, deaths, and line crossings; no uniform grid is used.",
  "- Essential H0 intervals are excluded before landscape construction.",
  "- Labels, outcomes, fusion, exact-500 recovery, and deletion remain closed.", "",
  paste0("Private input files: ", length(files), " per-unit diagram receipts (not copied to public output).")),
  file.path(output_dir, paste0("MV08I_LANDSCAPES_", format(Sys.Date()), ".md")))
cat("completed MV8-I landscapes: ", nrow(execution), " profiles; ", nrow(levels), " active levels\n", sep = "")
