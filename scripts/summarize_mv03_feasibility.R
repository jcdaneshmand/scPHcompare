#!/usr/bin/env Rscript

options(warn = 2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "Usage: summarize_mv03_feasibility.R <stage-a-metrics> ",
    "<stage-b-metrics> <stage-c-metrics> <stage-c-matched-cells> ",
    "<diagram-summary-csv> <stability-csv> <resource-summary-csv>",
    call. = FALSE
  )
}

metric_paths <- args[1:3]
matched_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
diagram_output <- args[[5L]]
stability_output <- args[[6L]]
resource_output <- args[[7L]]
dir.create(dirname(diagram_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(stability_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(resource_output), recursive = TRUE, showWarnings = FALSE)

metric_tables <- lapply(metric_paths, function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
})
metric_columns <- unique(unlist(lapply(metric_tables, names)))
metric_tables <- lapply(metric_tables, function(value) {
  missing <- setdiff(metric_columns, names(value))
  for (column in missing) value[[column]] <- NA
  value[metric_columns]
})
metrics <- do.call(rbind, metric_tables)
if (any(metrics$disposition != "completed") ||
    any(!file.exists(metrics$result_file))) {
  stop("Feasibility summary requires every staged job to have completed.",
       call. = FALSE)
}

diagram_rows <- list()
for (index in seq_len(nrow(metrics))) {
  metric <- metrics[index, , drop = FALSE]
  result <- readRDS(metric$result_file)
  for (dimension in 0:1) {
    intervals <- result$diagram[
      result$diagram[, "dimension"] == dimension &
        is.finite(result$diagram[, "death"]), , drop = FALSE
    ]
    persistence <- intervals[, "death"] - intervals[, "birth"]
    diagram_rows[[length(diagram_rows) + 1L]] <- data.frame(
      stage = metric$stage,
      cohort = metric$cohort,
      representation = metric$representation,
      sample_id = metric$sample_id,
      seed = metric$seed,
      view_id = metric$view_id,
      homology_dimension = paste0("H", dimension),
      finite_interval_count = length(persistence),
      total_persistence = sum(persistence),
      mean_persistence = if (length(persistence)) mean(persistence) else 0,
      median_persistence = if (length(persistence))
        stats::median(persistence) else 0,
      maximum_persistence = if (length(persistence)) max(persistence) else 0,
      squared_persistence_sum = sum(persistence ^ 2),
      essential_interval_count = sum(
        result$diagram[, "dimension"] == dimension &
          is.infinite(result$diagram[, "death"])
      ),
      diagram_sha256 = result$provenance$diagram_sha256,
      stringsAsFactors = FALSE
    )
  }
}
diagram_summary <- do.call(rbind, diagram_rows)
utils::write.csv(diagram_summary, diagram_output, row.names = FALSE)

coefficient_of_variation <- function(value) {
  value <- as.numeric(value)
  if (length(value) < 2L || !is.finite(mean(value)) || mean(value) == 0) {
    return(NA_real_)
  }
  stats::sd(value) / abs(mean(value))
}
relative_range <- function(value) {
  value <- as.numeric(value)
  if (!length(value) || !is.finite(mean(value)) || mean(value) == 0) {
    return(NA_real_)
  }
  diff(range(value)) / abs(mean(value))
}

stage_c <- diagram_summary[diagram_summary$stage == "C", , drop = FALSE]
matched <- utils::read.csv(
  matched_path, stringsAsFactors = FALSE, check.names = FALSE
)
matched_groups <- split(
  matched,
  interaction(matched$cohort, matched$sample_id, drop = TRUE, lex.order = TRUE)
)
cell_overlap <- lapply(matched_groups, function(value) {
  sets <- split(value$cell_id, value$seed)
  pairs <- utils::combn(names(sets), 2L, simplify = FALSE)
  jaccard <- vapply(pairs, function(pair) {
    left <- unique(sets[[pair[[1L]]]])
    right <- unique(sets[[pair[[2L]]]])
    length(intersect(left, right)) / length(union(left, right))
  }, numeric(1L))
  data.frame(
    cohort = value$cohort[[1L]], sample_id = value$sample_id[[1L]],
    selected_cell_pairwise_jaccard_min = min(jaccard),
    selected_cell_pairwise_jaccard_mean = mean(jaccard),
    selected_cell_pairwise_jaccard_max = max(jaccard),
    stringsAsFactors = FALSE
  )
})
cell_overlap <- do.call(rbind, cell_overlap)
groups <- interaction(
  stage_c$cohort, stage_c$representation, stage_c$sample_id,
  stage_c$view_id, stage_c$homology_dimension, drop = TRUE, lex.order = TRUE
)
stability_rows <- lapply(split(stage_c, groups), function(value) {
  overlap <- cell_overlap[
    cell_overlap$cohort == value$cohort[[1L]] &
      cell_overlap$sample_id == value$sample_id[[1L]], , drop = FALSE
  ]
  data.frame(
    cohort = value$cohort[[1L]],
    representation = value$representation[[1L]],
    sample_id = value$sample_id[[1L]],
    view_id = value$view_id[[1L]],
    homology_dimension = value$homology_dimension[[1L]],
    seeds = nrow(value),
    unique_diagram_hashes = length(unique(value$diagram_sha256)),
    selected_cell_pairwise_jaccard_min =
      overlap$selected_cell_pairwise_jaccard_min[[1L]],
    selected_cell_pairwise_jaccard_mean =
      overlap$selected_cell_pairwise_jaccard_mean[[1L]],
    selected_cell_pairwise_jaccard_max =
      overlap$selected_cell_pairwise_jaccard_max[[1L]],
    interval_count_min = min(value$finite_interval_count),
    interval_count_max = max(value$finite_interval_count),
    interval_count_cv = coefficient_of_variation(value$finite_interval_count),
    interval_count_relative_range = relative_range(value$finite_interval_count),
    total_persistence_min = min(value$total_persistence),
    total_persistence_max = max(value$total_persistence),
    total_persistence_cv = coefficient_of_variation(value$total_persistence),
    total_persistence_relative_range = relative_range(value$total_persistence),
    maximum_persistence_min = min(value$maximum_persistence),
    maximum_persistence_max = max(value$maximum_persistence),
    maximum_persistence_cv = coefficient_of_variation(value$maximum_persistence),
    maximum_persistence_relative_range = relative_range(value$maximum_persistence),
    interpretation =
      "descriptive_seed_sensitivity_only_no_biological_independence",
    stringsAsFactors = FALSE
  )
})
stability <- do.call(rbind, stability_rows)
utils::write.csv(stability, stability_output, row.names = FALSE)

resource_groups <- interaction(
  metrics$stage, metrics$cohort, metrics$representation, metrics$view_id,
  drop = TRUE, lex.order = TRUE
)
resource_rows <- lapply(split(metrics, resource_groups), function(value) {
  data.frame(
    stage = value$stage[[1L]],
    cohort = value$cohort[[1L]],
    representation = value$representation[[1L]],
    view_id = value$view_id[[1L]],
    jobs = nrow(value),
    completed_jobs = sum(value$disposition == "completed"),
    failed_jobs = sum(value$disposition != "completed"),
    worker_seconds_sum = sum(value$elapsed_seconds),
    elapsed_seconds_median = stats::median(value$elapsed_seconds),
    elapsed_seconds_max = max(value$elapsed_seconds),
    peak_process_tree_rss_bytes_max = max(
      value$peak_process_tree_rss_bytes
    ),
    h1_intervals_min = min(value$h1_intervals),
    h1_intervals_max = max(value$h1_intervals),
    stringsAsFactors = FALSE
  )
})
utils::write.csv(do.call(rbind, resource_rows), resource_output,
                 row.names = FALSE)
message("Summarized ", nrow(metrics), " completed MV-03 PH jobs and ",
        nrow(stability), " five-seed stability strata.")
