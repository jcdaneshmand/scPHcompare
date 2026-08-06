#!/usr/bin/env Rscript

parse_args <- function(args) {
  values <- list(
    output = "",
    report = "",
    repetitions = 3L,
    seed = 20260805L
  )
  for (arg in args) {
    if (startsWith(arg, "--output=")) {
      values$output <- sub("^--output=", "", arg)
    } else if (startsWith(arg, "--report=")) {
      values$report <- sub("^--report=", "", arg)
    } else if (startsWith(arg, "--repetitions=")) {
      values$repetitions <- as.integer(sub("^--repetitions=", "", arg))
    } else if (startsWith(arg, "--seed=")) {
      values$seed <- as.integer(sub("^--seed=", "", arg))
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  values
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!nzchar(args$output)) {
  stop("--output must be supplied explicitly.", call. = FALSE)
}
if (is.na(args$repetitions) || args$repetitions < 2L) {
  stop("--repetitions must be at least 2.", call. = FALSE)
}
if (is.na(args$seed)) {
  stop("--seed must be integer-compatible.", call. = FALSE)
}
if (dir.exists(args$output) &&
    length(list.files(args$output, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("--output must be empty to prevent mixed benchmark runs.", call. = FALSE)
}
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args$output, winslash = "/", mustWork = TRUE)

runs <- vector("list", args$repetitions)
manifests <- vector("list", args$repetitions)
for (repetition in seq_len(args$repetitions)) {
  run_dir <- file.path(output, sprintf("repetition-%02d", repetition))
  started_at <- Sys.time()
  result <- scPHcompare::run_toy_baseline(
    output_dir = run_dir,
    seed = args$seed
  )
  finished_at <- Sys.time()
  timings <- result$timings
  timings$repetition <- repetition
  end_to_end_seconds <- as.numeric(difftime(finished_at, started_at, units = "secs"))
  measured_stage_seconds <- sum(timings$elapsed_seconds)
  timings <- rbind(
    timings,
    data.frame(
      stage = "orchestration_and_namespace_load",
      sample_id = "",
      started_at = format(started_at, "%Y-%m-%dT%H:%M:%OS3%z"),
      finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%OS3%z"),
      elapsed_seconds = max(0, end_to_end_seconds - measured_stage_seconds),
      status = "completed",
      error_message = "",
      repetition = repetition,
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "end_to_end",
      sample_id = "",
      started_at = format(started_at, "%Y-%m-%dT%H:%M:%OS3%z"),
      finished_at = format(finished_at, "%Y-%m-%dT%H:%M:%OS3%z"),
      elapsed_seconds = end_to_end_seconds,
      status = "completed",
      error_message = "",
      repetition = repetition,
      stringsAsFactors = FALSE
    )
  )
  runs[[repetition]] <- timings
  manifests[[repetition]] <- result$manifest
}

reference_manifest <- manifests[[1L]]
identical_manifests <- vapply(
  manifests,
  identical,
  logical(1),
  y = reference_manifest
)
if (!all(identical_manifests)) {
  stop("Scientific manifests differed across repeated baseline runs.", call. = FALSE)
}

timings <- do.call(rbind, runs)
group_key <- interaction(timings$stage, timings$sample_id, drop = TRUE, lex.order = TRUE)
groups <- split(timings$elapsed_seconds, group_key)
group_rows <- split(timings[c("stage", "sample_id")], group_key)
summary_rows <- lapply(names(groups), function(key) {
  values <- groups[[key]]
  labels <- group_rows[[key]][1L, , drop = FALSE]
  data.frame(
    stage = labels$stage,
    sample_id = labels$sample_id,
    repetitions = length(values),
    mean_seconds = mean(values),
    median_seconds = stats::median(values),
    sd_seconds = if (length(values) > 1L) stats::sd(values) else NA_real_,
    min_seconds = min(values),
    max_seconds = max(values),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary <- summary[order(summary$stage, summary$sample_id), , drop = FALSE]

timings_file <- file.path(output, "profile_timings.csv")
summary_file <- file.path(output, "profile_summary.csv")
environment_file <- file.path(output, "profile_environment.csv")
utils::write.csv(timings, timings_file, row.names = FALSE, na = "")
utils::write.csv(summary, summary_file, row.names = FALSE, na = "")

environment <- data.frame(
  key = c(
    "timestamp", "r_version", "platform", "os", "machine", "logical_cores",
    "scPHcompare_version", "ripserr_version", "seed", "repetitions",
    "scientific_manifests_identical"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
    R.version.string,
    R.version$platform,
    unname(Sys.info()[["sysname"]]),
    unname(Sys.info()[["machine"]]),
    as.character(parallel::detectCores(logical = TRUE)),
    as.character(utils::packageVersion("scPHcompare")),
    as.character(utils::packageVersion("ripserr")),
    as.character(args$seed),
    as.character(args$repetitions),
    as.character(all(identical_manifests))
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(environment, environment_file, row.names = FALSE, na = "")

if (nzchar(args$report)) {
  dir.create(dirname(args$report), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, args$report, row.names = FALSE, na = "")
}

print(summary, row.names = FALSE)
cat("scientific_manifests_identical=TRUE\n")
cat("profile_output=", output, "\n", sep = "")
