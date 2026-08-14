#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1L) args[[1]] else {
  file.path("tmp", "realistic-ph-dimension-profile")
}
repeats <- if (length(args) >= 2L) as.integer(args[[2]]) else 2L

if (length(repeats) != 1L || is.na(repeats) || repeats < 1L) {
  stop("repeats must be one positive integer.", call. = FALSE)
}
if (dir.exists(output_root) && length(list.files(
  output_root, all.files = TRUE, no.. = TRUE
)) > 0L) {
  stop("output_root must be empty to prevent stale profile evidence.", call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

rows <- list()
index <- 0L
for (dimension in 0:1) {
  for (repeat_id in seq_len(repeats)) {
    run_dir <- file.path(
      output_root,
      paste0("h", dimension, "-repeat-", repeat_id)
    )
    elapsed <- system.time({
      result <- scPHcompare::run_realistic_fixture(
        run_dir,
        seed = 20260805L,
        num_cores = 1L,
        max_dimension = dimension,
        ph_poll_interval = 0.25,
        ph_max_time_per_sample = 20 * 60
      )
    })[["elapsed"]]
    index <- index + 1L
    rows[[index]] <- cbind(
      repeat_id = repeat_id,
      total_elapsed_seconds = as.numeric(elapsed),
      result$manifest
    )
  }
}

profile <- do.call(rbind, rows)
if (!all(profile$ph_completed == profile$samples_eligible)) {
  stop("At least one profile run did not complete every eligible sample.",
       call. = FALSE)
}
if (length(unique(profile$input_sha256)) != 1L ||
    length(unique(profile$harmony_sha256_rounded_8)) != 1L) {
  stop("Deterministic fixture inputs or Harmony matrices changed across runs.",
       call. = FALSE)
}
h0 <- profile[profile$ph_max_dimension == 0L, , drop = FALSE]
h1 <- profile[profile$ph_max_dimension == 1L, , drop = FALSE]
if (length(unique(h0$ph_sha256_rounded_8)) != 1L ||
    length(unique(h1$ph_sha256_rounded_8)) != 1L) {
  stop("Persistence diagrams were not repeatable within a dimension.",
       call. = FALSE)
}
if (!identical(unique(h0$h0_feature_counts), unique(h1$h0_feature_counts))) {
  stop("H0 features changed when H1 was requested.", call. = FALSE)
}
h1_counts <- as.integer(unlist(strsplit(h1$h1_feature_counts, ";", fixed = TRUE)))
if (anyNA(h1_counts) || any(h1_counts <= 0L)) {
  stop("The H1 fixture did not produce positive H1 feature counts.",
       call. = FALSE)
}

profile_file <- file.path(output_root, "ph_dimension_profile.csv")
utils::write.csv(profile, profile_file, row.names = FALSE, na = "")
summary <- do.call(rbind, lapply(split(profile, profile$ph_max_dimension), function(x) {
  data.frame(
    ph_max_dimension = x$ph_max_dimension[[1]],
    repeats = nrow(x),
    elapsed_median_seconds = stats::median(x$total_elapsed_seconds),
    elapsed_max_seconds = max(x$total_elapsed_seconds),
    ph_elapsed_median_seconds = stats::median(as.numeric(unlist(strsplit(
      x$ph_elapsed_seconds, ";", fixed = TRUE
    )))),
    process_tree_peak_rss_max_bytes = max(as.numeric(unlist(strsplit(
      x$ph_process_tree_peak_rss_bytes, ";", fixed = TRUE
    )))),
    h0_feature_counts = x$h0_feature_counts[[1]],
    h1_feature_counts = x$h1_feature_counts[[1]],
    ph_sha256_rounded_8 = x$ph_sha256_rounded_8[[1]],
    stringsAsFactors = FALSE
  )
}))
summary_file <- file.path(output_root, "ph_dimension_profile_summary.csv")
utils::write.csv(summary, summary_file, row.names = FALSE, na = "")
print(summary)
cat("Profile:", normalizePath(profile_file, winslash = "/"), "\n")
cat("Summary:", normalizePath(summary_file, winslash = "/"), "\n")
