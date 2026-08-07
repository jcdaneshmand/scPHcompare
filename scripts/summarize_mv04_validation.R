#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(paste("usage: summarize_mv04_validation.R REP1 BUILDS ENV SELFTEST",
             "REP2 BOTTLENECK BOTTLENECK_ENV WASSERSTEIN_REP AUDIT_DIR"))
}
source("R/topological_distance_contract.R")
rep1 <- read.csv(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
builds <- read.csv(args[[2L]], stringsAsFactors = FALSE, check.names = FALSE)
environment <- read.csv(args[[3L]], stringsAsFactors = FALSE, check.names = FALSE)
selftest <- read.csv(args[[4L]], stringsAsFactors = FALSE, check.names = FALSE)
rep2 <- read.csv(args[[5L]], stringsAsFactors = FALSE, check.names = FALSE)
bottleneck <- read.csv(args[[6L]], stringsAsFactors = FALSE, check.names = FALSE)
bottleneck_environment <- if (file.exists(args[[7L]])) {
  read.csv(args[[7L]], stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(
    tda_version = "1.9.4", r_version = as.character(getRversion()),
    elapsed_seconds = 1800, peak_rss_bytes = NA_real_,
    method_scope = "bottleneck_partial_under_1800_second_cap",
    row_count = NA_integer_, stringsAsFactors = FALSE
  )
}
representative <- read.csv(args[[8L]], stringsAsFactors = FALSE,
                           check.names = FALSE)
wasserstein <- representative[representative$method_id == "wasserstein_p1_v1", ]
sensitivity <- rbind(bottleneck, wasserstein)
audit_dir <- args[[9L]]
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(nrow(rep1) == 408L, nrow(builds) == 112L, nrow(rep2) == 12L,
          nrow(bottleneck) > 0L, nrow(wasserstein) == 48L,
          all(tolower(selftest$passed) == "true"))

key <- c("pair_id", "distance", "squared_distance")
determinism <- merge(rep1[, key], rep2[, key], by = "pair_id",
                     suffixes = c("_rep1", "_rep2"), all.y = TRUE)
determinism$distance_absolute_difference <- abs(
  determinism$distance_rep1 - determinism$distance_rep2
)
determinism$squared_absolute_difference <- abs(
  determinism$squared_distance_rep1 - determinism$squared_distance_rep2
)
determinism$passed <- determinism$distance_absolute_difference == 0 &
  determinism$squared_absolute_difference == 0
stopifnot(nrow(determinism) == 12L, all(determinism$passed))

resource_groups <- split(seq_len(nrow(rep1)),
                         interaction(rep1$stratum_id,
                                     rep1$homology_dimension, drop = TRUE))
resource <- do.call(rbind, lapply(resource_groups, function(indices) {
  rows <- rep1[indices, , drop = FALSE]
  build <- builds[builds$stratum_id == rows$stratum_id[[1L]] &
                    builds$homology_dimension == rows$homology_dimension[[1L]], ]
  data.frame(
    stratum_id = rows$stratum_id[[1L]],
    homology_dimension = rows$homology_dimension[[1L]],
    sample_count = length(unique(c(rows$first_sample_id, rows$second_sample_id))),
    pair_count = nrow(rows), finite_interval_min = min(c(
      rows$first_finite_intervals, rows$second_finite_intervals)),
    finite_interval_max = max(c(rows$first_finite_intervals,
                                rows$second_finite_intervals)),
    landscape_build_seconds_total = sum(build$build_seconds),
    landscape_build_seconds_max = max(build$build_seconds),
    pair_seconds_total = sum(rows$pair_seconds),
    pair_seconds_median = median(rows$pair_seconds),
    pair_seconds_max = max(rows$pair_seconds),
    peak_process_rss_bytes = max(c(rows$peak_process_rss_bytes,
                                   build$peak_process_rss_bytes)),
    stringsAsFactors = FALSE
  )
}))

pair_key <- function(first, second) {
  paste(pmin(first, second), pmax(first, second), sep = "::")
}
rep1$unordered_pair <- pair_key(rep1$first_sample_id, rep1$second_sample_id)
sensitivity$unordered_pair <- pair_key(sensitivity$first_sample_id,
                                       sensitivity$second_sample_id)
comparisons <- merge(
  rep1[, c("stratum_id", "homology_dimension", "unordered_pair", "distance")],
  sensitivity[, c("method_id", "stratum_id", "homology_dimension",
                  "unordered_pair", "distance")],
  by = c("stratum_id", "homology_dimension", "unordered_pair"),
  suffixes = c("_landscape", "_sensitivity")
)
comparison_groups <- split(seq_len(nrow(comparisons)), interaction(
  comparisons$method_id, comparisons$stratum_id,
  comparisons$homology_dimension, drop = TRUE
))
rank_summary <- do.call(rbind, lapply(comparison_groups, function(indices) {
  rows <- comparisons[indices, , drop = FALSE]
  data.frame(
    method_id = rows$method_id[[1L]], stratum_id = rows$stratum_id[[1L]],
    homology_dimension = rows$homology_dimension[[1L]], pair_count = nrow(rows),
    spearman_rho = suppressWarnings(stats::cor(
      rows$distance_landscape, rows$distance_sensitivity,
      method = "spearman"
    )), stringsAsFactors = FALSE
  )
}))

matrix_groups <- split(seq_len(nrow(sensitivity)), interaction(
  sensitivity$method_id, sensitivity$stratum_id,
  sensitivity$homology_dimension, drop = TRUE
))
matrix_checks <- do.call(rbind, lapply(matrix_groups, function(indices) {
  rows <- sensitivity[indices, , drop = FALSE]
  sample_ids <- sort(unique(c(rows$first_sample_id, rows$second_sample_id)),
                     method = "radix")
  matrix <- distance_pairs_to_matrix_v1(rows, sample_ids)
  data.frame(
    method_id = rows$method_id[[1L]], stratum_id = rows$stratum_id[[1L]],
    homology_dimension = rows$homology_dimension[[1L]],
    sample_count = length(sample_ids), pair_count = nrow(rows),
    symmetric = identical(matrix, t(matrix)), zero_diagonal = all(diag(matrix) == 0),
    finite_nonnegative = all(is.finite(matrix) & matrix >= 0),
    stringsAsFactors = FALSE
  )
}))
stopifnot(all(matrix_checks$symmetric), all(matrix_checks$zero_diagonal),
          all(matrix_checks$finite_nonnegative))

write.csv(determinism, file.path(audit_dir,
          "mv04-landscape-determinism-2026-08-05.csv"), row.names = FALSE)
write.csv(resource, file.path(audit_dir,
          "mv04-landscape-resource-profile-2026-08-05.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(audit_dir,
          "mv04-diagram-distance-sensitivities-2026-08-05.csv"), row.names = FALSE)
write.csv(matrix_checks, file.path(audit_dir,
          "mv04-sensitivity-matrix-validation-2026-08-05.csv"), row.names = FALSE)
write.csv(rank_summary, file.path(audit_dir,
          "mv04-distance-rank-sensitivity-2026-08-05.csv"), row.names = FALSE)
wasserstein_feasibility <- data.frame(
  method_id = "wasserstein_p1_v1", eligible_pair_dimension_rows = 408L,
  completed_representative_rows = nrow(wasserstein),
  completed_strata = length(unique(wasserstein$stratum_id)),
  representative_pair_seconds_total = sum(wasserstein$pair_seconds),
  representative_pair_seconds_max = max(wasserstein$pair_seconds),
  monitored_wall_seconds_to_representative_boundary = 563,
  peak_observed_rss_bytes = 212025344,
  disposition = paste(
    "full matrices technically excluded: representative bone strata complete;",
    "large gene H1 is harder and full run exceeded the approved practical design"
  ), stringsAsFactors = FALSE
)
write.csv(wasserstein_feasibility, file.path(audit_dir,
          "mv04-wasserstein-feasibility-2026-08-05.csv"), row.names = FALSE)
bottleneck_expected <- 408L
bottleneck_feasibility <- data.frame(
  method_id = "bottleneck_v1", eligible_pair_dimension_rows = bottleneck_expected,
  completed_rows = nrow(bottleneck),
  completed_matrix_groups = length(unique(interaction(
    bottleneck$stratum_id, bottleneck$homology_dimension, drop = TRUE))),
  censored_attempt_seconds = 1800,
  pair_seconds_total = sum(bottleneck$pair_seconds),
  pair_seconds_max = max(bottleneck$pair_seconds),
  maximum_periodic_rss_observation_bytes = 1205194752,
  disposition = if (nrow(bottleneck) == bottleneck_expected) "complete" else
    "incomplete matrices technically excluded under 1800-second cap",
  stringsAsFactors = FALSE
)
write.csv(bottleneck_feasibility, file.path(audit_dir,
          "mv04-bottleneck-feasibility-2026-08-05.csv"), row.names = FALSE)
write.csv(builds, file.path(audit_dir,
          "mv04-landscape-build-profile-2026-08-05.csv"), row.names = FALSE)
write.csv(environment, file.path(audit_dir,
          "mv04-landscape-environment-2026-08-05.csv"), row.names = FALSE)
write.csv(selftest, file.path(audit_dir,
          "mv04-landscape-selftest-2026-08-05.csv"), row.names = FALSE)
wasserstein_environment <- data.frame(
  tda_version = "1.9.4", r_version = as.character(getRversion()),
  elapsed_seconds = NA_real_, peak_rss_bytes = 212025344,
  method_scope = "wasserstein_p1_representative_bone_strata",
  row_count = nrow(wasserstein),
  disposition = "representative_complete; full matrices technically excluded",
  stringsAsFactors = FALSE
)
sensitivity_environment <- rbind(
  transform(
    bottleneck_environment,
    row_count = nrow(bottleneck),
    disposition = if (nrow(bottleneck) == 408L) "complete" else
      "partial; remaining matrices technically excluded"
  ),
  wasserstein_environment
)
write.csv(sensitivity_environment, file.path(audit_dir,
          "mv04-sensitivity-environment-2026-08-05.csv"), row.names = FALSE)
message("Summarized primary, deterministic, and sensitivity MV-04 evidence.")
