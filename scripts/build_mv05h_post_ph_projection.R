#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: build_mv05h_post_ph_projection.R MV5G_METRICS MV5H_METRICS ",
    "D4_METRICS D5_METRICS COMPONENT_OUT SUMMARY_OUT DECISION_OUT",
    call. = FALSE
  )
}
read_public <- function(path) utils::read.csv(
  normalizePath(path, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE, check.names = FALSE
)
source("R/provenance_utils.R")
g <- read_public(args[[1L]])
h <- read_public(args[[2L]])
d4 <- read_public(args[[3L]])
d5 <- read_public(args[[4L]])
for (value in list(g, h, d4, d5)) {
  if (any(value$outcome_label_state != "closed") ||
      any(as.logical(value$biological_outcomes_computed))) {
    stop("Post-PH projection inputs must remain label closed.", call. = FALSE)
  }
}
if (nrow(g) != 75L || nrow(h) != 75L ||
    any(g$disposition != "completed") || any(h$disposition != "completed")) {
  stop("Measured coordinate and PH stages must be complete.", call. = FALSE)
}
reserve <- 1.25
components <- data.frame(
  contract_id = "mv05h_post_ph_projection_components_v1",
  component = c(
    "integrated_coordinates_measured", "integrated_cell_ph_measured",
    "integrated_landscape_distances_projected",
    "integrated_baseline_retrieval_projected"
  ),
  evidence_basis = c(
    "complete_mv05g_measured", "complete_mv05h_measured",
    "accepted_sct_D4_with_25_percent_reserve",
    "accepted_sct_D5_with_25_percent_reserve"
  ),
  worker_seconds = c(
    sum(g$elapsed_seconds), sum(h$elapsed_seconds),
    sum(d4$elapsed_seconds) * reserve, sum(d5$elapsed_seconds) * reserve
  ),
  stage_complete = c(TRUE, TRUE, FALSE, FALSE),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
components$worker_hours <- components$worker_seconds / 3600
measured_storage <- sum(g$private_result_bytes) +
  sum(h$private_result_bytes)
future_storage <- sum(c(
  d4$private_result_bytes, d5$private_result_bytes
)) * reserve
peak_rss <- max(c(
  g$peak_process_tree_rss_bytes, h$peak_process_tree_rss_bytes
))
max_group_seconds <- max(c(g$elapsed_seconds, h$elapsed_seconds))
cap <- 21.6
total_hours <- sum(components$worker_hours)
total_storage <- measured_storage + future_storage
decision <- if (
  total_hours <= cap && total_storage <= 10 * 1000^3 &&
    peak_rss <= 8 * 1024^3 && max_group_seconds <= 1800
) "authorize_separate_integrated_landscape_distance_sprint" else
  "do_not_authorize_or_narrow_scope"
summary <- data.frame(
  contract_id = "mv05h_post_ph_projection_summary_v1",
  coordinate_groups_complete = nrow(g), ph_groups_complete = nrow(h),
  coordinate_views_complete = 6750L, ph_views_complete = 6750L,
  measured_coordinate_worker_hours = components$worker_hours[[1L]],
  measured_ph_worker_hours = components$worker_hours[[2L]],
  projected_future_worker_hours = sum(components$worker_hours[3:4]),
  projected_total_worker_hours = total_hours, worker_hour_cap = cap,
  measured_coordinate_ph_storage_bytes = measured_storage,
  projected_future_storage_bytes = future_storage,
  projected_total_storage_bytes = total_storage,
  storage_cap_bytes = 10 * 1000^3,
  measured_peak_rss_bytes = peak_rss, rss_cap_bytes = 8 * 1024^3,
  measured_max_group_seconds = max_group_seconds,
  group_cap_seconds = 1800, decision = decision,
  next_stage_authorized = decision ==
    "authorize_separate_integrated_landscape_distance_sprint",
  current_sprint_landscape_jobs = 0L,
  current_sprint_distance_jobs = 0L,
  biological_outcomes_computed = FALSE, outcome_label_state = "closed",
  stringsAsFactors = FALSE
)
authorization <- summary[c(
  "decision", "next_stage_authorized", "projected_total_worker_hours",
  "worker_hour_cap", "projected_total_storage_bytes", "storage_cap_bytes",
  "measured_peak_rss_bytes", "rss_cap_bytes", "measured_max_group_seconds",
  "group_cap_seconds", "current_sprint_landscape_jobs",
  "current_sprint_distance_jobs", "biological_outcomes_computed",
  "outcome_label_state"
)]
authorization$contract_id <-
  "mv05h_integrated_landscape_authorization_gate_v1"
authorization <- authorization[c("contract_id", setdiff(
  names(authorization), "contract_id"
))]
for (path in args[5:7]) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
write_provenance_csv(components, args[[5L]])
write_provenance_csv(summary, args[[6L]])
write_provenance_csv(authorization, args[[7L]])
message("Wrote mandatory MV5-H post-PH projection and stop decision.")
