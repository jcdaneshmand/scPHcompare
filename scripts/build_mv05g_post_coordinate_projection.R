#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "usage: build_mv05g_post_coordinate_projection.R MV5G_METRICS D3_METRICS ",
    "D4_METRICS D5_METRICS COMPONENT_OUTPUT SUMMARY_OUTPUT DECISION_OUTPUT",
    call. = FALSE
  )
}
read_public <- function(path) utils::read.csv(
  normalizePath(path, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE, check.names = FALSE
)
source("R/provenance_utils.R")
source("R/mv05f_integration_gate.R")
source("R/mv05g_coordinate_production.R")
g <- read_public(args[[1L]])
d3 <- read_public(args[[2L]])
d4 <- read_public(args[[3L]])
d5 <- read_public(args[[4L]])
mv05g_validate_group_metrics_v1(g)
for (value in list(d3, d4, d5)) {
  if (any(value$outcome_label_state != "closed") ||
      any(as.logical(value$biological_outcomes_computed))) {
    stop("Post-coordinate projection inputs must be label closed.", call. = FALSE)
  }
}
reserve <- 1.25
components <- data.frame(
  contract_id = "mv05g_post_coordinate_projection_components_v1",
  component = c(
    "integrated_coordinates_measured", "integrated_cell_ph_projected",
    "integrated_landscape_distances_projected",
    "integrated_baseline_retrieval_projected"
  ),
  evidence_basis = c(
    "complete_mv05g_measured", "accepted_sct_D3_with_25_percent_reserve",
    "accepted_sct_D4_with_25_percent_reserve",
    "accepted_sct_D5_with_25_percent_reserve"
  ),
  worker_seconds = c(
    sum(g$elapsed_seconds), sum(d3$elapsed_seconds) * reserve,
    sum(d4$elapsed_seconds) * reserve, sum(d5$elapsed_seconds) * reserve
  ),
  stage_complete = c(TRUE, FALSE, FALSE, FALSE),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
components$worker_hours <- components$worker_seconds / 3600
future_hours <- sum(components$worker_hours[-1L])
total_hours <- sum(components$worker_hours)
future_storage <- sum(c(
  d3$private_result_bytes, d4$private_result_bytes, d5$private_result_bytes
)) * reserve
measured_storage <- sum(g$private_result_bytes)
cap <- 21.6
decision <- if (
  total_hours <= cap && measured_storage + future_storage <= 10 * 1000^3 &&
    max(g$peak_process_tree_rss_bytes) <= 8 * 1024^3 &&
    max(g$elapsed_seconds) <= 1800
) "authorize_separate_integrated_cell_ph_sprint" else
  "do_not_authorize_or_narrow_scope"
summary <- data.frame(
  contract_id = "mv05g_post_coordinate_projection_summary_v1",
  coordinate_groups_complete = nrow(g), coordinate_views_complete = 6750L,
  query_mappings_complete = sum(g$completed_query_mappings),
  measured_coordinate_worker_hours = components$worker_hours[[1L]],
  projected_future_worker_hours = future_hours,
  projected_total_worker_hours = total_hours, worker_hour_cap = cap,
  measured_coordinate_storage_bytes = measured_storage,
  projected_future_storage_bytes = future_storage,
  projected_total_storage_bytes = measured_storage + future_storage,
  storage_cap_bytes = 10 * 1000^3,
  measured_peak_rss_bytes = max(g$peak_process_tree_rss_bytes),
  rss_cap_bytes = 8 * 1024^3,
  measured_max_group_seconds = max(g$elapsed_seconds),
  group_cap_seconds = 1800, decision = decision,
  next_stage_authorized = decision ==
    "authorize_separate_integrated_cell_ph_sprint",
  current_sprint_ph_jobs = 0L, biological_outcomes_computed = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
authorization <- summary[c(
  "contract_id", "decision", "next_stage_authorized",
  "projected_total_worker_hours", "worker_hour_cap",
  "projected_total_storage_bytes", "storage_cap_bytes",
  "measured_peak_rss_bytes", "rss_cap_bytes", "measured_max_group_seconds",
  "group_cap_seconds", "current_sprint_ph_jobs",
  "biological_outcomes_computed", "outcome_label_state"
)]
authorization$contract_id <- "mv05g_integrated_ph_authorization_gate_v1"
for (path in args[5:7]) dir.create(dirname(path), recursive = TRUE,
                                   showWarnings = FALSE)
write_provenance_csv(components, args[[5L]])
write_provenance_csv(summary, args[[6L]])
write_provenance_csv(authorization, args[[7L]])
message("Wrote mandatory MV5-G post-coordinate projection and stop decision.")
