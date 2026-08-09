#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(
    "usage: build_mv05f_resource_projection.R PILOT_METRICS D1_RESOURCES ",
    "D3_METRICS D4_METRICS D5_METRICS MV5C_STATUS COMPONENT_OUTPUT ",
    "SUMMARY_OUTPUT DECISION_OUTPUT", call. = FALSE
  )
}
read_public <- function(path) utils::read.csv(
  normalizePath(path, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE, check.names = FALSE
)
source("R/provenance_utils.R")
source("R/mv05f_integration_gate.R")
projection <- mv05f_project_full_workload_v2(
  read_public(args[[1L]]), read_public(args[[2L]]), read_public(args[[3L]]),
  read_public(args[[4L]]), read_public(args[[5L]]), read_public(args[[6L]])
)
for (path in args[7:9]) dir.create(dirname(path), recursive = TRUE,
                                   showWarnings = FALSE)
write_provenance_csv(projection$components, args[[7L]])
write_provenance_csv(projection$summary, args[[8L]])
decision <- projection$summary[c(
  "contract_id", "decision", "full_execution_authorized",
  "projected_worker_hours", "worker_hour_cap", "projected_storage_bytes",
  "storage_cap_bytes", "projected_peak_rss_bytes", "rss_cap_bytes",
  "projected_max_group_seconds", "group_cap_seconds",
  "biological_outcomes_computed", "outcome_label_state"
)]
decision$contract_id <- "mv05f_full_execution_authorization_gate_v1"
decision$next_action <- if (
  decision$decision ==
    "go_separately_authorized_full_label_closed_integrated_execution"
) "owner_may_authorize_separate_full_execution_goal" else
  "revise_scope_before_any_full_execution"
write_provenance_csv(decision, args[[9L]])
message("Wrote MV5-F full integrated workload projection and stop decision.")
