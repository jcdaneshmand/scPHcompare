#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: build_mv05i_post_projection.R PREVIOUS COMPLETION ",
    "RETRIEVAL_REFERENCE SUMMARY_OUTPUT DECISION_OUTPUT",
    call. = FALSE
  )
}
source("R/provenance_utils.R")
source("R/mv05i_post_projection.R")
previous <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
completion <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                              check.names = FALSE)
retrieval <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
result <- mv05i_measured_primary_projection_v1(previous, completion, retrieval)
decision <- result[c(
  "decision", "next_stage_authorized", "projected_total_worker_hours",
  "worker_hour_cap", "projected_total_storage_bytes", "storage_cap_bytes",
  "measured_peak_rss_bytes", "rss_cap_bytes", "measured_max_group_seconds",
  "group_cap_seconds", "current_sprint_retrieval_jobs",
  "current_sprint_clustering_jobs", "current_sprint_gene_view_jobs",
  "current_sprint_fusion_jobs", "current_sprint_new_data_jobs",
  "biological_outcomes_computed", "outcome_label_state"
)]
decision$contract_id <- "mv05i_integrated_retrieval_authorization_gate_v1"
decision <- decision[c("contract_id", setdiff(names(decision), "contract_id"))]
write_provenance_csv(result, args[[4L]])
write_provenance_csv(decision, args[[5L]])
message("Wrote mandatory MV5-I post-distance projection and stop decision.")
