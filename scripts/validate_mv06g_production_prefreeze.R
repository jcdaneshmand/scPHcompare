#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: validate_mv06g_production_prefreeze.R PARENT WORKLOAD ",
       "RUST_LIBRARY EVIDENCE_DIR OUTPUT", call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06g_production.R")
parent <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
workload <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                            check.names = FALSE)
policy <- utils::read.csv(file.path(args[[4L]], "mv06g-production-policy.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
sources <- utils::read.csv(file.path(args[[4L]], "mv06g-production-sources.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
production_queue <- utils::read.csv(
  file.path(args[[4L]], "mv06g-production-queue.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                       algo = "sha256", serialize = TRUE)
checks <- data.frame(
  category = c("parent_identity", "remaining_group_axis", "training_workload",
               "query_workload", "resource_policy", "implementation_identity",
               "stage1_equivalence_gate", "label_firewall",
               "production_closed"),
  passed = c(
    policy$parent_contract_sha256 == .mv06f_sha256(args[[1L]]) &&
      policy$queue_root_sha256 == parent$queue_root_sha256 &&
      policy$rust_library_sha256 == .mv06f_sha256(args[[3L]]),
    nrow(production_queue) == 74L &&
      identical(as.integer(production_queue$execution_order), 2:75) &&
      !anyDuplicated(production_queue$group_id),
    sum(production_queue$training_biological_pairs) == 260595L &&
      sum(production_queue$training_component_rows) == 1042380L,
    sum(production_queue$query_biological_pairs) == 33725L &&
      sum(production_queue$query_ranking_rows) == 303525L,
    policy$maximum_workers == 1L && !as.logical(policy$automatic_retry) &&
      policy$elapsed_cap_seconds_per_group == 1800L &&
      policy$rss_cap_bytes_per_group == 12 * 1024^3 &&
      policy$private_storage_cap_bytes == 5 * 1024^3 &&
      policy$total_worker_cap_seconds == 12 * 3600,
    identical(sources$path, mv06g_production_source_paths_v1()) &&
      all(file.exists(sources$path)) &&
      identical(unname(vapply(sources$path, .mv06f_sha256, character(1L))),
                unname(sources$sha256)) &&
      policy$production_implementation_root_sha256 == root,
    as.logical(policy$rebind_equivalence_authorized) &&
      policy$required_equivalence_artifacts == 3L,
    policy$outcome_label_state == "closed" &&
      !as.logical(policy$biological_outcomes_computed) &&
      policy$fusion_evaluations == 0L && policy$outcome_jobs == 0L,
    !as.logical(policy$remaining_production_authorized) &&
      policy$disposition == "prefreeze_pass_rebind_equivalence_only"
  ), stringsAsFactors = FALSE
)
checks$contract_id <- "mv06g_production_prefreeze_validation_v1"
checks$outcome_label_state <- "closed"
checks$biological_outcomes_computed <- FALSE
utils::write.csv(checks, args[[5L]], row.names = FALSE, na = "")
if (any(!checks$passed)) stop("MV6-G production prefreeze validation failed.",
                              call. = FALSE)
message("Validated MV6-G production prefreeze: 9/9 pass.")
