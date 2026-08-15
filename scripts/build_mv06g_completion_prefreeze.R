#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop(
  "usage: build_mv06g_completion_prefreeze.R QUEUE PARENT REBIND_POLICY ",
  "REBIND_EQUIVALENCE RUST_LIBRARY PRODUCTION_QUEUE OUTPUT_DIR",
  call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_completion.R")
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
queue <- read_csv(args[[1L]]); parent <- read_csv(args[[2L]])
rebind <- read_csv(args[[3L]]); equivalence <- read_csv(args[[4L]])
production_queue <- read_csv(args[[6L]])
rust <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
paths <- mv06g_completion_source_paths_v1()
sources <- data.frame(
  contract_id = "mv06g_completion_source_inventory_v1", path = paths,
  sha256 = unname(vapply(paths, .mv06f_sha256, character(1L))),
  stringsAsFactors = FALSE
)
root <- digest::digest(stats::setNames(sources$sha256, sources$path),
                       algo = "sha256", serialize = TRUE)
remaining <- queue[queue$execution_order != 1L, , drop = FALSE]
remaining <- remaining[order(remaining$execution_order), , drop = FALSE]
rownames(remaining) <- NULL
remaining <- mv06g_add_runner_schema_v1(remaining)
if (!identical(remaining$group_id, production_queue$group_id) ||
    nrow(equivalence) != 3L ||
    any(!as.logical(equivalence$sha256_identical)) ||
    any(!as.logical(equivalence$bytes_identical)) ||
    any(!as.logical(equivalence$resource_pass)) ||
    !as.logical(rebind$rebind_equivalence_authorized) ||
    as.logical(rebind$remaining_production_authorized) ||
    rebind$production_implementation_root_sha256 !=
      "9bf8614d8e2dbdfce43792e74f08620712674c8830770c7c8d70b1fea432a71c" ||
    rebind$rust_library_sha256 != .mv06f_sha256(rust)) {
  stop("MV6-G completion prefreeze inputs are stale.", call. = FALSE)
}
policy <- data.frame(
  contract_id = "mv06g_serial_completion_prefreeze_v1",
  parent_rebind_policy_sha256 = .mv06f_sha256(args[[3L]]),
  rebind_equivalence_sha256 = .mv06f_sha256(args[[4L]]),
  parent_contract_sha256 = .mv06f_sha256(args[[2L]]),
  queue_root_sha256 = parent$queue_root_sha256,
  scientific_implementation_root_sha256 =
    rebind$production_implementation_root_sha256,
  execution_implementation_root_sha256 = root,
  rust_library_sha256 = rebind$rust_library_sha256,
  remaining_groups = 74L,
  training_biological_pairs = sum(remaining$training_biological_pairs),
  training_component_rows = sum(remaining$training_component_rows),
  query_biological_pairs = sum(remaining$query_biological_pairs),
  query_ranking_rows = sum(remaining$query_ranking_rows),
  elapsed_cap_seconds_per_group = 1800L,
  rss_cap_bytes_per_group = 12 * 1024^3,
  private_storage_cap_bytes = 5 * 1024^3,
  total_worker_cap_seconds = 12 * 3600,
  maximum_workers = 1L, automatic_retry = FALSE,
  remaining_production_authorized = TRUE,
  complete_validation_required = TRUE,
  immutable_resume_required = TRUE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  fusion_evaluations = 0L, outcome_jobs = 0L,
  disposition = "prefreeze_pass_serial_label_closed_only",
  stringsAsFactors = FALSE
)
output <- args[[7L]]; dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(value, name) utils::write.csv(
  value, file.path(output, name), row.names = FALSE, na = ""
)
write_csv(policy, "mv06g-completion-policy.csv")
write_csv(sources, "mv06g-completion-sources.csv")
write_csv(remaining, "mv06g-completion-queue.csv")
message("Prefroze serial completion for 74 MV6-G groups.")
