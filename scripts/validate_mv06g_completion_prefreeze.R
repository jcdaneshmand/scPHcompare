#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop(
  "usage: validate_mv06g_completion_prefreeze.R PARENT REBIND_POLICY ",
  "REBIND_EQUIVALENCE RUST_LIBRARY EVIDENCE_DIR EXPECTED_QUEUE ",
  "PRIVATE_ROOT OUTPUT", call. = FALSE
)
source("R/mv06f_production.R")
source("R/mv06g_completion.R")
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
parent <- read_csv(args[[1L]]); parent$.file_path <- args[[1L]]
rebind <- read_csv(args[[2L]]); rebind$.file_sha256 <- .mv06f_sha256(args[[2L]])
equivalence <- read_csv(args[[3L]])
equivalence$.file_sha256 <- .mv06f_sha256(args[[3L]])
policy <- read_csv(file.path(args[[5L]], "mv06g-completion-policy.csv"))
sources <- read_csv(file.path(args[[5L]], "mv06g-completion-sources.csv"))
queue <- read_csv(file.path(args[[5L]], "mv06g-completion-queue.csv"))
expected <- read_csv(args[[6L]])
mv06g_validate_completion_policy_v1(
  policy, queue, parent, rebind, equivalence, sources, args[[4L]]
)
private_files <- if (dir.exists(args[[7L]])) list.files(
  args[[7L]], recursive = TRUE, all.files = TRUE, no.. = TRUE
) else character()
checks <- data.frame(
  contract_id = "mv06g_completion_prefreeze_validation_v1",
  category = c(
    "rebind_admission", "queue_axis", "training_workload", "query_workload",
    "resource_contract", "scientific_identity", "execution_identity",
    "serial_fail_closed", "label_firewall", "zero_production"
  ),
  passed = c(
    nrow(equivalence) == 3L && all(equivalence$sha256_identical) &&
      all(equivalence$bytes_identical) && all(equivalence$resource_pass),
    identical(queue$group_id, expected$group_id) &&
      identical(as.integer(queue$execution_order), 2:75),
    sum(queue$training_biological_pairs) == 260595L &&
      sum(queue$training_component_rows) == 1042380L,
    sum(queue$query_biological_pairs) == 33725L &&
      sum(queue$query_ranking_rows) == 303525L,
    policy$elapsed_cap_seconds_per_group == 1800L &&
      policy$rss_cap_bytes_per_group == 12 * 1024^3 &&
      policy$private_storage_cap_bytes == 5 * 1024^3 &&
      policy$total_worker_cap_seconds == 12 * 3600,
    policy$scientific_implementation_root_sha256 ==
      rebind$production_implementation_root_sha256,
    policy$execution_implementation_root_sha256 ==
      mv06g_completion_root_v1(sources),
    policy$maximum_workers == 1L && !as.logical(policy$automatic_retry) &&
      as.logical(policy$complete_validation_required) &&
      as.logical(policy$immutable_resume_required),
    policy$outcome_label_state == "closed" &&
      !as.logical(policy$biological_outcomes_computed) &&
      policy$fusion_evaluations == 0L && policy$outcome_jobs == 0L,
    length(private_files) == 0L
  ),
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(checks, args[[8L]], row.names = FALSE, na = "")
if (any(!checks$passed)) stop("MV6-G completion prefreeze failed.",
                              call. = FALSE)
message("Validated MV6-G completion prefreeze: 10/10 pass.")
