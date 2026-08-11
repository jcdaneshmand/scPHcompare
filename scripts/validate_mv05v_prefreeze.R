#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv05v_prefreeze.R AUDIT_DIR EXPECTED_HEAD",
       call. = FALSE)
}
source("R/provenance_utils.R")
source("R/mv05v_robustness_prefreeze.R")
audit_dir <- args[[1L]]
expected_head <- args[[2L]]
read_csv <- function(name) utils::read.csv(
  file.path(audit_dir, name), stringsAsFactors = FALSE, check.names = FALSE
)
queue <- read_csv("mv05v-full-group-queue-2026-08-10.csv")
scope <- read_csv("mv05v-base-pair-scope-2026-08-10.csv")
projection <- read_csv("mv05v-resource-projection-2026-08-10.csv")
decision <- read_csv("mv05v-prefreeze-decision-2026-08-10.csv")
summary <- read_csv("mv05v-prefreeze-summary-2026-08-10.csv")
sources <- read_csv("mv05v-source-freeze-2026-08-10.csv")
mv05v_validate_group_queue_v1(queue)
if (nrow(scope) != 150L || sum(scope$biological_pairs) != 70700L ||
    sum(scope$landscape_request_rows) != 141400L ||
    sum(scope$landscape_subchunks) != 720L ||
    nrow(sources[sources$source_class == "private_coordinate", ]) != 150L ||
    nrow(projection) != 3L ||
    !isTRUE(as.logical(decision$gate_complete)) ||
    isTRUE(as.logical(decision$full_execution_authorized)) ||
    isTRUE(as.logical(summary$orchestration_engine_bound)) ||
    isTRUE(as.logical(summary$execution_authorized)) ||
    any(as.logical(queue$outcomes_computed)) ||
    any(queue$outcome_label_state != "closed") ||
    length(unique(queue$prospective_head)) != 1L ||
    unique(queue$prospective_head) != expected_head) {
  stop("MV5-V independent prefreeze validation failed.", call. = FALSE)
}
message(
  "MV5-V prefreeze validation passed: groups=600 views=54000 ",
  "landscape_rows=565600 outcomes=0 execution_auth=0"
)
