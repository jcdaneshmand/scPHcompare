#!/usr/bin/env Rscript

options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 7L) {
  stop("usage: validate_mv06f_stage2_prefreeze.R QUEUE CONTRACT SOURCES ",
       "RESOURCE_PLAN OLD_CONTRACT RUST_LIBRARY OUTPUT_CSV", call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06f_stage2_execution.R")
read_csv <- function(index) utils::read.csv(
  args[[index]], stringsAsFactors = FALSE, check.names = FALSE
)
queue <- read_csv(1L); contract <- read_csv(2L); sources <- read_csv(3L)
plan <- read_csv(4L); old <- read_csv(5L)
checks <- list()
add <- function(category, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv06f_stage2_prefreeze_validation_v1",
    category = category, passed = isTRUE(passed), detail = detail,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
contract_ok <- tryCatch({
  mv06f_validate_stage2_rebind_contract_v1(
    queue, contract, sources, plan, args[[6L]]
  ); TRUE
}, error = function(error) FALSE)
add("rebind_contract", contract_ok,
    "queue/new implementation/Rust/resource roots validate")
add("parent_root",
    nrow(old) == 1L && contract$parent_implementation_root_sha256 ==
      old$implementation_root_sha256 &&
      old$implementation_root_sha256 ==
      "599074b3cd078cf27eb4a85148eb1df2ce3f84a5bdfd3160617b80a78f78c05e",
    "the accepted stage-one implementation root is preserved as parent")
stage2 <- queue[queue$stage == "stage_2", , drop = FALSE]
add("stage2_axis", nrow(stage2) == 74L &&
      identical(as.integer(stage2$execution_order), 2:75),
    "exactly the frozen remaining 74 groups retain execution order 2-75")
add("stage1_rebind_gate", as.logical(contract$stage1_reexecution_required) &&
      !as.logical(contract$stage2_authorized),
    "stage two remains closed until remediated stage one passes")
add("serial_admission", plan$value[plan$guard == "maximum_workers"] == 1,
    "one worker avoids unvalidated concurrency and retains live caps")
firewall <- all(sources$outcome_label_state == "closed") &&
  !any(as.logical(sources$biological_outcomes_computed)) &&
  contract$outcome_label_state == "closed" &&
  !as.logical(contract$biological_outcomes_computed) &&
  !any(tolower(names(contract)) %in% .mv06f_forbidden_fields)
add("label_firewall", firewall,
    "rebind and orchestration artifacts remain technical and label closed")
validation <- do.call(rbind, checks)
if (nrow(validation) != 6L || any(!validation$passed)) {
  print(validation)
  stop("MV6-F stage-two prefreeze validation failed.", call. = FALSE)
}
dir.create(dirname(args[[7L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(validation, args[[7L]], row.names = FALSE, na = "")
message("Validated MV6-F stage-two prefreeze: 6/6 categories pass.")
