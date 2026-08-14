#!/usr/bin/env Rscript

options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 9L) {
  stop("usage: validate_mv06f_stage2_complete.R QUEUE CONTRACT SOURCES ",
       "RESOURCE_PLAN ADMISSION RUST_LIBRARY PRIVATE_ROOT METRICS OUTPUT_DIR",
       call. = FALSE)
}
source("R/mv06f_production.R")
source("R/mv06f_stage2_execution.R")
read_csv <- function(index) utils::read.csv(
  args[[index]], stringsAsFactors = FALSE, check.names = FALSE
)
queue <- read_csv(1L); contract <- read_csv(2L); sources <- read_csv(3L)
plan <- read_csv(4L); admission <- read_csv(5L); metrics <- read_csv(8L)
mv06f_validate_stage2_rebind_contract_v1(
  queue, contract, sources, plan, args[[6L]]
)
mv06f_validate_stage2_admission_v1(admission, contract)
mv06f_validate_stage2_metrics_v1(metrics, queue, contract, complete = TRUE)
group_root <- file.path(args[[7L]], "groups")
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
rows <- vector("list", nrow(queue))
all_pair_ids <- character()
component_counts <- integer(4L)
names(component_counts) <- c("cell_H0", "cell_H1", "gene_H0", "gene_H1")
for (index in seq_len(nrow(queue))) {
  unit <- queue[index, , drop = FALSE]
  path <- file.path(group_root, safe_name(unit$group_id))
  status <- mv06f_validate_group_directory_v1(
    path, unit, contract$queue_root_sha256,
    contract$implementation_root_sha256, contract$rust_library_sha256
  )
  distances <- utils::read.csv(file.path(path, "distances.csv"),
                               stringsAsFactors = FALSE, check.names = FALSE)
  all_pair_ids <- c(all_pair_ids, distances$pair_id)
  keys <- paste(
    ifelse(distances$view_id == "cell_topology_v1", "cell", "gene"),
    distances$homology_dimension, sep = "_"
  )
  component_counts <- component_counts + table(factor(
    keys, levels = names(component_counts)
  ))
  files <- file.path(path, c(
    "diagrams.rds", "diagram-manifest.csv", "distances.csv", "metrics.csv",
    "status.csv"
  ))
  rows[[index]] <- data.frame(
    contract_id = "mv06f_complete_group_inventory_v1",
    group_id = unit$group_id, fold_id = unit$fold_id, seed = unit$seed,
    execution_order = unit$execution_order,
    biological_pairs = unit$biological_pairs,
    landscape_component_rows = unit$landscape_component_rows,
    diagrams_sha256 = status$diagrams_sha256,
    manifest_sha256 = status$diagram_manifest_sha256,
    distances_sha256 = status$distances_sha256,
    metrics_sha256 = status$metrics_sha256,
    status_sha256 = .mv06f_sha256(files[[5L]]),
    group_bytes = sum(file.info(files)$size),
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
inventory <- do.call(rbind, rows)
checks <- data.frame(
  contract_id = "mv06f_complete_validation_v1",
  category = c(
    "all_group_directories", "complete_workload", "unique_pair_axis",
    "balanced_component_totals", "stage2_resource_metrics",
    "label_and_downstream_firewall"
  ),
  passed = c(
    nrow(inventory) == 75L,
    sum(inventory$biological_pairs) == 35350L &&
      sum(inventory$landscape_component_rows) == 141400L,
    length(all_pair_ids) == 141400L && !anyDuplicated(all_pair_ids),
    identical(as.integer(component_counts), rep(35350L, 4L)),
    nrow(metrics) == 74L &&
      all(metrics$disposition %in% .mv06f_stage2_success),
    all(inventory$outcome_label_state == "closed") &&
      !any(as.logical(inventory$biological_outcomes_computed)) &&
      all(metrics$fusion_jobs == 0L) && all(metrics$clustering_jobs == 0L) &&
      all(metrics$outcome_jobs == 0L)
  ),
  detail = c(
    "75/75 remediated-root group directories validate",
    "35,350 biological pairs and 141,400 component rows reconstruct",
    "all component pair identities are globally unique",
    "35,350 rows exist for each cell/gene by H0/H1 component",
    "74/74 stage-two resource rows complete or validated reuse",
    "labels remain closed and fusion/clustering/outcome counts remain zero"
  ), outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (any(!checks$passed)) {
  print(checks)
  stop("MV6-F complete stage-two validation failed.", call. = FALSE)
}
dir.create(args[[9L]], recursive = TRUE, showWarnings = FALSE)
utils::write.csv(inventory, file.path(
  args[[9L]], "mv06f-complete-group-inventory.csv"
), row.names = FALSE, na = "")
utils::write.csv(checks, file.path(
  args[[9L]], "mv06f-complete-validation.csv"
), row.names = FALSE, na = "")
message("Validated complete MV6-F production: 6/6 categories pass.")
