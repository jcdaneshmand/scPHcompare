#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: validate_mv06g_fusion_prefreeze.R QUEUE COMPLETE_INVENTORY ",
       "COMPLETE_RESUME GROUP_ROOT RUST_LIBRARY EVIDENCE_DIR OUTPUT",
       call. = FALSE)
}
source("R/mv06f_production.R")
read_csv <- function(name) utils::read.csv(
  file.path(args[[6L]], name), stringsAsFactors = FALSE, check.names = FALSE
)
queue <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE,
                         check.names = FALSE)
inventory <- utils::read.csv(args[[2L]], stringsAsFactors = FALSE,
                             check.names = FALSE)
resume <- utils::read.csv(args[[3L]], stringsAsFactors = FALSE,
                          check.names = FALSE)
contract <- read_csv("mv06g-contract.csv")
sources <- read_csv("mv06g-source-groups.csv")
workload <- read_csv("mv06g-training-workload.csv")
methods <- read_csv("mv06g-method-panel.csv")
endpoints <- read_csv("mv06g-endpoint-plan.csv")
contrasts <- read_csv("mv06g-contrast-plan.csv")
inference <- read_csv("mv06g-inference-plan.csv")
firewall <- read_csv("mv06g-label-firewall.csv")
resources <- read_csv("mv06g-resource-plan.csv")
implementation <- read_csv("mv06g-implementation-sources.csv")
mv06f_validate_queue_v1(queue)
queue <- queue[order(queue$execution_order), , drop = FALSE]
rownames(queue) <- NULL
training_pairs <- queue$training_samples * (queue$training_samples - 1) / 2
expected_methods <- c(
  "cell_H0", "cell_H1", "gene_H0", "gene_H1", "cell_composite",
  "fusion_gene_weight_025", "fusion_gene_weight_050",
  "fusion_gene_weight_075", "gene_composite"
)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
group_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
source_ok <- nrow(sources) == 75L &&
  identical(sources$group_id, queue$group_id) &&
  all(vapply(seq_len(nrow(queue)), function(index) {
    directory <- file.path(group_root, safe_name(queue$group_id[[index]]))
    diagrams <- file.path(directory, "diagrams.rds")
    distances <- file.path(directory, "distances.csv")
    file.exists(diagrams) && file.exists(distances) &&
      .mv06f_sha256(diagrams) == sources$diagrams_sha256[[index]] &&
      .mv06f_sha256(distances) == sources$distances_sha256[[index]] &&
      sources$diagrams_sha256[[index]] == inventory$diagrams_sha256[
        inventory$group_id == queue$group_id[[index]]
      ] && sources$distances_sha256[[index]] == inventory$distances_sha256[
        inventory$group_id == queue$group_id[[index]]
      ]
  }, logical(1L)))
implementation_root <- digest::digest(
  stats::setNames(implementation$sha256, implementation$path),
  algo = "sha256", serialize = TRUE
)
implementation_ok <- nrow(implementation) == 5L &&
  !anyDuplicated(implementation$path) &&
  all(file.exists(implementation$path)) &&
  identical(
    unname(vapply(implementation$path, .mv06f_sha256, character(1L))),
    unname(implementation$sha256)
  ) && contract$implementation_root_sha256 == implementation_root
checks <- data.frame(
  category = c(
    "mv06f_completion_identity", "source_group_identity",
    "training_workload", "method_panel", "endpoint_and_contrasts",
    "blocked_inference", "label_firewall", "resource_staging",
    "implementation_identity", "contract_summary",
    "label_and_downstream_firewall", "stage1_only_admission"
  ),
  passed = c(
    nrow(inventory) == 75L && nrow(resume) == 376L &&
      all(as.logical(resume$sha256_unchanged)) &&
      contract$complete_inventory_sha256 == .mv06f_sha256(args[[2L]]) &&
      contract$complete_resume_sha256 == .mv06f_sha256(args[[3L]]),
    source_ok,
    nrow(workload) == 75L &&
      identical(as.integer(workload$training_biological_pairs),
                as.integer(training_pairs)) &&
      sum(workload$training_biological_pairs) == 262675L &&
      sum(workload$training_component_rows) == 1050700L &&
      sum(workload$query_biological_pairs) == 35350L &&
      sum(workload$query_ranking_rows) == 318150L,
    nrow(methods) == 9L && identical(methods$method_id, expected_methods) &&
      methods$method_id[methods$method_role == "fusion_primary"] ==
        "fusion_gene_weight_050" &&
      !any(as.logical(methods$selected_from_outcomes)),
    nrow(endpoints) == 2L && nrow(contrasts) == 2L &&
      all(contrasts$endpoint_id == "cross_study_tissue_mrr_v1") &&
      setequal(contrasts$comparator_id,
               c("cell_composite", "gene_composite")) &&
      all(as.logical(contrasts$fusion_benefit_requires_both)),
    nrow(inference) == 1L && inference$bootstrap_replicates == 2000L &&
      inference$randomization_replicates == 9999L &&
      grepl("Holm", inference$primary_adjustment),
    nrow(firewall) == 10L && all(as.logical(firewall$required)) &&
      sum(firewall$current_state == "passed") == 2L,
    nrow(resources) == 2L && resources$groups[[1L]] == 1L &&
      resources$groups[[2L]] == 74L &&
      all(resources$maximum_workers == 1L) &&
      all(resources$peak_process_tree_rss_cap_bytes == 12 * 1024^3) &&
      !as.logical(resources$authorized[[2L]]),
    implementation_ok,
    nrow(contract) == 1L && contract$groups == 75L &&
      contract$queue_root_sha256 == mv06f_queue_root_v1(queue) &&
      contract$training_biological_pairs == 262675L &&
      contract$training_component_rows == 1050700L &&
      contract$query_ranking_rows == 318150L &&
      contract$rust_library_sha256 == .mv06f_sha256(args[[5L]]) &&
      contract$metadata_sha256 ==
        "e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0",
    all(sources$outcome_label_state == "closed") &&
      all(workload$outcome_label_state == "closed") &&
      !any(as.logical(sources$biological_outcomes_computed)) &&
      !any(as.logical(workload$biological_outcomes_computed)) &&
      !as.logical(contract$outcome_labels_opened) &&
      !as.logical(contract$biological_outcomes_computed) &&
      contract$fusion_evaluations == 0L && contract$outcome_jobs == 0L,
    as.logical(contract$stage1_authorized) &&
      !as.logical(contract$stage2_authorized) &&
      contract$disposition ==
        "prefreeze_pass_stage1_training_scale_sentinel_only"
  ),
  stringsAsFactors = FALSE
)
checks$contract_id <- "mv06g_prefreeze_validation_v1"
checks$outcome_label_state <- "closed"
checks$biological_outcomes_computed <- FALSE
utils::write.csv(checks, args[[7L]], row.names = FALSE, na = "")
if (any(!checks$passed)) {
  stop("MV6-G fusion prefreeze validation failed: ",
       paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
}
message("Validated MV6-G fusion prefreeze: 12/12 categories pass.")
