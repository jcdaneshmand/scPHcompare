#!/usr/bin/env Rscript

options(warn = 2)
args <- getOption("mv06f.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 10L) {
  stop("usage: validate_mv06f_production_prefreeze.R CANDIDATE FOLDS ",
       "RESOURCES PANEL RUST_LIBRARY QUEUE SOURCES RESOURCE_PLAN ABORTS ",
       "CONTRACT", call. = FALSE)
}
source("R/mv06f_production.R")
read_csv <- function(index) utils::read.csv(
  args[[index]], stringsAsFactors = FALSE, check.names = FALSE
)
candidate <- read_csv(1L); folds <- read_csv(2L); resources <- read_csv(3L)
panel <- read_csv(4L); queue <- read_csv(6L); sources <- read_csv(7L)
resource_plan <- read_csv(8L); aborts <- read_csv(9L); contract <- read_csv(10L)

expected_hashes <- c(
  candidate = "842c047ba821f8eca317da52504910733509fb4fddd11d6f54f7e79d9f29d0b7",
  folds = "50379f98cd4927c5c8cb19dbd9ca8ecc7b7b3a9af2e04eb9c8358ecb0b722c6d",
  resources = "73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308",
  panel = "b3a5aff1a0bc01e871751fb9db0b3babfaf18835e68c5699346d8476d903d0ab",
  rust = "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d"
)
observed_hashes <- vapply(args[1:5], .mv06f_sha256, character(1L))
names(observed_hashes) <- names(expected_hashes)
checks <- list()
add <- function(category, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv06f_prefreeze_validation_v1", category = category,
    passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE
  )
}
add("frozen_input_hashes", identical(tolower(observed_hashes), expected_hashes),
    "candidate/fold/resource/panel/Rust hashes match the prospective freeze")
input_ok <- tryCatch({
  mv06f_validate_prefreeze_inputs_v1(candidate, folds, resources, panel); TRUE
}, error = function(error) FALSE)
add("input_contracts", input_ok,
    "90 samples, 75 folds, 450 caches, and 500 ordered genes validate")

# Reconstruct the queue independently from fold rows rather than trusting the
# builder's group IDs or execution ordering.
folds <- folds[order(folds$fold_id, as.integer(folds$seed), method = "radix"),
               , drop = FALSE]
independent <- data.frame(
  fold_id = folds$fold_id, seed = as.integer(folds$seed),
  training_samples = as.integer(folds$training_samples),
  held_out_samples = as.integer(folds$held_out_samples), stringsAsFactors = FALSE
)
independent$biological_pairs <- independent$training_samples *
  independent$held_out_samples
independent$landscape_component_rows <- 4L * independent$biological_pairs
queue_map <- match(paste(independent$fold_id, independent$seed, sep = "\r"),
                   paste(queue$fold_id, queue$seed, sep = "\r"))
queue_ok <- length(queue_map) == 75L && !anyNA(queue_map) &&
  all(independent$training_samples == queue$training_samples[queue_map]) &&
  all(independent$held_out_samples == queue$held_out_samples[queue_map]) &&
  all(independent$biological_pairs == queue$biological_pairs[queue_map]) &&
  all(independent$landscape_component_rows ==
        queue$landscape_component_rows[queue_map])
queue_contract_ok <- tryCatch({mv06f_validate_queue_v1(queue); TRUE},
                              error = function(error) FALSE)
add("independent_group_axis", queue_ok && queue_contract_ok,
    "all 75 fold-seed groups and per-group counts reconstruct independently")
add("complete_workload",
    sum(queue$cell_ph_jobs + queue$gene_ph_jobs) == 13500L &&
      sum(queue$diagram_dimension_records) == 27000L &&
      sum(queue$biological_pairs) == 35350L &&
      sum(queue$landscape_component_rows) == 141400L,
    "13,500 PH jobs, 27,000 dimensions, 35,350 pairs, and 141,400 rows")

maximum <- max(independent$biological_pairs)
expected_fold <- sort(unique(independent$fold_id[
  independent$biological_pairs == maximum
]), method = "radix")[[1L]]
stage <- queue[queue$stage == "stage_1_maximum", , drop = FALSE]
add("stage1_selection", nrow(stage) == 1L &&
      stage$fold_id == expected_fold && stage$seed == 20260807L &&
      stage$execution_order == 1L,
    "stage 1 is the label-free maximum-pair group with the frozen tie break")

source_ok <- nrow(sources) == 14L && !anyDuplicated(sources$path) &&
  all(file.exists(sources$path)) &&
  identical(tolower(unname(vapply(
    sources$path, .mv06f_sha256, character(1L)
  ))), tolower(unname(sources$sha256))) &&
  all(sources$outcome_label_state == "closed") &&
  !any(as.logical(sources$biological_outcomes_computed))
source_root <- if (source_ok) digest::digest(
  stats::setNames(sources$sha256, sources$path),
  algo = "sha256", serialize = TRUE
) else NA_character_
add("source_inventory", source_ok,
    "all accepted and new prefreeze sources match their frozen hashes")
add("resource_guards", nrow(resource_plan) == 7L &&
      setequal(resource_plan$guard, c(
        "group_elapsed_seconds", "group_process_tree_rss_bytes",
        "stage1_elapsed_seconds", "stage1_process_tree_rss_bytes",
        "concurrent_process_tree_rss_bytes", "production_worker_seconds",
        "private_root_bytes"
      )) && all(resource_plan$value > 0) &&
      all(resource_plan$outcome_label_state == "closed"),
    "all seven prospective time/RSS/work/storage guards are positive")
add("abort_rules", nrow(aborts) == 10L &&
      identical(as.integer(aborts$rule_order), 1:10) &&
      !any(as.logical(aborts$automatic_retry)) &&
      all(aborts$action == "stop_new_launches_preserve_validated_groups"),
    "ten fail-closed rules prohibit automatic retry")

queue_root <- if (queue_contract_ok) mv06f_queue_root_v1(queue) else NA_character_
contract_ok <- nrow(contract) == 1L &&
  contract$contract_id == "mv06f_production_prefreeze_v1" &&
  identical(contract$queue_root_sha256, queue_root) &&
  identical(contract$implementation_root_sha256, source_root) &&
  identical(contract$rust_library_sha256, expected_hashes[["rust"]]) &&
  !as.logical(contract$production_executed) &&
  contract$outcome_label_state == "closed" &&
  !as.logical(contract$biological_outcomes_computed) &&
  all(unlist(contract[c("fusion_jobs", "clustering_jobs", "outcome_jobs")]) == 0)
add("contract_roots_and_stop", contract_ok,
    "queue/source/Rust roots bind the prefreeze and production remains unexecuted")

public_frames <- list(queue, sources, resource_plan, aborts, contract)
firewall <- all(vapply(public_frames, function(value) {
  !any(tolower(names(value)) %in% .mv06f_forbidden_fields) &&
    all(value$outcome_label_state == "closed") &&
    !any(as.logical(value$biological_outcomes_computed))
}, logical(1L)))
add("public_label_firewall", firewall,
    "public prefreeze artifacts contain technical identities only")

validation <- do.call(rbind, checks)
if (nrow(validation) != 10L || any(!validation$passed)) {
  print(validation)
  stop("MV6-F prefreeze validation failed.", call. = FALSE)
}
output <- file.path(dirname(args[[10L]]), "mv06f-validation.csv")
utils::write.csv(validation, output, row.names = FALSE, na = "")
message("Validated MV6-F prefreeze: 10/10 categories pass.")
