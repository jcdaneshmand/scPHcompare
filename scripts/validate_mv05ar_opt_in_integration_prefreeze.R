#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: validate_mv05ar_opt_in_integration_prefreeze.R BUILD_A BUILD_B OUTPUT_CSV")
}
builds <- normalizePath(args[1:2], mustWork = TRUE)
output_csv <- normalizePath(args[[3L]], mustWork = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)

read_build <- function(build, id) utils::read.csv(file.path(
  build, paste0("mv05ar-", id, "-2026-08-12.csv")
), stringsAsFactors = FALSE, check.names = FALSE)
checks <- list()
record <- function(id, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.frame(
    contract_id = "mv05ar_independent_validation_v1",
    validation_id = id, passed = isTRUE(passed), evidence = evidence,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(passed)) stop("MV5-AR validation failed: ", id)
}

ledger_ids <- function(build) sub(
  "-[0-9]{4}-[0-9]{2}-[0-9]{2}\\.csv$", "",
  sub("^mv05ar-", "", basename(list.files(
    build, pattern = "^mv05ar-.*csv$"
  )))
)
ids <- ledger_ids(builds[[1L]])
record("artifact_set", length(ids) == 14L &&
         setequal(ids, ledger_ids(builds[[2L]])),
       "same 14 deterministic prefreeze ledgers")
repeat_ok <- vapply(ids, function(id) {
  identical(readLines(file.path(builds[[1L]], paste0("mv05ar-", id,
    "-2026-08-12.csv")), warn = FALSE),
    readLines(file.path(builds[[2L]], paste0("mv05ar-", id,
      "-2026-08-12.csv")), warn = FALSE))
}, logical(1L))
record("byte_repeat", all(repeat_ok), "all 14 ledgers byte-identical")

pathways <- read_build(builds[[1L]], "pathway-inventory")
namespace <- readLines("NAMESPACE", warn = FALSE)
actual_formals <- vapply(pathways$entrypoint, function(name) paste(
  names(formals(get(name, mode = "function"))), collapse = ";"), character(1))
record("pathway_rescan", nrow(pathways) == 6L &&
         identical(pathways$current_formals, unname(actual_formals)) &&
         all(paste0("export(", pathways$entrypoint, ")") %in% namespace) &&
         !any(pathways$corrected_control_currently_present),
       "six exported workflows rescanned; corrected control not implemented")

boundaries <- read_build(builds[[1L]], "boundaries")
record("additive_boundary", nrow(boundaries) == 6L &&
         !any(boundaries$default_behavior_changed) &&
         !any(boundaries$corrected_downstream_consumption) &&
         sum(grepl("add_null_default", boundaries$later_change)) == 2L,
       "only unified pass-through and postprocessing orchestration may change later")

controls <- read_build(builds[[1L]], "controls")
record("control_contract", nrow(controls) == 12L &&
         controls$value[controls$field == "downstream_use"] == "artifacts_only" &&
         controls$value[controls$field == "method"] == "auto" &&
         controls$value[controls$field == "exact_max_intervals"] == "500" &&
         controls$value[controls$field == "abs_tol"] == "1e-8" &&
         controls$value[controls$field == "rel_tol"] == "1e-8",
       "default-off artifacts-only strict auto/500 control")

artifacts <- read_build(builds[[1L]], "artifacts")
record("artifact_contract", nrow(artifacts) == 7L &&
         identical(artifacts$order, 1:7) && !any(artifacts$overwrite_allowed) &&
         !any(artifacts$legacy_filename_reused) &&
         tail(artifacts$artifact, 1L) == "completion-v1.csv",
       "seven non-colliding versioned artifacts; completion written last")

resources <- read_build(builds[[1L]], "resources")
value <- function(id) resources$value[resources$policy == id]
record("resource_policy", value("exact_exact_pair_seconds") == 30 &&
         value("adaptive_pair_seconds") == 240 &&
         value("iteration_overhead_seconds") == 30 &&
         value("minimum_rss_bytes") == 1610612736 &&
         value("workers") == 1 &&
         all(resources$outside_envelope_action[resources$policy %in% c(
           "observed_h0_min", "observed_h0_max", "observed_h1_min",
           "observed_h1_max")] == "profiling_required"),
       "one-worker budget admission with out-of-envelope refusal")

binding <- read_build(builds[[1L]], "accepted-evidence-binding")
record("accepted_evidence", binding$value[binding$evidence ==
         "mv05apr1_max_wall_seconds"] == "567.94" &&
         binding$value[binding$evidence == "mv05apr1_max_rss_bytes"] ==
           "990363648" &&
         binding$value[binding$evidence == "mv05aq_method"] == "auto" &&
         binding$value[binding$evidence == "mv05aq_exact_guard"] == "500" &&
         all(binding$accepted),
       "MV5-AP-R1 resources and MV5-AQ numerical policy bound exactly")

legacy <- read_build(builds[[1L]], "legacy-coexistence")
record("legacy_coexistence", nrow(legacy) == 8L &&
         !any(legacy$silent_conversion_allowed) && !any(legacy$overwrite_allowed),
       "legacy generator, keys, fields, modular and cross-iteration paths preserved")

atomic <- read_build(builds[[1L]], "atomic-resume-contract")
record("atomic_resume", nrow(atomic) == 10L &&
         sum(atomic$completion_visible) == 1L &&
         atomic$action[atomic$completion_visible] ==
           "write_hash_bound_completion_last",
       "pair shards precede cross-checked matrix and last completion marker")

validation <- read_build(builds[[1L]], "implementation-validation-plan")
aborts <- read_build(builds[[1L]], "implementation-abort-rules")
record("validation_and_aborts", nrow(validation) == 15L &&
         all(validation$required) && nrow(aborts) == 14L &&
         !any(aborts$partial_scientific_result_accepted),
       "15 required validation classes and 14 hard abort rules")

stages <- read_build(builds[[1L]], "stages")
record("migration_stages", nrow(stages) == 5L &&
         sum(stages$authorized_now) == 1L &&
         stages$stage_id[stages$authorized_now] == "prefreeze_additive_artifacts",
       "implementation, smoke, consumption, and default stages remain separate")

source_freeze <- read_build(builds[[1L]], "source-freeze")
actual_sha <- vapply(source_freeze$path, function(path) sub(" .*", "", system2(
  "sha256sum", path, stdout = TRUE)), character(1))
record("source_freeze", identical(unname(source_freeze$sha256),
                                   unname(actual_sha)),
       "11 source and namespace hashes reproduce")

legacy_sources <- c("R/unified_pipeline.R", "R/PH_PostProcessing_andAnalysis.R",
                    "R/cross_iteration_functions.R",
                    "R/custom_iteration_inputs_template.R", "NAMESPACE")
base_hash <- vapply(legacy_sources, function(path) system2(
  "git", c("rev-parse", paste0("f6df036:", path)), stdout = TRUE), character(1))
current_hash <- vapply(legacy_sources, function(path) system2(
  "git", c("hash-object", path), stdout = TRUE), character(1))
record("workflow_immutability", identical(unname(base_hash), unname(current_hash)),
       "four workflow sources and NAMESPACE equal MV5-AP-R1 completion")

prohibited <- read_build(builds[[1L]], "prohibited-change-counters")
record("prohibited_changes", nrow(prohibited) == 14L && all(prohibited$value == 0L),
       "all behavior/export/artifact/calculation/shortcut counters zero")

decision <- read_build(builds[[1L]], "continuation-decision")
record("bounded_decision", nrow(decision) == 1L &&
         decision$additive_implementation_authorized &&
         !decision$corrected_downstream_consumption_authorized &&
         !decision$workflow_default_change_authorized &&
         !decision$legacy_artifact_rewrite_authorized &&
         !decision$optimization_authorized && decision$next_sprint == "MV5-AS",
       "only additive opt-in artifact implementation authorized")

result <- do.call(rbind, checks)
utils::write.csv(result, output_csv, row.names = FALSE, na = "", quote = TRUE)
cat("MV5-AR independent validation passed:", nrow(result), "categories\n")
