#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: validate_mv05an_landscape_public_contract_prefreeze.R AUDIT_DIR OUTPUT_CSV")
audit_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE); output <- args[[2L]]
readc <- function(name) read.csv(file.path(audit_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
truth <- function(x) tolower(as.character(x)) == "true"
checks <- list(); add <- function(id, passed, evidence) checks[[length(checks)+1L]] <<- data.frame(
  contract_id="mv05an_independent_validation_v1", validation_id=id,
  passed=isTRUE(passed), evidence=evidence, production_helper_called=FALSE,
  behavior_or_export_changed=FALSE, project_calculation_executed=FALSE,
  stringsAsFactors=FALSE)

source <- readc("mv05an-source-freeze-2026-08-12.csv")
source_ok <- nrow(source)==20L && all(file.exists(source$artifact_locator)) &&
  all(vapply(source$artifact_locator, sha, character(1L))==source$sha256) &&
  !any(truth(source$behavior_changed))
add("source_hashes", source_ok, "20_of_20_exact")

inventory <- readc("mv05an-function-inventory-2026-08-12.csv")
function_files <- c("R/cross_iteration_functions.R", "R/landscape_contract.R",
  "R/landscape_convergence.R", "R/landscape_reference.R",
  "R/mv05am_four_panel_synthesis_gate.R", "R/PH_PostProcessing_andAnalysis.R")
observed <- list()
for(path in function_files){
  lines <- readLines(path, warn=FALSE)
  hit <- grep("^[A-Za-z0-9_.]*[Ll]andscape[A-Za-z0-9_.]*[[:space:]]*<-[[:space:]]*function|^ComputePersistenceLandscapes[[:space:]]*<-[[:space:]]*function", lines)
  for(i in hit) observed[[length(observed)+1L]] <- data.frame(file=basename(path),
    line=i,function_name=sub("[[:space:]]*<-[[:space:]]*function.*$","",trimws(lines[[i]])))
}
observed <- do.call(rbind, observed)
actual_r <- inventory[inventory$file != "mv05w_exact_landscape_group.py", c("file","line","function_name")]
rescan_ok <- nrow(observed)==45L && nrow(actual_r)==45L &&
  identical(observed$file,actual_r$file) && identical(observed$line,actual_r$line) &&
  identical(observed$function_name,actual_r$function_name) && nrow(inventory)==46L &&
  !any(truth(inventory$ambiguous)) && !any(truth(inventory$direct_export))
add("complete_function_rescan", rescan_ok, "45_R_functions_plus_1_exact_python_engine")

classes <- table(inventory$pathway_class)
class_ok <- all(c("corrected_exact_or_error_controlled_engine",
  "accepted_exact_production_engine", "legacy_grid_artifact_workflow",
  "legacy_curve_approximation", "internal_diagnostic",
  "versioned_fixed_grid_contract") %in% names(classes))
add("semantic_classification", class_ok, paste(names(classes), classes, collapse=";"))

entry <- readc("mv05an-exported-entrypoint-inventory-2026-08-12.csv")
namespace <- readLines("NAMESPACE",warn=FALSE)
entry_ok <- nrow(entry)==6L && all(paste0("export(",entry$entrypoint,")") %in% namespace) &&
  all(truth(entry$exported)) && !any(truth(entry$corrected_exact_default)) &&
  !any(truth(entry$behavior_changed_in_mv05an))
add("exported_entrypoint_rescan",entry_ok,"6_exported_workflows_legacy_exposure_no_change")

artifacts <- readc("mv05an-artifact-schema-2026-08-12.csv")
artifact_ok <- nrow(artifacts)==6L && !any(truth(artifacts$safe_to_overwrite)) &&
  sum(artifacts$semantic_class=="legacy_k1_unit_grid_unversioned")==4L &&
  all(artifacts$migration_action[1:3]=="detect_and_read_only")
add("artifact_schema_boundary",artifact_ok,"4_legacy_2_corrected_no_overwrite")

semantics <- readc("mv05an-existing-semantics-2026-08-12.csv")
semantics_ok <- nrow(semantics)==4L &&
  identical(semantics$specification,c("legacy_k1_unit_grid_v0","paper_k1_common_grid_v1",
    "full_l2_common_grid_v1","full_l2_error_controlled_v1")) &&
  identical(which(truth(semantics$currently_reached_by_exported_workflow)),1L) &&
  identical(which(truth(semantics$accepted_production_semantics)),4L)
add("existing_semantic_registry",semantics_ok,"legacy_exported_path_target_exact_path_distinct")

target <- readc("mv05an-target-public-contract-2026-08-12.csv")
states <- c("all_positive_persistence_finite_intervals","exclude_infinite_h0",
  "all_consecutive_active_levels_zero_pad_missing_depth",
  "compute_and_return_h0_h1_separately",
  "exact_linear_critical_pair_segments_when_within_explicit_guard",
  "no_universal_uniform_grid_or_level_cap")
target_ok <- nrow(target)==12L && all(truth(target$immutable_scientific_semantics)) &&
  all(states %in% target$required_state)
add("target_scientific_contract",target_ok,"12_items_all_levels_exact_or_error_controlled_h0_h1_separate")

api <- readc("mv05an-public-api-decision-2026-08-12.csv")
api_ok <- nrow(api)==2L && all(truth(api$export_in_later_sprint)) &&
  !any(truth(api$legacy_function_redirected)) && !any(truth(api$legacy_artifact_overwritten)) &&
  !any(truth(api$workflow_default_changed))
add("versioned_api_boundary",api_ok,"2_new_APIs_zero_redirect_overwrite_or_default_change")

compat <- readc("mv05an-backward-compatibility-2026-08-12.csv")
migration <- readc("mv05an-migration-plan-2026-08-12.csv")
migration_ok <- nrow(compat)==8L && !any(truth(compat$silent_conversion_allowed)) &&
  nrow(migration)==6L && identical(which(truth(migration$authorized_in_mv05ao)),1:4) &&
  !any(truth(migration$behavior_change[1:4])) && all(truth(migration$behavior_change[5:6]))
add("compatibility_and_staged_migration",migration_ok,"read_only_legacy_first_four_nonchanging_stages_only")

scope <- readc("mv05an-mv05ao-implementation-scope-2026-08-12.csv")
abort <- readc("mv05an-mv05ao-abort-rules-2026-08-12.csv")
guards_ok <- nrow(scope)==12L && all(truth(scope$required)) &&
  !any(truth(scope$default_change_authorized)) && !any(truth(scope$legacy_artifact_rewrite_authorized)) &&
  nrow(abort)==10L && all(abort$required_action=="abort_without_fallback_or_behavior_change")
add("implementation_scope_and_aborts",guards_ok,"12_required_items_10_abort_rules_no_default_or_rewrite")

decision <- readc("mv05an-continuation-decision-2026-08-12.csv")
decision_ok <- nrow(decision)==1L && decision$authorized_next_sprint=="MV5-AO" &&
  truth(decision$new_versioned_api_authorized) && !truth(decision$current_workflow_default_change_authorized) &&
  !truth(decision$existing_function_redirection_authorized) &&
  !truth(decision$legacy_artifact_overwrite_authorized) &&
  !truth(decision$new_scientific_calculation_on_project_data_authorized)
add("continuation_decision",decision_ok,"MV5_AO_new_versioned_API_only")

execution <- readc("mv05an-prohibited-change-counters-2026-08-12.csv")
zero_ok <- nrow(execution)==11L && all(execution$count==0L) && !any(truth(execution$executed))
add("zero_behavior_and_execution_changes",zero_ok,"11_change_categories_zero")

result <- do.call(rbind,checks)
if(nrow(result)!=12L || !all(result$passed)) stop("MV5-AN validation failed: ",paste(result$validation_id[!result$passed],collapse=", "))
if(file.exists(output)) stop("Refusing overwrite validation output")
write.csv(result,output,row.names=FALSE,na="")
message("MV5-AN independent validation passed: 12 categories; 46 pathways")
