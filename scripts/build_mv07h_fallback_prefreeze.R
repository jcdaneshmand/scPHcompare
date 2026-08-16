#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv07h_fallback_prefreeze.R BASE_PREFREEZE GATE_EVIDENCE PRIVATE_ROOT OUTPUT EXPECTED_HEAD")
}
base_dir <- args[[1L]]
gate_dir <- args[[2L]]
private_root <- args[[3L]]
output <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head) ||
    !grepl("^[0-9a-f]{40}$", expected_head)) {
  stop("MV7-H fallback prefreeze exact HEAD mismatch.")
}
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-H fallback prefreeze output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/mv07h_full_topology.R")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")

base_files <- sort(list.files(base_dir, pattern = "[.]csv$", full.names = TRUE),
                   method = "radix")
gate_files <- sort(list.files(gate_dir, pattern = "[.]csv$", full.names = TRUE),
                   method = "radix")
if (length(base_files) != 19L || length(gate_files) != 4L) {
  stop("MV7-H fallback source evidence file set differs.")
}
base_acceptance <- read.csv(file.path(base_dir, "mv07h-acceptance.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
base_independent <- read.csv(file.path(base_dir,
  "mv07h-independent-validation.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
base_decision <- read.csv(file.path(base_dir, "mv07h-decision.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
gate_decision <- read.csv(file.path(gate_dir,
  "mv07h-engine-policy-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
gate_resource <- read.csv(file.path(gate_dir, "mv07h-resource-gate.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
gate_gudhi <- read.csv(file.path(gate_dir, "mv07h-gudhi-feasibility.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
gate_equivalence <- read.csv(file.path(gate_dir,
  "mv07h-full-engine-equivalence.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
if (!all(truth(base_acceptance$passed)) ||
    !all(truth(base_independent$passed)) ||
    base_decision$decision !=
      "authorize_full_source_PH_then_one_landscape_stress_group" ||
    sum(truth(gate_decision$recommended)) != 1L ||
    gate_decision$option[truth(gate_decision$recommended)] !=
      "ripserr_primary_gudhi_exact_resource_fallback" ||
    !all(truth(gate_gudhi$passed)) || !all(truth(gate_equivalence$passed))) {
  stop("MV7-H fallback source evidence does not authorize amendment freeze.")
}

metrics_path <- file.path(private_root, "metrics.csv")
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(metrics) != 21L ||
    sum(metrics$stage == "source_views" & metrics$disposition == "completed") != 5L ||
    sum(metrics$stage == "cell_ph" & metrics$disposition == "completed") != 8L ||
    sum(metrics$stage == "gene_ph" & metrics$disposition == "completed") != 7L ||
    sum(metrics$stage == "gene_ph" &
          metrics$disposition == "rss_cap_exceeded") != 1L) {
  stop("MV7-H partial production ledger differs from the accepted resource gate.")
}
completed <- metrics[metrics$disposition == "completed",, drop = FALSE]
artifact_paths <- file.path(private_root, completed$output_file)
if (!all(file.exists(artifact_paths)) ||
    !identical(unname(vapply(artifact_paths, sha, character(1L))),
               unname(completed$output_sha256)) ||
    !identical(as.numeric(file.info(artifact_paths)$size),
               as.numeric(completed$output_bytes))) {
  stop("MV7-H partial production artifacts differ from their receipts.")
}
failed <- metrics[metrics$disposition != "completed",, drop = FALSE]
if (nrow(failed) != 1L ||
    file.exists(file.path(private_root, failed$output_file)) ||
    failed$job_id != "ph__20260805__SRA628554_SRS2664364__gene_topology_v1") {
  stop("MV7-H failed output boundary differs from the accepted gate.")
}

policy <- mv07h_ph_fallback_policy_v1()
authorization <- data.frame(
  contract_id = "mv07h_fallback_owner_authorization_v1",
  authorization_date = "2026-08-16",
  authorized_policy =
    "ripserr_primary_gudhi_exact_resource_fallback_12GiB",
  scientific_definition_change = FALSE,
  landscape_definition_change = FALSE,
  label_access_authorized = FALSE,
  outcome_access_authorized = FALSE,
  authorized_by = "project_owner",
  stringsAsFactors = FALSE
)
partial_state <- metrics[c(
  "job_id", "stage", "seed", "sample_id", "view_id", "disposition",
  "exit_status", "elapsed_seconds", "peak_process_tree_rss_bytes",
  "output_file", "output_bytes", "output_sha256", "elapsed_cap_seconds",
  "rss_cap_bytes", "outcome_label_state", "biological_outcomes_computed"
)]
partial_state$contract_id <- "mv07h_fallback_partial_state_v1"
partial_state <- partial_state[c("contract_id", setdiff(names(partial_state),
                                                         "contract_id"))]
partial_artifacts <- data.frame(
  contract_id = "mv07h_fallback_partial_artifact_v1",
  job_id = completed$job_id, output_file = completed$output_file,
  bytes = completed$output_bytes, sha256 = completed$output_sha256,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

implementation_paths <- c(
  "R/mv07h_full_topology.R", "scripts/run_mv07h_ph_fallback_entry.R",
  "scripts/run_mv07h_full_ph.R", "scripts/validate_mv07h_full_ph.R",
  "scripts/build_mv07h_fallback_prefreeze.R",
  "scripts/validate_mv07h_fallback_prefreeze.R",
  "scripts/validate_mv07h_prefreeze_repeat.R",
  "tests/testthat/test-mv07h-full-topology.R",
  "docs/specifications/MV07H_EXACT_PH_RESOURCE_FALLBACK_PREFREEZE_V1.md"
)
if (any(!file.exists(implementation_paths))) {
  stop("MV7-H fallback implementation freeze is incomplete.")
}
implementation_hashes <- vapply(implementation_paths, sha, character(1L))
implementation_root <- digest::digest(data.frame(
  path = implementation_paths, sha256 = unname(implementation_hashes),
  accepted_head = expected_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)
base_root <- digest::digest(data.frame(
  file = basename(base_files), sha256 = vapply(base_files, sha, character(1L)),
  stringsAsFactors = FALSE), algo = "sha256", serialize = TRUE)
gate_root <- digest::digest(data.frame(
  file = basename(gate_files), sha256 = vapply(gate_files, sha, character(1L)),
  stringsAsFactors = FALSE), algo = "sha256", serialize = TRUE)
partial_root <- digest::digest(partial_artifacts, algo = "sha256",
                               serialize = TRUE)
contract <- data.frame(
  contract_id = "mv07h_exact_ph_fallback_prefreeze_v1",
  accepted_head = expected_head,
  implementation_root_sha256 = implementation_root,
  base_prefreeze_root_sha256 = base_root,
  resource_gate_root_sha256 = gate_root,
  primary_ledger_sha256 = sha(metrics_path),
  partial_artifact_root_sha256 = partial_root,
  partial_source_jobs = 5L, partial_ph_jobs = 15L,
  partial_failed_primary_jobs = 1L,
  resume_source_jobs = 0L, resume_ph_jobs = 1225L,
  primary_engine = policy$primary_engine,
  fallback_engine = policy$fallback_engine,
  fallback_rss_cap_bytes = policy$fallback_rss_cap_bytes,
  fallback_elapsed_cap_seconds = policy$fallback_elapsed_cap_seconds,
  landscape_definition =
    "finite_positive;essential_H0_excluded;all_active_levels;exact_squared_L2;no_grid;no_level_cap;streamed_groups",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
freeze_paths <- c(base_files, gate_files, implementation_paths)
freeze <- data.frame(
  contract_id = "mv07h_fallback_source_freeze_v1",
  source_id = c(paste0("base_", basename(base_files)),
                paste0("gate_", basename(gate_files)),
                paste0("implementation_", basename(implementation_paths))),
  artifact_locator = gsub("\\\\", "/", freeze_paths),
  sha256 = vapply(freeze_paths, sha, character(1L)),
  bytes = as.numeric(file.info(freeze_paths)$size),
  accepted_head = expected_head, private_source = FALSE,
  stringsAsFactors = FALSE
)
firewall <- data.frame(
  contract_id = "mv07h_fallback_label_firewall_v1",
  allowed_trigger = "primary_ripserr_rss_cap_exceeded",
  forbidden_trigger = "diagram_content;landscape;cluster;label;outcome;claim",
  labels_may_select_engine = FALSE, labels_may_stop = FALSE,
  landscapes_open = FALSE, clustering_open = FALSE, labels_open = FALSE,
  outcomes_open = FALSE, claims_open = FALSE,
  outcome_label_state = "closed", stringsAsFactors = FALSE
)
acceptance <- data.frame(
  contract_id = "mv07h_fallback_prefreeze_acceptance_v1",
  category = c("owner_authorization", "base_prefreeze", "resource_gate",
               "full_engine_equivalence", "fallback_policy",
               "implementation_freeze", "partial_ledger",
               "partial_artifacts", "failed_output_atomicity",
               "landscape_invariance", "label_firewall"),
  passed = TRUE,
  detail = c("project owner approved exact fallback and 12 GiB cap",
             "19 exact v4 prefreeze artifacts", "one Ripserr RSS failure",
             "full 500-gene H0 H1 maximum error zero",
             "gene only after RSS cap; one GUDHI attempt",
             "nine exact implementation sources", "21 exact attempt receipts",
             "20 of 20 completed artifacts match", "failed output absent",
             "all active exact separate H0 H1 unchanged",
             "labels outcomes clustering claims closed"),
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv07h_fallback_prefreeze_decision_v1",
  decision = "authorize_resume_full_PH_with_exact_resource_fallback",
  existing_source_artifacts_reused = 5L,
  existing_ph_artifacts_reused = 15L,
  ph_jobs_remaining = 1225L,
  fallback_attempts_per_eligible_failure = 1L,
  fallback_repeat_required = TRUE,
  landscape_jobs_authorized = 0L, clustering_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  next_gate = "complete_full_PH_repeat_independent_validation",
  stringsAsFactors = FALSE
)

for (path in c(base_files, gate_files)) {
  if (!file.copy(path, file.path(output, basename(path)), overwrite = FALSE)) {
    stop("Failed to copy exact MV7-H fallback source evidence.")
  }
}
outputs <- list(
  "mv07h-ph-fallback-policy.csv" = policy,
  "mv07h-fallback-authorization.csv" = authorization,
  "mv07h-fallback-partial-state.csv" = partial_state,
  "mv07h-fallback-partial-artifacts.csv" = partial_artifacts,
  "mv07h-fallback-contract.csv" = contract,
  "mv07h-fallback-source-freeze.csv" = freeze,
  "mv07h-fallback-label-firewall.csv" = firewall,
  "mv07h-fallback-acceptance.csv" = acceptance,
  "mv07h-fallback-decision.csv" = decision
)
for (name in names(outputs)) write.csv(outputs[[name]], file.path(output, name),
                                       row.names = FALSE, na = "")
message("MV7-H exact PH fallback prefreeze complete; labels closed")
