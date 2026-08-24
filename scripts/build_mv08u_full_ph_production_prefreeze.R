#!/usr/bin/env Rscript

# Build the metadata-only MV8-U prefreeze for the 1,257 PH records that were
# not already closed by MV8-T. This builder reads only public receipts. It
# never opens a private source cache and never computes PH or landscapes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv08u_full_ph_production_prefreeze.R <output-dir>",
       call. = FALSE)
}
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
if (dir.exists(output_dir)) stop("refusing to overwrite MV8-U output", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)
dir.create(output_dir, recursive = TRUE)

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(value, path) {
  partial <- paste0(path, ".partial")
  writeLines(value, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
key4 <- function(x) paste(
  x$unit_id, x$seed, x$representation_id, x$panel_id, sep = "\r"
)

accepted_parent_head <- tolower(trimws(Sys.getenv(
  "MV08U_PARENT_HEAD", unset = ""
)))
if (!nzchar(accepted_parent_head)) {
  accepted_parent_head <- tolower(trimws(system2(
    "git", c("rev-parse", "HEAD"), stdout = TRUE
  )))
}
if (!grepl("^[0-9a-f]{40}$", accepted_parent_head)) {
  stop("cannot bind MV8-U parent HEAD", call. = FALSE)
}

r_root <- "docs/audits/mv08r-full-topology-production-prefreeze-v1"
s_root <- "docs/audits/mv08s-ph-sentinel-prefreeze-v1"
t_root <- "docs/audits/mv08t-ph-sentinel-closure-v1"
q_root <- "docs/audits/mv08q-full-source-production-closure-v1"
p_root <- "docs/audits/mv08p-full-source-production-prefreeze-v2"
inputs <- c(
  file.path(r_root, "mv08r-contract.csv"),
  file.path(r_root, "mv08r-source-cache-bindings.csv"),
  file.path(r_root, "mv08r-source-gene-view-bindings.csv"),
  file.path(r_root, "mv08r-ph-queue.csv"),
  file.path(r_root, "mv08r-backend-policy.csv"),
  file.path(r_root, "mv08r-landscape-contract.csv"),
  file.path(s_root, "mv08s-ph-sentinel-queue.csv"),
  file.path(t_root, "mv08t-private-artifact-rehash.csv"),
  file.path(t_root, "mv08t-resource-ledger.csv"),
  file.path(t_root, "mv08t-validation.csv"),
  file.path(t_root, "mv08t-decision.csv"),
  file.path(t_root, "mv08t-artifact-manifest.csv"),
  file.path(q_root, "mv08q-source-cache-rehash.csv"),
  file.path(q_root, "mv08q-validation.csv"),
  file.path(p_root, "mv08p-remaining-source-queue.csv"),
  file.path(s_root, "mv08s-external-input-bindings.csv")
)
if (!all(file.exists(inputs))) stop("MV8-U public prerequisite absent", call. = FALSE)

r_contract <- read_csv(inputs[[1L]])
source_bindings <- read_csv(inputs[[2L]])
views <- read_csv(inputs[[3L]])
r_queue <- read_csv(inputs[[4L]])
backend <- read_csv(inputs[[5L]])
landscape <- read_csv(inputs[[6L]])
s_queue <- read_csv(inputs[[7L]])
t_artifacts <- read_csv(inputs[[8L]])
t_resources <- read_csv(inputs[[9L]])
t_validation <- read_csv(inputs[[10L]])
t_decision <- read_csv(inputs[[11L]])
t_manifest <- read_csv(inputs[[12L]])
q_sources <- read_csv(inputs[[13L]])
q_validation <- read_csv(inputs[[14L]])
p_queue <- read_csv(inputs[[15L]])
s_external <- read_csv(inputs[[16L]])

if (nrow(r_contract) != 1L || nrow(source_bindings) != 132L ||
    nrow(views) != 1272L || nrow(r_queue) != 1280L ||
    nrow(s_queue) != 23L || nrow(q_sources) != 129L ||
    nrow(p_queue) != 129L || !all(t_validation$passed) ||
    !all(q_validation$passed) || t_decision$ph_records != 23L ||
    t_decision$full_ph_jobs_authorized != 0L || nrow(s_external) != 8L ||
    length(unique(s_external$panel_file_sha256)) != 1L) {
  stop("MV8-U prerequisite closure drift", call. = FALSE)
}
manifest_observed <- vapply(file.path(t_root, t_manifest$artifact), sha_file,
                            character(1L))
if (!identical(unname(manifest_observed), t_manifest$sha256)) {
  stop("MV8-T closure manifest drift", call. = FALSE)
}

closed <- t_artifacts[t_artifacts$artifact_type == "ph", , drop = FALSE]
if (nrow(closed) != 23L || anyDuplicated(closed$artifact_id) ||
    !all(closed$byte_identical & closed$independently_validated)) {
  stop("MV8-T closed PH identity drift", call. = FALSE)
}
remaining <- r_queue[!(r_queue$job_id %in% closed$artifact_id), , drop = FALSE]
if (nrow(remaining) != 1257L || any(remaining$view_kind != "gene_topology_v1") ||
    any(remaining$execution_role != "source_produced_gene_ph") ||
    sum(remaining$dataset_scope == "internal124") != 1236L ||
    sum(remaining$dataset_scope == "external8") != 21L) {
  stop("MV8-U remaining PH queue drift", call. = FALSE)
}

# Attach immutable view geometry and the accepted source-cache receipt.
view_index <- match(key4(remaining), key4(views))
if (anyNA(view_index) || anyDuplicated(key4(views))) {
  stop("MV8-U source-view join drift", call. = FALSE)
}
remaining$matrix_sha256 <- views$matrix_sha256[view_index]
remaining$distance_sha256 <- views$distance_sha256[view_index]
source_index <- match(remaining$unit_id, source_bindings$unit_id)
if (anyNA(source_index)) stop("MV8-U source binding absent", call. = FALSE)
remaining$source_origin <- source_bindings$source_origin[source_index]
remaining$fit_cells <- source_bindings$fit_cells[source_index]
remaining$source_cache_sha256 <- source_bindings$cache_sha256[source_index]

q_index <- match(remaining$unit_id, q_sources$unit_id)
root_role <- rep(NA_character_, nrow(remaining))
relative_file <- rep(NA_character_, nrow(remaining))
production_source <- !is.na(q_index)
q_order <- q_sources$job_order[q_index[production_source]]
p_index <- match(remaining$unit_id[production_source], p_queue$unit_id)
if (anyNA(p_index) || !identical(
      q_sources$cache_sha256[q_index[production_source]],
      remaining$source_cache_sha256[production_source])) {
  stop("MV8-U MV8-Q source locator drift", call. = FALSE)
}
root_role[production_source] <- ifelse(
  q_order <= 123L, "mv08p_original_v1",
  ifelse(q_order == 124L, "mv08pr_overlay_v1", "mv08ps_overlay_v1")
)
relative_file[production_source] <- file.path(
  "cache", p_queue$output_file[p_index]
)
primary_source <- !production_source
if (any(!remaining$unit_id[primary_source] %in%
        c("SRA701877_SRS3279688", "SRA742961_SRS3565197"))) {
  stop("MV8-U unexpected MV8-O primary source in remaining queue", call. = FALSE)
}
root_role[primary_source] <- "mv08o_internal_primary_v8"
relative_file[primary_source] <- file.path(
  "cache", paste0(remaining$unit_id[primary_source], "__primary.rds")
)

# Risk ordering is label-blind and derived only from measured topology burden:
# 475-gene jobs first, then exact500 residuals, then exact500 SCT data.
risk_stage <- ifelse(
  remaining$panel_id == "common475", 1L,
  ifelse(remaining$representation_id ==
           "sct_pearson_residual_all_qc_fit_selected384", 2L, 3L)
)
ord <- order(
  risk_stage, remaining$dataset_scope, remaining$unit_id, remaining$seed,
  remaining$representation_id, method = "radix"
)
remaining <- remaining[ord, , drop = FALSE]
root_role <- root_role[ord]
relative_file <- relative_file[ord]
risk_stage <- risk_stage[ord]
rownames(remaining) <- NULL

queue <- remaining
queue$contract_id <- "mv08u_full_ph_queue_v1"
queue$production_order <- seq_len(nrow(queue))
queue$mv08r_job_order <- queue$job_order
queue$risk_stage <- risk_stage
queue$risk_stratum <- paste(queue$representation_id, queue$panel_id, sep = ":")
queue$source_cache_root_role <- root_role
queue$source_cache_relative_file <- gsub("\\\\", "/", relative_file)
queue$output_file <- sprintf(
  "ph/%04d__%s.rds", queue$production_order,
  substr(vapply(queue$job_id, digest::digest, character(1L),
                algo = "sha256", serialize = FALSE), 1L, 20L)
)
queue$authorization_state <- "authorized_after_mv08u_commit"
queue$resume_policy <- "strict_completed_prefix_only"
queue$fallback_is_retry <- FALSE
queue$landscape_authorized <- FALSE
queue$comparison_authorized <- FALSE
queue$clustering_authorized <- FALSE
queue$fusion_authorized <- FALSE
queue$labels_authorized <- FALSE
queue$outcomes_authorized <- FALSE
queue <- queue[c(
  "contract_id", "production_order", "mv08r_job_order", "job_id",
  "risk_stage", "risk_stratum", "dataset_scope", "unit_id", "fit_cells",
  "seed", "representation_id", "panel_id", "panel_genes", "panel_sha256",
  "selected_cell_sha256", "matrix_sha256", "distance_sha256", "view_kind",
  "source_contract", "execution_role", "matrix_state",
  "source_origin", "source_cache_root_role", "source_cache_relative_file",
  "source_cache_sha256", "homology_dimensions", "filtration", "threshold",
  "field", "primary_engine", "primary_elapsed_cap_seconds",
  "primary_rss_cap_bytes", "fallback_engine", "fallback_trigger",
  "fallback_elapsed_cap_seconds", "fallback_rss_cap_bytes", "workers",
  "retries", "atomic_write", "resume_policy", "fallback_is_retry",
  "output_file", "authorization_state", "landscape_authorized",
  "comparison_authorized", "clustering_authorized", "fusion_authorized",
  "labels_authorized", "outcomes_authorized", "outcome_label_state",
  "biological_outcomes_computed"
)]
source_units <- queue[!duplicated(queue$unit_id), , drop = FALSE]
root_inventory <- do.call(rbind, lapply(
  sort(unique(source_units$source_cache_root_role)), function(role) {
    x <- source_units[source_units$source_cache_root_role == role, , drop = FALSE]
    data.frame(
      contract_id = "mv08u_source_root_inventory_v1",
      source_cache_root_role = role, source_units = nrow(x),
      expected_cache_bytes = sum(source_bindings$cache_bytes[
        match(x$unit_id, source_bindings$unit_id)
      ]),
      runtime_locator_policy = "root_argument_plus_frozen_relative_file",
      private_absolute_paths_published = FALSE, stringsAsFactors = FALSE
    )
  }
))
runtime_inputs <- data.frame(
  contract_id = "mv08u_runtime_input_binding_v1",
  input_role = "ordered_common475_panel_csv", expected_rows = 475L,
  file_sha256 = unique(s_external$panel_file_sha256),
  runtime_locator_policy = "explicit_launch_argument_no_public_absolute_path",
  stringsAsFactors = FALSE
)

# Reconstruct the seven source-produced sentinel measurements from committed
# MV8-S/MV8-T ledgers. No ignored runtime file is required by this prefreeze.
ph_resource <- t_resources[t_resources$stage == "ph_primary", , drop = FALSE]
ph_resource$sentinel_job_order <- as.integer(sub(
  "^ph_primary__", "", ph_resource$attempt_id
))
sentinel <- s_queue[s_queue$execution_role == "source_produced_gene_ph", , drop = FALSE]
resource_index <- match(sentinel$job_order, ph_resource$sentinel_job_order)
if (nrow(sentinel) != 7L || anyNA(resource_index)) {
  stop("MV8-U source sentinel resource mapping drift", call. = FALSE)
}
sentinel$elapsed_seconds <- ph_resource$elapsed_seconds[resource_index]
sentinel$peak_rss_bytes <- ph_resource$peak_process_tree_rss_bytes[resource_index]
sentinel$output_bytes <- ph_resource$output_bytes[resource_index]
sentinel$selected_engine <- "ripserr"
sentinel$fallback_used <- FALSE
sentinel$stratum <- paste(sentinel$representation_id, sentinel$panel_id, sep = ":")

strata <- sort(unique(sentinel$stratum))
measured <- do.call(rbind, lapply(strata, function(id) {
  x <- sentinel[sentinel$stratum == id, , drop = FALSE]
  data.frame(
    contract_id = "mv08u_sentinel_resource_summary_v1",
    risk_stratum = id, sentinel_records = nrow(x),
    elapsed_min_seconds = min(x$elapsed_seconds),
    elapsed_median_seconds = stats::median(x$elapsed_seconds),
    elapsed_max_seconds = max(x$elapsed_seconds),
    peak_rss_max_bytes = max(x$peak_rss_bytes),
    output_max_bytes = max(x$output_bytes), fallback_records = 0L,
    source_evidence = "MV8-T independently closed MV8-S primary PH",
    stringsAsFactors = FALSE
  )
}))

queue_counts <- as.data.frame(table(queue$risk_stratum), stringsAsFactors = FALSE)
names(queue_counts) <- c("risk_stratum", "production_records")
projection <- merge(queue_counts, measured, by = "risk_stratum", sort = TRUE)
if (nrow(projection) != 4L || sum(projection$production_records) != 1257L) {
  stop("MV8-U projection stratum coverage drift", call. = FALSE)
}
projection$median_projection_seconds <-
  projection$production_records * projection$elapsed_median_seconds
projection$measured_max_projection_seconds <-
  projection$production_records * projection$elapsed_max_seconds
projection$measured_max_storage_bytes <-
  projection$production_records * projection$output_max_bytes
projection$planning_safety_factor <- 2
projection$planning_projection_seconds <-
  2 * projection$measured_max_projection_seconds
projection$workers <- 1L
projection$outcome_label_state <- "closed"
projection$biological_outcomes_computed <- FALSE

aggregate_max_seconds <- sum(projection$measured_max_projection_seconds)
aggregate_planning_seconds <- sum(projection$planning_projection_seconds)
aggregate_storage_bytes <- sum(projection$measured_max_storage_bytes)
aggregate_elapsed_cap_seconds <- 72L * 3600L
private_storage_cap_bytes <- 1024^3

mem_available <- NA_real_
if (file.exists("/proc/meminfo")) {
  line <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
  if (length(line) == 1L) {
    mem_available <- as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1L]][[2L]]) * 1024
  }
}
disk_available <- suppressWarnings(as.numeric(system2(
  "df", c("-B1", "--output=avail", "."), stdout = TRUE
)[[2L]]))
launch_gate <- data.frame(
  contract_id = "mv08u_launch_resource_gate_v1",
  resource = c("available_memory", "workspace_filesystem_available"),
  observed_bytes_at_prefreeze = c(mem_available, disk_available),
  minimum_bytes_at_launch = c(16 * 1024^3, 6 * 1024^3),
  rationale = c(
    "12_GiB_fallback_cap_plus_4_GiB_headroom",
    "1_GiB_private_cap_plus_5_GiB_workspace_headroom"
  ),
  must_recheck_at_launch = TRUE, stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv08u_full_ph_production_prefreeze_v1",
  accepted_parent_head = accepted_parent_head,
  execution_head_state = "bind_exact_committed_head_at_launch_and_publish",
  mv08r_total_records = 1280L, mv08t_closed_records = 23L,
  production_records = 1257L, production_cell_records = 0L,
  production_gene_records = 1257L, production_internal_records = 1236L,
  production_external_records = 21L, source_cache_units = 131L,
  common_panel_file_sha256 = runtime_inputs$file_sha256,
  workers = 1L, within_run_retries = 0L,
  fallback_policy = "GUDHI_only_after_Ripserr_RSS_cap_exceeded",
  primary_elapsed_cap_seconds = 1800L,
  primary_rss_cap_bytes = 8 * 1024^3,
  fallback_elapsed_cap_seconds = 1800L,
  fallback_rss_cap_bytes = 12 * 1024^3,
  measured_max_projection_seconds = aggregate_max_seconds,
  planning_projection_seconds = aggregate_planning_seconds,
  aggregate_elapsed_cap_seconds = aggregate_elapsed_cap_seconds,
  measured_max_storage_projection_bytes = aggregate_storage_bytes,
  private_storage_cap_bytes = private_storage_cap_bytes,
  resume_policy = "atomic_outputs_strict_completed_prefix_no_automatic_retry",
  closure_policy = "MV8-W_reopen_reconstruct_MST_oracle_rehash_1257_and_bind_MV8T_23",
  homology_dimensions = "H0;H1_separate",
  landscape_state = "closed_pending_full_PH_closure_and_Rust_admission",
  comparison_state = "closed", clustering_state = "closed",
  fusion_state = "closed", outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

recovery <- data.frame(
  contract_id = "mv08u_resume_and_recovery_policy_v1",
  condition = c(
    "clean_interruption_between_jobs", "runner_or_child_interrupted_mid_job",
    "primary_rss_cap_exceeded", "primary_timeout", "child_failure",
    "aggregate_elapsed_cap_reached", "private_storage_cap_reached"
  ),
  action = c(
    "resume_after_rehashing_completed_strict_prefix",
    "preserve_all_artifacts_and_require_separate_recovery_prefreeze",
    "permit_one_exact_GUDHI_resource_fallback_for_same_job",
    "stop_preserve_no_fallback_no_retry", "stop_preserve_no_retry",
    "stop_before_next_job_and_preserve", "stop_before_next_job_and_preserve"
  ),
  automatic_retry = FALSE,
  scientific_contract_change = FALSE,
  deletion_allowed = FALSE, stringsAsFactors = FALSE
)

implementation_paths <- c(
  "R/toy_baseline.R", "R/dual_view_topology.R", "R/mv07g_sentinel.R",
  "R/mv07h_full_topology.R", "R/mv08s_ph_sentinel.R",
  "scripts/run_mv08v_ph_entry.R",
  "scripts/run_mv08v_full_ph_production.R",
  "scripts/build_mv08w_full_ph_production_closure.R"
)
if (!all(file.exists(implementation_paths))) {
  stop("MV8-U implementation file absent", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv08u_implementation_binding_v1",
  file = implementation_paths,
  sha256 = vapply(implementation_paths, sha_file, character(1L)),
  execution_state = "bound_not_executed", stringsAsFactors = FALSE
)
input_manifest <- data.frame(
  contract_id = "mv08u_input_manifest_v1", input = inputs,
  bytes = as.numeric(file.info(inputs)$size),
  sha256 = vapply(inputs, sha_file, character(1L)), stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "mv08t_closed", "queue_subtraction", "gene_only", "scope_counts",
    "representation_counts", "stratum_counts", "source_units",
    "source_hashes", "source_locators", "view_hashes", "serial_policy",
    "fallback_policy", "resource_projection", "aggregate_allowance",
    "storage_allowance", "launch_headroom", "resume_policy", "h0_h1",
    "landscape_definition_preserved", "downstream_firewalls",
    "labels_outcomes_closed", "implementation_bound", "no_execution"
  ),
  passed = c(
    nrow(closed) == 23L && all(closed$independently_validated),
    nrow(queue) == 1257L && !any(queue$job_id %in% closed$artifact_id),
    all(queue$view_kind == "gene_topology_v1"),
    sum(queue$dataset_scope == "internal124") == 1236L &&
      sum(queue$dataset_scope == "external8") == 21L,
    sum(queue$representation_id == "sct_data_all_qc_fit_selected384") == 625L &&
      sum(queue$representation_id == "sct_pearson_residual_all_qc_fit_selected384") == 632L,
    identical(as.integer(table(queue$risk_stage)), c(14L, 625L, 618L)),
    length(unique(queue$unit_id)) == 131L,
    all(grepl("^[0-9a-f]{64}$", queue$source_cache_sha256)),
    !anyNA(queue$source_cache_root_role) && !anyNA(queue$source_cache_relative_file) &&
      nrow(root_inventory) == 4L && sum(root_inventory$source_units) == 131L &&
      nrow(runtime_inputs) == 1L && runtime_inputs$expected_rows == 475L &&
      grepl("^[0-9a-f]{64}$", runtime_inputs$file_sha256),
    all(grepl("^[0-9a-f]{64}$", queue$matrix_sha256)) &&
      all(grepl("^[0-9a-f]{64}$", queue$distance_sha256)),
    all(queue$workers == 1L & queue$retries == 0L),
    all(queue$fallback_trigger == "rss_cap_exceeded_only") &&
      all(queue$fallback_rss_cap_bytes == 12 * 1024^3),
    nrow(projection) == 4L && all(projection$sentinel_records >= 1L),
    aggregate_planning_seconds < aggregate_elapsed_cap_seconds,
    aggregate_storage_bytes < private_storage_cap_bytes,
    all(launch_gate$observed_bytes_at_prefreeze >= launch_gate$minimum_bytes_at_launch),
    all(!recovery$automatic_retry & !recovery$deletion_allowed),
    all(queue$homology_dimensions == "H0;H1_separate"),
    all(landscape$owner_approved) &&
      any(landscape$required_state == "all_consecutive_active_levels") &&
      any(landscape$required_state == "no_universal_fixed_grid") &&
      any(landscape$required_state == "no_universal_level_cap"),
    !any(queue$landscape_authorized | queue$comparison_authorized |
         queue$clustering_authorized | queue$fusion_authorized),
    all(!queue$labels_authorized & !queue$outcomes_authorized &
        queue$outcome_label_state == "closed") &&
      !any(queue$biological_outcomes_computed),
    nrow(implementation) == length(implementation_paths), TRUE
  ),
  evidence = c(
    "23 independently reconstructed, repeated, and rehashed MV8-T PH records",
    "exact anti-join from the immutable 1,280-row MV8-R queue",
    "all eight cell records were already closed by MV8-T",
    "1,236 internal plus 21 external records",
    "625 SCT-data plus 632 Pearson-residual records",
    "14 common475, 625 exact500 residual, 618 exact500 data",
    "all remaining views resolve through 131 accepted source caches",
    "every source cache retains its accepted public SHA-256",
    "four root roles plus the ordered common475 runtime CSV are hash/count bound without public private paths",
    "every queue row carries immutable matrix and distance SHA-256",
    "one monitored child; zero automatic retries",
    "Ripserr primary; exact GUDHI only after an 8-GiB RSS stop",
    "all four remaining strata have closed sentinel measurements",
    "two-times measured-max planning projection fits inside 72 hours",
    "measured-max PH output projection fits inside one GiB",
    "current WSL memory and workspace disk exceed prospective minima",
    "strict-prefix resume and preservation-only failure policy",
    "H0 and H1 remain separate in every record",
    "all active levels, no fixed grid/cap, exact or error-controlled L2 preserved",
    "landscapes, comparisons, clustering and fusion remain closed",
    "labels and biological outcomes remain unopened",
    "entry, runner, closure and mathematical helpers SHA-bound",
    "builder reads public receipts only and performs no PH or landscape"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-U validation failed", call. = FALSE)

decision <- data.frame(
  contract_id = "mv08u_decision_v1",
  decision = "authorize_1257_record_serial_MV8V_after_commit_and_launch_recheck",
  ph_records_authorized = 1257L, workers = 1L, retries = 0L,
  landscape_groups_authorized = 0L, comparison_strata_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  next_gate = "MV8V_execution_then_MV8W_independent_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

atomic_csv(contract, file.path(output_dir, "mv08u-contract.csv"))
atomic_csv(queue, file.path(output_dir, "mv08u-full-ph-queue.csv"))
atomic_csv(root_inventory, file.path(output_dir, "mv08u-source-root-inventory.csv"))
atomic_csv(runtime_inputs, file.path(output_dir, "mv08u-runtime-input-bindings.csv"))
atomic_csv(measured, file.path(output_dir, "mv08u-sentinel-resource-summary.csv"))
atomic_csv(projection, file.path(output_dir, "mv08u-resource-projection.csv"))
atomic_csv(launch_gate, file.path(output_dir, "mv08u-launch-resource-gate.csv"))
atomic_csv(recovery, file.path(output_dir, "mv08u-resume-recovery-policy.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08u-implementation-bindings.csv"))
atomic_csv(input_manifest, file.path(output_dir, "mv08u-input-manifest.csv"))
atomic_csv(validation, file.path(output_dir, "mv08u-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08u-decision.csv"))

report <- c(
  "# MV8-U full-PH production prefreeze", "",
  "**Date:** 2026-08-23", "",
  "**Result:** 23/23 prospective checks pass; no PH or landscape was executed.", "",
  "## Exact scope", "",
  paste0(
    "MV8-T already closed 23 of MV8-R's 1,280 PH records. The exact remaining ",
    "queue is 1,257 gene-view records: 1,236 internal and 21 external. All eight ",
    "cell-view records are already closed. The remaining queue contains 625 ",
    "SCT-data and 632 Pearson-residual records across 131 accepted source caches."
  ), "",
  "## Evidence-based execution order and resources", "",
  paste0(
    "Execution is serial and label-blind: 14 common-475 records first, then 625 ",
    "exact-500 residual records, then 618 exact-500 SCT-data records. Every ",
    "remaining stratum was exercised by the closed sentinel. The measured-max ",
    "projection is ", round(aggregate_max_seconds / 3600, 2),
    " hours; the two-times planning projection is ",
    round(aggregate_planning_seconds / 3600, 2),
    " hours inside a 72-hour aggregate stop. The measured-max output projection ",
    "is ", round(aggregate_storage_bytes / 1024^2, 1),
    " MiB inside a 1-GiB private-output cap."
  ), "",
  paste0(
    "Ripserr remains primary with an 1,800-second and 8-GiB child cap. One exact ",
    "GUDHI attempt is permitted only after a Ripserr RSS-cap stop, with an ",
    "1,800-second and 12-GiB cap. A timeout or ordinary child failure stops the ",
    "run without retry."
  ), "",
  "## Resume and closure", "",
  paste0(
    "Outputs are atomic and resume accepts only a fully rehashed completed ",
    "prefix. Ambiguous mid-job evidence is preserved for a separate recovery ",
    "prefreeze. MV8-W must reopen and reconstruct all 1,257 new records, rerun ",
    "every H0 MST oracle, independently rehash them, and bind the 23 closed MV8-T ",
    "records for exact 1,280/1,280 coverage."
  ), "",
  "## Scientific firewall", "",
  paste0(
    "H0 and H1 remain separate. Landscapes retain the dissertation-aligned ",
    "all-active-level, no-fixed-grid, no-level-cap, exact/error-controlled ",
    "streamed definition, but no landscape is authorized here. Comparisons, ",
    "clustering, fusion, labels, outcomes, adoption, and manuscript claims ",
    "remain closed."
  )
)
atomic_text(report, file.path(
  output_dir, "MV08U_FULL_PH_PRODUCTION_PREFREEZE_2026-08-23.md"
))

artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08u-artifact-manifest.csv"]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08u-artifact-manifest.csv"))
cat("MV8-U checks=", sum(validation$passed), "/", nrow(validation),
    "; remaining PH=1257; no execution\n", sep = "")
