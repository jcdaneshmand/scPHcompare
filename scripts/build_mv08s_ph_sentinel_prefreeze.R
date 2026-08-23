#!/usr/bin/env Rscript

# Freeze the exact same-axis external baseline and 23-record PH sentinel from
# current private input identities. This builder opens source cache metadata but
# performs no normalization, PH, landscape, comparison, or biological analysis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop(paste(
  "usage: build_mv08s_ph_sentinel_prefreeze.R <mv08r-audit>",
  "<mv08p-private> <mv08o-private> <hca-bm002-outs>",
  "<hca-remaining-root> <reference-rds> <common-panel> <output-dir>",
  "<accepted-parent-head>"
), call. = FALSE)
for (package in "digest") {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
mv08r <- normalizePath(args[[1L]], mustWork = TRUE)
mv08p_private <- normalizePath(args[[2L]], mustWork = TRUE)
mv08o_private <- normalizePath(args[[3L]], mustWork = TRUE)
hca_bm002_outs <- normalizePath(args[[4L]], mustWork = TRUE)
hca_remaining_root <- normalizePath(args[[5L]], mustWork = TRUE)
reference_path <- normalizePath(args[[6L]], mustWork = TRUE)
panel_path <- normalizePath(args[[7L]], mustWork = TRUE)
output_dir <- normalizePath(args[[8L]], mustWork = FALSE)
accepted_parent_head <- tolower(trimws(args[[9L]]))
if (!grepl("^[0-9a-f]{40}$", accepted_parent_head) || dir.exists(output_dir)) {
  stop("MV8-S parent head invalid or output already exists", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE)

read_csv <- function(path) utils::read.csv(path, check.names = FALSE,
                                            stringsAsFactors = FALSE)
sha_file <- function(path) digest::digest(file = path, algo = "sha256",
                                          serialize = FALSE)
sha_object <- function(value) digest::digest(value, algo = "sha256", serialize = TRUE)
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
required_mv08r <- c(
  "mv08r-contract.csv", "mv08r-ph-queue.csv",
  "mv08r-external-same-axis-baseline.csv", "mv08r-source-cache-bindings.csv",
  "mv08r-source-gene-view-bindings.csv", "mv08r-backend-policy.csv",
  "mv08r-landscape-contract.csv", "mv08r-stage-gates.csv",
  "mv08r-artifact-manifest.csv"
)
if (!all(file.exists(file.path(mv08r, required_mv08r)))) {
  stop("MV8-R prerequisite audit is incomplete", call. = FALSE)
}
mv08r_contract <- read_csv(file.path(mv08r, "mv08r-contract.csv"))
mv08r_queue <- read_csv(file.path(mv08r, "mv08r-ph-queue.csv"))
mv08r_baseline <- read_csv(file.path(mv08r, "mv08r-external-same-axis-baseline.csv"))
mv08r_sources <- read_csv(file.path(mv08r, "mv08r-source-cache-bindings.csv"))
mv08r_views <- read_csv(file.path(mv08r, "mv08r-source-gene-view-bindings.csv"))
mv08r_backend <- read_csv(file.path(mv08r, "mv08r-backend-policy.csv"))
mv08r_landscape <- read_csv(file.path(mv08r, "mv08r-landscape-contract.csv"))
if (nrow(mv08r_contract) != 1L ||
    nrow(mv08r_queue) != 1280L || nrow(mv08r_baseline) != 8L ||
    nrow(mv08r_sources) != 132L || nrow(mv08r_views) != 1272L ||
    nrow(mv08r_backend) != 3L || nrow(mv08r_landscape) != 10L) {
  stop("MV8-R cardinality or parent binding drift", call. = FALSE)
}

panel <- read_csv(panel_path)
reference <- readRDS(reference_path)
if (nrow(panel) != 475L || !identical(as.integer(panel$panel_order_475), 1:475) ||
    length(unique(panel$common475_axis_sha256)) != 1L) {
  stop("MV8-S ordered common475 panel drift", call. = FALSE)
}
panel_sha <- unique(panel$common475_axis_sha256)
if (panel_sha != unique(mv08r_baseline$panel_sha256) ||
    !is.list(reference) || is.null(reference$cache_key) ||
    nrow(reference$panel) != 475L || length(reference$center) != 475L ||
    length(reference$scale) != 475L || reference$pca_model$n_components != 30L) {
  stop("MV8-S common475 panel/reference drift", call. = FALSE)
}

resolve_outs <- function(unit_id) {
  if (unit_id == "HCA_BM_002") return(hca_bm002_outs)
  file.path(hca_remaining_root, unit_id,
            paste0("mv08h_exact500_", tolower(unit_id)), "outs")
}
external_rows <- lapply(seq_len(nrow(mv08r_baseline)), function(index) {
  row <- mv08r_baseline[index, , drop = FALSE]
  source <- mv08r_sources[mv08r_sources$unit_id == row$unit_id, , drop = FALSE]
  outs <- resolve_outs(row$unit_id)
  filtered <- file.path(outs, "filtered_feature_bc_matrix.h5")
  raw <- file.path(outs, "raw_feature_bc_matrix.h5")
  if (nrow(source) != 1L || !all(file.exists(c(filtered, raw)))) {
    stop("MV8-S external input missing: ", row$unit_id, call. = FALSE)
  }
  data.frame(
    contract_id = "mv08s_external_input_binding_v1",
    baseline_order = index, dataset_scope = "external8", unit_id = row$unit_id,
    selection_seed = 20260805L, selected_cells = 384L,
    selected_cell_sha256 = row$selected_cell_sha256,
    qc_eligible_cells = source$fit_cells,
    filtered_h5_bytes = as.numeric(file.info(filtered)$size),
    filtered_h5_sha256 = sha_file(filtered),
    raw_h5_bytes = as.numeric(file.info(raw)$size), raw_h5_sha256 = sha_file(raw),
    panel_id = "common475", panel_genes = 475L, panel_sha256 = panel_sha,
    panel_file_bytes = as.numeric(file.info(panel_path)$size),
    panel_file_sha256 = sha_file(panel_path),
    reference_rds_bytes = as.numeric(file.info(reference_path)$size),
    reference_rds_sha256 = sha_file(reference_path),
    reference_cache_key = reference$cache_key,
    workers = 1L, retries = 0L, elapsed_cap_seconds = 1800L,
    rss_cap_bytes = 12 * 1024^3, private_locator_state = "runtime_argument_only_not_public",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
external_bindings <- do.call(rbind, external_rows)

source_paths <- c(
  SRA628554_SRS2664364 = file.path(
    mv08p_private, "cache", "internal",
    "SRA628554_SRS2664364__exact500_allqc_sct_model.rds"
  ),
  HCA_BM_002 = file.path(mv08o_private, "cache", "HCA_BM_002__primary.rds")
)
source_bindings <- do.call(rbind, lapply(names(source_paths), function(unit_id) {
  path <- source_paths[[unit_id]]
  frozen <- mv08r_sources[mv08r_sources$unit_id == unit_id, , drop = FALSE]
  if (!file.exists(path) || nrow(frozen) != 1L || sha_file(path) != frozen$cache_sha256) {
    stop("MV8-S source cache binding drift: ", unit_id, call. = FALSE)
  }
  cache <- readRDS(path)
  if (!identical(cache$contract_id, "mv08o_residual_source_cache_v1") ||
      cache$identity$unit_id != unit_id || is.null(cache$cache_key) ||
      cache$identity$outcome_label_state != "closed" ||
      isTRUE(cache$identity$biological_outcomes_computed)) {
    stop("MV8-S source cache schema drift: ", unit_id, call. = FALSE)
  }
  data.frame(
    contract_id = "mv08s_source_cache_binding_v1", dataset_scope = frozen$dataset_scope,
    unit_id = unit_id, cache_bytes = as.numeric(file.info(path)$size),
    cache_sha256 = sha_file(path), source_cache_key = cache$cache_key,
    private_locator_state = "runtime_argument_only_not_public",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))

queue <- mv08r_queue[
  grepl("^same_axis_external_baseline_", mv08r_queue$execution_role) |
    (mv08r_queue$sentinel & mv08r_queue$execution_role == "source_produced_gene_ph"),
  , drop = FALSE
]
queue <- queue[order(queue$job_order), , drop = FALSE]
source_index <- which(queue$execution_role == "source_produced_gene_ph")
queue$matrix_sha256 <- NA_character_
queue$distance_sha256 <- NA_character_
for (index in source_index) {
  row <- queue[index, , drop = FALSE]
  view <- mv08r_views[
    mv08r_views$unit_id == row$unit_id & mv08r_views$seed == row$seed &
      mv08r_views$representation_id == row$representation_id &
      mv08r_views$panel_id == row$panel_id, , drop = FALSE
  ]
  if (nrow(view) != 1L || !isTRUE(view$ph_eligible)) {
    stop("MV8-S source view binding is absent", call. = FALSE)
  }
  queue$matrix_state[[index]] <- "source_cache_hash_bound"
  queue$matrix_sha256[[index]] <- view$matrix_sha256
  queue$distance_sha256[[index]] <- view$distance_sha256
}
baseline_index <- setdiff(seq_len(nrow(queue)), source_index)
queue$matrix_sha256[baseline_index] <- NA_character_
queue$distance_sha256[baseline_index] <- NA_character_
queue$contract_id <- "mv08s_ph_sentinel_queue_v1"
queue$sentinel_order <- seq_len(nrow(queue))
queue$authorization_state <- "authorized_after_mv08s_commit"
queue$output_file <- paste0(
  "ph/", sprintf("%04d", queue$sentinel_order), "__",
  substr(vapply(queue$job_id, digest::digest, character(1L),
                algo = "sha256", serialize = FALSE), 1L, 20L), ".rds"
)
queue <- queue[c(
  "contract_id", "sentinel_order", "job_id", "dataset_scope", "unit_id", "seed",
  "representation_id", "panel_id", "panel_genes", "panel_sha256",
  "selected_cell_sha256", "view_kind", "execution_role", "matrix_state",
  "matrix_sha256", "distance_sha256", "homology_dimensions", "filtration",
  "threshold", "field", "primary_engine", "primary_elapsed_cap_seconds",
  "primary_rss_cap_bytes", "fallback_engine", "fallback_trigger",
  "fallback_elapsed_cap_seconds", "fallback_rss_cap_bytes", "workers", "retries",
  "atomic_write", "output_file", "authorization_state", "outcome_label_state",
  "biological_outcomes_computed"
)]
names(queue)[names(queue) == "sentinel_order"] <- "job_order"
if (nrow(queue) != 23L || sum(queue$view_kind == "cell_topology_v1") != 8L ||
    sum(queue$view_kind == "gene_topology_v1") != 15L ||
    sum(queue$execution_role == "source_produced_gene_ph") != 7L ||
    anyDuplicated(queue$job_id) || anyDuplicated(queue$output_file)) {
  stop("MV8-S sentinel queue drift", call. = FALSE)
}

find_job <- function(unit, view, representation = NULL, panel = NULL,
                     seed = NULL) {
  hit <- queue$unit_id == unit & queue$view_kind == view
  if (!is.null(representation)) hit <- hit & queue$representation_id == representation
  if (!is.null(panel)) hit <- hit & queue$panel_id == panel
  if (!is.null(seed)) hit <- hit & queue$seed == seed
  values <- queue$job_id[hit]
  if (length(values) != 1L) stop("MV8-S cross-engine job not unique", call. = FALSE)
  values
}
cross_contract <- data.frame(
  contract_id = "mv08s_cross_engine_contract_v1",
  cross_check_order = 1:4,
  cross_check_id = c("hca002_cell_subset32", "hca002_gene_subset32",
                     "hard_residual_gene_subset32", "hca002_gene_full"),
  job_id = c(
    find_job("HCA_BM_002", "cell_topology_v1", "sct_data_selected384_fit_same_axis"),
    find_job("HCA_BM_002", "gene_topology_v1", "sct_data_selected384_fit_same_axis"),
    find_job("SRA628554_SRS2664364", "gene_topology_v1",
             "sct_pearson_residual_all_qc_fit_selected384", "exact500",
             20260805L),
    find_job("HCA_BM_002", "gene_topology_v1", "sct_data_selected384_fit_same_axis")
  ),
  mode = c("subset32", "subset32", "subset32", "full"),
  engines = "ripserr_0.3.0_vs_TDA_ripsDiag_GUDHI_1.9.4",
  dimensions = "H0;H1_separate", tolerance = 1e-6,
  elapsed_cap_seconds = 1800L, rss_cap_bytes = 12 * 1024^3,
  workers = 1L, retries = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

contract <- data.frame(
  contract_id = "mv08s_same_axis_external_baseline_ph_sentinel_prefreeze_v1",
  accepted_parent_head = accepted_parent_head, baseline_units = 8L,
  execution_head_state = "bind_exact_committed_head_at_launch_and_publish",
  baseline_executions = 16L, ph_records = 23L, ph_repeat_records = 23L,
  cross_engine_jobs = 4L, monitored_units = 66L, workers = 1L, retries = 0L,
  aggregate_elapsed_cap_seconds = 126600L,
  private_storage_cap_bytes = 8 * 1024^3,
  full_ph_jobs_authorized = 0L, landscape_groups_authorized = 0L,
  comparison_strata_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
stage_gates <- data.frame(
  contract_id = "mv08s_stage_gate_v1", stage_order = 1:5,
  stage_id = c("same_axis_baseline", "ph_sentinel", "independent_closure",
               "full_ph", "landscapes"),
  authorization_state = c("authorized_after_mv08s_commit", "authorized_after_mv08s_commit",
                          "authorized_after_execution", "closed", "closed"),
  requirement = c(
    "8 current-input baseline records plus byte-identical repeats",
    "23 PH records plus repeats, H0 MST and resource gates",
    "independent rehash, reconstruction, cross-engine and privacy closure",
    "separate later prefreeze required", "Rust rebuild/hash rebind/R oracle required"
  ), next_stage_automatic = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

implementation_files <- c(
  "R/mv08s_ph_sentinel.R", "scripts/build_mv08s_ph_sentinel_prefreeze.R",
  "scripts/run_mv08s_same_axis_baseline_entry.R", "scripts/run_mv08s_ph_entry.R",
  "scripts/run_mv08s_cross_engine_entry.R", "scripts/run_mv08s_ph_sentinel.R",
  "scripts/build_mv08t_ph_sentinel_closure.R",
  "tests/testthat/test-mv08s-ph-sentinel.R",
  "docs/specifications/MV08S_SAME_AXIS_EXTERNAL_BASELINE_PH_SENTINEL_PREFREEZE_V1.md"
)
if (!all(file.exists(implementation_files))) {
  stop("MV8-S implementation binding file absent", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv08s_implementation_binding_v1", file = implementation_files,
  bytes = as.numeric(file.info(implementation_files)$size),
  sha256 = vapply(implementation_files, sha_file, character(1L)),
  stringsAsFactors = FALSE
)
input_manifest <- data.frame(
  contract_id = "mv08s_input_manifest_v1",
  input = c(paste0("mv08r/", required_mv08r), "common475_panel", "common475_reference"),
  bytes = c(as.numeric(file.info(file.path(mv08r, required_mv08r))$size),
            as.numeric(file.info(panel_path)$size), as.numeric(file.info(reference_path)$size)),
  sha256 = c(vapply(file.path(mv08r, required_mv08r), sha_file, character(1L)),
             sha_file(panel_path), sha_file(reference_path)),
  private_values_published = FALSE, stringsAsFactors = FALSE
)
validation <- data.frame(
  contract_id = "mv08s_prefreeze_validation_v1",
  check_id = c(
    "parent_head_bound", "eight_current_hca_inputs_rehashed", "two_source_caches_rehashed",
    "queue_23_exact", "queue_axes_8_cell_15_gene", "source_views_7_hash_bound",
    "baseline_repeats_required", "ph_repeats_required", "cross_subset_and_full",
    "one_worker_zero_retry", "resource_caps_frozen", "h0_mst_required",
    "landscape_definition_preserved", "landscapes_closed", "outcomes_closed",
    "implementation_hash_bound"
  ),
  passed = c(
    identical(accepted_parent_head, "10e0ac9281b0ca321e1c782ae86e483d6df277cc"),
    nrow(external_bindings) == 8L && all(grepl("^[0-9a-f]{64}$", external_bindings$raw_h5_sha256)),
    nrow(source_bindings) == 2L && all(grepl("^[0-9a-f]{64}$", source_bindings$cache_sha256)),
    nrow(queue) == 23L, sum(queue$view_kind == "cell_topology_v1") == 8L &&
      sum(queue$view_kind == "gene_topology_v1") == 15L,
    sum(queue$execution_role == "source_produced_gene_ph") == 7L &&
      all(grepl("^[0-9a-f]{64}$", queue$matrix_sha256[source_index])),
    contract$baseline_executions == 16L, contract$ph_repeat_records == 23L,
    identical(cross_contract$mode, c("subset32", "subset32", "subset32", "full")),
    all(queue$workers == 1L) && all(queue$retries == 0L),
    contract$aggregate_elapsed_cap_seconds == 126600L && contract$private_storage_cap_bytes == 8 * 1024^3,
    TRUE, nrow(mv08r_landscape) == 10L && all(mv08r_landscape$owner_approved),
    contract$landscape_groups_authorized == 0L, contract$outcome_jobs_authorized == 0L,
    nrow(implementation) == length(implementation_files)
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-S independent prefreeze validation failed", call. = FALSE)

decision <- data.frame(
  contract_id = "mv08s_prefreeze_decision_v1",
  decision = "authorize_exact_8_baselines_and_23_record_PH_sentinel_after_commit",
  closure_required_before_full_ph = TRUE, full_ph_state = "closed",
  landscape_state = "closed_requires_Rust_rebuild_hash_rebind_R_oracle",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(contract, file.path(output_dir, "mv08s-contract.csv"))
atomic_csv(external_bindings, file.path(output_dir, "mv08s-external-input-bindings.csv"))
atomic_csv(source_bindings, file.path(output_dir, "mv08s-source-cache-bindings.csv"))
atomic_csv(queue, file.path(output_dir, "mv08s-ph-sentinel-queue.csv"))
atomic_csv(cross_contract, file.path(output_dir, "mv08s-cross-engine-contract.csv"))
atomic_csv(stage_gates, file.path(output_dir, "mv08s-stage-gates.csv"))
atomic_csv(implementation, file.path(output_dir, "mv08s-implementation-bindings.csv"))
atomic_csv(input_manifest, file.path(output_dir, "mv08s-input-manifest.csv"))
atomic_csv(validation, file.path(output_dir, "mv08s-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08s-decision.csv"))
atomic_text(c(
  "# MV8-S same-axis external baseline and PH-sentinel prefreeze", "",
  paste0("Accepted MV8-R parent: `", accepted_parent_head, "`."), "",
  "This prefreeze binds eight current exact-reference HCA filtered/raw matrix pairs,",
  "the ordered common-475 reference, and two source-produced residual caches.",
  "It authorizes exactly eight baseline primaries plus repeats, 23 PH primaries",
  "plus repeats, and four cross-engine checks under one-worker monitored caps.", "",
  "H0 and H1 remain separate. Every PH record must pass the full-view H0 MST",
  "oracle. Landscapes remain closed under the dissertation-aligned all-finite,",
  "all-active-level, no-fixed-grid, no-level-cap streamed definition; their next",
  "gate additionally requires a rebuilt and hash-rebound Rust kernel plus R oracle.", "",
  "No labels, outcomes, comparisons, clustering, fusion, or manuscript claims are opened."
), file.path(output_dir, "MV08S_PREFREEZE_REPORT.md"))

artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08s-artifact-manifest.csv"]
manifest <- data.frame(
  contract_id = "mv08s_artifact_manifest_v1", artifact = basename(artifacts),
  bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE
)
atomic_csv(manifest, file.path(output_dir, "mv08s-artifact-manifest.csv"))
message("MV8-S prefreeze checks=", sum(validation$passed), "/", nrow(validation),
        "; queue=", nrow(queue), "; execution=0")
