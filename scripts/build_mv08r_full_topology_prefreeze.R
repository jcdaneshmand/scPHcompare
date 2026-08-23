#!/usr/bin/env Rscript

# Build the metadata-only MV8-R topology-production prefreeze. This script
# reconciles public source/view receipts and prior topology evidence. It never
# opens a private expression cache and never calls PH or landscape code.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: build_mv08r_full_topology_prefreeze.R <output-dir>",
       call. = FALSE)
}

output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
if (dir.exists(output_dir)) {
  stop("refusing to overwrite MV8-R prefreeze output", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE)

read_csv <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
sha_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) {
    stop("failed to publish ", basename(path), call. = FALSE)
  }
}
atomic_text <- function(x, path) {
  partial <- paste0(path, ".partial")
  writeLines(x, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) {
    stop("failed to publish ", basename(path), call. = FALSE)
  }
}

accepted_head <- tolower(trimws(Sys.getenv("MV08R_GIT_HEAD", unset = "")))
if (!nzchar(accepted_head)) {
  accepted_head <- tolower(trimws(system2(
    "git", c("rev-parse", "HEAD"), stdout = TRUE
  )))
}
if (!grepl("^[0-9a-f]{40}$", accepted_head)) {
  stop("cannot bind accepted git HEAD", call. = FALSE)
}

mv08n_root <- "docs/audits/mv08n-pearson-residual-migration-prefreeze-v1"
mv08o_root <- "docs/audits/mv08o-residual-source-sentinel-closure-v1"
mv08q_root <- "docs/audits/mv08q-full-source-production-closure-v1"
mv07h_root <- "docs/audits/mv07h-full-ph-evidence"
mv07hl_root <- "docs/audits/mv07h-landscape-complete-validation"
mv08i_root <- "docs/audits/mv08i-hca-validation-v1"
mv08h_root <- "docs/audits/mv08h-common475-topology-review-v1"

required_inputs <- c(
  file.path(mv08n_root, "mv08n-contract.csv"),
  file.path(mv08n_root, "mv08n-gene-view-queue.csv"),
  file.path(mv08n_root, "mv08n-landscape-queue.csv"),
  file.path(mv08n_root, "mv08n-comparison-contract.csv"),
  file.path(mv08o_root, "mv08o-source-sentinel-resource.csv"),
  file.path(mv08o_root, "mv08o-source-sentinel-view-summary.csv"),
  file.path(mv08q_root, "mv08q-source-cache-rehash.csv"),
  file.path(mv08q_root, "mv08q-source-view-summary.csv"),
  file.path(mv08q_root, "mv08q-validation.csv"),
  file.path(mv07h_root, "mv07h-ph-metrics.csv"),
  file.path(mv07h_root, "mv07h-ph-engine-attempts.csv"),
  file.path(mv07hl_root, "mv07h-landscape-complete-group-inventory.csv"),
  file.path(mv08i_root, "mv08i-input-identity.csv"),
  file.path(mv08i_root, "mv08i-topology-summary.csv"),
  file.path(mv08h_root, "mv08h-topology-review-identity.csv")
)
if (!all(file.exists(required_inputs))) {
  stop("one or more frozen MV8-R public inputs are absent", call. = FALSE)
}

mv08n_contract <- read_csv(required_inputs[[1L]])
mv08n_views <- read_csv(required_inputs[[2L]])
mv08n_landscapes <- read_csv(required_inputs[[3L]])
mv08n_comparisons <- read_csv(required_inputs[[4L]])
mv08o_resource <- read_csv(required_inputs[[5L]])
mv08o_views <- read_csv(required_inputs[[6L]])
mv08q_sources <- read_csv(required_inputs[[7L]])
mv08q_views <- read_csv(required_inputs[[8L]])
mv08q_validation <- read_csv(required_inputs[[9L]])
mv07h_ph <- read_csv(required_inputs[[10L]])
mv07h_attempts <- read_csv(required_inputs[[11L]])
mv07h_landscape_inventory <- read_csv(required_inputs[[12L]])
mv08i_inputs <- read_csv(required_inputs[[13L]])
mv08i_topology <- read_csv(required_inputs[[14L]])
mv08h_identity <- read_csv(required_inputs[[15L]])

if (nrow(mv08n_contract) != 1L || nrow(mv08n_views) != 1280L ||
    nrow(mv08n_landscapes) != 28L || nrow(mv08n_comparisons) != 40L ||
    nrow(mv08q_sources) != 129L || nrow(mv08q_views) != 1248L ||
    !all(mv08q_validation$passed) || nrow(mv07h_ph) != 1240L ||
    nrow(mv07h_landscape_inventory) != 20L || nrow(mv08i_inputs) != 8L ||
    nrow(mv08i_topology) != 32L || nrow(mv08h_identity) != 1L) {
  stop("MV8-R prerequisite cardinality or closure drift", call. = FALSE)
}

exact_sha <- "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e"
common_sha <- "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"
seeds <- 20260805:20260809
hard_unit <- "SRA628554_SRS2664364"

# Bind the 129 production caches plus the three accepted MV8-O primary caches.
primary <- mv08o_resource[
  mv08o_resource$run_role == "primary",
  c("unit_id", "dataset_scope", "cache_bytes", "cache_sha256",
    "worker_audit_sha256", "outcome_label_state",
    "biological_outcomes_computed"),
  drop = FALSE
]
primary$fit_cells <- vapply(primary$unit_id, function(id) {
  as.numeric(unique(mv08o_views$fit_cells[
    mv08o_views$run_role == "primary" & mv08o_views$unit_id == id
  ]))
}, numeric(1L))
primary$source_origin <- "mv08o_primary"
production <- mv08q_sources[c(
  "unit_id", "dataset_scope", "cache_bytes", "cache_sha256",
  "worker_audit_sha256", "outcome_label_state",
  "biological_outcomes_computed", "fit_cells"
)]
production$source_origin <- "mv08q_production"
source_bindings <- rbind(production, primary[names(production)])
source_bindings <- source_bindings[
  order(source_bindings$dataset_scope, source_bindings$unit_id,
        method = "radix"), , drop = FALSE
]
rownames(source_bindings) <- NULL
if (nrow(source_bindings) != 132L || anyDuplicated(source_bindings$unit_id) ||
    sum(source_bindings$dataset_scope == "internal124") != 124L ||
    sum(source_bindings$dataset_scope == "external8") != 8L ||
    any(source_bindings$outcome_label_state != "closed") ||
    any(source_bindings$biological_outcomes_computed) ||
    any(!grepl("^[0-9a-f]{64}$", source_bindings$cache_sha256))) {
  stop("full 132-source cache binding drift", call. = FALSE)
}
source_bindings$source_order <- seq_len(nrow(source_bindings))
source_bindings$contract_id <- "mv08r_source_cache_binding_v1"
source_bindings$private_locator_state <- "runtime_argument_only_not_public"
source_bindings <- source_bindings[c(
  "contract_id", "source_order", "dataset_scope", "unit_id", "fit_cells",
  "source_origin", "cache_bytes", "cache_sha256", "worker_audit_sha256",
  "private_locator_state", "outcome_label_state",
  "biological_outcomes_computed"
)]

# Reconcile source-produced geometry receipts. Common475 rows were emitted with
# the exact500 hash by the worker; panel identity is corrected from the frozen
# ordered panel, while matrix/distance hashes remain immutable evidence.
sentinel_views <- mv08o_views[mv08o_views$run_role == "primary", , drop = FALSE]
source_views <- rbind(mv08q_views, sentinel_views)
source_views <- source_views[
  order(source_views$dataset_scope, source_views$unit_id, source_views$seed,
        source_views$representation_id, source_views$panel_id,
        method = "radix"), , drop = FALSE
]
rownames(source_views) <- NULL
if (nrow(source_views) != 1272L ||
    anyDuplicated(paste(source_views$dataset_scope, source_views$unit_id,
                        source_views$seed, source_views$representation_id,
                        source_views$panel_id, sep = "\r")) ||
    sum(source_views$dataset_scope == "internal124") != 1240L ||
    sum(source_views$dataset_scope == "external8") != 32L) {
  stop("source-produced view axis drift", call. = FALSE)
}
source_views$observed_panel_sha256 <- source_views$panel_sha256
source_views$panel_sha256 <- ifelse(
  source_views$panel_id == "common475", common_sha, exact_sha
)
source_views$panel_metadata_reconciled <-
  source_views$observed_panel_sha256 != source_views$panel_sha256
source_views$ph_eligible <- source_views$dataset_scope == "internal124" |
  source_views$panel_id == "common475" |
  (source_views$panel_id == "exact500" &
     source_views$representation_id ==
       "sct_pearson_residual_all_qc_fit_selected384")
source_views$authorization_state <- ifelse(
  source_views$ph_eligible, "closed_until_mv08r_commit", "diagnostic_no_ph"
)
source_views$contract_id <- "mv08r_source_gene_view_binding_v1"
source_views$view_order <- seq_len(nrow(source_views))
source_views$outcome_label_state <- "closed"
source_views$biological_outcomes_computed <- FALSE
source_views <- source_views[c(
  "contract_id", "view_order", "dataset_scope", "unit_id", "seed",
  "representation_id", "panel_id", "panel_genes", "panel_sha256",
  "observed_panel_sha256", "panel_metadata_reconciled", "selected_cells",
  "selected_cell_sha256", "matrix_sha256", "distance_sha256",
  "values_finite", "zero_variance_gene_count", "correlation_chord_valid",
  "ph_eligible", "authorization_state", "outcome_label_state",
  "biological_outcomes_computed"
)]

# The old MV8-I common475 baseline used a different processed matrix and cell
# selection. Reconstruct selected-384-fit SCT data from each current exact-
# reference raw-read unit so external fit/layer effects share one cell axis.
external_axis <- unique(source_views[
  source_views$dataset_scope == "external8",
  c("unit_id", "seed", "selected_cells", "selected_cell_sha256"),
  drop = FALSE
])
external_axis <- external_axis[order(external_axis$unit_id), , drop = FALSE]
old_axis <- mv08i_inputs[c("unit_id", "selected_cell_sha256")]
names(old_axis)[[2L]] <- "mv08i_selected_cell_sha256"
external_axis <- merge(external_axis, old_axis, by = "unit_id", sort = FALSE)
external_axis <- external_axis[order(external_axis$unit_id), , drop = FALSE]
external_axis$old_axis_equal <-
  external_axis$selected_cell_sha256 == external_axis$mv08i_selected_cell_sha256
if (nrow(external_axis) != 8L || any(external_axis$old_axis_equal) ||
    any(external_axis$selected_cells != 384L)) {
  stop("external same-axis reconciliation drift", call. = FALSE)
}
external_baseline <- data.frame(
  contract_id = "mv08r_external_same_axis_baseline_v1",
  baseline_order = seq_len(8L), dataset_scope = "external8",
  unit_id = external_axis$unit_id, seed = as.integer(external_axis$seed),
  selected_cells = as.integer(external_axis$selected_cells),
  selected_cell_sha256 = external_axis$selected_cell_sha256,
  prior_mv08i_selected_cell_sha256 = external_axis$mv08i_selected_cell_sha256,
  prior_axis_eligible_for_paired_effect = FALSE,
  source_input = "current_exact_reference_raw_read_filtered_plus_raw_h5",
  fit_scope = "same_frozen_selected384_only",
  representation_id = "sct_data_selected384_fit_same_axis",
  panel_id = "common475", panel_genes = 475L,
  panel_sha256 = common_sha,
  cell_view = "frozen_internal_common475_projection",
  gene_view = "pearson_correlation_chord",
  authorization_state = "closed_until_mv08r_commit",
  workers = 1L, retries = 0L, elapsed_cap_seconds = 1800L,
  rss_cap_bytes = 12 * 1024^3, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

# Build 1,272 gene PH records: 1,264 eligible source-produced records and
# eight corrected same-axis selected-fit external baselines. Add eight external
# cell records from that same baseline substage, for 1,280 total PH jobs.
gene_source <- source_views[source_views$ph_eligible, , drop = FALSE]
gene_source$view_kind <- "gene_topology_v1"
gene_source$source_contract <- "mv08r_source_gene_view_binding_v1"
gene_source$execution_role <- "source_produced_gene_ph"
gene_source$matrix_state <- "immutable_source_cache"
gene_source <- gene_source[c(
  "dataset_scope", "unit_id", "seed", "representation_id", "panel_id",
  "panel_genes", "panel_sha256", "selected_cell_sha256", "view_kind",
  "source_contract", "execution_role", "matrix_state"
)]
gene_baseline <- data.frame(
  dataset_scope = "external8", unit_id = external_baseline$unit_id,
  seed = external_baseline$seed,
  representation_id = external_baseline$representation_id,
  panel_id = "common475", panel_genes = 475L, panel_sha256 = common_sha,
  selected_cell_sha256 = external_baseline$selected_cell_sha256,
  view_kind = "gene_topology_v1",
  source_contract = "mv08r_external_same_axis_baseline_v1",
  execution_role = "same_axis_external_baseline_gene_ph",
  matrix_state = "construct_before_ph", stringsAsFactors = FALSE
)
cell_baseline <- gene_baseline
cell_baseline$view_kind <- "cell_topology_v1"
cell_baseline$execution_role <- "same_axis_external_baseline_cell_ph"
ph_queue <- rbind(gene_baseline, cell_baseline, gene_source)
ph_queue$sentinel <-
  (ph_queue$unit_id == hard_unit &
     ph_queue$seed %in% c(20260805L, 20260806L)) |
  ph_queue$unit_id == "HCA_BM_002"
ph_queue$stage_order <- ifelse(
  grepl("same_axis_external_baseline", ph_queue$execution_role), 1L,
  ifelse(ph_queue$sentinel, 2L, 3L)
)
ph_queue <- ph_queue[
  order(ph_queue$stage_order, ph_queue$dataset_scope, ph_queue$unit_id,
        ph_queue$seed, ph_queue$representation_id, ph_queue$view_kind,
        method = "radix"), , drop = FALSE
]
rownames(ph_queue) <- NULL
ph_queue$job_order <- seq_len(nrow(ph_queue))
ph_queue$job_id <- paste(
  "mv08r_ph", ph_queue$dataset_scope, ph_queue$unit_id, ph_queue$seed,
  ph_queue$representation_id, ph_queue$panel_id, ph_queue$view_kind,
  sep = ":"
)
ph_queue$homology_dimensions <- "H0;H1_separate"
ph_queue$filtration <- "complete_vietoris_rips"
ph_queue$threshold <- -1
ph_queue$field <- 2L
ph_queue$primary_engine <- "ripserr_0.3.0"
ph_queue$primary_elapsed_cap_seconds <- ifelse(
  ph_queue$view_kind == "cell_topology_v1", 600L, 1800L
)
ph_queue$primary_rss_cap_bytes <- ifelse(
  ph_queue$view_kind == "cell_topology_v1", 4 * 1024^3, 8 * 1024^3
)
ph_queue$fallback_engine <- ifelse(
  ph_queue$view_kind == "gene_topology_v1", "TDA_ripsDiag_GUDHI_1.9.4", "none"
)
ph_queue$fallback_trigger <- ifelse(
  ph_queue$view_kind == "gene_topology_v1", "rss_cap_exceeded_only", "none"
)
ph_queue$fallback_elapsed_cap_seconds <- ifelse(
  ph_queue$view_kind == "gene_topology_v1", 1800L, NA_integer_
)
ph_queue$fallback_rss_cap_bytes <- ifelse(
  ph_queue$view_kind == "gene_topology_v1", 12 * 1024^3, NA_real_
)
ph_queue$workers <- 1L
ph_queue$retries <- 0L
ph_queue$atomic_write <- TRUE
ph_queue$authorization_state <- ifelse(
  ph_queue$stage_order <= 2L, "authorized_after_mv08r_commit",
  "closed_pending_sentinel_closure"
)
ph_queue$outcome_label_state <- "closed"
ph_queue$biological_outcomes_computed <- FALSE
ph_queue$contract_id <- "mv08r_ph_queue_v1"
ph_queue <- ph_queue[c(
  "contract_id", "job_order", "job_id", "stage_order", "sentinel",
  "dataset_scope", "unit_id", "seed", "representation_id", "panel_id",
  "panel_genes", "panel_sha256", "selected_cell_sha256", "view_kind",
  "source_contract", "execution_role", "matrix_state",
  "homology_dimensions", "filtration", "threshold", "field",
  "primary_engine", "primary_elapsed_cap_seconds", "primary_rss_cap_bytes",
  "fallback_engine", "fallback_trigger", "fallback_elapsed_cap_seconds",
  "fallback_rss_cap_bytes", "workers", "retries", "atomic_write",
  "authorization_state", "outcome_label_state",
  "biological_outcomes_computed"
)]
if (nrow(ph_queue) != 1280L || anyDuplicated(ph_queue$job_id) ||
    sum(ph_queue$view_kind == "gene_topology_v1") != 1272L ||
    sum(ph_queue$view_kind == "cell_topology_v1") != 8L ||
    sum(ph_queue$execution_role == "source_produced_gene_ph") != 1264L ||
    sum(ph_queue$execution_role ==
          "same_axis_external_baseline_gene_ph") != 8L ||
    sum(ph_queue$authorization_state ==
          "authorized_after_mv08r_commit") != 23L) {
  stop("corrected MV8-R PH queue drift", call. = FALSE)
}

# Bind existing internal selected-fit gene and cell PH as immutable comparator
# evidence; these records are not recomputed by MV8-R.
comparator_bindings <- mv07h_ph[c(
  "seed", "sample_id", "view_id", "point_count", "diagram_sha256",
  "ph_cache_key", "output_sha256", "outcome_label_state",
  "biological_outcomes_computed"
)]
names(comparator_bindings)[names(comparator_bindings) == "sample_id"] <- "unit_id"
comparator_bindings$dataset_scope <- "internal124"
comparator_bindings$representation_id <- ifelse(
  comparator_bindings$view_id == "gene_topology_v1",
  "existing_selectedfit_data_exact500", "immutable_cell_topology_exact500"
)
comparator_bindings$panel_id <- "exact500"
comparator_bindings$reuse_state <- "immutable_reuse_no_recompute"
comparator_bindings$contract_id <- "mv08r_internal_comparator_binding_v1"
comparator_bindings <- comparator_bindings[c(
  "contract_id", "dataset_scope", "unit_id", "seed", "representation_id",
  "panel_id", "view_id", "point_count", "diagram_sha256", "ph_cache_key",
  "output_sha256", "reuse_state", "outcome_label_state",
  "biological_outcomes_computed"
)]

# Preserve the revised landscape definition verbatim and keep the maximum-
# burden stress group selection label-blind until PH metrics exist.
landscape_queue <- mv08n_landscapes
landscape_queue$representation_id[
  landscape_queue$representation_id == "sct_data_selected384_fit"
] <- "sct_data_selected384_fit_same_axis"
landscape_queue$integration <- "exact_streamed_squared_L2"
landscape_queue$level_policy <- "all_consecutive_active_levels"
landscape_queue$grid_policy <- "none"
landscape_queue$level_cap <- "none"
landscape_queue$essential_h0_policy <- "exclude_infinite_interval"
landscape_queue$dimension_policy <- "H0_H1_separate"
landscape_queue$stress_order <- "defer_label_blind_maximum_interval_burden"
landscape_queue$engine_policy <-
  "accepted_rust_kernel_after_rebuild_rebind;R_exact_or_error_controlled_oracle"
landscape_queue$rust_binary_state <- "absent_rebuild_required_before_execution"
landscape_queue$authorization_state <- "closed_pending_full_ph_and_rust_rebind"
landscape_queue$contract_id <- "mv08r_landscape_queue_v1"

comparison_firewall <- mv08n_comparisons
comparison_firewall$left_stack[
  comparison_firewall$left_stack == "selectedfit_data_common475"
] <- "same_axis_selectedfit_data_common475"
comparison_firewall$right_stack[
  comparison_firewall$right_stack == "selectedfit_data_common475"
] <- "same_axis_selectedfit_data_common475"
comparison_firewall$axis_policy <- ifelse(
  comparison_firewall$dataset_scope == "external8",
  "same_exact_reference_input_and_selected384_digest_required",
  "same_internal_sample_seed_axis_required"
)
comparison_firewall$authorization_state <-
  "closed_pending_complete_immutable_landscapes"
comparison_firewall$contract_id <- "mv08r_comparison_firewall_v1"

backend_policy <- data.frame(
  contract_id = "mv08r_backend_policy_v1",
  view_kind = c("cell_topology_v1", "gene_topology_v1",
                "gene_topology_v1"),
  attempt_scope = c("primary", "primary", "resource_fallback"),
  engine = c("ripserr_0.3.0", "ripserr_0.3.0",
             "TDA_ripsDiag_GUDHI_1.9.4"),
  trigger = c("always", "always", "primary_rss_cap_exceeded_only"),
  elapsed_cap_seconds = c(600L, 1800L, 1800L),
  rss_cap_bytes = c(4, 8, 12) * 1024^3,
  workers = 1L, retries = 0L,
  mathematical_estimand =
    "complete_vietoris_rips_H0_H1_field_2_threshold_minus_1",
  essential_h0_normalization = c(
    "native_infinite", "native_infinite", "replace_capped_class_with_infinite"
  ),
  substitution_allowed = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

landscape_contract <- data.frame(
  contract_id = "mv08r_landscape_contract_v1",
  item = c(
    "finite_intervals", "essential_h0", "level_policy", "integration",
    "dimension_policy", "grid_policy", "level_cap_policy", "streaming",
    "combined_summary", "engine_gate"
  ),
  required_state = c(
    "all_finite_positive_persistence_intervals",
    "exclude_infinite_interval",
    "all_consecutive_active_levels",
    "exact_or_error_controlled_squared_l2_on_dimension_support",
    "h0_h1_separate_primary_outputs", "no_universal_fixed_grid",
    "no_universal_level_cap",
    "stream_or_chunk_without_dense_landscape_materialization",
    "secondary_only_after_separate_h0_h1",
    "rust_rebuild_hash_rebind_and_R_oracle_before_production"
  ),
  owner_approved = TRUE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

reconciliation <- data.frame(
  contract_id = "mv08r_metadata_reconciliation_v1",
  discrepancy = c(
    "common475_worker_panel_sha256",
    "external_selectedfit_baseline_cell_axis",
    "current_rust_binary_availability"
  ),
  observed_state = c(
    "common475_rows_recorded_exact500_sha256",
    "MV8I_selected384_digests_differ_from_current_exact_reference_all_8_units",
    "accepted_binary_not_present_in_current_WSL_workspace"
  ),
  resolution = c(
    "bind_common475_to_frozen_ordered_common475_sha_without_altering_matrix_distance_hashes",
    "reconstruct_selected384_fit_baseline_from_current_exact_reference_raw_reads_on_identical_axes",
    "rebuild_from_locked_source_then_require_hash_equivalence_and_R_oracles_before_landscapes"
  ),
  changes_estimand = FALSE,
  historical_artifacts_mutated = FALSE,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

stage_gates <- data.frame(
  contract_id = "mv08r_stage_gate_v1", stage_order = 1:7,
  stage_id = c(
    "source_closure", "same_axis_external_baseline", "ph_sentinel",
    "full_ph", "landscape_stress", "landscape_completion",
    "paired_comparison"
  ),
  authorization_state = c(
    "complete", "authorized_after_commit", "authorized_after_commit",
    rep("closed", 4L)
  ),
  pass_requirement = c(
    "132 source caches and required geometry hash-bound",
    "8 exact-reference selectedfit baselines match frozen selected384 digests",
    "23 jobs complete; repeats, MST oracle, subset/full cross-engine and caps pass",
    "1280/1280 atomic PH records complete with H0/H1 separate",
    "accepted Rust binary rebound; max interval-burden group plus R oracle and repeat pass",
    "28/28 groups complete; immutable resume and component balance pass",
    "all prerequisite distance hashes complete; labels and outcomes still closed"
  ),
  next_stage_automatic = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)

implementation_paths <- c(
  "R/dual_view_topology.R", "R/landscape_reference.R",
  "R/landscape_public_api.R", "R/landscape_rust_prototype.R",
  "R/mv07h_full_topology.R", "scripts/run_mv07h_full_ph.R",
  "scripts/run_mv07h_landscape_group.R",
  "scripts/run_mv08h_common475_topology_review.R",
  "rust/scph_landscape_kernel/src/lib.rs",
  "rust/scph_landscape_kernel/Cargo.toml",
  "rust/scph_landscape_kernel/Cargo.lock"
)
if (!all(file.exists(implementation_paths))) {
  stop("MV8-R implementation binding file absent", call. = FALSE)
}
implementation <- data.frame(
  contract_id = "mv08r_implementation_binding_v1",
  role = c(
    "typed_view_and_ripserr", "canonical_R_landscape_oracle",
    "public_landscape_contract", "Rust_FFI_shim", "GUDHI_fallback_pattern",
    "resource_monitored_PH_pattern", "streamed_landscape_pattern",
    "same_axis_external_dual_view_pattern", "Rust_kernel_source",
    "Rust_build_manifest", "Rust_lockfile"
  ),
  file = implementation_paths,
  sha256 = vapply(implementation_paths, sha_file, character(1L)),
  execution_state = "bound_not_executed", stringsAsFactors = FALSE
)

package_version_or_absent <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    "absent"
  }
}
runtime <- data.frame(
  contract_id = "mv08r_runtime_audit_v1",
  component = c("R", "ripserr", "TDA", "digest", "cargo", "rustc",
                "accepted_Rust_binary"),
  observed_version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    package_version_or_absent("ripserr"), package_version_or_absent("TDA"),
    package_version_or_absent("digest"),
    if (nzchar(Sys.which("cargo"))) "present" else "absent",
    if (nzchar(Sys.which("rustc"))) "present" else "absent", "absent"
  ),
  required_before_stage = c(
    "ph_sentinel", "ph_sentinel", "ph_sentinel", "prefreeze",
    "landscape_stress", "landscape_stress", "landscape_stress"
  ),
  current_gate_passed = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  resolution = c(
    "none", "none", "none", "none",
    "reinstall_or_restore_locked_Rust_toolchain", "same_as_cargo",
    "rebuild_and_hash_rebind_from_locked_source"
  ), stringsAsFactors = FALSE
)
if (!identical(runtime$observed_version[2:4], c("0.3.0", "1.9.4", "0.6.37"))) {
  stop("current WSL R backend versions differ from the accepted PH environment",
       call. = FALSE)
}

contract <- data.frame(
  contract_id = "mv08r_full_topology_production_prefreeze_v1",
  accepted_head = accepted_head, source_fits = 132L,
  source_produced_views = 1272L, diagnostic_no_ph_views = 8L,
  same_axis_external_baselines = 8L, gene_ph_records = 1272L,
  external_cell_ph_records = 8L, total_ph_jobs = 1280L,
  ph_sentinel_jobs = 23L, internal_immutable_comparator_records = 1240L,
  landscape_groups = 28L, internal_landscape_groups = 20L,
  external_landscape_groups = 8L, paired_comparison_strata = 40L,
  homology_dimensions = "H0;H1_separate",
  filtration = "complete_vietoris_rips", threshold = -1, field = 2L,
  landscape_definition = paste(
    "finite_positive;essential_H0_excluded;all_active_levels;",
    "exact_or_error_controlled_squared_L2;no_grid;no_level_cap;streamed",
    sep = ""
  ),
  external_axis_policy = "same_raw_read_input_same_selected384_digest",
  comparison_state = "closed", clustering_state = "closed",
  fusion_state = "closed", default_adoption_state = "closed",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check_id = c(
    "source_closure", "source_scope", "source_views", "view_geometry",
    "diagnostic_exclusion", "panel_reconciliation", "external_axis_drift",
    "same_axis_baseline", "ph_queue", "ph_dimensions", "backend_policy",
    "internal_comparators", "landscape_contract", "landscape_queue",
    "comparison_firewall", "runtime_ph_ready", "rust_gate_closed",
    "implementation_bound", "label_outcome_firewall", "no_execution"
  ),
  passed = c(
    nrow(source_bindings) == 132L,
    sum(source_bindings$dataset_scope == "internal124") == 124L &&
      sum(source_bindings$dataset_scope == "external8") == 8L,
    nrow(source_views) == 1272L,
    all(source_views$values_finite) &&
      all(source_views$correlation_chord_valid[source_views$ph_eligible]),
    sum(!source_views$ph_eligible) == 8L,
    sum(source_views$panel_metadata_reconciled) == 16L,
    nrow(external_axis) == 8L && !any(external_axis$old_axis_equal),
    nrow(external_baseline) == 8L &&
      all(grepl("^[0-9a-f]{64}$", external_baseline$selected_cell_sha256)),
    nrow(ph_queue) == 1280L && !anyDuplicated(ph_queue$job_id),
    all(ph_queue$homology_dimensions == "H0;H1_separate"),
    nrow(backend_policy) == 3L &&
      backend_policy$rss_cap_bytes[[3L]] == 12 * 1024^3,
    nrow(comparator_bindings) == 1240L,
    nrow(landscape_contract) == 10L &&
      all(landscape_contract$owner_approved),
    nrow(landscape_queue) == 28L &&
      all(landscape_queue$grid_policy == "none") &&
      all(landscape_queue$level_cap == "none"),
    nrow(comparison_firewall) == 40L &&
      all(comparison_firewall$authorization_state ==
            "closed_pending_complete_immutable_landscapes"),
    all(runtime$current_gate_passed[runtime$required_before_stage %in%
                                      c("prefreeze", "ph_sentinel")]),
    !any(runtime$current_gate_passed[
      runtime$required_before_stage == "landscape_stress"
    ]),
    nrow(implementation) == length(implementation_paths),
    all(source_bindings$outcome_label_state == "closed") &&
      !any(source_bindings$biological_outcomes_computed),
    !any(c(source_views$authorization_state == "executed",
           ph_queue$authorization_state == "executed"))
  ),
  evidence = c(
    "129 MV8-Q plus three MV8-O primary caches",
    "124 internal and eight external fits", "1,272 cache-backed view receipts",
    "all PH-eligible geometries finite and correlation-chord valid",
    "eight external SCT-data exact500 views remain diagnostic-only",
    "16 common475 source rows rebound to the canonical common475 panel SHA",
    "all eight MV8-I axes differ from the current exact-reference axes",
    "eight selectedfit baselines require current selected384 digests",
    "1,272 gene plus eight external cell PH records",
    "H0/H1 retained separately in every PH record",
    "Ripserr primary and 12-GiB exact GUDHI resource fallback",
    "620 gene plus 620 cell internal comparator records immutable",
    "all finite intervals, all active levels, exact/error-controlled L2",
    "20 internal plus eight external streamed groups",
    "40 descriptive strata closed until all distance hashes complete",
    "R/Ripserr/TDA environment ready for sentinel PH",
    "Rust toolchain and accepted binary must be restored before landscapes",
    "11 implementation files hash-bound",
    "labels, outcomes, clustering, fusion and adoption closed",
    "builder opened no private matrix and ran no PH or landscapes"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("MV8-R independent metadata validation failed", call. = FALSE)
}

decision <- data.frame(
  contract_id = "mv08r_decision_v1",
  decision = "authorize_same_axis_external_baseline_and_23_job_PH_sentinel_after_commit",
  source_fits_reused = 132L, external_baseline_jobs_authorized = 8L,
  ph_jobs_authorized = 23L, full_ph_jobs_closed = 1257L,
  landscape_groups_authorized = 0L, comparison_strata_authorized = 0L,
  clustering_jobs_authorized = 0L, fusion_jobs_authorized = 0L,
  label_jobs_authorized = 0L, outcome_jobs_authorized = 0L,
  next_gate = "same_axis_external_baseline_plus_PH_sentinel_closure",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)

atomic_csv(contract, file.path(output_dir, "mv08r-contract.csv"))
atomic_csv(source_bindings,
           file.path(output_dir, "mv08r-source-cache-bindings.csv"))
atomic_csv(source_views,
           file.path(output_dir, "mv08r-source-gene-view-bindings.csv"))
atomic_csv(external_baseline,
           file.path(output_dir, "mv08r-external-same-axis-baseline.csv"))
atomic_csv(ph_queue, file.path(output_dir, "mv08r-ph-queue.csv"))
atomic_csv(comparator_bindings,
           file.path(output_dir, "mv08r-internal-comparator-bindings.csv"))
atomic_csv(backend_policy,
           file.path(output_dir, "mv08r-backend-policy.csv"))
atomic_csv(landscape_contract,
           file.path(output_dir, "mv08r-landscape-contract.csv"))
atomic_csv(landscape_queue,
           file.path(output_dir, "mv08r-landscape-queue.csv"))
atomic_csv(comparison_firewall,
           file.path(output_dir, "mv08r-comparison-firewall.csv"))
atomic_csv(reconciliation,
           file.path(output_dir, "mv08r-metadata-reconciliation.csv"))
atomic_csv(stage_gates, file.path(output_dir, "mv08r-stage-gates.csv"))
atomic_csv(implementation,
           file.path(output_dir, "mv08r-implementation-bindings.csv"))
atomic_csv(runtime, file.path(output_dir, "mv08r-runtime-audit.csv"))
atomic_csv(validation,
           file.path(output_dir, "mv08r-independent-validation.csv"))
atomic_csv(decision, file.path(output_dir, "mv08r-decision.csv"))

report <- c(
  "# MV8-R full topology-production prefreeze", "",
  "**Date:** 2026-08-23", "",
  "**Result:** 20/20 metadata-only checks pass; no topology was executed.", "",
  "## Corrected execution scope", "",
  paste0(
    "The 132 completed all-QC source fits are hash-bound to 1,272 produced ",
    "gene-view receipts. Eight external SCT-data exact-500 views remain ",
    "diagnostic-only. The executable topology axis contains 1,272 gene PH ",
    "records plus eight same-axis external cell records."
  ), "",
  "## Two prospective corrections", "",
  paste0(
    "The source worker recorded the exact-500 SHA on common-475 rows. MV8-R ",
    "rebinds those rows to the frozen ordered common-475 SHA without changing ",
    "their matrix or distance hashes."
  ), "",
  paste0(
    "The older MV8-I external baseline used different selected-cell digests ",
    "from all eight current raw-read exact-reference units. It is valid ",
    "historical technical evidence but is not a paired fit-scope comparator. ",
    "MV8-R therefore requires a selected-384-fit baseline reconstructed from ",
    "the same current raw-read inputs and the exact same selected-cell digest."
  ), "",
  "## PH and landscape policy", "",
  paste0(
    "PH is complete Vietoris-Rips over field 2 through H1. Ripserr is primary; ",
    "a 12-GiB TDA/GUDHI exact fallback is permitted for gene views only after ",
    "an 8-GiB Ripserr RSS stop. H0 and H1 remain separate."
  ), "",
  paste0(
    "Landscapes use all finite positive intervals, exclude essential H0, retain ",
    "every consecutive active level, and use exact or error-controlled ",
    "squared-L2 integration without a universal grid or level cap. Work is ",
    "streamed. The accepted Rust binary is currently absent, so every ",
    "landscape remains closed until the locked kernel is rebuilt, hash-bound, ",
    "and checked against the canonical R oracle."
  ), "",
  "## Authorized next gate", "",
  paste0(
    "After commit, only the eight same-axis external baseline jobs and the ",
    "23-job PH sentinel are authorized. Full PH, all 28 landscape groups, all ",
    "40 descriptive comparisons, clustering, fusion, labels, outcomes, ",
    "default adoption, and manuscript claims remain closed."
  )
)
atomic_text(report, file.path(
  output_dir, "MV08R_FULL_TOPOLOGY_PRODUCTION_PREFREEZE_2026-08-23.md"
))

artifacts <- list.files(output_dir, full.names = TRUE)
artifacts <- artifacts[
  basename(artifacts) != "mv08r-artifact-manifest.csv"
]
manifest <- data.frame(
  artifact = basename(artifacts), bytes = as.numeric(file.info(artifacts)$size),
  sha256 = vapply(artifacts, sha_file, character(1L)),
  stringsAsFactors = FALSE
)
atomic_csv(manifest,
           file.path(output_dir, "mv08r-artifact-manifest.csv"))

cat(
  "MV8-R prefreeze checks=", sum(validation$passed), "/", nrow(validation),
  "; sources=132; PH queue=1280; topology not executed\n", sep = ""
)
