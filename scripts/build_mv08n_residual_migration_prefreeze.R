#!/usr/bin/env Rscript

# Build the label-closed MV8-N migration prefreeze.  This script hashes and
# inventories inputs only; it never opens expression matrices or runs SCT/PH.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: build_mv08n_residual_migration_prefreeze.R <primary-raw-root> <added-raw-root> <hca-bm002-filtered.h5> <hca-remaining-root> <output-dir>",
    call. = FALSE
  )
}

primary_root <- normalizePath(args[[1L]], mustWork = TRUE)
added_root <- normalizePath(args[[2L]], mustWork = TRUE)
bm002_path <- normalizePath(args[[3L]], mustWork = TRUE)
hca_root <- normalizePath(args[[4L]], mustWork = TRUE)
output_dir <- normalizePath(args[[5L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
sha_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
safe_id <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

head_commit <- tolower(trimws(Sys.getenv("MV08N_GIT_HEAD", unset = "")))
if (!nzchar(head_commit)) {
  head_commit <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
}
if (!grepl("^[0-9a-f]{40}$", head_commit)) stop("cannot bind git HEAD", call. = FALSE)

primary <- read_csv("docs/audits/mv08f-cache-recovery-evidence/mv08f-raw-recovery.csv")
added <- read_csv("docs/audits/mv07f-production-evidence/mv07f-raw-production.csv")
selection <- read_csv("docs/audits/mv07fp-prefreeze-evidence-v4/mv07fp-cache-manifest.csv")
exact <- read_csv("docs/audits/mv07h-prefreeze-evidence-v4/mv07h-panel.csv")
common <- read_csv("docs/audits/mv08e-reference-reconciliation-evidence/mv08e-common475-panel.csv")
remaining_summary <- read_csv("docs/audits/mv08h-exact500-remaining-v2/mv08h-exact500-unit-summary.csv")
remaining_binding <- read_csv("docs/audits/mv08h-exact500-remaining-v2/mv08h-exact500-remaining-input-binding.csv")
mv08m_identity <- read_csv("docs/audits/mv08m-exact500-gene-representation-v1/mv08m-identity.csv")
mv08k_identity <- read_csv("docs/audits/mv08k-exact500-transform-contract-v1/mv08k-exact500-transform-identity.csv")

if (nrow(primary) != 90L || nrow(added) != 34L || nrow(selection) != 620L ||
    nrow(exact) != 500L || nrow(common) != 475L || nrow(remaining_summary) != 7L ||
    nrow(remaining_binding) != 7L || nrow(mv08m_identity) != 1L || nrow(mv08k_identity) != 1L) {
  stop("frozen input cardinality drift", call. = FALSE)
}
if (anyDuplicated(c(primary$sample_id, added$sample_id)) ||
    !setequal(unique(selection$sample_id), c(primary$sample_id, added$sample_id)) ||
    !all(selection$selected_cells == 384L) ||
    !identical(sort(unique(selection$seed)), 20260805:20260809) ||
    any(table(selection$sample_id) != 5L) ||
    any(selection$outcome_label_state != "closed") ||
    any(as.logical(selection$biological_outcomes_computed))) {
  stop("internal sample/selection axis drift", call. = FALSE)
}

exact <- exact[order(exact$panel_order), , drop = FALSE]
common <- common[order(common$panel_order_475), , drop = FALSE]
exact_sha <- unique(tolower(exact$panel_sha256))
common_sha <- unique(tolower(common$common475_axis_sha256))
if (length(exact_sha) != 1L || length(common_sha) != 1L ||
    exact_sha != "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e" ||
    common_sha != "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba" ||
    !all(common$feature_id %in% exact$feature_id)) {
  stop("panel identity drift", call. = FALSE)
}

primary_axis <- data.frame(
  source_tier = "primary90", sample_id = primary$sample_id,
  qc_cells = as.integer(primary$cells), raw_genes = as.integer(primary$genes),
  counts_sha256 = tolower(primary$counts_sha256),
  raw_file = primary$private_raw_cache_file,
  expected_bytes = as.numeric(primary$private_raw_cache_size_bytes),
  expected_sha256 = tolower(primary$private_raw_cache_sha256), stringsAsFactors = FALSE
)
added_axis <- data.frame(
  source_tier = "added34", sample_id = added$sample_id,
  qc_cells = as.integer(added$observed_cells), raw_genes = as.integer(added$observed_genes),
  counts_sha256 = tolower(added$counts_sha256), raw_file = added$private_cache_file,
  expected_bytes = as.numeric(added$private_cache_bytes),
  expected_sha256 = tolower(added$private_cache_sha256), stringsAsFactors = FALSE
)
internal_private <- rbind(primary_axis, added_axis)
internal_private$root <- ifelse(internal_private$source_tier == "primary90", primary_root, added_root)
internal_private$path <- file.path(internal_private$root, internal_private$raw_file)
internal_private$exists <- file.exists(internal_private$path)
internal_private$live_bytes <- ifelse(internal_private$exists, file.info(internal_private$path)$size, NA_real_)
internal_private$live_sha256 <- vapply(internal_private$path, function(path) {
  if (!file.exists(path)) return(NA_character_)
  sha_file(path)
}, character(1L))
internal_private$identity_passed <- internal_private$exists &
  internal_private$live_bytes == internal_private$expected_bytes &
  internal_private$live_sha256 == internal_private$expected_sha256
if (!all(internal_private$identity_passed)) stop("internal raw cache identity drift", call. = FALSE)

selection_groups <- split(selection, selection$sample_id)
selection_summary <- do.call(rbind, lapply(names(selection_groups), function(id) {
  x <- selection_groups[[id]][order(selection_groups[[id]]$seed), , drop = FALSE]
  data.frame(
    sample_id = id, selection_axes = nrow(x), selected_cells_per_axis = unique(x$selected_cells),
    selection_axis_sha256 = sha_object(data.frame(
      seed = as.integer(x$seed), selected_cell_sha256 = as.character(x$selected_cell_sha256),
      stringsAsFactors = FALSE
    )), stringsAsFactors = FALSE
  )
}))
internal_private <- merge(internal_private, selection_summary, by = "sample_id", sort = FALSE)
internal_private <- internal_private[order(internal_private$source_tier, internal_private$sample_id), , drop = FALSE]

internal_axis <- data.frame(
  contract_id = "mv08n_internal_source_axis_v1", source_tier = internal_private$source_tier,
  sample_id = internal_private$sample_id, qc_cells = internal_private$qc_cells,
  raw_genes = internal_private$raw_genes, counts_sha256 = internal_private$counts_sha256,
  raw_cache_bytes = internal_private$live_bytes, raw_cache_sha256 = internal_private$live_sha256,
  selection_axes = internal_private$selection_axes,
  selected_cells_per_axis = internal_private$selected_cells_per_axis,
  selection_axis_sha256 = internal_private$selection_axis_sha256,
  normalization_fit_scope = "all_frozen_qc_cells_per_sample",
  topology_observation_scope = "frozen_selected384_per_seed",
  exact500_panel_sha256 = exact_sha, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
atomic_csv(internal_axis, file.path(output_dir, "mv08n-internal-source-axis.csv"))

remaining_units <- sort(as.character(remaining_summary$unit_id))
remaining_paths <- file.path(
  hca_root, remaining_units, paste0("mv08h_exact500_", tolower(remaining_units)),
  "outs", "filtered_feature_bc_matrix.h5"
)
external_paths <- c(HCA_BM_002 = bm002_path, stats::setNames(remaining_paths, remaining_units))
external_paths <- external_paths[order(names(external_paths))]
if (length(external_paths) != 8L || any(!file.exists(external_paths))) {
  stop("external exact-reference H5 inventory incomplete", call. = FALSE)
}
external_bytes <- file.info(external_paths)$size
external_hashes <- vapply(external_paths, sha_file, character(1L))
qc_cells <- c(
  HCA_BM_002 = as.integer(mv08m_identity$eligible_cells),
  stats::setNames(as.integer(remaining_summary$post_qc_cells), remaining_summary$unit_id)
)[names(external_paths)]
selected_sha <- rep(NA_character_, 8L); names(selected_sha) <- names(external_paths)
selected_sha[["HCA_BM_002"]] <- as.character(mv08k_identity$selected_cell_sha256)
external_axis <- data.frame(
  contract_id = "mv08n_external_source_axis_v1", unit_id = names(external_paths),
  filtered_h5_bytes = as.numeric(external_bytes), filtered_h5_sha256 = external_hashes,
  qc_eligible_cells = as.integer(qc_cells), selection_seed = 20260805L,
  selected_cells = 384L, selected_cell_sha256 = unname(selected_sha),
  selected_axis_state = ifelse(names(external_paths) == "HCA_BM_002", "frozen_mv08k", "freeze_in_source_preflight"),
  exact500_present = 500L, common475_present = 475L,
  exact500_panel_sha256 = exact_sha, common475_panel_sha256 = common_sha,
  normalization_fit_scope = "all_frozen_qc_cells_per_unit",
  topology_observation_scope = "frozen_selected384_seed20260805",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(external_axis, file.path(output_dir, "mv08n-external-source-axis.csv"))

min_id <- internal_axis$sample_id[order(internal_axis$qc_cells, internal_axis$sample_id)][[1L]]
max_id <- internal_axis$sample_id[order(-internal_axis$qc_cells, internal_axis$sample_id)][[1L]]
source_queue <- rbind(
  data.frame(
    dataset_scope = "internal124", source_tier = internal_axis$source_tier,
    unit_id = internal_axis$sample_id, fit_cells = internal_axis$qc_cells,
    selected_axes = 5L, selected_cells_per_axis = 384L,
    source_sha256 = internal_axis$raw_cache_sha256, stringsAsFactors = FALSE
  ),
  data.frame(
    dataset_scope = "external8", source_tier = "hca_exact_reference",
    unit_id = external_axis$unit_id, fit_cells = external_axis$qc_eligible_cells,
    selected_axes = 1L, selected_cells_per_axis = 384L,
    source_sha256 = external_axis$filtered_h5_sha256, stringsAsFactors = FALSE
  )
)
source_queue <- source_queue[order(source_queue$dataset_scope, source_queue$unit_id), , drop = FALSE]
source_queue$job_order <- seq_len(nrow(source_queue))
source_queue$sentinel_role <- ifelse(
  source_queue$unit_id == min_id, "internal_minimum_cell_sentinel",
  ifelse(source_queue$unit_id == max_id, "internal_maximum_cell_sentinel",
         ifelse(source_queue$unit_id == "HCA_BM_002", "external_mv08m_bridge_sentinel", "none"))
)
source_queue$authorization_state <- ifelse(source_queue$sentinel_role == "none",
                                           "closed_pending_sentinel", "source_view_sentinel_authorized")
source_queue$repeat_required <- source_queue$unit_id %in% c(max_id, "HCA_BM_002")
source_queue$output_file <- paste0(
  ifelse(source_queue$dataset_scope == "internal124", "internal/", "external/"),
  safe_id(source_queue$unit_id), "__exact500_allqc_sct_model.rds"
)
source_queue$elapsed_cap_seconds <- 1800L
source_queue$rss_cap_bytes <- 12 * 1024^3
source_queue$workers <- 1L; source_queue$retries <- 0L
source_queue$outcome_label_state <- "closed"; source_queue$biological_outcomes_computed <- FALSE
source_queue <- source_queue[, c(
  "job_order", "dataset_scope", "source_tier", "unit_id", "fit_cells", "selected_axes",
  "selected_cells_per_axis", "source_sha256", "sentinel_role", "authorization_state",
  "repeat_required", "output_file", "elapsed_cap_seconds", "rss_cap_bytes", "workers",
  "retries", "outcome_label_state", "biological_outcomes_computed"
)]
source_queue$contract_id <- "mv08n_residual_source_queue_v1"
source_queue <- source_queue[, c("contract_id", setdiff(names(source_queue), "contract_id"))]
atomic_csv(source_queue, file.path(output_dir, "mv08n-residual-source-queue.csv"))

internal_views <- do.call(rbind, lapply(seq_len(nrow(selection)), function(i) {
  row <- selection[i, , drop = FALSE]
  data.frame(
    dataset_scope = "internal124", unit_id = row$sample_id, seed = as.integer(row$seed),
    selected_cell_sha256 = row$selected_cell_sha256,
    representation_id = c("sct_data_all_qc_fit_selected384", "sct_pearson_residual_all_qc_fit_selected384"),
    panel_id = "exact500", panel_genes = 500L, panel_sha256 = exact_sha,
    view_role = c("fit_scope_control", "migration_candidate"), ph_planned = TRUE,
    stringsAsFactors = FALSE
  )
}))
external_representation <- data.frame(
  representation_id = c(
    "sct_data_selected384_fit", "sct_data_all_qc_fit_selected384",
    "sct_pearson_residual_all_qc_fit_selected384",
    "sct_data_all_qc_fit_selected384", "sct_pearson_residual_all_qc_fit_selected384"
  ),
  panel_id = c("common475", "common475", "common475", "exact500", "exact500"),
  panel_genes = c(475L, 475L, 475L, 500L, 500L),
  view_role = c("current_h5_baseline", "fit_scope_control", "layer_candidate",
                "exact500_feasibility_diagnostic", "migration_candidate"),
  ph_planned = c(TRUE, TRUE, TRUE, FALSE, TRUE), stringsAsFactors = FALSE
)
external_views <- do.call(rbind, lapply(seq_len(nrow(external_axis)), function(i) {
  row <- external_axis[i, , drop = FALSE]
  x <- external_representation
  data.frame(
    dataset_scope = "external8", unit_id = row$unit_id, seed = 20260805L,
    selected_cell_sha256 = row$selected_cell_sha256,
    representation_id = x$representation_id, panel_id = x$panel_id,
    panel_genes = x$panel_genes,
    panel_sha256 = ifelse(x$panel_id == "exact500", exact_sha, common_sha),
    view_role = x$view_role, ph_planned = x$ph_planned, stringsAsFactors = FALSE
  )
}))
view_queue <- rbind(internal_views, external_views)
view_queue$job_order <- seq_len(nrow(view_queue))
view_queue$view_id <- paste(
  "mv08n_gene_topology_v1", view_queue$dataset_scope, view_queue$unit_id,
  view_queue$seed, view_queue$representation_id, view_queue$panel_id, sep = ":"
)
view_queue$selected_cells <- 384L
view_queue$metric <- "pearson_correlation_chord"
view_queue$exact_panel_required <- TRUE
view_queue$finite_required <- TRUE
view_queue$nonzero_variance_required <- TRUE
view_queue$ph_authorized <- FALSE
view_queue$outcome_label_state <- "closed"; view_queue$biological_outcomes_computed <- FALSE
view_queue$contract_id <- "mv08n_gene_view_queue_v1"
view_queue <- view_queue[, c(
  "contract_id", "job_order", "view_id", "dataset_scope", "unit_id", "seed",
  "selected_cells", "selected_cell_sha256", "representation_id", "panel_id", "panel_genes",
  "panel_sha256", "view_role", "metric", "exact_panel_required", "finite_required",
  "nonzero_variance_required", "ph_planned", "ph_authorized", "outcome_label_state",
  "biological_outcomes_computed"
)]
atomic_csv(view_queue, file.path(output_dir, "mv08n-gene-view-queue.csv"))

ph_queue <- view_queue[view_queue$ph_planned, c(
  "dataset_scope", "unit_id", "seed", "representation_id", "panel_id", "panel_genes"
), drop = FALSE]
ph_queue$job_order <- seq_len(nrow(ph_queue))
ph_queue$job_id <- paste("mv08n_ph", ph_queue$dataset_scope, ph_queue$unit_id,
                         ph_queue$seed, ph_queue$representation_id, ph_queue$panel_id, sep = ":")
ph_queue$homology_dimensions <- "H0;H1_separate"
ph_queue$filtration <- "complete_vietoris_rips"
ph_queue$field <- 2L
ph_queue$engine_policy <- "Ripserr_primary;GUDHI_exact_resource_fallback"
ph_queue$elapsed_cap_seconds <- 1800L; ph_queue$rss_cap_bytes <- 12 * 1024^3
ph_queue$workers <- 1L; ph_queue$retries <- 0L; ph_queue$authorization_state <- "closed"
ph_queue$outcome_label_state <- "closed"; ph_queue$biological_outcomes_computed <- FALSE
ph_queue$contract_id <- "mv08n_ph_queue_v1"
ph_queue <- ph_queue[, c("contract_id", "job_order", "job_id", setdiff(names(ph_queue), c("contract_id", "job_order", "job_id")))]
atomic_csv(ph_queue, file.path(output_dir, "mv08n-ph-queue.csv"))

dims <- c("H0", "H1")
internal_landscape <- expand.grid(
  representation_id = c("sct_data_all_qc_fit_selected384", "sct_pearson_residual_all_qc_fit_selected384"),
  seed = 20260805:20260809, homology_dimension = dims, stringsAsFactors = FALSE
)
internal_landscape$dataset_scope <- "internal124"; internal_landscape$panel_id <- "exact500"
internal_landscape$units <- 124L; internal_landscape$unordered_pairs <- choose(124L, 2L)
external_landscape <- expand.grid(
  stack_id = c("selectedfit_data_common475", "allqc_data_common475",
               "allqc_residual_common475", "allqc_residual_exact500"),
  homology_dimension = dims, stringsAsFactors = FALSE
)
external_landscape$representation_id <- ifelse(
  external_landscape$stack_id == "selectedfit_data_common475", "sct_data_selected384_fit",
  ifelse(external_landscape$stack_id == "allqc_data_common475", "sct_data_all_qc_fit_selected384",
         "sct_pearson_residual_all_qc_fit_selected384")
)
external_landscape$seed <- 20260805L; external_landscape$dataset_scope <- "external8"
external_landscape$panel_id <- ifelse(grepl("exact500$", external_landscape$stack_id), "exact500", "common475")
external_landscape$units <- 8L; external_landscape$unordered_pairs <- choose(8L, 2L)
external_landscape$stack_id <- NULL
landscape_queue <- rbind(
  internal_landscape[, c("dataset_scope", "representation_id", "panel_id", "seed", "homology_dimension", "units", "unordered_pairs")],
  external_landscape[, c("dataset_scope", "representation_id", "panel_id", "seed", "homology_dimension", "units", "unordered_pairs")]
)
landscape_queue$group_order <- seq_len(nrow(landscape_queue))
landscape_queue$group_id <- paste(
  "mv08n_landscape", landscape_queue$dataset_scope, landscape_queue$representation_id,
  landscape_queue$panel_id, landscape_queue$seed, landscape_queue$homology_dimension, sep = ":"
)
landscape_queue$integration <- "exact_critical_breakpoint_squared_L2"
landscape_queue$level_policy <- "all_active_consecutive_levels"
landscape_queue$grid_policy <- "none"; landscape_queue$level_cap <- "none"
landscape_queue$elapsed_cap_seconds <- 3600L; landscape_queue$rss_cap_bytes <- 12 * 1024^3
landscape_queue$workers <- 1L; landscape_queue$retries <- 0L
landscape_queue$authorization_state <- "closed"; landscape_queue$outcome_label_state <- "closed"
landscape_queue$biological_outcomes_computed <- FALSE; landscape_queue$contract_id <- "mv08n_landscape_queue_v1"
landscape_queue <- landscape_queue[, c("contract_id", "group_order", "group_id", setdiff(names(landscape_queue), c("contract_id", "group_order", "group_id")))]
atomic_csv(landscape_queue, file.path(output_dir, "mv08n-landscape-queue.csv"))

internal_contrasts <- data.frame(
  contrast_id = c("fit_scope_effect", "layer_effect", "net_migration_effect"),
  left_stack = c("existing_selectedfit_data_exact500", "allqc_data_exact500", "existing_selectedfit_data_exact500"),
  right_stack = c("allqc_data_exact500", "allqc_residual_exact500", "allqc_residual_exact500"),
  stringsAsFactors = FALSE
)
internal_comparison <- merge(
  internal_contrasts,
  expand.grid(seed = 20260805:20260809, homology_dimension = dims, stringsAsFactors = FALSE),
  all = TRUE
)
internal_comparison$dataset_scope <- "internal124"; internal_comparison$unordered_pairs <- choose(124L, 2L)
external_contrasts <- data.frame(
  contrast_id = c("fit_scope_effect", "layer_effect", "common475_net_migration_effect",
                  "residual_panel_effect", "exact500_total_migration_effect"),
  left_stack = c("selectedfit_data_common475", "allqc_data_common475", "selectedfit_data_common475",
                 "allqc_residual_common475", "selectedfit_data_common475"),
  right_stack = c("allqc_data_common475", "allqc_residual_common475", "allqc_residual_common475",
                  "allqc_residual_exact500", "allqc_residual_exact500"),
  stringsAsFactors = FALSE
)
external_comparison <- merge(
  external_contrasts,
  data.frame(seed = 20260805L, homology_dimension = dims, stringsAsFactors = FALSE),
  all = TRUE
)
external_comparison$dataset_scope <- "external8"; external_comparison$unordered_pairs <- choose(8L, 2L)
comparison <- rbind(internal_comparison, external_comparison)
comparison$comparison_order <- seq_len(nrow(comparison))
comparison$metrics <- "pearson;spearman;median_abs_change;p95_abs_change;top10_or_all7_neighbor_overlap;hashes"
comparison$interpretation <- "descriptive_no_equivalence_threshold"
comparison$authorization_state <- "closed"; comparison$outcome_label_state <- "closed"
comparison$biological_outcomes_computed <- FALSE; comparison$contract_id <- "mv08n_comparison_contract_v1"
comparison <- comparison[, c("contract_id", "comparison_order", "dataset_scope", "contrast_id",
                             "seed", "homology_dimension", "left_stack", "right_stack",
                             "unordered_pairs", "metrics", "interpretation", "authorization_state",
                             "outcome_label_state", "biological_outcomes_computed")]
atomic_csv(comparison, file.path(output_dir, "mv08n-comparison-contract.csv"))

stage_gates <- data.frame(
  stage_order = 1:6,
  stage_id = c("input_closure", "source_view_sentinel", "full_source_production",
               "ph_sentinel_and_full_ph", "landscape_stress_and_completion", "paired_comparison"),
  authorization_state = c("complete_in_prefreeze", "authorized_after_commit", rep("closed", 4L)),
  pass_requirement = c(
    "124 raw + 620 selections + 8 H5 + panels + comparator roots exact",
    "internal min/max five-seed views + HCA_BM_002; max/HCA repeat; no PH",
    "separate sentinel closure and storage projection at or below 40 GiB",
    "separate exact-engine/cross-engine/repeat/resource prefreeze",
    "maximum interval-burden stress + R oracle + repeat + immutable resume",
    "all sources/PH/landscapes complete and immutable"
  ),
  next_stage_automatic = FALSE, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, contract_id = "mv08n_stage_gate_v1",
  stringsAsFactors = FALSE
)
stage_gates <- stage_gates[, c("contract_id", "stage_order", "stage_id", "authorization_state",
                               "pass_requirement", "next_stage_automatic", "outcome_label_state",
                               "biological_outcomes_computed")]
atomic_csv(stage_gates, file.path(output_dir, "mv08n-stage-gates.csv"))

source_files <- c(
  "docs/audits/mv08f-cache-recovery-evidence/mv08f-raw-recovery.csv",
  "docs/audits/mv07f-production-evidence/mv07f-raw-production.csv",
  "docs/audits/mv07fp-prefreeze-evidence-v4/mv07fp-cache-manifest.csv",
  "docs/audits/mv07h-prefreeze-evidence-v4/mv07h-panel.csv",
  "docs/audits/mv08e-reference-reconciliation-evidence/mv08e-common475-panel.csv",
  "docs/audits/mv08h-exact500-remaining-v2/mv08h-exact500-unit-summary.csv",
  "docs/audits/mv08h-exact500-remaining-v2/mv08h-exact500-remaining-input-binding.csv",
  "docs/audits/mv08m-exact500-gene-representation-v1/mv08m-identity.csv",
  "docs/audits/mv08k-exact500-transform-contract-v1/mv08k-exact500-transform-identity.csv",
  "docs/audits/mv07h-full-ph-evidence/mv07h-independent-validation.csv",
  "docs/audits/mv07h-landscape-complete-validation/mv07h-landscape-complete-independent-validation.csv",
  "docs/audits/mv08i-hca-validation-v1/mv08i-independent-validation.csv"
)
source_freeze <- data.frame(
  contract_id = "mv08n_source_freeze_v1", artifact = source_files,
  bytes = file.info(source_files)$size,
  sha256 = vapply(source_files, sha_file, character(1L)),
  role = c("primary90_raw_identity", "added34_raw_identity", "620_selection_axes",
           "exact500_panel", "common475_panel", "remaining7_hca_qc", "remaining7_input_binding",
           "mv08m_representation_evidence", "mv08k_selected_axis_bridge",
           "internal_existing_ph_validation", "internal_existing_landscape_validation",
           "historical_external_common475_validation"),
  private_content = FALSE, stringsAsFactors = FALSE
)
atomic_csv(source_freeze, file.path(output_dir, "mv08n-source-freeze.csv"))

contract <- data.frame(
  contract_id = "mv08n_pearson_residual_migration_prefreeze_v1",
  accepted_head = head_commit, internal_samples = 124L, internal_seeds = 5L,
  internal_selection_axes = 620L, external_units = 8L, selected_cells_per_view = 384L,
  all_qc_model_fits = 132L, proposed_new_gene_views = nrow(view_queue),
  proposed_ph_records = nrow(ph_queue), proposed_landscape_groups = nrow(landscape_queue),
  paired_comparison_strata = nrow(comparison), exact500_panel_sha256 = exact_sha,
  common475_panel_sha256 = common_sha,
  candidate_representation = "sct_pearson_residual_all_qc_fit_selected384",
  topology_metric = "pearson_correlation_chord",
  landscape_definition = "finite_positive;essential_H0_excluded;all_active_levels;exact_squared_L2;no_grid;no_level_cap",
  cell_topology_state = "immutable_reuse_no_recompute",
  default_adoption_state = "candidate_approved_not_default",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
atomic_csv(contract, file.path(output_dir, "mv08n-contract.csv"))

validation <- data.frame(
  check_id = c(
    "internal_raw_count", "internal_raw_hashes", "internal_tiers", "selection_axis",
    "panel_identity", "external_h5_count", "external_h5_hashes", "external_qc_minimum",
    "source_queue", "sentinel_scope", "gene_view_queue", "ph_queue",
    "landscape_queue", "comparison_contract", "stage_firewall", "cell_topology_immutable",
    "label_outcome_firewall", "public_privacy"
  ),
  passed = c(
    nrow(internal_axis) == 124L, all(internal_private$identity_passed),
    sum(internal_axis$source_tier == "primary90") == 90L && sum(internal_axis$source_tier == "added34") == 34L,
    nrow(selection) == 620L && all(internal_axis$selection_axes == 5L) &&
      all(internal_axis$selected_cells_per_axis == 384L),
    exact_sha == "48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e" &&
      common_sha == "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba",
    nrow(external_axis) == 8L, all(nchar(external_axis$filtered_h5_sha256) == 64L),
    min(external_axis$qc_eligible_cells) >= 384L,
    nrow(source_queue) == 132L && sum(source_queue$authorization_state == "source_view_sentinel_authorized") == 3L,
    setequal(source_queue$sentinel_role[source_queue$sentinel_role != "none"],
             c("internal_minimum_cell_sentinel", "internal_maximum_cell_sentinel", "external_mv08m_bridge_sentinel")),
    nrow(view_queue) == 1280L && sum(view_queue$dataset_scope == "internal124") == 1240L &&
      sum(view_queue$dataset_scope == "external8") == 40L,
    nrow(ph_queue) == 1272L && all(ph_queue$authorization_state == "closed"),
    nrow(landscape_queue) == 28L && all(landscape_queue$authorization_state == "closed"),
    nrow(comparison) == 40L && all(comparison$authorization_state == "closed"),
    stage_gates$authorization_state[[2L]] == "authorized_after_commit" &&
      all(stage_gates$authorization_state[3:6] == "closed") && !any(stage_gates$next_stage_automatic),
    contract$cell_topology_state == "immutable_reuse_no_recompute",
    all(internal_axis$outcome_label_state == "closed") && all(external_axis$outcome_label_state == "closed") &&
      !any(internal_axis$biological_outcomes_computed) && !any(external_axis$biological_outcomes_computed),
    !any(grepl("[/\\\\]tmp[/\\\\]|private", unlist(lapply(list(internal_axis, external_axis), names)), ignore.case = TRUE))
  ),
  evidence = c(
    "124 hash-bound raw-count shards", "all live bytes and SHA-256 match accepted audits",
    "90 primary + 34 added", "620 frozen sample-seed axes; 384 cells each",
    "ordered exact500 and common475 hashes match", "8 exact-reference filtered H5 files",
    "all eight live H5 SHA-256 identities frozen", paste0("minimum external QC cells=", min(external_axis$qc_eligible_cells)),
    "132 source fits; only three sentinel units authorized", paste0("internal min=", min_id, "; internal max=", max_id, "; external=HCA_BM_002"),
    "1,240 internal + 40 external proposed views", "1,272 planned PH records remain closed",
    "20 internal + 8 external landscape groups remain closed", "30 internal + 10 external paired strata remain closed",
    "no automatic stage progression", "cell topology and PCA are immutable comparators",
    "labels/outcomes closed", "public axes contain identities and aggregate metadata only"
  ), stringsAsFactors = FALSE
)
atomic_csv(validation, file.path(output_dir, "mv08n-independent-validation.csv"))

decision <- data.frame(
  contract_id = "mv08n_decision_v1",
  decision = if (all(validation$passed)) "prefreeze_pass_source_view_sentinel_authorized_after_commit" else "prefreeze_fail",
  candidate_representation = "sct_pearson_residual_all_qc_fit_selected384",
  default_adopted = FALSE, source_view_sentinel_authorized = all(validation$passed),
  full_source_production_authorized = FALSE, ph_authorized = FALSE,
  landscapes_authorized = FALSE, comparisons_authorized = FALSE,
  clustering_authorized = FALSE, fusion_authorized = FALSE,
  labels_authorized = FALSE, outcomes_authorized = FALSE,
  next_gate = "MV8-O_internal_min_max_plus_HCA_BM_002_source_view_sentinel",
  stringsAsFactors = FALSE
)
atomic_csv(decision, file.path(output_dir, "mv08n-decision.csv"))

cat("MV8-N prefreeze checks=", sum(validation$passed), "/", nrow(validation),
    "; internal=", nrow(internal_axis), "; external=", nrow(external_axis),
    "; views=", nrow(view_queue), "; ph=", nrow(ph_queue), "\n", sep = "")
quit(status = if (all(validation$passed)) 0L else 2L)
