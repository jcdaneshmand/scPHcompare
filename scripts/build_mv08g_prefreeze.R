#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "ripserr", "TDA", "cluster", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop("usage: build_mv08g_prefreeze.R MV08F_RECOVERY MV08F_PREFREEZE MV08E_EVIDENCE MV07H_PREFREEZE MV07H_PH_EVIDENCE MV07H_LANDSCAPE_EVIDENCE RUST_LIBRARY OUTPUT EXPECTED_HEAD")
}
recovery_dir <- args[[1L]]; recovery_prefreeze <- args[[2L]]
mv08e_dir <- args[[3L]]; mv07h_prefreeze <- args[[4L]]
mv07h_ph <- args[[5L]]; mv07h_landscape <- args[[6L]]
rust <- normalizePath(args[[7L]], winslash = "/", mustWork = TRUE)
output <- args[[8L]]; expected_head <- tolower(trimws(args[[9L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head || !grepl("^[0-9a-f]{40}$", expected_head)) {
  stop("MV8-G exact HEAD mismatch.")
}
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV8-G prefreeze output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv07h_full_topology.R")
source("R/mv08g_panel_sensitivity.R")
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
sha <- function(path) digest::digest(file = path, algo = "sha256",
                                     serialize = FALSE)

recovery <- read.csv(file.path(recovery_dir, "mv08f-recovery-decision.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
resource <- read.csv(file.path(recovery_dir,
  "mv08f-recovery-resource-summary.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
recovered <- read.csv(file.path(recovery_dir,
  "mv08f-recovered-cache-validation.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
comparator <- read.csv(file.path(recovery_prefreeze,
  "mv08f-existing-comparator-status.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
panel <- read.csv(file.path(mv08e_dir, "mv08e-common475-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
manifest <- read.csv(file.path(mv07h_prefreeze, "mv07h-cache-manifest.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
ph_decision <- read.csv(file.path(mv07h_ph, "mv07h-full-ph-decision.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
landscape_decision <- read.csv(file.path(mv07h_landscape,
  "mv07h-landscape-complete-validation-decision.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
landscape_inventory <- read.csv(file.path(mv07h_landscape,
  "mv07h-landscape-complete-group-inventory.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(recovery) != 1L ||
    recovery$decision != "recovery_exact_authorize_475_source_prefreeze" ||
    recovery$raw_exact != 90L || recovery$sct_exact != 450L ||
    recovery$unexpected_cache_files != 0L ||
    nrow(resource) != 1L || !truth(resource$all_resource_caps_passed) ||
    nrow(recovered) != 450L || !all(truth(recovered$exact_identity_passed)) ||
    nrow(comparator) != 3L || !all(truth(comparator$exact_hash_identity)) ||
    !identical(comparator$artifact_kind,
      c("source_bundle", "ph_record", "landscape_group")) ||
    !identical(as.integer(comparator$live_count), c(5L, 1240L, 20L))) {
  stop("MV8-F cache recovery is not independently closed.")
}
mv08g_validate_common475_panel_v1(panel)
mv07h_validate_cache_manifest_v1(manifest)
if (ph_decision$source_jobs != 5L || ph_decision$ph_jobs != 1240L ||
    ph_decision$typed_views != 1240L ||
    ph_decision$outcome_label_state != "closed" ||
    landscape_decision$landscape_groups_complete != 20L ||
    landscape_decision$component_rows_complete != 152520L ||
    landscape_decision$outcome_label_state != "closed" ||
    nrow(landscape_inventory) != 20L) {
  stop("Accepted 500-gene comparator ledgers are incomplete.")
}
axis <- mv08g_sample_seed_axis_v1(manifest)
source_queue <- mv08g_source_queue_v1(axis)
ph_queue <- mv08g_ph_queue_v1(axis)
landscape_queue <- mv08g_landscape_queue_v1()
shift_queue <- mv08g_matched_shift_queue_v1()
caps <- mv08g_resource_caps_v1()
rust_sha <- sha(rust)
accepted_rust_sha <-
  "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d"
if (rust_sha != accepted_rust_sha) stop("MV8-G Rust binary identity drift.")

implementation_paths <- c(
  "R/mv08g_panel_sensitivity.R", "scripts/run_mv08g_source_entry.R",
  "scripts/run_mv08g_source_monitor.R",
  "scripts/validate_mv08g_sources.R",
  "scripts/run_mv08g_ph_entry.R", "scripts/run_mv08g_ph_fallback_entry.R",
  "scripts/run_mv08g_landscape_group.R",
  "scripts/run_mv08g_matched_shift_group.R",
  "scripts/run_mv08g_comparison.R",
  "scripts/build_mv08g_prefreeze.R", "scripts/validate_mv08g_prefreeze.R",
  "tests/testthat/test-mv08g-panel-sensitivity.R",
  "docs/specifications/MV08G_COMMON475_PAIRED_REFERENCE_SENSITIVITY_PREFREEZE_V2.md")
if (any(!file.exists(implementation_paths))) stop("MV8-G implementation is incomplete.")
implementation_hashes <- vapply(implementation_paths, sha, character(1L))
implementation_root <- digest::digest(data.frame(
  path = implementation_paths, sha256 = unname(implementation_hashes),
  accepted_head = expected_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)

contract <- data.frame(
  contract_id = "mv08g_common475_paired_sensitivity_prefreeze_v2",
  accepted_head = expected_head,
  implementation_root_sha256 = implementation_root,
  rust_library_sha256 = rust_sha, samples = 124L, seeds = 5L,
  panel_genes = 475L, cells_per_sample = 384L, pca_components = 30L,
  source_jobs = 5L, source_repeat_jobs = 1L, typed_views = 1240L,
  ph_jobs = 1240L,
  landscape_groups = 20L, landscape_rows = 152520L,
  matched_shift_groups = 20L, matched_shift_rows = 2480L,
  comparison_components = 4L, candidate_k = "2:10",
  homology_dimensions = "H0;H1_separate",
  filtration = "complete_vietoris_rips", threshold = -1, field = 2L,
  landscape_definition =
    "finite_positive;essential_H0_excluded;all_active_levels;exact_squared_L2;no_grid;no_level_cap;streamed_groups",
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
comparison <- data.frame(
  contract_id = "mv08g_comparison_contract_v1",
  comparison_id = c("distance_rank", "scaled_stress", "top10_neighbors",
    "matched_landscape_shift", "pam_k_2_10", "panel_selected_k",
    "fixed_500_selected_k_pam", "fixed_500_selected_k_average"),
  metric = c("spearman_and_kendall", "nonnegative_LS_normalized_stress",
    "mean_per_sample_overlap", "shift_over_median_nonzero_500_distance",
    "adjusted_rand_index", "exact_k_agreement", "adjusted_rand_index",
    "adjusted_rand_index"),
  components = 4L, seeds = 5L,
  label_free = TRUE, primary = TRUE, stringsAsFactors = FALSE)
thresholds <- data.frame(
  contract_id = "mv08g_interpretation_threshold_v1",
  metric = c("median_spearman", "median_top10_overlap",
             "median_fixed_k_pam_ari"), threshold = c(0.95, 0.80, 0.80),
  direction = "greater_than_or_equal",
  high_rule = "all_four_components_pass_all_three",
  material_rule = "any_component_misses_at_least_two_of_three",
  equivalence_test = FALSE, stringsAsFactors = FALSE)
firewall <- data.frame(
  contract_id = "mv08g_label_firewall_v1",
  allowed_execution_inputs =
    "sample_id;seed;cache_identity;selected_cell_identity;panel_identity;transform_identity;PH_identity;distance_identity",
  forbidden_execution_inputs =
    "tissue;approach;study;class;label;outcome;prior_result",
  labels_may_select = FALSE, labels_may_stop = FALSE,
  hca_expression_jobs = 0L, hca_fastq_download_jobs = 0L,
  raw_reprocessing_jobs = 0L, label_jobs = 0L, outcome_jobs = 0L,
  label_access_state = "closed", stringsAsFactors = FALSE)
projection <- data.frame(
  contract_id = "mv08g_resource_projection_v1",
  component = c("source_and_PH", "within475_landscapes",
                "matched_shifts", "comparisons", "total"),
  projected_worker_hours = c(6.0, 1.25, 0.25, 0.5, 8.0),
  basis = c("accepted_MV7H_measured_plus_margin",
            "accepted_MV7H_landscape_measured_plus_margin",
            "row_scaled_landscape_plus_margin", "bounded_4_component_summary",
            "sum_of_conservative_components"), stringsAsFactors = FALSE)
runtime <- data.frame(
  contract_id = "mv08g_runtime_v1", r_version = as.character(getRversion()),
  ripserr_version = as.character(utils::packageVersion("ripserr")),
  tda_version = as.character(utils::packageVersion("TDA")),
  cluster_version = as.character(utils::packageVersion("cluster")),
  mclust_version = as.character(utils::packageVersion("mclust")),
  rust_library_sha256 = rust_sha, rust_engine_version = 1L,
  omp_threads = 1L, openblas_threads = 1L, mkl_threads = 1L,
  stringsAsFactors = FALSE)
source_files <- c(
  recovery_decision = file.path(recovery_dir, "mv08f-recovery-decision.csv"),
  recovery_validation = file.path(recovery_dir,
    "mv08f-recovered-cache-validation.csv"),
  accepted_500_live_identity = file.path(recovery_prefreeze,
    "mv08f-existing-comparator-status.csv"),
  common475_panel = file.path(mv08e_dir, "mv08e-common475-panel.csv"),
  accepted_cache_manifest = file.path(mv07h_prefreeze,
    "mv07h-cache-manifest.csv"),
  accepted_500_ph_decision = file.path(mv07h_ph, "mv07h-full-ph-decision.csv"),
  accepted_500_landscape_decision = file.path(mv07h_landscape,
    "mv07h-landscape-complete-validation-decision.csv"),
  accepted_500_landscape_inventory = file.path(mv07h_landscape,
    "mv07h-landscape-complete-group-inventory.csv"),
  implementation_paths)
freeze <- data.frame(
  contract_id = "mv08g_source_freeze_v2", source_id = names(source_files),
  artifact_locator = unname(source_files),
  sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size), accepted_head = expected_head,
  private_source = FALSE, stringsAsFactors = FALSE)
acceptance <- data.frame(
  contract_id = "mv08g_prefreeze_acceptance_v2",
  category = c("recovery", "panel", "cache_axis", "source_queue", "ph_queue",
    "landscape_queue", "matched_shift_queue", "accepted_500", "definition",
    "comparison", "resources", "rust", "label_firewall", "stop_boundary"),
  passed = c(
    nrow(recovered) == 450L && all(truth(recovered$exact_identity_passed)),
    nrow(panel) == 475L, nrow(axis) == 620L && all(table(axis$seed) == 124L),
    nrow(source_queue) == 5L && sum(source_queue$typed_view_count) == 1240L,
    nrow(ph_queue) == 1240L && all(table(ph_queue$view_id) == 620L),
    nrow(landscape_queue) == 20L &&
      sum(landscape_queue$component_rows) == 152520L,
    nrow(shift_queue) == 20L && sum(shift_queue$component_rows) == 2480L,
    all(truth(comparator$exact_hash_identity)) && ph_decision$ph_jobs == 1240L &&
      landscape_decision$component_rows_complete == 152520L,
    grepl("all_active_levels", contract$landscape_definition, fixed = TRUE) &&
      grepl("no_grid", contract$landscape_definition, fixed = TRUE),
    nrow(comparison) == 8L && nrow(thresholds) == 3L,
    projection$projected_worker_hours[projection$component == "total"] <= 12,
    rust_sha == accepted_rust_sha,
    firewall$label_access_state == "closed" &&
      sum(firewall[c("hca_expression_jobs", "hca_fastq_download_jobs",
        "raw_reprocessing_jobs", "label_jobs", "outcome_jobs")]) == 0,
    TRUE),
  detail = c("90 raw and 450 SCT exact", "ordered exact common475_v1",
    "124 biological samples by five seeds", "five global refits",
    "620 cell plus 620 gene PH", "20 groups and 152,520 rows",
    "20 groups and 2,480 matched rows", "immutable comparator complete",
    "dissertation-aligned exact landscape", "eight frozen comparison families",
    "conservative projection at most 12 worker-hours", "accepted Rust hash",
    "HCA and biological labels closed", "prefreeze authorizes source only"),
  stringsAsFactors = FALSE)
if (!all(acceptance$passed)) stop("MV8-G prefreeze acceptance failed: ",
  paste(acceptance$category[!acceptance$passed], collapse = ", "))
decision <- data.frame(
  contract_id = "mv08g_prefreeze_decision_v2",
  decision = "authorize_five_common475_source_bundles_and_one_repeat",
  source_jobs_authorized = 5L, source_repeat_jobs_authorized = 1L,
  ph_jobs_authorized = 0L,
  landscape_groups_authorized = 0L, matched_shift_groups_authorized = 0L,
  comparison_jobs_authorized = 0L, hca_fastq_download_authorized = FALSE,
  raw_reprocessing_authorized = FALSE, label_access_authorized = FALSE,
  next_gate = "independent_source_identity_and_repeat_validation",
  stringsAsFactors = FALSE)

write_provenance_csv(panel, file.path(output, "mv08g-panel.csv"))
write_provenance_csv(manifest, file.path(output, "mv08g-cache-manifest.csv"))
write_provenance_csv(landscape_inventory,
  file.path(output, "mv08g-accepted500-landscape-inventory.csv"))
outputs <- list(
  "mv08g-contract.csv" = contract,
  "mv08g-sample-seed-axis.csv" = axis,
  "mv08g-source-queue.csv" = source_queue,
  "mv08g-ph-queue.csv" = ph_queue,
  "mv08g-landscape-queue.csv" = landscape_queue,
  "mv08g-matched-shift-queue.csv" = shift_queue,
  "mv08g-resource-caps.csv" = caps,
  "mv08g-resource-projection.csv" = projection,
  "mv08g-comparison-contract.csv" = comparison,
  "mv08g-interpretation-thresholds.csv" = thresholds,
  "mv08g-label-firewall.csv" = firewall,
  "mv08g-runtime.csv" = runtime,
  "mv08g-source-freeze.csv" = freeze,
  "mv08g-acceptance.csv" = acceptance,
  "mv08g-decision.csv" = decision)
for (name in names(outputs)) write_provenance_csv(
  outputs[[name]], file.path(output, name))
message("MV8-G prefreeze complete: source-only gate; HCA raw reads remain closed")
