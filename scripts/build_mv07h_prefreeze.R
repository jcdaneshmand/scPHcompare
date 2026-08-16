#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "ripserr", "TDA")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: build_mv07h_prefreeze.R MV07G_PREFREEZE MV07G_EVIDENCE RUST_LIBRARY OUTPUT EXPECTED_HEAD")
}
mv07g_prefreeze <- args[[1L]]
mv07g_evidence <- args[[2L]]
rust <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
output <- args[[4L]]
expected_head <- tolower(trimws(args[[5L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head || !grepl("^[0-9a-f]{40}$", expected_head)) {
  stop("MV7-H exact HEAD mismatch.")
}
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-H prefreeze output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
source("R/mv07h_full_topology.R")
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")

panel <- read.csv(file.path(mv07g_prefreeze, "mv07g-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
manifest <- read.csv(file.path(mv07g_prefreeze, "mv07g-cache-manifest.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
sentinel_axis <- read.csv(file.path(mv07g_prefreeze,
  "mv07g-sentinel-axis.csv"), stringsAsFactors = FALSE, check.names = FALSE)
landscape_contract <- read.csv(file.path(mv07g_prefreeze,
  "mv07g-landscape-contract.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
mv07g_decision <- read.csv(file.path(mv07g_evidence,
  "mv07g-validation-decision.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
mv07g_checks <- read.csv(file.path(mv07g_evidence,
  "mv07g-independent-validation.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
mv07g_ph <- read.csv(file.path(mv07g_evidence, "mv07g-ph-metrics.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
mv07g_projection <- read.csv(file.path(mv07g_evidence,
  "mv07g-full-ph-projection.csv"), stringsAsFactors = FALSE,
  check.names = FALSE)
if (mv07g_decision$decision != "authorize_MV7H_full_PH_landscape_prefreeze" ||
    !all(truth(mv07g_checks$passed)) || nrow(mv07g_ph) != 60L ||
    nrow(panel) != 500L || length(unique(panel$panel_sha256)) != 1L) {
  stop("MV7-G is not durably closed for MV7-H.")
}
mv07h_validate_cache_manifest_v1(manifest)
axis <- mv07h_sample_seed_axis_v1(manifest)
source_queue <- mv07h_source_queue_v1(axis)
ph_queue <- mv07h_ph_queue_v1(axis)
landscape_queue <- mv07h_landscape_queue_v1(mv07g_ph)
caps <- mv07h_resource_caps_v1()
rust_sha <- .mv07h_sha256(rust)
accepted_rust_sha <-
  "51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d"
if (rust_sha != accepted_rust_sha) {
  stop("MV7-H Rust library differs from the accepted accelerator.")
}
implementation_paths <- c(
  "R/mv07h_full_topology.R", "scripts/run_mv07h_source_entry.R",
  "scripts/run_mv07h_ph_entry.R", "scripts/run_mv07h_full_ph.R",
  "scripts/validate_mv07h_full_ph.R",
  "scripts/run_mv07h_landscape_group.R",
  "scripts/run_mv07h_landscape_stress.R",
  "scripts/validate_mv07h_landscape_stress.R"
)
implementation_hashes <- vapply(implementation_paths, .mv07h_sha256,
                                 character(1L))
implementation_root <- digest::digest(data.frame(
  path = implementation_paths, sha256 = unname(implementation_hashes),
  accepted_head = expected_head, stringsAsFactors = FALSE),
  algo = "sha256", serialize = TRUE)

mv07g_hours <- sum(mv07g_projection$projected_worker_hours)
historical_seconds <- 21538.531
historical_rows <- 141400L
projected_landscape_seconds <- historical_seconds * 152520 / historical_rows
resource_projection <- data.frame(
  contract_id = "mv07h_resource_projection_v1",
  component = c("full_source_and_PH", "all_landscape_groups", "total"),
  projected_worker_seconds = c(
    mv07g_hours * 3600, projected_landscape_seconds,
    mv07g_hours * 3600 + projected_landscape_seconds),
  projected_worker_hours = c(
    mv07g_hours, projected_landscape_seconds / 3600,
    mv07g_hours + projected_landscape_seconds / 3600),
  basis = c(
    "MV7-G conservative measured projection",
    "MV6-F combined source_PH_landscape time scaled by component rows",
    "sum of conservative upper-bound components"),
  stage1_stress_required = c(FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
contract <- data.frame(
  contract_id = "mv07h_full124_dual_view_topology_landscape_prefreeze_v1",
  accepted_head = expected_head, implementation_root_sha256 = implementation_root,
  rust_library_sha256 = rust_sha, samples = 124L, seeds = 5L,
  panel_genes = 500L, cells_per_sample = 384L, pca_components = 30L,
  source_jobs = 5L, typed_views = 1240L, ph_jobs = 1240L,
  landscape_groups = 20L, unordered_pairs_per_seed = 7626L,
  view_specific_pairs = 76260L, landscape_component_rows = 152520L,
  homology_dimensions = "H0;H1_separate", filtration = "complete_vietoris_rips",
  threshold = -1, field = 2L,
  landscape_definition =
    "finite_positive;essential_H0_excluded;all_active_levels;exact_squared_L2;no_grid;no_level_cap;streamed_groups",
  full_PH_cross_engine_checks = 20L, source_PH_repeat_artifacts = 13L,
  stress_R_oracle_checks = 3L, stress_analytic_oracle_checks = 1L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
firewall <- data.frame(
  contract_id = "mv07h_label_firewall_v1",
  allowed_execution_inputs =
    "sample_id;seed;cache_identity;selected_cell_identity;panel_identity;transform_identity;PH_identity",
  forbidden_execution_inputs =
    "tissue;approach;study;class;label;outcome;prior_benchmark_result",
  labels_may_select = FALSE, labels_may_stop = FALSE,
  cross_seed_pairs = 0L, combined_dimension_jobs = 0L, clustering_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, label_access_state = "closed",
  stringsAsFactors = FALSE
)
runtime <- data.frame(
  contract_id = "mv07h_runtime_v1", r_version = as.character(getRversion()),
  ripserr_version = as.character(utils::packageVersion("ripserr")),
  tda_version = as.character(utils::packageVersion("TDA")),
  rust_library_sha256 = rust_sha, rust_engine_version = 1L,
  omp_threads = 1L, openblas_threads = 1L, mkl_threads = 1L,
  stringsAsFactors = FALSE
)
source_files <- c(
  panel = file.path(mv07g_prefreeze, "mv07g-panel.csv"),
  cache_manifest = file.path(mv07g_prefreeze, "mv07g-cache-manifest.csv"),
  sentinel_axis = file.path(mv07g_prefreeze, "mv07g-sentinel-axis.csv"),
  landscape_contract = file.path(mv07g_prefreeze,
                                 "mv07g-landscape-contract.csv"),
  mv07g_validation_decision = file.path(mv07g_evidence,
                                        "mv07g-validation-decision.csv"),
  mv07g_independent_validation = file.path(mv07g_evidence,
                                           "mv07g-independent-validation.csv"),
  mv07g_ph_metrics = file.path(mv07g_evidence, "mv07g-ph-metrics.csv"),
  mv07g_projection = file.path(mv07g_evidence,
                               "mv07g-full-ph-projection.csv"),
  topology_runtime = "R/dual_view_topology.R",
  standardization_runtime = "R/mv03_pilot.R",
  cache_runtime = "R/mv05_resource_safe_execution.R",
  sentinel_runtime = "R/mv07g_sentinel.R",
  landscape_rust_shim = "R/landscape_rust_prototype.R",
  landscape_contract_runtime = "R/landscape_contract.R",
  landscape_reference_runtime = "R/landscape_reference.R",
  mv07h_runtime = "R/mv07h_full_topology.R",
  source_runner = "scripts/run_mv07h_source_entry.R",
  ph_runner = "scripts/run_mv07h_ph_entry.R",
  full_ph_monitor = "scripts/run_mv07h_full_ph.R",
  full_ph_validator = "scripts/validate_mv07h_full_ph.R",
  landscape_group_runner = "scripts/run_mv07h_landscape_group.R",
  stress_monitor = "scripts/run_mv07h_landscape_stress.R",
  stress_validator = "scripts/validate_mv07h_landscape_stress.R",
  builder = "scripts/build_mv07h_prefreeze.R",
  prefreeze_validator = "scripts/validate_mv07h_prefreeze.R",
  repeat_validator = "scripts/validate_mv07h_prefreeze_repeat.R",
  tests = "tests/testthat/test-mv07h-full-topology.R",
  specification =
    "docs/specifications/MV07H_FULL124_DUAL_VIEW_TOPOLOGY_LANDSCAPE_PREFREEZE_V1.md",
  mv06f_resource_audit =
    "docs/audits/MV06F_STAGE2_COMPLETE_PREFREEZE_2026-08-14.md",
  rust_manifest = "rust/scph_landscape_kernel/Cargo.toml",
  rust_lock = "rust/scph_landscape_kernel/Cargo.lock",
  rust_source = "rust/scph_landscape_kernel/src/lib.rs"
)
if (any(!file.exists(source_files))) stop("MV7-H source freeze is incomplete.")
freeze <- data.frame(
  contract_id = "mv07h_source_freeze_v1", source_id = names(source_files),
  artifact_locator = unname(source_files),
  sha256 = vapply(source_files, .mv07h_sha256, character(1L)),
  bytes = as.numeric(file.info(source_files)$size),
  accepted_head = expected_head, private_source = FALSE,
  stringsAsFactors = FALSE
)
acceptance <- data.frame(
  contract_id = "mv07h_prefreeze_acceptance_v1",
  category = c("mv07g_closure", "panel", "cache_axis", "source_queue",
               "ph_queue", "landscape_queue", "component_balance",
               "landscape_definition", "rust_identity", "resource_caps",
               "resource_projection", "repeat_oracles", "label_firewall",
               "no_execution"),
  passed = c(
    all(truth(mv07g_checks$passed)),
    nrow(panel) == 500L && length(unique(panel$panel_sha256)) == 1L,
    nrow(axis) == 620L && all(table(axis$seed) == 124L),
    nrow(source_queue) == 5L && sum(source_queue$typed_view_count) == 1240L,
    nrow(ph_queue) == 1240L && all(table(ph_queue$view_id) == 620L),
    nrow(landscape_queue) == 20L &&
      sum(landscape_queue$stage == "stage_1_stress") == 1L,
    all(table(landscape_queue$view_id) == 10L) &&
      all(table(landscape_queue$homology_dimension) == 10L) &&
      sum(landscape_queue$component_rows) == 152520L,
    nrow(landscape_contract) == 8L &&
      identical(landscape_contract$item, c("finite_intervals", "essential_h0",
        "level_policy", "integration", "dimension_policy", "grid_policy",
        "level_cap_policy", "streaming")),
    rust_sha == accepted_rust_sha,
    nrow(caps) == 5L && all(caps$workers == 1L) && all(caps$retries == 0L) &&
      caps$rss_cap_bytes[caps$stage == "landscape_group"] == 12 * 1024^3,
    resource_projection$projected_worker_hours[resource_projection$component ==
      "total"] <= 48,
    contract$source_PH_repeat_artifacts == 13L &&
      contract$stress_R_oracle_checks == 3L,
    firewall$label_access_state == "closed" &&
      sum(firewall[c("cross_seed_pairs", "combined_dimension_jobs",
                     "clustering_jobs", "label_jobs", "outcome_jobs")]) == 0,
    TRUE
  ),
  detail = c("11 accepted MV7-G categories", "exact MV7-FP 500",
             "124 samples by five seeds", "five reused-transform bundles",
             "620 cell plus 620 gene PH", "20 groups; one stress first",
             "76,260 pairs; 152,520 separate components",
             "eight dissertation-aligned requirements", "accepted Rust binary",
             "serial zero-retry bounded groups", "conservative total under 48h",
             "source/PH repeat plus stress R oracles", "labels outcomes closed",
             "zero MV7-H production jobs"),
  stringsAsFactors = FALSE
)
if (!all(acceptance$passed)) stop("MV7-H prefreeze acceptance failed: ",
  paste(acceptance$category[!acceptance$passed], collapse = ", "))
decision <- data.frame(
  contract_id = "mv07h_prefreeze_decision_v1",
  decision = "authorize_full_source_PH_then_one_landscape_stress_group",
  source_jobs_authorized = 5L, ph_jobs_authorized = 1240L,
  landscape_groups_authorized = 1L, landscape_groups_closed = 19L,
  clustering_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L,
  next_gate = "independent_full_PH_then_stress_repeat_R_oracle_resume",
  stringsAsFactors = FALSE
)
write_provenance_csv(panel, file.path(output, "mv07h-panel.csv"))
write_provenance_csv(manifest, file.path(output, "mv07h-cache-manifest.csv"))
write_provenance_csv(sentinel_axis, file.path(output, "mv07h-sentinel-axis.csv"))
outputs <- list(
  "mv07h-contract.csv" = contract,
  "mv07h-sample-seed-axis.csv" = axis,
  "mv07h-source-queue.csv" = source_queue,
  "mv07h-ph-queue.csv" = ph_queue,
  "mv07h-landscape-queue.csv" = landscape_queue,
  "mv07h-resource-caps.csv" = caps,
  "mv07h-resource-projection.csv" = resource_projection,
  "mv07h-landscape-contract.csv" = landscape_contract,
  "mv07h-label-firewall.csv" = firewall,
  "mv07h-runtime.csv" = runtime,
  "mv07h-source-freeze.csv" = freeze,
  "mv07h-acceptance.csv" = acceptance,
  "mv07h-decision.csv" = decision
)
for (name in names(outputs)) write_provenance_csv(
  outputs[[name]], file.path(output, name))
message("MV7-H prefreeze complete: 1,240 PH; 20 landscape groups; labels closed")
