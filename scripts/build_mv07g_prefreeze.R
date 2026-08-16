#!/usr/bin/env Rscript

options(warn = 2)
for (package in c("digest", "pkgload", "ripserr", "TDA")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: build_mv07g_prefreeze.R MV07FP_DIR MV07FP_PREFREEZE MV07D_SENTINEL MV07E_DIR OUTPUT EXPECTED_HEAD")
}
panel_dir <- args[[1]]
panel_prefreeze <- args[[2]]
sentinel_path <- args[[3]]
mv07e <- args[[4]]
output <- args[[5]]
expected_head <- tolower(trimws(args[[6]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head) ||
    !grepl("^[0-9a-f]{40}$", expected_head)) {
  stop("MV7-G exact HEAD mismatch.")
}
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE))) {
  stop("MV7-G prefreeze output must be empty.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source("R/provenance_utils.R")
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
truth <- function(value) tolower(as.character(value)) %in% c("true", "t", "1")
panel <- read.csv(file.path(panel_dir, "mv07fp-panel.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
panel_decision <- read.csv(file.path(panel_dir, "mv07fp-decision.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
panel_validation <- read.csv(
  file.path(panel_dir, "mv07fp-independent-validation.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
panel_repeat <- read.csv(file.path(panel_dir, "mv07fp-repeat-validation.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
seed_validation <- read.csv(
  file.path(panel_dir, "mv07fp-seed-stability-validation.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(panel) != 500L ||
    panel_decision$decision != "lock_exact_124_panel_and_authorize_MV7G_prefreeze" ||
    !all(truth(panel_validation$passed)) || !all(truth(panel_repeat$sha256_equal)) ||
    !all(truth(seed_validation$passed))) {
  stop("MV7-FP is not durably closed for MV7-G.")
}
manifest_path <- file.path(panel_prefreeze, "mv07fp-cache-manifest.csv")
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE,
                     check.names = FALSE)
sentinels <- read.csv(sentinel_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
axis <- mv07g_sentinel_axis_v1(sentinels, manifest)
queue <- mv07g_queue_v1(axis)
caps <- mv07g_resource_caps_v1()
fit_scopes <- read.csv(file.path(mv07e, "mv07e-fit-scopes.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
typed <- read.csv(file.path(mv07e, "mv07e-typed-topology.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
landscape_path <- file.path(mv07e, "mv07e-landscape-contract.csv")
landscape <- read.csv(landscape_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
label <- read.csv(file.path(mv07e, "mv07e-label-firewall.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
contract <- data.frame(
  contract_id = "mv07g_typed_view_ph_sentinel_prefreeze_v1",
  samples = 6L, seeds = 5L, global_fit_jobs = 5L,
  global_fit_samples_per_seed = 124L,
  global_fit_cells_per_seed = 47616L, panel_genes = 500L,
  pca_components = 30L, typed_views = 60L, ph_jobs = 60L,
  homology_dimensions = "H0;H1", filtration = "complete_vietoris_rips",
  threshold = -1, field = 2L, cross_engine_checks = 24L,
  representative_repeat_artifacts = 13L,
  full_projection_worker_hours_cap = 24,
  full_projection_private_bytes_cap = 8 * 1024^3,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
cross <- data.frame(
  contract_id = "mv07g_cross_engine_contract_v1",
  seed = min(axis$seed), sentinel_samples = 6L,
  views = "cell_topology_v1;gene_topology_v1",
  deterministic_subset = "first_32_ordered_points_per_view",
  homology_dimensions = "H0;H1", expected_checks = 24L,
  primary_engine = "ripserr", corroborating_engine = "TDA_ripsDiag_GUDHI",
  maximum_absolute_error_tolerance = 1e-6,
  biological_interpretation = FALSE, stringsAsFactors = FALSE
)
firewall <- data.frame(
  contract_id = "mv07g_label_firewall_v1",
  allowed_execution_inputs =
    "frozen_sample_id;seed;cache_identity;feature_identity;panel_identity",
  sentinel_design_labels_used_at_execution = FALSE,
  labels_may_select = FALSE, labels_may_stop = FALSE,
  landscape_jobs = 0L, distance_jobs = 0L, clustering_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L, label_access_state = "closed",
  stringsAsFactors = FALSE
)
runtime <- data.frame(
  contract_id = "mv07g_runtime_v1",
  r_version = as.character(getRversion()),
  ripserr_version = as.character(utils::packageVersion("ripserr")),
  tda_version = as.character(utils::packageVersion("TDA")),
  omp_threads = 1L, openblas_threads = 1L, mkl_threads = 1L,
  stringsAsFactors = FALSE
)
source_files <- c(
  panel = file.path(panel_dir, "mv07fp-panel.csv"),
  panel_decision = file.path(panel_dir, "mv07fp-decision.csv"),
  panel_validation = file.path(panel_dir, "mv07fp-independent-validation.csv"),
  panel_repeat = file.path(panel_dir, "mv07fp-repeat-validation.csv"),
  seed_validation = file.path(panel_dir,
                              "mv07fp-seed-stability-validation.csv"),
  cache_manifest = manifest_path, sentinel = sentinel_path,
  fit_scopes = file.path(mv07e, "mv07e-fit-scopes.csv"),
  typed_topology = file.path(mv07e, "mv07e-typed-topology.csv"),
  landscape_contract = landscape_path,
  label_firewall = file.path(mv07e, "mv07e-label-firewall.csv"),
  topology_runtime = "R/dual_view_topology.R",
  standardization_runtime = "R/mv03_pilot.R",
  cache_runtime = "R/mv05_resource_safe_execution.R",
  sentinel_helper = "R/mv07g_sentinel.R",
  source_runner = "scripts/run_mv07g_source_entry.R",
  ph_runner = "scripts/run_mv07g_ph_entry.R",
  monitor = "scripts/run_mv07g_sentinel.R",
  validator = "scripts/validate_mv07g_sentinel.R",
  builder = "scripts/build_mv07g_prefreeze.R",
  prefreeze_validator = "scripts/validate_mv07g_prefreeze.R",
  tests = "tests/testthat/test-mv07g-sentinel.R",
  specification =
    "docs/specifications/MV07G_TYPED_VIEW_PH_SENTINEL_PREFREEZE_V1.md"
)
if (any(!file.exists(source_files))) stop("MV7-G source freeze incomplete.")
freeze <- data.frame(
  contract_id = "mv07g_source_freeze_v1", source_id = names(source_files),
  artifact_locator = unname(source_files),
  sha256 = vapply(source_files, sha, character(1L)),
  bytes = as.numeric(file.info(source_files)$size),
  accepted_head = expected_head, private_source = FALSE,
  stringsAsFactors = FALSE
)
acceptance <- data.frame(
  contract_id = "mv07g_prefreeze_acceptance_v1",
  category = c("mv07fp_closure", "cache_axis", "sentinel_axis", "fit_scopes",
               "typed_topology", "queue", "resource_caps", "cross_engine",
               "repeat", "landscape_unchanged", "label_firewall",
               "no_execution"),
  passed = c(
    TRUE, nrow(manifest) == 620L && length(unique(manifest$sample_id)) == 124L,
    nrow(axis) == 30L && all(table(axis$seed) == 6L),
    nrow(fit_scopes) == 5L && all(fit_scopes$fit_samples == 124L),
    nrow(typed) == 2L && setequal(typed$view,
      c("cell_topology_v1", "gene_topology_v1")),
    nrow(queue) == 65L && sum(queue$stage == "global_fit_views") == 5L &&
      sum(queue$stage %in% c("cell_ph", "gene_ph")) == 60L,
    nrow(caps) == 4L && all(caps$workers == 1L) && all(caps$retries == 0L),
    cross$expected_checks == 24L, contract$representative_repeat_artifacts == 13L,
    nrow(landscape) == 8L && !any(truth(landscape$changed_by_mv07e)),
    firewall$label_access_state == "closed" &&
      sum(firewall[c("landscape_jobs", "distance_jobs", "clustering_jobs",
                     "label_jobs", "outcome_jobs")]) == 0,
    TRUE
  ),
  detail = c("500-gene canonical panel locked", "124 by five immutable caches",
             "six by five frozen identities", "five global descriptive fits",
             "384-cell and 500-gene views", "five fits plus sixty PH jobs",
             "serial no retry atomic", "first seed twelve views H0/H1",
             "one seed thirteen artifacts", "eight requirements unchanged",
             "labels outcomes landscapes closed", "zero PCA and PH jobs"),
  stringsAsFactors = FALSE
)
if (!all(acceptance$passed)) stop("MV7-G prefreeze acceptance failed.")
decision <- data.frame(
  contract_id = "mv07g_prefreeze_decision_v1",
  decision = "authorize_serial_six_sample_five_seed_typed_view_PH_sentinel",
  global_fit_jobs_authorized = 5L, ph_jobs_authorized = 60L,
  landscape_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L,
  next_gate = "MV7H_full_PH_landscape_prefreeze_after_measured_validation",
  stringsAsFactors = FALSE
)
write_provenance_csv(panel, file.path(output, "mv07g-panel.csv"))
write_provenance_csv(manifest, file.path(output, "mv07g-cache-manifest.csv"))
outputs <- list(
  "mv07g-contract.csv" = contract,
  "mv07g-sentinel-axis.csv" = axis,
  "mv07g-queue.csv" = queue,
  "mv07g-resource-caps.csv" = caps,
  "mv07g-cross-engine-contract.csv" = cross,
  "mv07g-label-firewall.csv" = firewall,
  "mv07g-runtime.csv" = runtime,
  "mv07g-landscape-contract.csv" = landscape,
  "mv07g-source-freeze.csv" = freeze,
  "mv07g-acceptance.csv" = acceptance,
  "mv07g-decision.csv" = decision
)
for (name in names(outputs)) {
  write_provenance_csv(outputs[[name]], file.path(output, name))
}
message("MV7-G prefreeze complete: 5 fits, 60 PH jobs, labels closed")
