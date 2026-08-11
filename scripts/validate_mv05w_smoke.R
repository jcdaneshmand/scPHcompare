#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 15L) {
  stop(paste(
    "usage: validate_mv05w_smoke.R QUEUE INVENTORY D1_ROOT G_ROOT",
    "RESULT_ROOT RESOURCE_CSV REPEAT_ROOT REPEAT_RESOURCE_CSV",
    "RESUME_BEFORE RESUME_AFTER RESUME_RESOURCE OUTPUT_DIR EXECUTION_HEAD",
    "PYTHON_EXECUTABLE PYTHON_SCRIPT_SHA256"
  ), call. = FALSE)
}
source("R/provenance_utils.R")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv03_pilot.R")
source("R/mv05_resource_safe_execution.R")
source("R/mv05_benchmark_execution.R")
source("R/mv05_inductive_mapping.R")
source("R/mv05d2_ph_profiling.R")
source("R/mv05d5_retrieval_inputs.R")
source("R/mv05f_integration_gate.R")
source("R/mv05h_integrated_ph_production.R")
source("R/mv05n_clustering_gate.R")
source("R/mv05t_robustness_gate.R")
source("R/mv05u_robustness_admission.R")
source("R/mv05v_robustness_prefreeze.R")
source("R/mv05w_launch_readiness.R")
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
queue <- read_csv(args[[1L]])
inventory <- read_csv(args[[2L]])
d1_root <- normalizePath(args[[3L]], mustWork = TRUE)
g_root <- normalizePath(args[[4L]], mustWork = TRUE)
result_root <- normalizePath(args[[5L]], mustWork = TRUE)
resources <- read_csv(args[[6L]])
repeat_root <- normalizePath(args[[7L]], mustWork = TRUE)
repeat_resources <- read_csv(args[[8L]])
before <- read_csv(args[[9L]])
after <- read_csv(args[[10L]])
resume_resources <- read_csv(args[[11L]])
output_dir <- args[[12L]]
execution_head <- args[[13L]]
python <- normalizePath(args[[14L]], winslash = "/", mustWork = TRUE)
python_sha <- args[[15L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
mv05w_validate_smoke_queue_v1(queue)
if (trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)) !=
      execution_head) stop("MV5-W validation HEAD drifted.", call. = FALSE)
unit <- queue[queue$execution_authorized, , drop = FALSE]
safe <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
first_dir <- file.path(result_root, safe(unit$robustness_group_id))
second_dir <- file.path(repeat_root, safe(unit$robustness_group_id))
validate_dir <- function(path) {
  manifest <- read_csv(file.path(path, "artifact_manifest.csv"))
  status <- read_csv(file.path(path, "status.csv"))
  if (nrow(manifest) != 9L || nrow(status) != 1L ||
      status$status != "completed" ||
      status$robustness_group_id != unit$robustness_group_id ||
      any(!vapply(seq_len(nrow(manifest)), function(index) {
        artifact <- file.path(path, manifest$artifact_file[[index]])
        file.exists(artifact) && file_sha(artifact) == manifest$sha256[[index]] &&
          file.info(artifact)$size == manifest$bytes[[index]]
      }, logical(1L)))) stop("MV5-W artifact validation failed.", call. = FALSE)
  manifest
}
first_manifest <- validate_dir(first_dir)
second_manifest <- validate_dir(second_dir)
deterministic <- first_manifest[first_manifest$deterministic_repeat_required, ]
repeat_match <- merge(
  deterministic[c("artifact_file", "sha256", "bytes")],
  second_manifest[c("artifact_file", "sha256", "bytes")],
  by = "artifact_file", suffixes = c("_first", "_repeat"), sort = TRUE
)
if (nrow(repeat_match) != 8L ||
    any(repeat_match$sha256_first != repeat_match$sha256_repeat) ||
    any(repeat_match$bytes_first != repeat_match$bytes_repeat)) {
  stop("MV5-W clean repeat differs.", call. = FALSE)
}
if (nrow(before) != 11L || nrow(after) != 11L ||
    !identical(before, after) || nrow(resume_resources) != 1L ||
    resume_resources$disposition != "reused_validated") {
  stop("MV5-W immutable resume failed.", call. = FALSE)
}

views <- read_csv(file.path(first_dir, "view_metrics.csv"))
pairs <- read_csv(file.path(first_dir, "pair_scope.csv"))
intervals <- read_csv(file.path(first_dir, "finite_intervals.csv"))
landscape <- read_csv(file.path(first_dir, "landscape_pairs.csv"))
energy <- read_csv(file.path(first_dir, "energy_pairs.csv"))
methods <- read_csv(file.path(first_dir, "method_rows.csv"))
if (nrow(views) != 90L || nrow(pairs) != unit$biological_pairs ||
    nrow(landscape) != unit$landscape_request_rows ||
    nrow(energy) != unit$energy_request_rows ||
    nrow(methods) != unit$assembled_method_rows ||
    length(unique(landscape$subchunk_id)) != unit$landscape_subchunks ||
    any(table(landscape$subchunk_id) > 250L)) {
  stop("MV5-W smoke cardinality failed.", call. = FALSE)
}

fold_study <- sub("^large_loso_v1:", "", unit$fold_id)
d1_path <- file.path(d1_root, paste0(fold_study, "__", unit$seed,
                                      "__sct_cell_fold.rds"))
g_stem <- paste0("mv05g_group__", fold_study, "__", unit$seed)
g_path <- file.path(g_root, "results", g_stem, paste0(g_stem, ".rds"))
d1_expected <- inventory[inventory$source_type == "sct" &
                           inventory$fold_study == fold_study &
                           inventory$seed == unit$seed, ]
g_expected <- inventory[inventory$source_type == "integrated" &
                          inventory$fold_study == fold_study &
                          inventory$seed == unit$seed, ]
if (file_sha(d1_path) != d1_expected$sha256 ||
    file_sha(g_path) != g_expected$sha256) {
  stop("MV5-W selected source hash failed.", call. = FALSE)
}
d1 <- readRDS(d1_path)
g <- readRDS(g_path)
ids <- sort(names(g$payload$coordinates), method = "radix")
source_views <- lapply(ids, function(sample_id) {
  mv05h_new_integrated_cell_view_v1(
    g$payload$coordinates[[sample_id]], sample_id, d1$identity$fold_id,
    d1$identity$fit_scope_id, unit$seed, g$cache_key,
    g$payload$coordinate_set_sha256
  )
})
names(source_views) <- ids
rebuilt <- lapply(source_views, mv05u_transform_view_v1,
                  configuration_id = unit$configuration_id)
expected_pairs <- mv05w_full_pair_coverage_v1(
  d1, unit$robustness_group_id, unit$base_pair_axis_sha256
)
if (!identical(pairs, expected_pairs) ||
    !identical(sort(names(rebuilt)), sort(views$sample_id))) {
  stop("MV5-W independent pair/view axis failed.", call. = FALSE)
}
all_mst <- TRUE
for (sample_id in names(rebuilt)) {
  observed <- sort(intervals$death[
    intervals$sample_id == sample_id & intervals$homology_dimension == "H0"
  ])
  expected <- sort(stats::hclust(stats::dist(rebuilt[[sample_id]]$payload),
                                 method = "single")$height)
  tolerance <- max(1e-7, max(expected) * 1e-7)
  if (length(observed) != nrow(rebuilt[[sample_id]]$payload) - 1L ||
      max(abs(observed - expected)) > tolerance) all_mst <- FALSE
}
if (!all_mst) stop("MV5-W all-view MST oracle failed.", call. = FALSE)
direct_energy <- mv05n_training_energy_pairs_v1(rebuilt, expected_pairs)
direct_energy <- direct_energy[match(energy$pair_request_id,
                                     direct_energy$pair_request_id), ]
if (max(abs(energy$distance - direct_energy$distance)) > 1e-9) {
  stop("MV5-W independent energy oracle failed.", call. = FALSE)
}
assembled <- mv05w_assemble_method_rows_v1(
  landscape, energy, unit$robustness_group_id
)
non_distance <- setdiff(names(methods), "distance")
if (!identical(methods[non_distance], assembled[non_distance]) ||
    max(abs(methods$distance - assembled$distance)) > 1e-12) {
  stop("MV5-W method assembly reconstruction failed.", call. = FALSE)
}
oracle <- system2(
  python, c("scripts/mv05w_exact_landscape_group.py", "--oracle-only",
            "--script-sha256", python_sha), stdout = TRUE, stderr = TRUE
)
if (!any(grepl("passed", oracle, fixed = TRUE))) {
  stop("MV5-W exact landscape oracle failed.", call. = FALSE)
}
if (nrow(resources) != 1L || resources$elapsed_seconds > 600 ||
    resources$peak_process_tree_rss_bytes > 4 * 1024^3 ||
    resources$cumulative_private_bytes > 4 * 1024^3 ||
    resources$labels_opened || resources$outcomes_computed ||
    nrow(repeat_resources) != 1L) {
  stop("MV5-W resource or label boundary failed.", call. = FALSE)
}

validation <- data.frame(
  contract_id = "mv05w_independent_validation_v1",
  validation_id = c(
    "bound_identity", "exact_group_axis", "complete_pair_axis",
    "landscape_subchunks", "all_view_h0_mst", "exact_landscape_oracle",
    "independent_energy", "four_method_assembly", "artifact_hashes",
    "clean_repeat", "immutable_resume", "resource_caps", "label_safety"
  ),
  passed = TRUE, labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
repeat_public <- data.frame(
  contract_id = "mv05w_repeat_validation_v1",
  artifact_file = repeat_match$artifact_file,
  sha256 = repeat_match$sha256_first, bytes = repeat_match$bytes_first,
  exact_repeat = TRUE, labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
resume_public <- data.frame(
  contract_id = "mv05w_resume_validation_v1", files_before = nrow(before),
  files_after = nrow(after), differences = 0L, reused_validated = TRUE,
  rebuilt = FALSE, labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv05w_launch_decision_v1", launch_readiness_complete = TRUE,
  configuration_execution_authorized = TRUE,
  authorized_configuration_id = unit$configuration_id,
  authorized_group_count = 150L,
  all_four_configurations_authorized = FALSE,
  next_action = "prospectively_bind_and_execute_first_configuration_only",
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(validation, file.path(output_dir,
  "mv05w-independent-validation-2026-08-10.csv"))
write_provenance_csv(repeat_public, file.path(output_dir,
  "mv05w-deterministic-repeat-2026-08-10.csv"))
write_provenance_csv(resume_public, file.path(output_dir,
  "mv05w-resume-validation-2026-08-10.csv"))
write_provenance_csv(decision, file.path(output_dir,
  "mv05w-launch-decision-2026-08-10.csv"))
message("MV5-W independent smoke validation passed: views=90 pairs=425 outcomes=0")
