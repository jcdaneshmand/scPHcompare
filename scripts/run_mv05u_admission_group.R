#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(paste(
    "usage: run_mv05u_admission_group.R EXECUTION_QUEUE PRIVATE_INVENTORY",
    "D1_ROOT G_ROOT RESULT_ROOT ADMISSION_UNIT_ID EXPECTED_HEAD",
    "IMPLEMENTATION_SHA256 PYTHON_EXECUTABLE PYTHON_SCRIPT_SHA256 RUN_MODE"
  ), call. = FALSE)
}

for (package in c("digest", "ripserr")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for MV5-U admission.", call. = FALSE)
  }
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

queue_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
inventory_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
d1_root <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
g_root <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
result_root <- args[[5L]]
unit_id <- args[[6L]]
expected_head <- args[[7L]]
expected_implementation_sha <- args[[8L]]
python_executable <- normalizePath(args[[9L]], winslash = "/", mustWork = TRUE)
python_script_sha <- args[[10L]]
run_mode <- match.arg(args[[11L]], c("build_or_resume", "validate_resume"))

read_public <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
safe_name <- function(value) gsub("[^A-Za-z0-9_.-]", "_", value)
bool <- function(value) as.logical(value)

head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
if (!identical(head, expected_head)) {
  stop("MV5-U runner is not at its bound prospective HEAD.", call. = FALSE)
}
implementation_files <- c(
  dual_view = "R/dual_view_topology.R",
  ph_profile = "R/mv05d2_ph_profiling.R",
  energy = "R/mv05d5_retrieval_inputs.R",
  clustering_helpers = "R/mv05n_clustering_gate.R",
  mv05t = "R/mv05t_robustness_gate.R",
  mv05u = "R/mv05u_robustness_admission.R",
  landscape = "scripts/mv05u_exact_landscape_group.py",
  group_runner = "scripts/run_mv05u_admission_group.R"
)
if (any(!file.exists(implementation_files))) {
  stop("MV5-U implementation source is missing.", call. = FALSE)
}
implementation_sha <- digest::digest(
  stats::setNames(vapply(implementation_files, file_sha, character(1L)),
                  names(implementation_files)),
  algo = "sha256", serialize = TRUE
)
if (!identical(implementation_sha, expected_implementation_sha) ||
    file_sha("scripts/mv05u_exact_landscape_group.py") != python_script_sha) {
  stop("MV5-U implementation identity is stale.", call. = FALSE)
}

queue <- read_public(queue_path)
mv05u_validate_execution_queue_v1(queue)
unit <- queue[queue$admission_unit_id == unit_id, , drop = FALSE]
if (nrow(unit) != 1L || unit$implementation_sha256 != implementation_sha) {
  stop("MV5-U requested admission unit is absent or stale.", call. = FALSE)
}
environment_lines <- system2(
  python_executable,
  c(
    "scripts/mv05u_exact_landscape_group.py", "--environment-only",
    "--script-sha256", shQuote(python_script_sha)
  ), stdout = TRUE, stderr = TRUE
)
environment_status <- attr(environment_lines, "status")
parts <- strsplit(environment_lines, "\t", fixed = TRUE)
runtime <- if (length(parts) == 6L && all(lengths(parts) == 2L)) {
  stats::setNames(vapply(parts, `[[`, character(1L), 2L),
                  vapply(parts, `[[`, character(1L), 1L))
} else character()
if ((!is.null(environment_status) && environment_status != 0L) ||
    length(runtime) != 6L ||
    file_sha(python_executable) != unit$python_executable_sha256 ||
    runtime[["python_version"]] != unit$python_version ||
    runtime[["persim_version"]] != unit$persim_version ||
    runtime[["numpy_version"]] != unit$numpy_version ||
    runtime[["scipy_version"]] != unit$scipy_version) {
  stop("MV5-U Python runtime identity is stale.", call. = FALSE)
}
inventory <- read_public(inventory_path)
required_inventory <- c(
  "source_type", "fold_study", "seed", "sha256", "labels_opened",
  "outcomes_computed", "admission_executed"
)
if (nrow(inventory) != 150L ||
    !all(required_inventory %in% names(inventory)) ||
    any(bool(inventory$labels_opened)) || any(bool(inventory$outcomes_computed)) ||
    any(bool(inventory$admission_executed))) {
  stop("MV5-U private inventory violates its closed source contract.",
       call. = FALSE)
}

fold_study <- sub("^large_loso_v1:", "", unit$fold_id)
seed <- as.integer(unit$seed)
d1_path <- file.path(
  d1_root, paste0(fold_study, "__", seed, "__sct_cell_fold.rds")
)
g_stem <- paste0("mv05g_group__", fold_study, "__", seed)
g_path <- file.path(g_root, "results", g_stem, paste0(g_stem, ".rds"))
expected_d1 <- inventory[
  inventory$source_type == "sct" & inventory$fold_study == fold_study &
    inventory$seed == seed, , drop = FALSE
]
expected_g <- inventory[
  inventory$source_type == "integrated" & inventory$fold_study == fold_study &
    inventory$seed == seed, , drop = FALSE
]
if (nrow(expected_d1) != 1L || nrow(expected_g) != 1L ||
    !file.exists(d1_path) || !file.exists(g_path) ||
    file_sha(d1_path) != expected_d1$sha256 ||
    file_sha(g_path) != expected_g$sha256) {
  stop("MV5-U selected private coordinate source is missing or stale.",
       call. = FALSE)
}

final_dir <- file.path(result_root, safe_name(unit_id))
validate_published <- function(path) {
  status_path <- file.path(path, "status.csv")
  manifest_path <- file.path(path, "artifact_manifest.csv")
  if (!dir.exists(path) || !file.exists(status_path) ||
      !file.exists(manifest_path)) return(FALSE)
  status <- read_public(status_path)
  manifest <- read_public(manifest_path)
  if (nrow(status) != 1L || status$status != "completed" ||
      status$admission_unit_id != unit_id ||
      status$implementation_sha256 != implementation_sha ||
      status$python_executable_sha256 != unit$python_executable_sha256 ||
      status$execution_queue_sha256 != file_sha(queue_path) ||
      status$source_freeze_sha256 != unit$source_freeze_sha256 ||
      nrow(manifest) != 8L || anyDuplicated(manifest$artifact_file)) {
    return(FALSE)
  }
  all(vapply(seq_len(nrow(manifest)), function(index) {
    artifact <- file.path(path, manifest$artifact_file[[index]])
    file.exists(artifact) && file_sha(artifact) == manifest$sha256[[index]] &&
      unname(file.info(artifact)$size) == manifest$bytes[[index]]
  }, logical(1L)))
}
if (dir.exists(final_dir)) {
  if (!validate_published(final_dir)) {
    stop("Existing MV5-U unit is partial, stale, or hash-invalid.",
         call. = FALSE)
  }
  message("MV5-U unit reused: ", unit_id)
  quit(status = 0L, save = "no")
}
if (run_mode == "validate_resume") {
  stop("MV5-U resume validation found a missing unit.", call. = FALSE)
}

dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
partial <- tempfile(pattern = paste0(".mv05u_", safe_name(unit_id), "_"),
                    tmpdir = result_root)
dir.create(partial, recursive = FALSE, showWarnings = FALSE)
published <- FALSE
on.exit(if (!published && dir.exists(partial)) unlink(partial, recursive = TRUE),
        add = TRUE)

started_total <- proc.time()[["elapsed"]]
d1 <- readRDS(d1_path)
g <- readRDS(g_path)
mv05d1_validate_cell_fold_record_v1(d1)
mv05f_validate_group_record_v1(g)
sct_views <- d1$payload$cell_views
integrated_coordinates <- g$payload$coordinates
if (!identical(names(sct_views), names(integrated_coordinates)) ||
    length(sct_views) != 90L || d1$identity$fold_id != unit$fold_id ||
    d1$identity$seed != seed) {
  stop("MV5-U selected source axes drifted.", call. = FALSE)
}
for (sample_id in names(sct_views)) {
  mv05t_validate_coordinate_pair_v1(
    sct_views[[sample_id]], integrated_coordinates[[sample_id]]
  )
}

if (unit$representation == "sct_whole") {
  source_views <- sct_views
} else if (unit$representation == "inductive_integrated") {
  ids <- sort(names(integrated_coordinates), method = "radix")
  source_views <- lapply(ids, function(sample_id) {
    mv05h_new_integrated_cell_view_v1(
      integrated_coordinates[[sample_id]], sample_id,
      d1$identity$fold_id, d1$identity$fit_scope_id, seed,
      g$cache_key, g$payload$coordinate_set_sha256
    )
  })
  names(source_views) <- ids
} else {
  stop("MV5-U representation is not admitted.", call. = FALSE)
}

transform_started <- proc.time()[["elapsed"]]
views <- lapply(source_views, mv05u_transform_view_v1,
                configuration_id = unit$configuration_id)
transform_seconds <- proc.time()[["elapsed"]] - transform_started
if (length(views) != 90L || any(vapply(views, function(view) {
  !identical(dim(view$payload), c(as.integer(unit$cells),
                                  as.integer(unit$coordinates))) ||
    any(!is.finite(view$payload))
}, logical(1L)))) {
  stop("MV5-U transformed view axis is incomplete.", call. = FALSE)
}

pairs <- mv05u_pair_coverage_v1(d1)
pair_path <- file.path(partial, "pair_coverage.csv")
write_provenance_csv(pairs, pair_path)

energy_started <- proc.time()[["elapsed"]]
energy <- mv05n_training_energy_pairs_v1(views, pairs)
energy_seconds <- proc.time()[["elapsed"]] - energy_started
energy$contract_id <- "mv05u_energy_pair_admission_v1"
energy$admission_unit_id <- unit_id
energy$fold_id <- unit$fold_id
energy$seed <- seed
energy$representation <- unit$representation
energy$configuration_id <- unit$configuration_id
energy$outcome_label_state <- "closed"
energy$biological_outcomes_computed <- FALSE
energy <- energy[c(
  "contract_id", "admission_unit_id", "fold_id", "seed", "representation",
  "configuration_id", "pair_request_id", "pair_ordinal", "pair_scope",
  "first_sample_id", "second_sample_id", "distance",
  "outcome_label_state", "biological_outcomes_computed"
)]
write_provenance_csv(energy, file.path(partial, "energy_pairs.csv"))

view_rows <- vector("list", length(views))
interval_rows <- list()
interval_cursor <- 0L
ph_seconds <- 0
mst_seconds <- 0
for (index in seq_along(views)) {
  view <- views[[index]]
  ph_started <- proc.time()[["elapsed"]]
  result <- run_topology_view_ph(view, max_dim = 1L, threshold = -1, field = 2L)
  ph_seconds <- ph_seconds + proc.time()[["elapsed"]] - ph_started
  mst_started <- proc.time()[["elapsed"]]
  oracle <- mv05u_validate_ph_result_v1(result, view)
  mst_seconds <- mst_seconds + proc.time()[["elapsed"]] - mst_started
  diagram <- result$diagram
  sample_id <- view$sample_id
  h0 <- diagram[diagram[, "dimension"] == 0, , drop = FALSE]
  h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
  view_rows[[index]] <- data.frame(
    contract_id = "mv05u_view_ph_metric_v1", admission_unit_id = unit_id,
    fold_id = unit$fold_id, seed = seed,
    representation = unit$representation,
    configuration_id = unit$configuration_id, sample_id = sample_id,
    point_count = nrow(view$payload), coordinate_count = ncol(view$payload),
    point_metric = view$point_metric,
    minimum_row_norm = min(sqrt(rowSums(view$payload^2))),
    maximum_row_norm = max(sqrt(rowSums(view$payload^2))),
    h0_intervals = nrow(h0), h1_intervals = nrow(h1),
    finite_intervals = sum(is.finite(diagram[, "death"])),
    essential_h0_intervals = sum(
      diagram[, "dimension"] == 0 & is.infinite(diagram[, "death"])
    ),
    h0_mst_maximum_absolute_error = oracle$maximum_absolute_error,
    h0_mst_tolerance = oracle$tolerance, h0_mst_oracle_passed = oracle$passed,
    source_view_cache_key = source_views[[sample_id]]$cache_key,
    transformed_view_cache_key = view$cache_key,
    transformed_payload_sha256 = view$payload_sha256,
    diagram_sha256 = result$provenance$diagram_sha256,
    outcome_label_state = "closed", outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  finite <- diagram[is.finite(diagram[, "death"]), , drop = FALSE]
  if (nrow(finite)) {
    interval_cursor <- interval_cursor + 1L
    interval_rows[[interval_cursor]] <- data.frame(
      contract_id = "mv05u_finite_interval_staging_v1",
      admission_unit_id = unit_id, sample_id = sample_id,
      homology_dimension = paste0("H", as.integer(finite[, "dimension"])),
      birth = finite[, "birth"], death = finite[, "death"],
      diagram_sha256 = result$provenance$diagram_sha256,
      essential_h0_policy = "exclude", outcome_label_state = "closed",
      outcomes_computed = FALSE, stringsAsFactors = FALSE
    )
  }
}
view_metrics <- do.call(rbind, view_rows)
view_metrics <- view_metrics[order(view_metrics$sample_id, method = "radix"), ]
intervals <- do.call(rbind, interval_rows)
intervals <- intervals[order(
  intervals$sample_id, intervals$homology_dimension, intervals$birth,
  intervals$death, method = "radix"
), ]
rownames(view_metrics) <- rownames(intervals) <- NULL
if (nrow(view_metrics) != 90L || any(!view_metrics$h0_mst_oracle_passed) ||
    any(view_metrics$h0_intervals != as.integer(unit$cells)) ||
    any(view_metrics$essential_h0_intervals != 1L)) {
  stop("MV5-U PH/MST view completion failed.", call. = FALSE)
}
write_provenance_csv(view_metrics, file.path(partial, "view_metrics.csv"))
write_provenance_csv(intervals, file.path(partial, "finite_intervals.csv"))

source_identity <- data.frame(
  contract_id = "mv05u_unit_source_identity_v1",
  admission_unit_id = unit_id, fold_id = unit$fold_id, seed = seed,
  representation = unit$representation,
  configuration_id = unit$configuration_id,
  source_freeze_sha256 = unit$source_freeze_sha256,
  execution_queue_sha256 = file_sha(queue_path),
  d1_source_sha256 = expected_d1$sha256,
  integrated_source_sha256 = expected_g$sha256,
  implementation_sha256 = implementation_sha,
  python_executable_sha256 = unit$python_executable_sha256,
  python_version = unit$python_version, persim_version = unit$persim_version,
  numpy_version = unit$numpy_version, scipy_version = unit$scipy_version,
  prospective_head = expected_head, outcome_label_state = "closed",
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)
write_provenance_csv(source_identity, file.path(partial, "source_identity.csv"))

python_stdout <- system2(
  python_executable,
  c(
    "scripts/mv05u_exact_landscape_group.py",
    "--views", shQuote(file.path(partial, "view_metrics.csv")),
    "--intervals", shQuote(file.path(partial, "finite_intervals.csv")),
    "--pairs", shQuote(pair_path),
    "--summary-output", shQuote(file.path(partial, "landscape_summary.csv")),
    "--pair-output", shQuote(file.path(partial, "landscape_pairs.csv")),
    "--admission-unit-id", shQuote(unit_id),
    "--script-sha256", shQuote(python_script_sha)
  ), stdout = TRUE, stderr = TRUE
)
python_status <- attr(python_stdout, "status")
if (!is.null(python_status) && python_status != 0L) {
  stop("MV5-U exact landscape subprocess failed: ",
       paste(python_stdout, collapse = "\n"), call. = FALSE)
}
landscape_summary <- read_public(file.path(partial, "landscape_summary.csv"))
landscape_pairs <- read_public(file.path(partial, "landscape_pairs.csv"))
if (nrow(landscape_summary) != 180L || nrow(landscape_pairs) != 64L ||
    any(!is.finite(landscape_pairs$distance)) ||
    any(landscape_pairs$distance < 0) ||
    any(bool(landscape_pairs$outcomes_computed))) {
  stop("MV5-U exact landscape outputs are incomplete.", call. = FALSE)
}

resource <- data.frame(
  contract_id = "mv05u_internal_resource_metric_v1",
  admission_unit_id = unit_id, transformed_views = 90L,
  pair_coverage_rows = 32L, landscape_summary_rows = 180L,
  landscape_pair_rows = 64L, energy_pair_rows = 32L,
  transform_seconds = transform_seconds, ph_seconds = ph_seconds,
  mst_seconds = mst_seconds, energy_seconds = energy_seconds,
  runner_elapsed_seconds = proc.time()[["elapsed"]] - started_total,
  outcome_label_state = "closed", outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
write_provenance_csv(resource, file.path(partial, "internal_resources.csv"))

artifact_files <- c(
  "source_identity.csv", "pair_coverage.csv", "view_metrics.csv",
  "finite_intervals.csv", "energy_pairs.csv", "landscape_summary.csv",
  "landscape_pairs.csv", "internal_resources.csv"
)
manifest <- data.frame(
  contract_id = "mv05u_unit_artifact_manifest_v1",
  admission_unit_id = unit_id, artifact_file = artifact_files,
  sha256 = vapply(file.path(partial, artifact_files), file_sha, character(1L)),
  bytes = as.numeric(file.info(file.path(partial, artifact_files))$size),
  deterministic_repeat_required = artifact_files != "internal_resources.csv",
  private_artifact = TRUE, outcome_label_state = "closed",
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)
write_provenance_csv(manifest, file.path(partial, "artifact_manifest.csv"))
status <- data.frame(
  contract_id = "mv05u_unit_status_v1", admission_unit_id = unit_id,
  execution_order = unit$execution_order, status = "completed",
  completed_views = 90L, landscape_summary_rows = 180L,
  landscape_pair_rows = 64L, energy_pair_rows = 32L,
  artifact_manifest_sha256 = file_sha(file.path(partial, "artifact_manifest.csv")),
  execution_queue_sha256 = file_sha(queue_path),
  source_freeze_sha256 = unit$source_freeze_sha256,
  implementation_sha256 = implementation_sha,
  python_executable_sha256 = unit$python_executable_sha256,
  python_version = unit$python_version, persim_version = unit$persim_version,
  prospective_head = expected_head, outcome_label_state = "closed",
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)
write_provenance_csv(status, file.path(partial, "status.csv"))

if (dir.exists(final_dir) || !file.rename(partial, final_dir)) {
  stop("MV5-U could not atomically publish the completed unit.",
       call. = FALSE)
}
published <- TRUE
if (!validate_published(final_dir)) {
  stop("MV5-U published unit failed immediate immutable validation.",
       call. = FALSE)
}
message("MV5-U unit built: ", unit_id)
