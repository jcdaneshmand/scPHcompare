#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste("usage: build_mv05ab_cosine_execution.R OUTPUT_DIR",
             "PROSPECTIVE_HEAD PYTHON_EXECUTABLE"), call. = FALSE)
}
source("R/provenance_utils.R")
source("R/mv05ab_cosine_execution.R")
output_dir <- args[[1L]]
head <- args[[2L]]
python <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
if (!grepl("^[0-9a-f]{40}$", head)) stop("MV5-AB HEAD is invalid.", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE)
write_once <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

implementation_files <- c(
  provenance = "R/provenance_utils.R", toy = "R/toy_baseline.R",
  dual_view = "R/dual_view_topology.R", mv03 = "R/mv03_pilot.R",
  resource_safe = "R/mv05_resource_safe_execution.R",
  benchmark = "R/mv05_benchmark_execution.R",
  inductive = "R/mv05_inductive_mapping.R",
  ph_profile = "R/mv05d2_ph_profiling.R",
  energy = "R/mv05d5_retrieval_inputs.R",
  integration_gate = "R/mv05f_integration_gate.R",
  integrated_ph = "R/mv05h_integrated_ph_production.R",
  clustering_helpers = "R/mv05n_clustering_gate.R",
  mv05t = "R/mv05t_robustness_gate.R",
  mv05u = "R/mv05u_robustness_admission.R",
  mv05v = "R/mv05v_robustness_prefreeze.R",
  mv05w = "R/mv05w_launch_readiness.R",
  mv05ab = "R/mv05ab_cosine_execution.R",
  landscape = "scripts/mv05w_exact_landscape_group.py",
  group_runner = "scripts/run_mv05ab_cosine_group.R",
  monitor = "scripts/monitor_mv05ab_cosine_configuration.R")
support_files <- c(
  builder = "scripts/build_mv05ab_cosine_execution.R",
  specification = "docs/specifications/MV05AB_LABEL_CLOSED_COSINE_EXECUTION_SPECIFICATION_V1.md",
  tests = "tests/testthat/test-mv05ab-cosine-execution.R",
  aa_queue = "docs/audits/mv05aa-cosine-execution-queue-2026-08-11.csv",
  aa_decision = "docs/audits/mv05aa-continuation-decision-2026-08-11.csv",
  aa_source_freeze = "docs/audits/mv05aa-source-freeze-2026-08-11.csv",
  v_source_freeze = "docs/audits/mv05v-source-freeze-2026-08-10.csv",
  private_inventory = "docs/audits/mv05t-private-coordinate-inventory-2026-08-10.csv")
files <- c(implementation_files, support_files)
if (any(!file.exists(files))) {
  stop("MV5-AB source roster is incomplete: ",
       paste(names(files)[!file.exists(files)], collapse = ", "), call. = FALSE)
}
implementation_sha <- digest::digest(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)),
  names(implementation_files)), algo = "sha256", serialize = TRUE)
python_sha <- file_sha(python)
python_script_sha <- file_sha("scripts/mv05w_exact_landscape_group.py")
environment <- system2(python, c(
  "scripts/mv05w_exact_landscape_group.py", "--environment-only",
  "--script-sha256", python_script_sha), stdout = TRUE, stderr = TRUE)
parts <- strsplit(environment, "\t", fixed = TRUE)
runtime <- if (length(parts) == 6L && all(lengths(parts) == 2L)) {
  stats::setNames(vapply(parts, `[[`, character(1L), 2L),
                  vapply(parts, `[[`, character(1L), 1L))
} else character()
if (length(runtime) != 6L) stop("MV5-AB runtime identity failed.", call. = FALSE)

inventory <- read_csv(support_files[["private_inventory"]])
if (nrow(inventory) != 150L || anyDuplicated(inventory$sha256) ||
    any(as.logical(inventory$labels_opened)) ||
    any(as.logical(inventory$outcomes_computed)) ||
    any(as.logical(inventory$admission_executed))) {
  stop("MV5-AB private coordinate inventory drifted.", call. = FALSE)
}
private_sources <- data.frame(
  contract_id = "mv05ab_source_freeze_v1",
  source_class = "private_coordinate",
  source_id = paste(inventory$source_type, inventory$fold_study,
                    inventory$seed, sep = ":"),
  source_locator = inventory$private_locator,
  sha256 = inventory$sha256, bytes = as.numeric(inventory$bytes),
  labels_opened = FALSE, rankings_computed = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE)
engine_sources <- data.frame(
  contract_id = "mv05ab_source_freeze_v1",
  source_class = "committed_engine", source_id = names(files),
  source_locator = unname(files), sha256 = vapply(files, file_sha, character(1L)),
  bytes = as.numeric(file.info(files)$size), labels_opened = FALSE,
  rankings_computed = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
runtime_source <- data.frame(
  contract_id = "mv05ab_source_freeze_v1",
  source_class = "private_runtime", source_id = "python_executable",
  source_locator = "private_runtime:python", sha256 = python_sha,
  bytes = as.numeric(file.info(python)$size), labels_opened = FALSE,
  rankings_computed = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE)
sources <- rbind(engine_sources, private_sources, runtime_source)
sources <- sources[order(sources$source_class, sources$source_id,
                         method = "radix"), , drop = FALSE]
source_sha <- digest::digest(stats::setNames(
  sources$sha256, paste(sources$source_class, sources$source_id, sep = "\r")),
  algo = "sha256", serialize = TRUE)

queue <- read_csv(support_files[["aa_queue"]])
queue$contract_id <- "mv05ab_cosine_execution_queue_v1"
queue$source_freeze_sha256 <- source_sha
queue$implementation_sha256 <- implementation_sha
queue$prospective_head <- head
queue$python_executable_sha256 <- python_sha
queue$python_version <- runtime[["python_version"]]
queue$persim_version <- runtime[["persim_version"]]
queue$numpy_version <- runtime[["numpy_version"]]
queue$scipy_version <- runtime[["scipy_version"]]
queue$configuration_execution_order <- seq_len(nrow(queue))
queue$execution_authorized <- TRUE
queue$execution_completed <- FALSE
queue$labels_opened <- FALSE
queue$rankings_computed <- FALSE
queue$outcomes_computed <- FALSE
mv05ab_validate_cosine_queue_v1(queue)
scope <- mv05ab_cosine_scope_v1(queue)
scope$implementation_sha256 <- implementation_sha
scope$python_script_sha256 <- python_script_sha
scope$prospective_head <- head
scope$execution_authorized <- TRUE
scope$execution_completed <- FALSE

write_once(sources, file.path(output_dir, "mv05ab-source-freeze-2026-08-11.csv"))
write_once(queue, file.path(output_dir, "mv05ab-cosine-execution-queue-2026-08-11.csv"))
write_once(scope, file.path(output_dir, "mv05ab-cosine-execution-scope-2026-08-11.csv"))
message("MV5-AB bound 150 cosine groups: implementation_sha=", implementation_sha)
