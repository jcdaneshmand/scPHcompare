#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05x_configuration_execution.R OUTPUT_DIR PROSPECTIVE_HEAD PYTHON_EXECUTABLE",
    call. = FALSE
  )
}
source("R/provenance_utils.R")
source("R/mv05v_robustness_prefreeze.R")
source("R/mv05x_configuration_execution.R")
output_dir <- args[[1L]]
head <- args[[2L]]
python <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
if (!grepl("^[0-9a-f]{40}$", head)) stop("MV5-X HEAD is invalid.", call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE
)
file_sha <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
write_once <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  write_provenance_csv(value, path)
}

implementation_files <- c(
  dual_view = "R/dual_view_topology.R",
  ph_profile = "R/mv05d2_ph_profiling.R",
  energy = "R/mv05d5_retrieval_inputs.R",
  clustering_helpers = "R/mv05n_clustering_gate.R",
  mv05t = "R/mv05t_robustness_gate.R",
  mv05u = "R/mv05u_robustness_admission.R",
  mv05v = "R/mv05v_robustness_prefreeze.R",
  mv05w = "R/mv05w_launch_readiness.R",
  mv05x = "R/mv05x_configuration_execution.R",
  landscape = "scripts/mv05w_exact_landscape_group.py",
  group_runner = "scripts/run_mv05w_full_group.R",
  monitor = "scripts/monitor_mv05x_configuration.R"
)
support_files <- c(
  builder = "scripts/build_mv05x_configuration_execution.R",
  snapshot = "scripts/snapshot_mv05w_resume.R",
  tests = "tests/testthat/test-mv05x-configuration-execution.R",
  mv05w_queue = "docs/audits/mv05w-smoke-queue-2026-08-10.csv",
  mv05w_sources = "docs/audits/mv05w-source-freeze-2026-08-10.csv",
  mv05w_decision = "docs/audits/mv05w-launch-decision-2026-08-10.csv"
)
files <- c(implementation_files, support_files)
if (any(!file.exists(files))) stop("MV5-X source roster is incomplete.", call. = FALSE)
implementation_sha <- digest::digest(stats::setNames(
  vapply(implementation_files, file_sha, character(1L)),
  names(implementation_files)
), algo = "sha256", serialize = TRUE)
python_sha <- file_sha(python)
python_script_sha <- file_sha("scripts/mv05w_exact_landscape_group.py")
environment <- system2(
  python, c(
    "scripts/mv05w_exact_landscape_group.py", "--environment-only",
    "--script-sha256", python_script_sha
  ), stdout = TRUE, stderr = TRUE
)
parts <- strsplit(environment, "\t", fixed = TRUE)
runtime <- stats::setNames(
  vapply(parts, `[[`, character(1L), 2L),
  vapply(parts, `[[`, character(1L), 1L)
)
if (length(runtime) != 6L) stop("MV5-X runtime identity failed.", call. = FALSE)

old_sources <- read_csv(support_files[["mv05w_sources"]])
private_sources <- old_sources[old_sources$source_class == "private_coordinate", ]
private_sources$contract_id <- "mv05x_source_freeze_v1"
engine_sources <- data.frame(
  contract_id = "mv05x_source_freeze_v1",
  source_class = "committed_engine", source_id = names(files),
  source_locator = unname(files), sha256 = vapply(files, file_sha, character(1L)),
  bytes = as.numeric(file.info(files)$size), labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)
runtime_source <- data.frame(
  contract_id = "mv05x_source_freeze_v1", source_class = "private_runtime",
  source_id = "python_executable", source_locator = "private_runtime:python",
  sha256 = python_sha, bytes = as.numeric(file.info(python)$size),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
sources <- rbind(engine_sources, private_sources, runtime_source)
sources <- sources[order(sources$source_class, sources$source_id,
                         method = "radix"), , drop = FALSE]
source_sha <- digest::digest(stats::setNames(
  sources$sha256, paste(sources$source_class, sources$source_id, sep = "\r")
), algo = "sha256", serialize = TRUE)

queue <- read_csv(support_files[["mv05w_queue"]])
queue$contract_id <- "mv05x_pc20_execution_queue_v1"
queue$source_freeze_sha256 <- source_sha
queue$implementation_sha256 <- implementation_sha
queue$prospective_head <- head
queue$python_executable_sha256 <- python_sha
queue$python_version <- runtime[["python_version"]]
queue$persim_version <- runtime[["persim_version"]]
queue$numpy_version <- runtime[["numpy_version"]]
queue$scipy_version <- runtime[["scipy_version"]]
queue$execution_authorized <-
  queue$configuration_id == "cells384_pc20_euclidean_v1"
queue$execution_completed <- FALSE
queue$configuration_execution_order <- NA_integer_
authorized <- which(queue$execution_authorized)
queue$configuration_execution_order[authorized] <- seq_along(authorized)
mv05x_validate_configuration_queue_v1(queue)
scope <- mv05x_configuration_scope_v1(queue)
scope$implementation_sha256 <- implementation_sha
scope$python_script_sha256 <- python_script_sha
scope$prospective_head <- head
scope$execution_authorized <- TRUE

write_once(sources, file.path(output_dir, "mv05x-source-freeze-2026-08-10.csv"))
write_once(queue, file.path(output_dir, "mv05x-pc20-execution-queue-2026-08-10.csv"))
write_once(scope, file.path(output_dir, "mv05x-pc20-execution-scope-2026-08-10.csv"))
message("MV5-X bound 150 PC20 groups: implementation_sha=", implementation_sha)
