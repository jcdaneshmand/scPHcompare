#!/usr/bin/env Rscript

options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: build_mv05w_launch_readiness.R OUTPUT_DIR PROSPECTIVE_HEAD PYTHON_EXECUTABLE",
    call. = FALSE
  )
}
source("R/provenance_utils.R")
source("R/mv05v_robustness_prefreeze.R")
source("R/mv05w_launch_readiness.R")
output_dir <- args[[1L]]
head <- args[[2L]]
python <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
if (!grepl("^[0-9a-f]{40}$", head)) stop("MV5-W HEAD is invalid.", call. = FALSE)
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
  landscape = "scripts/mv05w_exact_landscape_group.py",
  group_runner = "scripts/run_mv05w_full_group.R",
  monitor = "scripts/monitor_mv05w_smoke.R"
)
support_files <- c(
  builder = "scripts/build_mv05w_launch_readiness.R",
  snapshot = "scripts/snapshot_mv05w_resume.R",
  tests = "tests/testthat/test-mv05w-launch-readiness.R",
  mv05v_queue = "docs/audits/mv05v-full-group-queue-2026-08-10.csv",
  mv05v_sources = "docs/audits/mv05v-source-freeze-2026-08-10.csv",
  mv05v_scope = "docs/audits/mv05v-base-pair-scope-2026-08-10.csv"
)
files <- c(implementation_files, support_files)
if (any(!file.exists(files))) stop("MV5-W source roster is incomplete.", call. = FALSE)
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
if (length(runtime) != 6L) stop("MV5-W runtime identity failed.", call. = FALSE)

old_sources <- read_csv(support_files[["mv05v_sources"]])
new_sources <- data.frame(
  contract_id = "mv05w_source_freeze_v1", source_class = "committed_engine",
  source_id = names(files), source_locator = unname(files),
  sha256 = vapply(files, file_sha, character(1L)),
  bytes = as.numeric(file.info(files)$size), labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)
private_sources <- old_sources[old_sources$source_class == "private_coordinate", ]
private_sources$contract_id <- "mv05w_source_freeze_v1"
runtime_source <- data.frame(
  contract_id = "mv05w_source_freeze_v1", source_class = "private_runtime",
  source_id = "python_executable", source_locator = "private_runtime:python",
  sha256 = python_sha, bytes = as.numeric(file.info(python)$size),
  labels_opened = FALSE, outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
sources <- rbind(new_sources, private_sources, runtime_source)
sources <- sources[order(sources$source_class, sources$source_id,
                         method = "radix"), ]
source_sha <- digest::digest(stats::setNames(
  sources$sha256, paste(sources$source_class, sources$source_id, sep = "\r")
), algo = "sha256", serialize = TRUE)

queue <- read_csv(support_files[["mv05v_queue"]])
queue$contract_id <- "mv05w_smoke_queue_v1"
queue$source_freeze_sha256 <- source_sha
queue$implementation_sha256 <- implementation_sha
queue$prospective_head <- head
queue$python_executable_sha256 <- python_sha
queue$python_version <- runtime[["python_version"]]
queue$persim_version <- runtime[["persim_version"]]
queue$numpy_version <- runtime[["numpy_version"]]
queue$scipy_version <- runtime[["scipy_version"]]
queue$execution_authorized <- FALSE
queue$execution_authorized[[1L]] <- TRUE
queue$execution_completed <- FALSE
mv05w_validate_smoke_queue_v1(queue)

summary <- data.frame(
  contract_id = "mv05w_launch_readiness_prefreeze_v1",
  sources = nrow(sources), queue_rows = nrow(queue),
  authorized_smoke_rows = sum(queue$execution_authorized),
  authorized_configuration_rows = 0L,
  selected_smoke_group_id = queue$robustness_group_id[[1L]],
  selected_smoke_pairs = queue$biological_pairs[[1L]],
  selected_smoke_landscape_rows = queue$landscape_request_rows[[1L]],
  selected_smoke_energy_rows = queue$energy_request_rows[[1L]],
  selected_smoke_method_rows = queue$assembled_method_rows[[1L]],
  implementation_sha256 = implementation_sha,
  python_script_sha256 = python_script_sha,
  execution_authorized = FALSE, labels_opened = FALSE,
  outcomes_computed = FALSE, stringsAsFactors = FALSE
)
write_once(sources, file.path(output_dir, "mv05w-source-freeze-2026-08-10.csv"))
write_once(queue, file.path(output_dir, "mv05w-smoke-queue-2026-08-10.csv"))
write_once(summary, file.path(output_dir, "mv05w-prefreeze-summary-2026-08-10.csv"))
message("MV5-W bound one smoke group: ", queue$robustness_group_id[[1L]],
        " implementation_sha=", implementation_sha)
