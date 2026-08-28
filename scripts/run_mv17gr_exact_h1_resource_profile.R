#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    paste(
      "usage: run_mv17gr_exact_h1_resource_profile.R",
      "<public-prefreeze> <private-prefreeze> <matrix-root>",
      "<private-output> <public-output>"
    ), call. = FALSE
  )
}
public_prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
private_prefreeze <- normalizePath(args[[2L]], mustWork = TRUE)
matrix_root <- normalizePath(args[[3L]], mustWork = TRUE)
private <- args[[4L]]
public <- args[[5L]]
public_partial <- paste0(public, ".partial")
if (dir.exists(private) || dir.exists(public) || dir.exists(public_partial)) {
  stop("MV17-GR output root exists", call. = FALSE)
}

source("R/mv08z_landscape_production.R")
source("R/mv17_calibration.R")
source("R/mv17g_parallel_recovery.R")
read_csv <- .mv08z_read_csv
write_csv <- .mv08z_atomic_csv
sha256 <- .mv08z_sha256_file
thread_environment <- mv17g_parallel_thread_environment_v1()
do.call(Sys.setenv, as.list(thread_environment))

manifest <- read_csv(file.path(public_prefreeze, "mv17gr-prefreeze-manifest.csv"))
manifest_paths <- file.path(public_prefreeze, manifest$artifact)
if (!all(file.exists(manifest_paths)) ||
    !identical(unname(as.numeric(file.info(manifest_paths)$size)),
               unname(as.numeric(manifest$bytes))) ||
    !identical(unname(vapply(manifest_paths, sha256, character(1L))),
               unname(tolower(manifest$sha256)))) {
  stop("MV17-GR public prefreeze drift", call. = FALSE)
}
private_binding <- read_csv(file.path(
  public_prefreeze, "mv17gr-private-binding.csv"
))
private_paths <- file.path(private_prefreeze, private_binding$artifact)
if (!all(file.exists(private_paths)) ||
    !identical(unname(as.numeric(file.info(private_paths)$size)),
               unname(as.numeric(private_binding$bytes))) ||
    !identical(unname(vapply(private_paths, sha256, character(1L))),
               unname(tolower(private_binding$sha256)))) {
  stop("MV17-GR private prefreeze drift", call. = FALSE)
}
contract <- read_csv(file.path(public_prefreeze, "mv17gr-contract.csv"))
queue <- read_csv(file.path(private_prefreeze, "mv17gr-profile-queue.csv"))
if (nrow(contract) != 1L || nrow(queue) != 10L ||
    contract$workers != 1L || contract$retries != 0L ||
    contract$attempt_timeout_seconds != 3600L ||
    contract$attempt_address_space_cap_bytes != 80 * 1024 ^ 3 ||
    !contract$execution_authorized_after_commit) {
  stop("MV17-GR execution contract drift", call. = FALSE)
}

dir.create(file.path(private, "results"), recursive = TRUE)
dir.create(file.path(private, "logs"), recursive = TRUE)
worker <- normalizePath(
  "scripts/run_mv17gr_exact_h1_profile_worker.R", mustWork = TRUE
)
rows <- vector("list", nrow(queue))
for (i in seq_len(nrow(queue))) {
  q <- queue[i, , drop = FALSE]
  stem <- sprintf("%02d__%s__%s", q$profile_order, q$case_role, q$engine)
  result_path <- file.path(private, "results", paste0(stem, ".rds"))
  time_path <- file.path(private, "logs", paste0(stem, ".time.txt"))
  stdout_path <- file.path(private, "logs", paste0(stem, ".stdout.txt"))
  stderr_path <- file.path(private, "logs", paste0(stem, ".stderr.txt"))
  matrix_path <- file.path(
    matrix_root, "matrices", sprintf("gene__%03d.rds", q$unit_order)
  )
  if (!file.exists(matrix_path)) stop("MV17-GR matrix missing", call. = FALSE)
  status <- suppressWarnings(system2(
    "/usr/bin/time",
    c(
      "-v", "-o", shQuote(time_path),
      "timeout", as.character(contract$attempt_timeout_seconds),
      "prlimit", paste0("--as=", contract$attempt_address_space_cap_bytes),
      "Rscript", "--vanilla", shQuote(worker), shQuote(matrix_path),
      q$null_family, as.character(q$seed), q$engine, shQuote(result_path)
    ),
    stdout = stdout_path, stderr = stderr_path
  ))
  resource <- mv17c_parse_gnu_time_v1(time_path)
  success <- identical(as.integer(status), 0L) &&
    resource$exit_status == 0L && file.exists(result_path)
  if (!success && file.exists(result_path)) {
    stop("MV17-GR failed attempt promoted a result", call. = FALSE)
  }
  rows[[i]] <- data.frame(
    contract_id = "mv17gr_profile_ledger_v1",
    profile_order = q$profile_order,
    case_role = q$case_role,
    engine = q$engine,
    exit_status = resource$exit_status,
    success = success,
    wall_seconds = resource$wall_seconds,
    maximum_RSS_bytes = resource$maximum_RSS_bytes,
    result_present = file.exists(result_path),
    stdout_bytes = as.numeric(file.info(stdout_path)$size),
    stderr_bytes = as.numeric(file.info(stderr_path)$size),
    labels_opened = FALSE,
    outcomes_opened = FALSE,
    stringsAsFactors = FALSE
  )
  write_csv(do.call(rbind, rows[seq_len(i)]),
            file.path(private, "mv17gr-private-progress.csv"))
}
ledger <- do.call(rbind, rows)
if (nrow(ledger) != 10L || anyNA(ledger) ||
    any(ledger$wall_seconds > contract$attempt_timeout_seconds + 10) ||
    any(ledger$maximum_RSS_bytes > contract$attempt_address_space_cap_bytes) ||
    length(list.files(private, pattern = "[.]partial$", recursive = TRUE))) {
  stop("MV17-GR profile evidence/cap drift", call. = FALSE)
}
write_csv(ledger, file.path(private, "mv17gr-private-ledger.csv"))
private_files <- sort(list.files(private, recursive = TRUE, full.names = TRUE))
binding <- data.frame(
  contract_id = "mv17gr_execution_private_binding_v1",
  files = length(private_files),
  bytes = sum(as.numeric(file.info(private_files)$size)),
  artifact_set_sha256 = digest::digest(
    paste(sort(vapply(private_files, sha256, character(1L))), collapse = "\n"),
    algo = "sha256", serialize = FALSE
  ),
  tracking_state = "private_not_tracked",
  stringsAsFactors = FALSE
)
aggregate <- do.call(rbind, lapply(
  split(ledger, interaction(ledger$case_role, ledger$engine, drop = TRUE)),
  function(z) data.frame(
    contract_id = "mv17gr_public_resource_v1",
    case_role = z$case_role[[1L]], engine = z$engine[[1L]],
    attempts = nrow(z), successes = sum(z$success),
    aggregate_wall_seconds = sum(z$wall_seconds),
    maximum_RSS_bytes = max(z$maximum_RSS_bytes),
    stringsAsFactors = FALSE
  )
))
status <- data.frame(
  contract_id = "mv17gr_execution_status_v1",
  attempts = nrow(ledger), successes = sum(ledger$success),
  failures = sum(!ledger$success), workers = 1L, retries = 0L,
  partials = 0L, labels_opened = FALSE, outcomes_opened = FALSE,
  downstream_surfaces = "closed", stringsAsFactors = FALSE
)
dir.create(public_partial, recursive = TRUE)
write_csv(aggregate, file.path(public_partial, "mv17gr-resource-summary.csv"))
write_csv(status, file.path(public_partial, "mv17gr-status.csv"))
write_csv(binding, file.path(public_partial, "mv17gr-private-result-binding.csv"))
files <- sort(list.files(public_partial))
write_csv(data.frame(
  contract_id = "mv17gr_execution_manifest_v1", artifact = files,
  bytes = as.numeric(file.info(file.path(public_partial, files))$size),
  sha256 = vapply(file.path(public_partial, files), sha256, character(1L)),
  stringsAsFactors = FALSE
), file.path(public_partial, "mv17gr-artifact-manifest.csv"))
if (!file.rename(public_partial, public)) {
  stop("MV17-GR public promotion failed", call. = FALSE)
}
message("MV17-GR profile complete; attempts=10; successes=", sum(ledger$success))
