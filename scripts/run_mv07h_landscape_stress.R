#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RAYON_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv07h_landscape_stress.R PREFREEZE PH_ROOT RUST_LIBRARY PRIVATE_ROOT PUBLIC_DIR EXPECTED_HEAD")
}
prefreeze <- args[[1L]]
ph_root <- args[[2L]]
rust <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
private_root <- args[[4L]]
public_dir <- args[[5L]]
expected_head <- tolower(trimws(args[[6L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (head != expected_head) stop("MV7-H stress exact HEAD mismatch.")
source("R/mv07h_full_topology.R")
queue <- read.csv(file.path(prefreeze, "mv07h-landscape-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
contract <- read.csv(file.path(prefreeze, "mv07h-contract.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
row <- queue[queue$stage == "stage_1_stress", , drop = FALSE]
if (nrow(row) != 1L || contract$rust_library_sha256 != .mv07h_sha256(rust)) {
  stop("MV7-H stress group or Rust identity is stale.")
}
public_names <- c("mv07h-stress-resource.csv",
                  "mv07h-stress-repeat-resource.csv",
                  "mv07h-stress-scientific-repeat.csv",
                  "mv07h-stress-decision.csv")
public_resume <- dir.exists(public_dir)
if (public_resume && !all(file.exists(file.path(public_dir, public_names)))) {
  stop("MV7-H stress public resume state is incomplete.")
}
for (subdir in c("landscape", "logs", "repeat/landscape", "repeat/logs")) {
  dir.create(file.path(private_root, subdir), recursive = TRUE,
             showWarnings = FALSE)
}
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(e) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(e) 0), numeric(1L)))
}
run_once <- function(repeat_mode = FALSE) {
  root <- file.path(private_root, if (repeat_mode) "repeat/landscape" else
    "landscape")
  log_root <- file.path(private_root, if (repeat_mode) "repeat/logs" else
    "logs")
  ledger_path <- file.path(private_root, if (repeat_mode)
    "stress-repeat-metric.csv" else "stress-metric.csv")
  safe <- gsub(":", "_", row$group_id, fixed = TRUE)
  output_dir <- file.path(root, safe)
  if (file.exists(ledger_path) || dir.exists(output_dir)) {
    if (!file.exists(ledger_path) || !dir.exists(output_dir)) {
      stop("MV7-H stress resume state is ambiguous.")
    }
    metric <- read.csv(ledger_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
    distance <- file.path(output_dir, "distances.csv")
    if (nrow(metric) != 1L || metric$disposition != "completed" ||
        metric$distances_sha256 != .mv07h_sha256(distance)) {
      stop("MV7-H stress resume receipt is stale.")
    }
    return(metric)
  }
  prefix <- if (repeat_mode) "repeat__" else ""
  stdout <- file.path(log_root, paste0(prefix, "stress__stdout.txt"))
  stderr <- file.path(log_root, paste0(prefix, "stress__stderr.txt"))
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"), args = c(
      "--vanilla", "scripts/run_mv07h_landscape_group.R", prefreeze,
      ph_root, rust, row$group_id, root, contract$implementation_root_sha256
    ), stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  peak <- 0; cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25); peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > row$elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > row$rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  disposition <- if (nzchar(cap_failure)) cap_failure else if (
    identical(status, 0L) && dir.exists(output_dir)
  ) "completed" else "failed"
  distance <- file.path(output_dir, "distances.csv")
  metric <- data.frame(
    contract_id = "mv07h_landscape_stress_resource_v1",
    group_id = row$group_id, repeat_mode = repeat_mode,
    disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = row$elapsed_cap_seconds,
    rss_cap_bytes = row$rss_cap_bytes,
    distances_bytes = if (file.exists(distance)) file.info(distance)$size else NA,
    distances_sha256 = if (file.exists(distance)) .mv07h_sha256(distance) else NA,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  partial <- tempfile(pattern = basename(ledger_path),
                      tmpdir = dirname(ledger_path))
  write.csv(metric, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, ledger_path)) stop("Failed to publish stress receipt.")
  if (disposition != "completed") stop("MV7-H landscape stress failed.")
  metric
}
primary <- run_once(FALSE)
repeat_metric <- run_once(TRUE)
scientific_repeat <- data.frame(
  contract_id = "mv07h_landscape_stress_repeat_v1", group_id = row$group_id,
  artifact = "distances.csv", primary_bytes = primary$distances_bytes,
  repeat_bytes = repeat_metric$distances_bytes,
  primary_sha256 = primary$distances_sha256,
  repeat_sha256 = repeat_metric$distances_sha256,
  bytes_equal = primary$distances_bytes == repeat_metric$distances_bytes,
  sha256_equal = primary$distances_sha256 == repeat_metric$distances_sha256,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
decision <- data.frame(
  contract_id = "mv07h_landscape_stress_decision_v1",
  decision = if (primary$disposition == "completed" &&
    repeat_metric$disposition == "completed" && scientific_repeat$sha256_equal &&
    primary$elapsed_seconds <= row$elapsed_cap_seconds &&
    primary$peak_process_tree_rss_bytes <= row$rss_cap_bytes) {
      "stress_complete_await_independent_validation"
    } else "stop_stress_failure",
  group_id = row$group_id, component_rows = row$component_rows,
  remaining_groups_authorized = 0L, clustering_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
if (decision$decision == "stop_stress_failure") stop("MV7-H stress gate failed.")
parent <- dirname(public_dir); dir.create(parent, recursive = TRUE,
                                          showWarnings = FALSE)
staging <- tempfile(pattern = "mv07h-stress-public-", tmpdir = parent)
dir.create(staging)
write.csv(primary, file.path(staging, public_names[[1L]]), row.names = FALSE,
          na = "")
write.csv(repeat_metric, file.path(staging, public_names[[2L]]), row.names = FALSE,
          na = "")
write.csv(scientific_repeat, file.path(staging, public_names[[3L]]),
          row.names = FALSE, na = "")
write.csv(decision, file.path(staging, public_names[[4L]]), row.names = FALSE,
          na = "")
if (public_resume) {
  same <- vapply(public_names, function(name) {
    a <- read.csv(file.path(public_dir, name), stringsAsFactors = FALSE,
                  check.names = FALSE)
    b <- read.csv(file.path(staging, name), stringsAsFactors = FALSE,
                  check.names = FALSE)
    identical(names(a), names(b)) && nrow(a) == nrow(b) &&
      isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE))
  }, logical(1L))
  unlink(staging, recursive = TRUE)
  if (!all(same)) stop("MV7-H stress public resume evidence differs.")
  message("MV7-H landscape stress immutable resume complete")
} else if (!file.rename(staging, public_dir)) {
  unlink(staging, recursive = TRUE); stop("Failed to publish stress evidence.")
} else message("MV7-H landscape stress complete; validation required")
