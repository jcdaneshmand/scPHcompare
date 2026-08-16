#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: run_mv07g_sentinel.R PREFREEZE PRIMARY_CACHE ADDED_CACHE PRIVATE_ROOT PUBLIC_DIR EXPECTED_HEAD")
}
prefreeze <- args[[1]]
primary_cache <- args[[2]]
added_cache <- args[[3]]
private_root <- args[[4]]
public_dir <- args[[5]]
expected_head <- tolower(trimws(args[[6]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) stop("MV7-G exact HEAD mismatch.")
public_names <- c(
  "mv07g-source-metrics.csv", "mv07g-ph-metrics.csv",
  "mv07g-repeat-validation.csv", "mv07g-full-ph-projection.csv",
  "mv07g-decision.csv"
)
public_resume <- dir.exists(public_dir)
if (public_resume &&
    !all(file.exists(file.path(public_dir, public_names)))) {
  stop("MV7-G public resume state is incomplete or ambiguous.")
}
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
sha <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
queue <- read.csv(file.path(prefreeze, "mv07g-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
caps <- read.csv(file.path(prefreeze, "mv07g-resource-caps.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(queue) != 65L || any(queue$workers != 1L) || any(queue$retries != 0L)) {
  stop("MV7-G queue violates the serial no-retry contract.")
}
for (subdir in c("source", "ph", "logs", "repeat/source", "repeat/ph",
                 "repeat/logs")) {
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
read_metrics <- function(path) {
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE,
                                  check.names = FALSE) else data.frame()
}
write_metrics <- function(value, path) {
  partial <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  write.csv(value, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) {
    unlink(partial)
    stop("Failed to atomically publish private MV7-G metrics.")
  }
}
run_unit <- function(row, repeat_mode = FALSE) {
  ledger_path <- file.path(private_root,
    if (repeat_mode) "repeat-metrics.csv" else "metrics.csv")
  ledger <- read_metrics(ledger_path)
  job_id <- if (repeat_mode) paste0("repeat__", row$job_id) else row$job_id
  relative <- row$output_file
  output <- file.path(private_root, if (repeat_mode) "repeat" else "", relative)
  hit <- if (nrow(ledger)) ledger[ledger$job_id == job_id, , drop = FALSE] else
    data.frame()
  if (file.exists(output) || nrow(hit)) {
    if (nrow(hit) != 1L || !file.exists(output) ||
        hit$output_sha256 != sha(output) ||
        hit$output_bytes != as.numeric(file.info(output)$size) ||
        hit$disposition != "completed") {
      stop("MV7-G ambiguous or stale resume state for ", job_id)
    }
    return(hit)
  }
  script <- if (row$stage == "global_fit_views") {
    "scripts/run_mv07g_source_entry.R"
  } else {
    "scripts/run_mv07g_ph_entry.R"
  }
  script_args <- if (row$stage == "global_fit_views") {
    c(prefreeze, primary_cache, added_cache, output, as.character(row$seed))
  } else {
    source_path <- file.path(private_root, if (repeat_mode) "repeat" else "",
      "source", paste0("mv07g__", row$seed, "__source.rds"))
    c(source_path, row$sample_id, row$view_id, output)
  }
  log_root <- file.path(private_root, if (repeat_mode) "repeat/logs" else "logs")
  log_stem <- gsub("[^A-Za-z0-9_.-]", "_", job_id)
  stdout <- file.path(log_root, paste0(log_stem, "__stdout.txt"))
  stderr <- file.path(log_root, paste0(log_stem, "__stderr.txt"))
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"), args = c("--vanilla", script, script_args),
    stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  peak <- 0
  cap_failure <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > row$elapsed_cap_seconds) {
      cap_failure <- "elapsed_cap_exceeded"
      process$kill_tree()
    } else if (peak > row$rss_cap_bytes) {
      cap_failure <- "rss_cap_exceeded"
      process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  disposition <- if (nzchar(cap_failure)) cap_failure else if (
    identical(status, 0L) && file.exists(output)
  ) "completed" else "failed"
  result <- data.frame(
    contract_id = "mv07g_private_resource_metric_v1", job_id = job_id,
    stage = row$stage, seed = row$seed, sample_id = row$sample_id,
    view_id = row$view_id, disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    output_file = relative,
    output_bytes = if (file.exists(output)) as.numeric(file.info(output)$size) else NA,
    output_sha256 = if (file.exists(output)) sha(output) else NA_character_,
    elapsed_cap_seconds = row$elapsed_cap_seconds,
    rss_cap_bytes = row$rss_cap_bytes,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  ledger <- if (nrow(ledger)) rbind(ledger, result) else result
  write_metrics(ledger, ledger_path)
  if (disposition != "completed") stop("MV7-G job failed: ", job_id)
  result
}
primary <- lapply(seq_len(nrow(queue)), function(index) run_unit(queue[index,]))
primary <- do.call(rbind, primary)
repeat_queue <- queue[queue$seed == min(queue$seed), , drop = FALSE]
repeated <- lapply(seq_len(nrow(repeat_queue)), function(index) {
  run_unit(repeat_queue[index,], repeat_mode = TRUE)
})
repeated <- do.call(rbind, repeated)
repeat_rows <- lapply(seq_len(nrow(repeat_queue)), function(index) {
  row <- repeat_queue[index,]
  first <- file.path(private_root, row$output_file)
  second <- file.path(private_root, "repeat", row$output_file)
  data.frame(
    contract_id = "mv07g_repeat_validation_v1", job_id = row$job_id,
    file = row$output_file, original_bytes = as.numeric(file.info(first)$size),
    repeat_bytes = as.numeric(file.info(second)$size),
    original_sha256 = sha(first), repeat_sha256 = sha(second),
    bytes_equal = as.numeric(file.info(first)$size) ==
      as.numeric(file.info(second)$size),
    sha256_equal = identical(sha(first), sha(second)),
    stringsAsFactors = FALSE
  )
})
repeat_validation <- do.call(rbind, repeat_rows)
if (!all(repeat_validation$bytes_equal & repeat_validation$sha256_equal)) {
  stop("MV7-G representative repeat is not byte-identical.")
}
source_rows <- lapply(primary$seed[primary$stage == "global_fit_views"],
                      function(seed) {
  metric <- primary[primary$stage == "global_fit_views" & primary$seed == seed,]
  path <- file.path(private_root, metric$output_file)
  record <- readRDS(path)
  mv07g_validate_source_record_v1(record)
  data.frame(
    contract_id = "mv07g_source_metric_v1", seed = seed,
    fit_scope_id = record$identity$fit_scope_id, fit_samples = 124L,
    fit_cells = 47616L, panel_size = 500L, pca_components = 30L,
    sentinel_samples = length(record$views), typed_views = 2L * length(record$views),
    source_cache_key = record$cache_key,
    pca_model_cache_key = record$pca_model$cache_key,
    elapsed_seconds = metric$elapsed_seconds,
    peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
    output_bytes = metric$output_bytes, output_sha256 = metric$output_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
source_metrics <- do.call(rbind, source_rows)
ph_rows <- lapply(which(primary$stage %in% c("cell_ph", "gene_ph")),
                  function(index) {
  metric <- primary[index,]
  record <- readRDS(file.path(private_root, metric$output_file))
  mv07g_validate_ph_record_v1(record)
  data.frame(
    contract_id = "mv07g_ph_metric_v1", job_id = metric$job_id,
    seed = metric$seed, sample_id = metric$sample_id,
    view_id = metric$view_id,
    point_count = record$h0_mst_oracle$point_count,
    finite_h0_intervals = record$h0_mst_oracle$finite_h0_intervals,
    finite_h1_intervals = record$h0_mst_oracle$finite_h1_intervals,
    h0_mst_maximum_absolute_error =
      record$h0_mst_oracle$maximum_absolute_error,
    h0_mst_tolerance = record$h0_mst_oracle$tolerance,
    h0_mst_passed = record$h0_mst_oracle$passed,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    ph_cache_key = record$cache_key,
    elapsed_seconds = metric$elapsed_seconds,
    peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
    output_bytes = metric$output_bytes, output_sha256 = metric$output_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
ph_metrics <- do.call(rbind, ph_rows)
source_factor <- 124 / 6
projection <- data.frame(
  contract_id = "mv07g_full_ph_projection_v1",
  component = c("five_global_fits_and_all_typed_views", "cell_ph", "gene_ph"),
  sentinel_units = c(5L, 30L, 30L), full_units = c(5L, 620L, 620L),
  projected_worker_seconds = c(
    sum(source_metrics$elapsed_seconds) * source_factor,
    mean(ph_metrics$elapsed_seconds[ph_metrics$view_id == "cell_topology_v1"]) * 620,
    mean(ph_metrics$elapsed_seconds[ph_metrics$view_id == "gene_topology_v1"]) * 620
  ),
  projection_basis = c("conservative_full_source_bundle_linear_from_six_views",
                       "mean_sentinel_view", "mean_sentinel_view"),
  stringsAsFactors = FALSE
)
projection$projected_worker_hours <- projection$projected_worker_seconds / 3600
mean_source_bytes <- mean(source_metrics$output_bytes)
mean_cell_bytes <- mean(ph_metrics$output_bytes[
  ph_metrics$view_id == "cell_topology_v1"])
mean_gene_bytes <- mean(ph_metrics$output_bytes[
  ph_metrics$view_id == "gene_topology_v1"])
projection$projected_private_bytes <- c(
  mean_source_bytes * 5 * source_factor,
  mean_cell_bytes * 620, mean_gene_bytes * 620)
aggregate_elapsed <- sum(primary$elapsed_seconds) + sum(repeated$elapsed_seconds)
private_files <- c(
  list.files(file.path(private_root, "source"), recursive = TRUE,
             full.names = TRUE, all.files = TRUE, no.. = TRUE),
  list.files(file.path(private_root, "ph"), recursive = TRUE,
             full.names = TRUE, all.files = TRUE, no.. = TRUE),
  list.files(file.path(private_root, "logs"), recursive = TRUE,
             full.names = TRUE, all.files = TRUE, no.. = TRUE),
  list.files(file.path(private_root, "repeat"), recursive = TRUE,
             full.names = TRUE, all.files = TRUE, no.. = TRUE),
  file.path(private_root, c("metrics.csv", "repeat-metrics.csv"))
)
private_files <- unique(private_files[file.exists(private_files)])
private_files <- private_files[file.info(private_files)$isdir %in% FALSE]
private_bytes <- sum(file.info(private_files)$size)
aggregate_cap <- caps[caps$stage == "aggregate",]
decision <- data.frame(
  contract_id = "mv07g_sentinel_decision_v1",
  decision = if (nrow(source_metrics) == 5L && nrow(ph_metrics) == 60L &&
    all(ph_metrics$h0_mst_passed) && nrow(repeat_validation) == 13L &&
    all(repeat_validation$sha256_equal) &&
    aggregate_elapsed <= aggregate_cap$elapsed_cap_seconds &&
    private_bytes <= aggregate_cap$storage_cap_bytes)
      "sentinel_complete_await_independent_validation" else
      "stop_sentinel_failure",
  source_jobs = nrow(source_metrics), ph_jobs = nrow(ph_metrics),
  repeat_artifacts = nrow(repeat_validation),
  aggregate_elapsed_seconds = aggregate_elapsed,
  aggregate_elapsed_cap_seconds = aggregate_cap$elapsed_cap_seconds,
  private_bytes = private_bytes,
  private_storage_cap_bytes = aggregate_cap$storage_cap_bytes,
  landscape_jobs = 0L, distance_jobs = 0L, clustering_jobs = 0L,
  label_jobs = 0L, outcome_jobs = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
if (decision$decision == "stop_sentinel_failure") stop("MV7-G sentinel failed.")
parent <- dirname(public_dir)
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
staging <- tempfile(pattern = "mv07g-public-", tmpdir = parent)
dir.create(staging)
write.csv(source_metrics, file.path(staging, "mv07g-source-metrics.csv"),
          row.names = FALSE, na = "")
write.csv(ph_metrics, file.path(staging, "mv07g-ph-metrics.csv"),
          row.names = FALSE, na = "")
write.csv(repeat_validation,
          file.path(staging, "mv07g-repeat-validation.csv"),
          row.names = FALSE, na = "")
write.csv(projection, file.path(staging, "mv07g-full-ph-projection.csv"),
          row.names = FALSE, na = "")
write.csv(decision, file.path(staging, "mv07g-decision.csv"),
          row.names = FALSE, na = "")
if (public_resume) {
  same <- vapply(public_names, function(name) {
    first <- file.path(public_dir, name)
    second <- file.path(staging, name)
    accepted <- read.csv(first, stringsAsFactors = FALSE, check.names = FALSE)
    regenerated <- read.csv(second, stringsAsFactors = FALSE,
                            check.names = FALSE)
    identical(names(accepted), names(regenerated)) &&
      nrow(accepted) == nrow(regenerated) &&
      isTRUE(all.equal(accepted, regenerated, tolerance = 1e-12,
                       check.attributes = FALSE))
  }, logical(1L))
  unlink(staging, recursive = TRUE)
  if (!all(same)) stop("MV7-G regenerated public evidence differs on resume.")
  message("MV7-G immutable resume complete; public evidence unchanged")
} else if (!file.rename(staging, public_dir)) {
  unlink(staging, recursive = TRUE)
  stop("Failed to atomically publish MV7-G public evidence.")
} else {
  message("MV7-G sentinel production complete; independent validation required")
}
