#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: run_mv07h_full_ph.R PREFREEZE PRIMARY_CACHE ADDED_CACHE MV07G_ROOT PRIVATE_ROOT PUBLIC_DIR EXPECTED_HEAD")
}
prefreeze <- args[[1L]]
primary_cache <- args[[2L]]
added_cache <- args[[3L]]
mv07g_root <- args[[4L]]
private_root <- args[[5L]]
public_dir <- args[[6L]]
expected_head <- tolower(trimws(args[[7L]]))
head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (!identical(head, expected_head)) stop("MV7-H exact HEAD mismatch.")
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv07h_full_topology.R")

source_queue <- read.csv(file.path(prefreeze, "mv07h-source-queue.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
ph_queue <- read.csv(file.path(prefreeze, "mv07h-ph-queue.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
sentinel_axis <- read.csv(file.path(prefreeze, "mv07h-sentinel-axis.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
caps <- read.csv(file.path(prefreeze, "mv07h-resource-caps.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
fallback_policy <- read.csv(file.path(prefreeze, "mv07h-ph-fallback-policy.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
mv07h_validate_ph_fallback_policy_v1(fallback_policy)
if (nrow(source_queue) != 5L || nrow(ph_queue) != 1240L ||
    any(source_queue$workers != 1L) || any(ph_queue$workers != 1L) ||
    any(source_queue$retries != 0L) || any(ph_queue$retries != 0L)) {
  stop("MV7-H full-PH queue violates its serial no-retry contract.")
}
public_names <- c(
  "mv07h-source-metrics.csv", "mv07h-ph-metrics.csv",
  "mv07h-ph-engine-attempts.csv",
  "mv07h-sentinel-equivalence.csv", "mv07h-repeat-validation.csv",
  "mv07h-full-ph-decision.csv"
)
public_resume <- dir.exists(public_dir)
if (public_resume && !all(file.exists(file.path(public_dir, public_names)))) {
  stop("MV7-H public resume state is incomplete or ambiguous.")
}
for (subdir in c("source", "ph", "logs", "fallback/logs",
                 "repeat/source", "repeat/ph", "repeat/logs",
                 "repeat/fallback/logs")) {
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
read_ledger <- function(path) {
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE,
                                  check.names = FALSE) else data.frame()
}
write_ledger <- function(value, path) {
  partial <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  write.csv(value, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) {
    unlink(partial)
    stop("Failed to atomically publish MV7-H private metrics.")
  }
}
augment_metric <- function(value, ph_engine, attempt_kind,
                           fallback_trigger = "") {
  value$ph_engine <- ph_engine
  value$attempt_kind <- attempt_kind
  value$fallback_trigger <- fallback_trigger
  value
}
run_attempt <- function(row, repeat_mode = FALSE, fallback = FALSE) {
  ledger_name <- if (fallback) {
    if (repeat_mode) "repeat-fallback-metrics.csv" else "fallback-metrics.csv"
  } else if (repeat_mode) "repeat-metrics.csv" else "metrics.csv"
  ledger_path <- file.path(private_root, ledger_name)
  ledger <- read_ledger(ledger_path)
  job_id <- if (repeat_mode) paste0("repeat__", row$job_id) else row$job_id
  prefix <- if (repeat_mode) "repeat" else ""
  output <- file.path(private_root, prefix, row$output_file)
  hit <- if (nrow(ledger)) ledger[ledger$job_id == job_id, , drop = FALSE] else
    data.frame()
  if (file.exists(output) || nrow(hit)) {
    completed <- if (nrow(hit)) hit[hit$disposition == "completed",,
                                    drop = FALSE] else data.frame()
    if (nrow(completed) == 1L && file.exists(output) &&
        completed$output_sha256 == .mv07h_sha256(output) &&
        completed$output_bytes == as.numeric(file.info(output)$size)) {
      engine <- if (fallback) "TDA_ripsDiag_GUDHI" else if (
        row$stage == "source_views") "source_reconstruction" else "ripserr"
      return(augment_metric(
        completed, engine,
        if (fallback) "exact_resource_fallback" else "primary",
        if (fallback) "rss_cap_exceeded" else ""))
    }
    if (!file.exists(output) && !fallback && nrow(hit) == 1L &&
        hit$disposition == "rss_cap_exceeded") {
      return(augment_metric(hit, if (row$stage == "source_views") {
        "source_reconstruction"
      } else "ripserr", "primary", ""))
    }
    if (!file.exists(output) && fallback && nrow(hit) == 1L &&
        hit$disposition != "completed") {
      return(augment_metric(hit, "TDA_ripsDiag_GUDHI",
                            "exact_resource_fallback", "rss_cap_exceeded"))
    }
    {
      stop("MV7-H ambiguous or stale resume state for ", job_id)
    }
  }
  is_source <- row$stage == "source_views"
  if (fallback && (is_source || row$view_id != "gene_topology_v1")) {
    stop("MV7-H GUDHI fallback is restricted to gene PH jobs.")
  }
  script <- if (is_source) {
    "scripts/run_mv07h_source_entry.R"
  } else if (fallback) {
    "scripts/run_mv07h_ph_fallback_entry.R"
  } else "scripts/run_mv07h_ph_entry.R"
  script_args <- if (is_source) {
    parent <- file.path(mv07g_root, "source",
                        paste0("mv07g__", row$seed, "__source.rds"))
    c(prefreeze, primary_cache, added_cache, parent, output,
      as.character(row$seed))
  } else {
    source_path <- file.path(private_root, prefix, "source",
                             paste0("mv07h__", row$seed, "__source.rds"))
    c(source_path, row$sample_id, row$view_id, output)
  }
  log_root <- if (fallback) file.path(private_root, prefix, "fallback", "logs")
    else file.path(private_root, prefix, "logs")
  stem <- gsub("[^A-Za-z0-9_.-]", "_", job_id)
  stdout <- file.path(log_root, paste0(stem, "__stdout.txt"))
  stderr <- file.path(log_root, paste0(stem, "__stderr.txt"))
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"), args = c("--vanilla", script, script_args),
    stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  peak <- 0
  cap_failure <- ""
  elapsed_cap <- if (fallback) {
    fallback_policy$fallback_elapsed_cap_seconds
  } else row$elapsed_cap_seconds
  rss_cap <- if (fallback) fallback_policy$fallback_rss_cap_bytes else
    row$rss_cap_bytes
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > elapsed_cap) {
      cap_failure <- "elapsed_cap_exceeded"; process$kill_tree()
    } else if (peak > rss_cap) {
      cap_failure <- "rss_cap_exceeded"; process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  disposition <- if (nzchar(cap_failure)) cap_failure else if (
    identical(status, 0L) && file.exists(output)
  ) "completed" else "failed"
  result <- data.frame(
    contract_id = "mv07h_private_resource_metric_v1", job_id = job_id,
    stage = row$stage, seed = row$seed,
    sample_id = if (is_source) NA_character_ else row$sample_id,
    view_id = if (is_source) NA_character_ else row$view_id,
    disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    output_file = row$output_file,
    output_bytes = if (file.exists(output)) as.numeric(file.info(output)$size)
      else NA_real_,
    output_sha256 = if (file.exists(output)) .mv07h_sha256(output)
      else NA_character_,
    elapsed_cap_seconds = elapsed_cap,
    rss_cap_bytes = rss_cap,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  ledger <- if (nrow(ledger)) rbind(ledger, result) else result
  write_ledger(ledger, ledger_path)
  engine <- if (fallback) "TDA_ripsDiag_GUDHI" else if (is_source) {
    "source_reconstruction"
  } else "ripserr"
  augment_metric(result, engine,
                 if (fallback) "exact_resource_fallback" else "primary",
                 if (fallback) "rss_cap_exceeded" else "")
}
run_fallback_unit <- function(row, repeat_mode = FALSE) {
  result <- run_attempt(row, repeat_mode = repeat_mode, fallback = TRUE)
  if (result$disposition != "completed") {
    stop("MV7-H exact GUDHI fallback failed: ", result$job_id)
  }
  result
}
run_unit <- function(row, repeat_mode = FALSE) {
  primary <- run_attempt(row, repeat_mode = repeat_mode, fallback = FALSE)
  if (primary$disposition == "completed") return(primary)
  eligible <- mv07h_ph_fallback_eligible_v1(
    row$stage, row$view_id, primary$disposition, fallback_policy)
  if (!eligible) stop("MV7-H job failed: ", primary$job_id)
  run_fallback_unit(row, repeat_mode = repeat_mode)
}

primary_source <- do.call(rbind, lapply(seq_len(nrow(source_queue)),
  function(index) run_unit(source_queue[index, , drop = FALSE])))
primary_ph <- do.call(rbind, lapply(seq_len(nrow(ph_queue)),
  function(index) run_unit(ph_queue[index, , drop = FALSE])))
primary <- rbind(primary_source, primary_ph)

repeat_seed <- min(source_queue$seed)
repeat_source_row <- source_queue[source_queue$seed == repeat_seed,, drop = FALSE]
repeat_source <- run_unit(repeat_source_row, repeat_mode = TRUE)
sentinel_ids <- sort(unique(sentinel_axis$sample_id[
  sentinel_axis$seed == repeat_seed]), method = "radix")
repeat_ph_queue <- ph_queue[ph_queue$seed == repeat_seed &
  ph_queue$sample_id %in% sentinel_ids, , drop = FALSE]
repeat_ph <- do.call(rbind, lapply(seq_len(nrow(repeat_ph_queue)),
  function(index) run_unit(repeat_ph_queue[index, , drop = FALSE], TRUE)))
default_repeat_queue <- rbind(
  repeat_source_row[c("job_id", "output_file")],
  repeat_ph_queue[c("job_id", "output_file")]
)
fallback_job_ids <- primary_ph$job_id[
  primary_ph$ph_engine == "TDA_ripsDiag_GUDHI"]
fallback_repeat_queue <- ph_queue[
  ph_queue$job_id %in% setdiff(fallback_job_ids,
                               default_repeat_queue$job_id),, drop = FALSE]
fallback_repeated <- if (nrow(fallback_repeat_queue)) do.call(rbind, lapply(
  seq_len(nrow(fallback_repeat_queue)), function(index) run_fallback_unit(
    fallback_repeat_queue[index, , drop = FALSE], repeat_mode = TRUE
  ))) else primary_ph[FALSE,, drop = FALSE]
repeated <- rbind(repeat_source, repeat_ph, fallback_repeated)
repeat_queue <- unique(rbind(
  default_repeat_queue,
  fallback_repeat_queue[c("job_id", "output_file")]
))
repeat_validation <- do.call(rbind, lapply(seq_len(nrow(repeat_queue)),
  function(index) {
    row <- repeat_queue[index, , drop = FALSE]
    first <- file.path(private_root, row$output_file)
    second <- file.path(private_root, "repeat", row$output_file)
    data.frame(
      contract_id = "mv07h_repeat_validation_v1", job_id = row$job_id,
      file = row$output_file, original_bytes = as.numeric(file.info(first)$size),
      repeat_bytes = as.numeric(file.info(second)$size),
      original_sha256 = .mv07h_sha256(first),
      repeat_sha256 = .mv07h_sha256(second),
      bytes_equal = file.info(first)$size == file.info(second)$size,
      sha256_equal = .mv07h_sha256(first) == .mv07h_sha256(second),
      stringsAsFactors = FALSE
    )
  }))
if (nrow(repeat_validation) != nrow(repeat_queue) ||
    !all(repeat_validation$sha256_equal)) {
  stop("MV7-H representative source/PH repeat is not byte-identical.")
}

source_rows <- lapply(seq_len(nrow(primary_source)), function(index) {
  metric <- primary_source[index, , drop = FALSE]
  record <- readRDS(file.path(private_root, metric$output_file))
  mv07h_validate_source_record_v1(record)
  data.frame(
    contract_id = "mv07h_source_metric_v1", seed = metric$seed,
    fit_scope_id = record$identity$fit_scope_id, fit_samples = 124L,
    fit_cells = 47616L, panel_size = 500L, pca_components = 30L,
    samples = length(record$views), typed_views = 2L * length(record$views),
    parent_mv07g_source_cache_key =
      record$identity$parent_mv07g_source_cache_key,
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

ph_rows <- lapply(seq_len(nrow(primary_ph)), function(index) {
  metric <- primary_ph[index, , drop = FALSE]
  record <- readRDS(file.path(private_root, metric$output_file))
  mv07h_validate_ph_record_v1(record)
  data.frame(
    contract_id = "mv07h_ph_metric_v1", job_id = metric$job_id,
    seed = metric$seed, sample_id = metric$sample_id,
    view_id = metric$view_id, point_count = record$h0_mst_oracle$point_count,
    finite_h0_intervals = record$h0_mst_oracle$finite_h0_intervals,
    finite_h1_intervals = record$h0_mst_oracle$finite_h1_intervals,
    h0_mst_maximum_absolute_error =
      record$h0_mst_oracle$maximum_absolute_error,
    h0_mst_tolerance = record$h0_mst_oracle$tolerance,
    h0_mst_passed = record$h0_mst_oracle$passed,
    ph_engine = record$topology_result$provenance$ph_engine,
    ph_engine_version = record$topology_result$provenance$ph_engine_version,
    attempt_kind = metric$attempt_kind,
    fallback_trigger = metric$fallback_trigger,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    ph_cache_key = record$cache_key, elapsed_seconds = metric$elapsed_seconds,
    peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
    output_bytes = metric$output_bytes, output_sha256 = metric$output_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
})
ph_metrics <- do.call(rbind, ph_rows)

equivalence <- list()
for (seed in source_queue$seed) {
  full <- readRDS(file.path(private_root, "source",
                            paste0("mv07h__", seed, "__source.rds")))
  parent <- readRDS(file.path(mv07g_root, "source",
                              paste0("mv07g__", seed, "__source.rds")))
  for (sample_id in names(parent$views)) for (view_id in .mv07h_views) {
    first <- parent$views[[sample_id]][[view_id]]
    second <- full$views[[sample_id]][[view_id]]
    equivalence[[length(equivalence) + 1L]] <- data.frame(
      contract_id = "mv07h_mv07g_sentinel_equivalence_v1", seed = seed,
      sample_id = sample_id, view_id = view_id,
      mv07g_cache_key = first$cache_key, mv07h_cache_key = second$cache_key,
      mv07g_payload_sha256 = first$payload_sha256,
      mv07h_payload_sha256 = second$payload_sha256,
      cache_key_equal = first$cache_key == second$cache_key,
      payload_sha256_equal = first$payload_sha256 == second$payload_sha256,
      passed = first$cache_key == second$cache_key &&
        first$payload_sha256 == second$payload_sha256,
      stringsAsFactors = FALSE
    )
  }
}
equivalence <- do.call(rbind, equivalence)
if (nrow(equivalence) != 60L || !all(equivalence$passed)) {
  stop("MV7-H did not reproduce every accepted MV7-G sentinel view.")
}

primary_attempts <- read_ledger(file.path(private_root, "metrics.csv"))
fallback_attempts <- read_ledger(file.path(private_root, "fallback-metrics.csv"))
repeat_attempts <- read_ledger(file.path(private_root, "repeat-metrics.csv"))
repeat_fallback_attempts <- read_ledger(file.path(
  private_root, "repeat-fallback-metrics.csv"))
attempts <- list(
  primary = primary_attempts[primary_attempts$stage != "source_views",,
                             drop = FALSE],
  fallback = fallback_attempts,
  repeat_primary = repeat_attempts[repeat_attempts$stage != "source_views",,
                                   drop = FALSE],
  repeat_fallback = repeat_fallback_attempts
)
attempt_rows <- lapply(names(attempts), function(scope) {
  value <- attempts[[scope]]
  if (!nrow(value)) return(NULL)
  value$attempt_scope <- scope
  value$ph_engine <- if (grepl("fallback", scope, fixed = TRUE)) {
    "TDA_ripsDiag_GUDHI"
  } else "ripserr"
  value
})
attempt_rows <- Filter(Negate(is.null), attempt_rows)
engine_attempts <- do.call(rbind, attempt_rows)
engine_attempts <- engine_attempts[c(
  "contract_id", "job_id", "attempt_scope", "stage", "seed", "sample_id",
  "view_id", "ph_engine", "disposition", "exit_status", "elapsed_seconds",
  "peak_process_tree_rss_bytes", "elapsed_cap_seconds", "rss_cap_bytes",
  "outcome_label_state", "biological_outcomes_computed"
)]
engine_attempts$contract_id <- "mv07h_ph_engine_attempt_v1"
aggregate_elapsed <- sum(primary_attempts$elapsed_seconds) +
  sum(fallback_attempts$elapsed_seconds) + sum(repeat_attempts$elapsed_seconds) +
  sum(repeat_fallback_attempts$elapsed_seconds)
private_files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                            all.files = TRUE, no.. = TRUE)
private_files <- private_files[file.info(private_files)$isdir %in% FALSE]
private_bytes <- sum(file.info(private_files)$size)
aggregate_cap <- caps[caps$stage == "aggregate", , drop = FALSE]
decision <- data.frame(
  contract_id = "mv07h_full_ph_decision_v1",
  decision = if (nrow(source_metrics) == 5L && nrow(ph_metrics) == 1240L &&
    all(ph_metrics$h0_mst_passed) && all(ph_metrics$finite_h1_intervals > 0L) &&
    nrow(equivalence) == 60L && all(equivalence$passed) &&
    nrow(repeat_validation) == nrow(repeat_queue) &&
    all(repeat_validation$sha256_equal) &&
    aggregate_elapsed <= aggregate_cap$elapsed_cap_seconds &&
    private_bytes <= aggregate_cap$storage_cap_bytes) {
      "full_PH_complete_await_independent_validation"
    } else "stop_full_PH_failure",
  source_jobs = nrow(source_metrics), typed_views = sum(source_metrics$typed_views),
  ph_jobs = nrow(ph_metrics), sentinel_equivalence_checks = nrow(equivalence),
  repeat_artifacts = nrow(repeat_validation),
  ripserr_selected_records = sum(ph_metrics$ph_engine == "ripserr"),
  gudhi_fallback_records = sum(
    ph_metrics$ph_engine == "TDA_ripsDiag_GUDHI"),
  rss_triggered_fallbacks = sum(
    engine_attempts$attempt_scope == "fallback" &
      engine_attempts$disposition == "completed"),
  aggregate_elapsed_seconds = aggregate_elapsed,
  aggregate_elapsed_cap_seconds = aggregate_cap$elapsed_cap_seconds,
  private_bytes = private_bytes,
  private_storage_cap_bytes = aggregate_cap$storage_cap_bytes,
  landscape_jobs = 0L, clustering_jobs = 0L, label_jobs = 0L,
  outcome_jobs = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
if (decision$decision == "stop_full_PH_failure") stop("MV7-H full PH failed.")
parent <- dirname(public_dir)
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
staging <- tempfile(pattern = "mv07h-full-ph-public-", tmpdir = parent)
dir.create(staging)
write.csv(source_metrics, file.path(staging, public_names[[1L]]),
          row.names = FALSE, na = "")
write.csv(ph_metrics, file.path(staging, public_names[[2L]]),
          row.names = FALSE, na = "")
write.csv(engine_attempts, file.path(staging, public_names[[3L]]),
          row.names = FALSE, na = "")
write.csv(equivalence, file.path(staging, public_names[[4L]]),
          row.names = FALSE, na = "")
write.csv(repeat_validation, file.path(staging, public_names[[5L]]),
          row.names = FALSE, na = "")
write.csv(decision, file.path(staging, public_names[[6L]]),
          row.names = FALSE, na = "")
if (public_resume) {
  same <- vapply(public_names, function(name) {
    accepted <- read.csv(file.path(public_dir, name), stringsAsFactors = FALSE,
                         check.names = FALSE)
    regenerated <- read.csv(file.path(staging, name), stringsAsFactors = FALSE,
                            check.names = FALSE)
    identical(names(accepted), names(regenerated)) &&
      nrow(accepted) == nrow(regenerated) && isTRUE(all.equal(
        accepted, regenerated, tolerance = 1e-12, check.attributes = FALSE))
  }, logical(1L))
  unlink(staging, recursive = TRUE)
  if (!all(same)) stop("MV7-H regenerated public evidence differs on resume.")
  message("MV7-H full-PH immutable resume complete; public evidence unchanged")
} else if (!file.rename(staging, public_dir)) {
  unlink(staging, recursive = TRUE)
  stop("Failed to atomically publish MV7-H full-PH public evidence.")
} else {
  message("MV7-H full PH complete; independent validation required")
}
