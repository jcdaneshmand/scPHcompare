#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(10L, 11L))) stop(paste(
  "usage: run_mv08s_ph_sentinel.R <prefreeze> <mv08p-private>",
  "<mv08o-private> <hca-bm002-outs> <hca-remaining-root> <reference-rds>",
  "<common-panel> <private-root> <public-root> <expected-head> [--resume]"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
mv08p_private <- normalizePath(args[[2L]], mustWork = TRUE)
mv08o_private <- normalizePath(args[[3L]], mustWork = TRUE)
hca_bm002_outs <- normalizePath(args[[4L]], mustWork = TRUE)
hca_remaining_root <- normalizePath(args[[5L]], mustWork = TRUE)
reference_path <- normalizePath(args[[6L]], mustWork = TRUE)
panel_path <- normalizePath(args[[7L]], mustWork = TRUE)
private_root <- normalizePath(args[[8L]], mustWork = FALSE)
public_root <- normalizePath(args[[9L]], mustWork = FALSE)
expected_head <- tolower(trimws(args[[10L]]))
resume <- length(args) == 11L && identical(args[[11L]], "--resume")
if (length(args) == 11L && !resume) stop("unknown MV8-S runner flag", call. = FALSE)
observed_head <- tolower(trimws(Sys.getenv("MV08S_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", expected_head) ||
    !identical(observed_head, expected_head)) {
  stop("MV8-S exact committed HEAD binding failed", call. = FALSE)
}
if ((dir.exists(private_root) || dir.exists(public_root)) && !resume) {
  stop("MV8-S output roots already exist; explicit --resume required", call. = FALSE)
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(public_root, recursive = TRUE, showWarnings = FALSE)
for (subdir in c(
  "baseline", "repeat/baseline", "ph", "repeat/ph", "cross", "logs"
)) dir.create(file.path(private_root, subdir), recursive = TRUE, showWarnings = FALSE)
source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/mv07g_sentinel.R")
source("R/mv08s_ph_sentinel.R")

read_csv <- function(path) utils::read.csv(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
sha_file <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)
atomic_csv <- function(value, path) {
  partial <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) {
    unlink(partial); stop("failed to atomically publish ", basename(path), call. = FALSE)
  }
}
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(
    ps::ps_children(root, recursive = TRUE), error = function(error) list()
  ))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(error) 0
  ), numeric(1L)))
}
classify_stderr <- function(path, stage) {
  text <- if (file.exists(path)) paste(readLines(path, warn = FALSE), collapse = "\n") else ""
  if (!nzchar(trimws(text))) return("empty")
  known <- stage == "baseline" &&
    grepl("could not find glmGamPoi installed", text, fixed = TRUE) &&
    grepl("Falling back to native (slower) implementation", text, fixed = TRUE)
  if (known) "known_glmGamPoi_native_fallback" else "unexpected"
}
resolve_hca_outs <- function(unit_id) {
  if (unit_id == "HCA_BM_002") return(hca_bm002_outs)
  file.path(
    hca_remaining_root, unit_id,
    paste0("mv08h_exact500_", tolower(unit_id)), "outs"
  )
}
resolve_source_cache <- function(unit_id) {
  if (unit_id == "SRA628554_SRS2664364") return(file.path(
    mv08p_private, "cache", "internal",
    paste0(unit_id, "__exact500_allqc_sct_model.rds")
  ))
  if (unit_id == "HCA_BM_002") return(file.path(
    mv08o_private, "cache", paste0(unit_id, "__primary.rds")
  ))
  stop("MV8-S source cache request is outside the sentinel", call. = FALSE)
}

contract <- read_csv(file.path(prefreeze, "mv08s-contract.csv"))
bindings <- read_csv(file.path(prefreeze, "mv08s-external-input-bindings.csv"))
source_bindings <- read_csv(file.path(prefreeze, "mv08s-source-cache-bindings.csv"))
queue <- read_csv(file.path(prefreeze, "mv08s-ph-sentinel-queue.csv"))
cross_contract <- read_csv(file.path(prefreeze, "mv08s-cross-engine-contract.csv"))
if (nrow(contract) != 1L || contract$accepted_parent_head !=
      "10e0ac9281b0ca321e1c782ae86e483d6df277cc" ||
    nrow(bindings) != 8L || nrow(source_bindings) != 2L ||
    nrow(queue) != 23L || nrow(cross_contract) != 4L ||
    any(queue$authorization_state != "authorized_after_mv08s_commit") ||
    any(queue$workers != 1L) || any(queue$retries != 0L)) {
  stop("MV8-S committed execution contract drift", call. = FALSE)
}
for (index in seq_len(nrow(bindings))) {
  row <- bindings[index, , drop = FALSE]
  outs <- resolve_hca_outs(row$unit_id)
  filtered <- file.path(outs, "filtered_feature_bc_matrix.h5")
  raw <- file.path(outs, "raw_feature_bc_matrix.h5")
  if (!all(file.exists(c(filtered, raw))) ||
      sha_file(filtered) != row$filtered_h5_sha256 ||
      sha_file(raw) != row$raw_h5_sha256) {
    stop("MV8-S HCA input rehash failed: ", row$unit_id, call. = FALSE)
  }
}
for (index in seq_len(nrow(source_bindings))) {
  row <- source_bindings[index, , drop = FALSE]
  path <- resolve_source_cache(row$unit_id)
  if (!file.exists(path) || sha_file(path) != row$cache_sha256) {
    stop("MV8-S source-cache rehash failed: ", row$unit_id, call. = FALSE)
  }
}
if (sha_file(reference_path) != unique(bindings$reference_rds_sha256) ||
    sha_file(panel_path) != unique(bindings$panel_file_sha256)) {
  stop("MV8-S reference or panel rehash failed", call. = FALSE)
}

ledger_path <- file.path(public_root, "mv08s-resource-ledger.csv")
progress_path <- file.path(public_root, "mv08s-progress.csv")
ledger <- if (file.exists(ledger_path)) read_csv(ledger_path) else data.frame()
progress <- if (file.exists(progress_path)) read_csv(progress_path) else
  data.frame(
    contract_id = "mv08s_progress_v1", stage = "prefight",
    execution_head = expected_head,
    completed_units = 0L, total_units = 66L,
    last_unit = "none", disposition = "ready",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
publish_progress <- function(stage, completed, last, disposition) {
  progress <<- data.frame(
    contract_id = "mv08s_progress_v1", stage = stage,
    execution_head = expected_head,
    completed_units = completed, total_units = 66L,
    last_unit = last, disposition = disposition,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  atomic_csv(progress, progress_path)
}
run_child <- function(attempt_id, stage, script, script_args, output,
                      elapsed_cap, rss_cap) {
  hit <- if (nrow(ledger)) ledger[ledger$attempt_id == attempt_id, , drop = FALSE]
    else data.frame()
  if (nrow(hit)) {
    valid <- nrow(hit) == 1L && hit$disposition == "completed" &&
      file.exists(output) && hit$output_sha256 == sha_file(output) &&
      hit$output_bytes == as.numeric(file.info(output)$size)
    resumable_rss <- nrow(hit) == 1L && stage == "ph_primary" &&
      hit$disposition == "rss_cap_exceeded" && !file.exists(output)
    if (valid || resumable_rss) return(hit)
    stop("MV8-S prior attempt is not resumable: ", attempt_id, call. = FALSE)
  }
  if (file.exists(output)) {
    stop("MV8-S unowned output exists: ", basename(output), call. = FALSE)
  }
  stem <- gsub("[^A-Za-z0-9_.-]", "_", attempt_id)
  stdout <- file.path(private_root, "logs", paste0(stem, "__stdout.txt"))
  stderr <- file.path(private_root, "logs", paste0(stem, "__stderr.txt"))
  started <- Sys.time()
  process <- processx::process$new(
    Sys.which("Rscript"), c("--vanilla", script, script_args),
    stdout = stdout, stderr = stderr, cleanup_tree = TRUE
  )
  peak <- 0
  cap_failure <- ""
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
  stderr_class <- classify_stderr(stderr, if (grepl("baseline", stage)) "baseline" else stage)
  valid <- identical(status, 0L) && file.exists(output) &&
    stderr_class != "unexpected"
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid)
    "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08s_resource_metric_v1",
    execution_head = expected_head,
    attempt_order = nrow(ledger) + 1L,
    attempt_id = attempt_id,
    stage = stage,
    disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed,
    peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = elapsed_cap,
    rss_cap_bytes = rss_cap,
    output_file = basename(output),
    output_bytes = if (file.exists(output)) as.numeric(file.info(output)$size) else NA_real_,
    output_sha256 = if (file.exists(output)) sha_file(output) else NA_character_,
    stdout_bytes = if (file.exists(stdout)) as.numeric(file.info(stdout)$size) else 0,
    stderr_bytes = if (file.exists(stderr)) as.numeric(file.info(stderr)$size) else 0,
    stderr_class = stderr_class,
    workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  ledger <<- if (nrow(ledger)) rbind(ledger, metric) else metric
  atomic_csv(ledger, ledger_path)
  metric
}

completed_units <- 0L
baseline_repeat_rows <- list()
for (index in seq_len(nrow(bindings))) {
  binding <- bindings[index, , drop = FALSE]
  outs <- resolve_hca_outs(binding$unit_id)
  script_args <- c(
    prefreeze, binding$unit_id,
    file.path(outs, "filtered_feature_bc_matrix.h5"),
    file.path(outs, "raw_feature_bc_matrix.h5"),
    reference_path, panel_path
  )
  primary <- file.path(private_root, "baseline", paste0(binding$unit_id, ".rds"))
  repeated <- file.path(private_root, "repeat", "baseline",
                        paste0(binding$unit_id, ".rds"))
  first <- run_child(
    paste0("baseline__", binding$unit_id), "baseline",
    "scripts/run_mv08s_same_axis_baseline_entry.R",
    c(script_args, primary), primary, 1800, 12 * 1024^3
  )
  if (first$disposition != "completed") {
    stop("MV8-S baseline stopped: ", binding$unit_id, call. = FALSE)
  }
  completed_units <- completed_units + 1L
  publish_progress("baseline", completed_units, binding$unit_id, "completed")
  second <- run_child(
    paste0("baseline_repeat__", binding$unit_id), "baseline_repeat",
    "scripts/run_mv08s_same_axis_baseline_entry.R",
    c(script_args, repeated), repeated, 1800, 12 * 1024^3
  )
  exact <- second$disposition == "completed" &&
    first$output_bytes == second$output_bytes &&
    first$output_sha256 == second$output_sha256
  if (!exact) stop("MV8-S baseline repeat differs: ", binding$unit_id, call. = FALSE)
  completed_units <- completed_units + 1L
  publish_progress("baseline_repeat", completed_units, binding$unit_id, "completed")
  baseline_repeat_rows[[index]] <- data.frame(
    contract_id = "mv08s_baseline_repeat_validation_v1",
    unit_id = binding$unit_id,
    primary_sha256 = first$output_sha256,
    repeat_sha256 = second$output_sha256,
    bytes_equal = first$output_bytes == second$output_bytes,
    sha256_equal = first$output_sha256 == second$output_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}
baseline_repeats <- do.call(rbind, baseline_repeat_rows)

selected_rows <- list()
ph_repeat_rows <- list()
for (index in seq_len(nrow(queue))) {
  row <- queue[index, , drop = FALSE]
  source_path <- if (row$execution_role == "source_produced_gene_ph") {
    resolve_source_cache(row$unit_id)
  } else file.path(private_root, "baseline", paste0(row$unit_id, ".rds"))
  output <- file.path(private_root, row$output_file)
  primary <- run_child(
    paste0("ph_primary__", row$job_order), "ph_primary",
    "scripts/run_mv08s_ph_entry.R",
    c(prefreeze, row$job_id, source_path, panel_path, "ripserr", output),
    output, row$primary_elapsed_cap_seconds, row$primary_rss_cap_bytes
  )
  chosen <- primary
  engine <- "ripserr"
  if (primary$disposition != "completed") {
    if (row$view_kind != "gene_topology_v1" ||
        primary$disposition != "rss_cap_exceeded") {
      stop("MV8-S primary PH stopped: ", row$job_id, call. = FALSE)
    }
    chosen <- run_child(
      paste0("ph_fallback__", row$job_order), "ph_fallback",
      "scripts/run_mv08s_ph_entry.R",
      c(prefreeze, row$job_id, source_path, panel_path, "gudhi", output),
      output, row$fallback_elapsed_cap_seconds, row$fallback_rss_cap_bytes
    )
    engine <- "TDA_ripsDiag_GUDHI"
  }
  if (chosen$disposition != "completed") {
    stop("MV8-S selected PH attempt failed: ", row$job_id, call. = FALSE)
  }
  selected_rows[[index]] <- data.frame(
    contract_id = "mv08s_selected_ph_metric_v1",
    execution_head = expected_head,
    job_order = row$job_order, job_id = row$job_id, unit_id = row$unit_id,
    seed = row$seed, representation_id = row$representation_id,
    panel_id = row$panel_id, view_kind = row$view_kind,
    selected_engine = engine,
    elapsed_seconds = chosen$elapsed_seconds,
    peak_process_tree_rss_bytes = chosen$peak_process_tree_rss_bytes,
    output_file = row$output_file, output_bytes = chosen$output_bytes,
    output_sha256 = chosen$output_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  completed_units <- completed_units + 1L
  publish_progress("ph", completed_units, row$job_id, "completed")

  repeat_output <- file.path(private_root, "repeat", row$output_file)
  repeated <- run_child(
    paste0("ph_repeat__", row$job_order, "__", engine), "ph_repeat",
    "scripts/run_mv08s_ph_entry.R",
    c(prefreeze, row$job_id, source_path, panel_path,
      if (engine == "ripserr") "ripserr" else "gudhi", repeat_output),
    repeat_output,
    if (engine == "ripserr") row$primary_elapsed_cap_seconds else
      row$fallback_elapsed_cap_seconds,
    if (engine == "ripserr") row$primary_rss_cap_bytes else
      row$fallback_rss_cap_bytes
  )
  exact <- repeated$disposition == "completed" &&
    chosen$output_bytes == repeated$output_bytes &&
    chosen$output_sha256 == repeated$output_sha256
  if (!exact) stop("MV8-S PH repeat differs: ", row$job_id, call. = FALSE)
  ph_repeat_rows[[index]] <- data.frame(
    contract_id = "mv08s_ph_repeat_validation_v1",
    job_order = row$job_order, job_id = row$job_id,
    selected_engine = engine,
    primary_sha256 = chosen$output_sha256,
    repeat_sha256 = repeated$output_sha256,
    bytes_equal = chosen$output_bytes == repeated$output_bytes,
    sha256_equal = chosen$output_sha256 == repeated$output_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  completed_units <- completed_units + 1L
  publish_progress("ph_repeat", completed_units, row$job_id, "completed")
}
selected <- do.call(rbind, selected_rows)
ph_repeats <- do.call(rbind, ph_repeat_rows)

cross_rows <- list()
for (index in seq_len(nrow(cross_contract))) {
  spec <- cross_contract[index, , drop = FALSE]
  row <- queue[queue$job_id == spec$job_id, , drop = FALSE]
  source_path <- if (row$execution_role == "source_produced_gene_ph") {
    resolve_source_cache(row$unit_id)
  } else file.path(private_root, "baseline", paste0(row$unit_id, ".rds"))
  output <- file.path(private_root, "cross", paste0(spec$cross_check_id, ".csv"))
  metric <- run_child(
    paste0("cross__", spec$cross_check_id), "cross_engine",
    "scripts/run_mv08s_cross_engine_entry.R",
    c(prefreeze, spec$cross_check_id, source_path, panel_path, spec$mode, output),
    output, spec$elapsed_cap_seconds, spec$rss_cap_bytes
  )
  if (metric$disposition != "completed") {
    stop("MV8-S cross-engine check failed: ", spec$cross_check_id, call. = FALSE)
  }
  cross_rows[[index]] <- read_csv(output)
  completed_units <- completed_units + 1L
  publish_progress("cross_engine", completed_units, spec$cross_check_id, "completed")
}
cross_results <- do.call(rbind, cross_rows)

baseline_metrics <- do.call(rbind, lapply(seq_len(nrow(bindings)), function(index) {
  binding <- bindings[index, , drop = FALSE]
  path <- file.path(private_root, "baseline", paste0(binding$unit_id, ".rds"))
  record <- readRDS(path); mv08s_validate_baseline_record_v1(record, binding)
  metric <- ledger[ledger$attempt_id == paste0("baseline__", binding$unit_id),,
                   drop = FALSE]
  data.frame(
    contract_id = "mv08s_baseline_metric_v1", unit_id = binding$unit_id,
    selected_cells = record$identity$selected_cells,
    selected_cell_sha256 = record$identity$selected_cell_sha256,
    panel_genes = 475L, panel_sha256 = record$identity$panel_sha256,
    cell_view_payload_sha256 = record$identity$cell_view_payload_sha256,
    gene_view_payload_sha256 = record$identity$gene_view_payload_sha256,
    elapsed_seconds = metric$elapsed_seconds,
    peak_process_tree_rss_bytes = metric$peak_process_tree_rss_bytes,
    output_bytes = metric$output_bytes, output_sha256 = metric$output_sha256,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))
ph_metrics <- do.call(rbind, lapply(seq_len(nrow(selected)), function(index) {
  metric <- selected[index, , drop = FALSE]
  path <- file.path(private_root, metric$output_file)
  record <- readRDS(path)
  mv08s_validate_ph_record_v1(record)
  data.frame(
    metric,
    point_count = record$h0_mst_oracle$point_count,
    finite_h0_intervals = record$h0_mst_oracle$finite_h0_intervals,
    finite_h1_intervals = record$h0_mst_oracle$finite_h1_intervals,
    h0_mst_maximum_absolute_error =
      record$h0_mst_oracle$maximum_absolute_error,
    h0_mst_tolerance = record$h0_mst_oracle$tolerance,
    h0_mst_passed = record$h0_mst_oracle$passed,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    ph_cache_key = record$cache_key,
    stringsAsFactors = FALSE
  )
}))
private_files <- list.files(private_root, recursive = TRUE, full.names = TRUE)
private_files <- private_files[file.info(private_files)$isdir %in% FALSE]
private_bytes <- sum(as.numeric(file.info(private_files)$size))
aggregate_elapsed <- sum(ledger$elapsed_seconds)
decision <- data.frame(
  contract_id = "mv08s_execution_decision_v1",
  execution_head = expected_head,
  decision = "23_record_PH_sentinel_complete_await_independent_closure",
  baseline_units = nrow(baseline_metrics), baseline_repeat_units = nrow(baseline_repeats),
  ph_records = nrow(ph_metrics), ph_repeat_records = nrow(ph_repeats),
  ripserr_selected_records = sum(ph_metrics$selected_engine == "ripserr"),
  gudhi_fallback_records = sum(ph_metrics$selected_engine == "TDA_ripsDiag_GUDHI"),
  cross_engine_checks = nrow(cross_results),
  all_mst_oracles_passed = all(ph_metrics$h0_mst_passed),
  all_repeats_byte_identical = all(baseline_repeats$sha256_equal) &&
    all(ph_repeats$sha256_equal),
  all_cross_engine_checks_passed = all(cross_results$passed),
  all_resource_caps_passed = all(
    (ledger$disposition == "completed" &
       ledger$peak_process_tree_rss_bytes <= ledger$rss_cap_bytes &
       ledger$elapsed_seconds <= ledger$elapsed_cap_seconds + 5) |
    (ledger$stage == "ph_primary" &
       ledger$disposition == "rss_cap_exceeded" &
       ledger$elapsed_seconds <= ledger$elapsed_cap_seconds + 5 &
       is.na(ledger$output_sha256))
  ),
  unexpected_stderr_records = sum(ledger$stderr_class == "unexpected"),
  aggregate_elapsed_seconds = aggregate_elapsed,
  aggregate_elapsed_cap_seconds = contract$aggregate_elapsed_cap_seconds,
  private_bytes = private_bytes,
  private_storage_cap_bytes = contract$private_storage_cap_bytes,
  full_ph_jobs_authorized = 0L, landscape_groups_authorized = 0L,
  comparison_strata_authorized = 0L, clustering_jobs_authorized = 0L,
  fusion_jobs_authorized = 0L, label_jobs_authorized = 0L,
  outcome_jobs_authorized = 0L, outcome_label_state = "closed",
  biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
)
if (completed_units != 66L || nrow(ph_metrics) != 23L ||
    nrow(cross_results) != 8L || !all(ph_metrics$h0_mst_passed) ||
    !all(baseline_repeats$sha256_equal) || !all(ph_repeats$sha256_equal) ||
    !all(cross_results$passed) || !decision$all_resource_caps_passed ||
    decision$unexpected_stderr_records != 0L ||
    aggregate_elapsed > contract$aggregate_elapsed_cap_seconds ||
    private_bytes > contract$private_storage_cap_bytes) {
  stop("MV8-S terminal execution gates failed", call. = FALSE)
}
atomic_csv(baseline_metrics, file.path(public_root, "mv08s-baseline-metrics.csv"))
atomic_csv(baseline_repeats, file.path(public_root, "mv08s-baseline-repeat-validation.csv"))
atomic_csv(ph_metrics, file.path(public_root, "mv08s-ph-metrics.csv"))
atomic_csv(ph_repeats, file.path(public_root, "mv08s-ph-repeat-validation.csv"))
atomic_csv(cross_results, file.path(public_root, "mv08s-cross-engine-results.csv"))
atomic_csv(decision, file.path(public_root, "mv08s-execution-decision.csv"))
publish_progress("complete", 66L, "mv08s", "complete_await_independent_closure")
message("MV8-S complete: eight baseline units, 23 PH records, 23 repeats, four cross-engine jobs")
