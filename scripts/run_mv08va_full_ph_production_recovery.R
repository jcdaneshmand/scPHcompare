#!/usr/bin/env Rscript

# Resume MV8-V from MV8-VA's independently validated, byte-preserved job-1
# prefix. The original stopped roots remain immutable; this runner writes only
# fresh recovery roots and starts new PH work at production order 2.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("digest", "processx", "ps")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " required", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(9L, 10L))) stop(paste(
  "usage: run_mv08va_full_ph_production_recovery.R <mv08u-prefreeze>",
  "<mv08p-private> <mv08pr-private> <mv08ps-private>",
  "<mv08o-internal-private> <common-panel> <private-output>",
  "<public-output> <expected-head> [--resume]"
), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
source_roots <- list(
  mv08p_original_v1 = normalizePath(args[[2L]], mustWork = TRUE),
  mv08pr_overlay_v1 = normalizePath(args[[3L]], mustWork = TRUE),
  mv08ps_overlay_v1 = normalizePath(args[[4L]], mustWork = TRUE),
  mv08o_internal_primary_v8 = normalizePath(args[[5L]], mustWork = TRUE)
)
panel_path <- normalizePath(args[[6L]], mustWork = TRUE)
private_root <- normalizePath(args[[7L]], mustWork = FALSE)
public_root <- normalizePath(args[[8L]], mustWork = FALSE)
expected_head <- tolower(trimws(args[[9L]]))
resume <- length(args) == 10L && identical(args[[10L]], "--resume")
if (length(args) == 10L && !resume) stop("unknown MV8-V runner flag", call. = FALSE)
environment_head <- tolower(trimws(Sys.getenv("MV08V_GIT_HEAD", unset = "")))
if (!grepl("^[0-9a-f]{40}$", expected_head) ||
    !identical(expected_head, environment_head)) {
  stop("MV8-V exact committed HEAD binding failed", call. = FALSE)
}
if ((dir.exists(private_root) || dir.exists(public_root)) && !resume) {
  stop("MV8-V output roots exist; explicit --resume required", call. = FALSE)
}
if (resume && (!dir.exists(private_root) || !dir.exists(public_root))) {
  stop("MV8-V resume requires both existing roots", call. = FALSE)
}
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(public_root, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("ph", "logs")) {
  dir.create(file.path(private_root, subdir), recursive = TRUE,
             showWarnings = FALSE)
}
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
  partial <- paste0(path, ".partial")
  utils::write.csv(value, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
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
classify_stderr <- function(path) {
  text <- if (file.exists(path)) paste(readLines(path, warn = FALSE), collapse = "\n") else ""
  if (!nzchar(trimws(text))) "empty" else "unexpected"
}
private_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) 0 else sum(as.numeric(file.info(files)$size))
}
memory_available <- function() {
  line <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
  as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1L]][[2L]]) * 1024
}
disk_available <- function(path) suppressWarnings(as.numeric(system2(
  "df", c("-B1", "--output=avail", path), stdout = TRUE
)[[2L]]))
resolve_source <- function(row) {
  role <- row$source_cache_root_role
  if (!(role %in% names(source_roots))) {
    stop("MV8-V unknown source root role: ", role, call. = FALSE)
  }
  file.path(source_roots[[role]], row$source_cache_relative_file)
}

contract <- read_csv(file.path(prefreeze, "mv08u-contract.csv"))
queue <- read_csv(file.path(prefreeze, "mv08u-full-ph-queue.csv"))
gate <- read_csv(file.path(prefreeze, "mv08u-launch-resource-gate.csv"))
runtime_inputs <- read_csv(file.path(prefreeze, "mv08u-runtime-input-bindings.csv"))
implementation <- read_csv(file.path(prefreeze, "mv08u-implementation-bindings.csv"))
input_manifest <- read_csv(file.path(prefreeze, "mv08u-input-manifest.csv"))
decision <- read_csv(file.path(prefreeze, "mv08u-decision.csv"))
recovery_prefreeze <- normalizePath(
  Sys.getenv("MV08VA_RECOVERY_PREFREEZE", unset = ""), mustWork = TRUE
)
recovery_decision <- read_csv(file.path(recovery_prefreeze, "mv08va-decision.csv"))
recovery_implementation <- read_csv(file.path(
  recovery_prefreeze, "mv08va-implementation-bindings.csv"
))
execution_amendment <- normalizePath(
  Sys.getenv("MV08VD_RECOVERY_PREFREEZE", unset = ""), mustWork = TRUE
)
amendment_decision <- read_csv(file.path(execution_amendment, "mv08vd-decision.csv"))
amendment_binding <- read_csv(file.path(
  execution_amendment, "mv08vd-implementation-bindings.csv"
))
amendment_manifest <- read_csv(file.path(
  execution_amendment, "mv08vd-artifact-manifest.csv"
))
amended_files <- c(
  "scripts/bootstrap_mv08va_full_ph_recovery.R",
  "scripts/run_mv08va_full_ph_production_recovery.R"
)
historical_unamended <- !(recovery_implementation$file %in% amended_files)
if (nrow(recovery_decision) != 1L ||
    recovery_decision$decision !=
      "authorize_no_retry_job1_bootstrap_and_resume_at_job2" ||
    recovery_decision$accepted_completed_records != 1L ||
    recovery_decision$retry_records_authorized != 0L ||
    !all(file.exists(recovery_implementation$file[historical_unamended])) ||
    !all(vapply(recovery_implementation$file[historical_unamended], sha_file,
                character(1L)) ==
           recovery_implementation$sha256[historical_unamended]) ||
    nrow(amendment_decision) != 1L ||
    amendment_decision$decision != "authorize_amendment_bound_resume_at_job2" ||
    !amendment_decision$runner_resume_authorized ||
    nrow(amendment_binding) != length(amended_files) ||
    !setequal(amendment_binding$file, amended_files) ||
    !all(amendment_binding$mv08va_sha256 == recovery_implementation$sha256[
      match(amendment_binding$file, recovery_implementation$file)
    ]) ||
    !all(vapply(amendment_binding$file, sha_file, character(1L)) ==
           amendment_binding$sha256) ||
    !all(vapply(file.path(execution_amendment, amendment_manifest$artifact),
                sha_file, character(1L)) == amendment_manifest$sha256)) {
  stop("MV8-VA committed recovery binding failed", call. = FALSE)
}
if (nrow(contract) != 1L || nrow(queue) != 1257L || nrow(decision) != 1L ||
    decision$ph_records_authorized != 1257L ||
    any(queue$authorization_state != "authorized_after_mv08u_commit") ||
    any(queue$workers != 1L | queue$retries != 0L) ||
    contract$execution_head_state != "bind_exact_committed_head_at_launch_and_publish") {
  stop("MV8-V committed authorization drift", call. = FALSE)
}
if (!all(file.exists(implementation$file)) ||
    !all(vapply(implementation$file, sha_file, character(1L)) ==
         implementation$sha256) ||
    !all(file.exists(input_manifest$input)) ||
    !all(vapply(input_manifest$input, sha_file, character(1L)) ==
         input_manifest$sha256)) {
  stop("MV8-V committed implementation or public input drift", call. = FALSE)
}
if (memory_available() < gate$minimum_bytes_at_launch[gate$resource == "available_memory"] ||
    disk_available(private_root) < gate$minimum_bytes_at_launch[
      gate$resource == "workspace_filesystem_available"]) {
  stop("MV8-V launch resource gate not met", call. = FALSE)
}
panel <- read_csv(panel_path)
if (nrow(runtime_inputs) != 1L ||
    sha_file(panel_path) != runtime_inputs$file_sha256 ||
    nrow(panel) != 475L ||
    !identical(as.integer(panel$panel_order_475), 1:475)) {
  stop("MV8-V common475 panel drift", call. = FALSE)
}

# Rehash every unique source once before starting or resuming any PH.
unit_index <- !duplicated(queue$unit_id)
for (i in which(unit_index)) {
  path <- resolve_source(queue[i, , drop = FALSE])
  if (!file.exists(path) || sha_file(path) != queue$source_cache_sha256[[i]]) {
    stop("MV8-V source cache rehash failed: ", queue$unit_id[[i]], call. = FALSE)
  }
}

ledger_path <- file.path(public_root, "mv08v-resource-ledger.csv")
metrics_path <- file.path(public_root, "mv08v-selected-ph-metrics.csv")
progress_path <- file.path(public_root, "mv08v-progress.csv")
ledger <- if (file.exists(ledger_path)) read_csv(ledger_path) else data.frame()
metrics <- if (file.exists(metrics_path)) read_csv(metrics_path) else data.frame()
if (resume && nrow(metrics)) {
  if (!identical(as.integer(metrics$production_order), seq_len(nrow(metrics))) ||
      any(metrics$disposition != "completed") ||
      any(metrics$job_id != queue$job_id[seq_len(nrow(metrics))])) {
    stop("MV8-V resume metrics are not a completed strict prefix", call. = FALSE)
  }
  ledger_orders <- as.integer(sub(
    "^ph_(primary|fallback)__", "", ledger$attempt_id
  ))
  if (!nrow(ledger) || anyNA(ledger_orders) ||
      any(ledger_orders > nrow(metrics)) ||
      any(!ledger$stage %in% c("ph_primary", "ph_fallback"))) {
    stop("MV8-V resume ledger extends beyond the completed strict prefix",
         call. = FALSE)
  }
  for (i in seq_len(nrow(metrics))) {
    output <- file.path(private_root, metrics$output_file[[i]])
    if (!file.exists(output) || sha_file(output) != metrics$output_sha256[[i]] ||
        as.numeric(file.info(output)$size) != metrics$output_bytes[[i]]) {
      stop("MV8-V completed-prefix artifact drift at ", i, call. = FALSE)
    }
    primary_row <- ledger[ledger$attempt_id == paste0("ph_primary__", i),,
                          drop = FALSE]
    fallback_row <- ledger[ledger$attempt_id == paste0("ph_fallback__", i),,
                           drop = FALSE]
    if (metrics$selected_engine[[i]] == "ripserr") {
      valid_attempt <- nrow(primary_row) == 1L &&
        primary_row$disposition == "completed" && !nrow(fallback_row) &&
        primary_row$output_sha256 == metrics$output_sha256[[i]]
    } else {
      valid_attempt <- nrow(primary_row) == 1L &&
        primary_row$disposition == "rss_cap_exceeded" &&
        nrow(fallback_row) == 1L && fallback_row$disposition == "completed" &&
        fallback_row$output_sha256 == metrics$output_sha256[[i]]
    }
    if (!valid_attempt) {
      stop("MV8-V completed-prefix attempt drift at ", i, call. = FALSE)
    }
  }
}
if (!resume && (file.exists(ledger_path) || file.exists(metrics_path) ||
                file.exists(progress_path))) {
  stop("MV8-V public output unexpectedly pre-exists", call. = FALSE)
}

publish_progress <- function(state, last_order, last_job) {
  value <- data.frame(
    contract_id = "mv08v_progress_v1", execution_head = expected_head,
    state = state, completed_records = nrow(metrics), total_records = 1257L,
    last_production_order = last_order, last_job_id = last_job,
    aggregate_attempt_seconds = if (nrow(ledger)) sum(ledger$elapsed_seconds) else 0,
    private_bytes = private_bytes(), workers = 1L, retries = 0L,
    landscape_records = 0L, comparison_records = 0L,
    clustering_records = 0L, fusion_records = 0L, label_records = 0L,
    biological_outcome_records = 0L, outcome_label_state = "closed",
    biological_outcomes_computed = FALSE, stringsAsFactors = FALSE
  )
  atomic_csv(value, progress_path)
}
run_child <- function(attempt_id, stage, script_args, output,
                      elapsed_cap, rss_cap) {
  if (file.exists(output) || file.exists(paste0(output, ".partial"))) {
    stop("MV8-V unowned or partial output exists: ", basename(output), call. = FALSE)
  }
  stem <- gsub("[^A-Za-z0-9_.-]", "_", attempt_id)
  stdout <- file.path(private_root, "logs", paste0(stem, "__stdout.txt"))
  stderr <- file.path(private_root, "logs", paste0(stem, "__stderr.txt"))
  if (any(file.exists(c(stdout, stderr)))) {
    stop("MV8-V ambiguous prior attempt logs require recovery audit", call. = FALSE)
  }
  started <- Sys.time()
  process <- processx::process$new(
    Sys.which("Rscript"), c("--vanilla", "scripts/run_mv08v_ph_entry.R", script_args),
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
  stderr_class <- classify_stderr(stderr)
  valid <- identical(status, 0L) && file.exists(output) &&
    !file.exists(paste0(output, ".partial")) && stderr_class == "empty"
  disposition <- if (nzchar(cap_failure)) cap_failure else if (valid) "completed" else "failed"
  metric <- data.frame(
    contract_id = "mv08v_resource_metric_v1", execution_head = expected_head,
    attempt_order = nrow(ledger) + 1L, attempt_id = attempt_id, stage = stage,
    disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    elapsed_cap_seconds = elapsed_cap, rss_cap_bytes = rss_cap,
    output_file = gsub("\\\\", "/", substring(
      output, nchar(private_root) + 2L
    )),
    output_bytes = if (file.exists(output)) as.numeric(file.info(output)$size) else NA_real_,
    output_sha256 = if (file.exists(output)) sha_file(output) else NA_character_,
    stdout_bytes = if (file.exists(stdout)) as.numeric(file.info(stdout)$size) else 0,
    stderr_bytes = if (file.exists(stderr)) as.numeric(file.info(stderr)$size) else 0,
    stderr_class = stderr_class, workers = 1L, retries = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  ledger <<- if (nrow(ledger)) rbind(ledger, metric) else metric
  atomic_csv(ledger, ledger_path)
  metric
}

start <- nrow(metrics) + 1L
publish_progress("running", nrow(metrics), if (nrow(metrics)) tail(metrics$job_id, 1L) else "none")
for (i in seq.int(start, nrow(queue))) {
  if ((nrow(ledger) && sum(ledger$elapsed_seconds) >= contract$aggregate_elapsed_cap_seconds) ||
      private_bytes() >= contract$private_storage_cap_bytes) {
    publish_progress("stopped_before_next_job_resource_cap", i - 1L,
                     if (i > 1L) queue$job_id[[i - 1L]] else "none")
    stop("MV8-V aggregate resource cap reached; evidence preserved", call. = FALSE)
  }
  row <- queue[i, , drop = FALSE]
  source_path <- resolve_source(row)
  output <- file.path(private_root, row$output_file)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  primary <- run_child(
    paste0("ph_primary__", i), "ph_primary",
    c(prefreeze, row$job_id, source_path, panel_path, "ripserr", output),
    output, row$primary_elapsed_cap_seconds, row$primary_rss_cap_bytes
  )
  chosen <- primary
  selected_engine <- "ripserr"
  if (primary$disposition != "completed") {
    if (primary$disposition != "rss_cap_exceeded" ||
        row$fallback_trigger != "rss_cap_exceeded_only") {
      publish_progress("stopped", i - 1L, row$job_id)
      stop("MV8-V primary stopped without authorized fallback: ", row$job_id,
           call. = FALSE)
    }
    if (file.exists(output) || file.exists(paste0(output, ".partial"))) {
      publish_progress("stopped_ambiguous_primary_output", i - 1L, row$job_id)
      stop("MV8-V RSS stop left an artifact; recovery audit required", call. = FALSE)
    }
    chosen <- run_child(
      paste0("ph_fallback__", i), "ph_fallback",
      c(prefreeze, row$job_id, source_path, panel_path, "gudhi", output),
      output, row$fallback_elapsed_cap_seconds, row$fallback_rss_cap_bytes
    )
    selected_engine <- "TDA_ripsDiag_GUDHI"
  }
  if (chosen$disposition != "completed") {
    publish_progress("stopped", i - 1L, row$job_id)
    stop("MV8-V selected PH attempt failed: ", row$job_id, call. = FALSE)
  }
  record <- readRDS(output)
  mv08s_validate_ph_record_v1(record, row)
  metric <- data.frame(
    contract_id = "mv08v_selected_ph_metric_v1", execution_head = expected_head,
    production_order = i, mv08r_job_order = row$mv08r_job_order,
    job_id = row$job_id, unit_id = row$unit_id, seed = row$seed,
    representation_id = row$representation_id, panel_id = row$panel_id,
    view_kind = row$view_kind, selected_engine = selected_engine,
    disposition = "completed", elapsed_seconds = chosen$elapsed_seconds,
    peak_process_tree_rss_bytes = chosen$peak_process_tree_rss_bytes,
    output_file = row$output_file, output_bytes = chosen$output_bytes,
    output_sha256 = chosen$output_sha256,
    point_count = record$h0_mst_oracle$point_count,
    finite_h0_intervals = record$h0_mst_oracle$finite_h0_intervals,
    finite_h1_intervals = record$h0_mst_oracle$finite_h1_intervals,
    h0_mst_passed = record$h0_mst_oracle$passed,
    diagram_sha256 = record$topology_result$provenance$diagram_sha256,
    ph_cache_key = record$cache_key, landscape_records = 0L,
    comparison_records = 0L, clustering_records = 0L, fusion_records = 0L,
    label_records = 0L, biological_outcome_records = 0L,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  metrics <- if (nrow(metrics)) rbind(metrics, metric) else metric
  atomic_csv(metrics, metrics_path)
  publish_progress("running", i, row$job_id)
}
publish_progress("ph_production_complete_closure_pending", 1257L,
                 tail(queue$job_id, 1L))
cat("MV8-V full PH production completed 1257/1257; MV8-W closure required\n")
