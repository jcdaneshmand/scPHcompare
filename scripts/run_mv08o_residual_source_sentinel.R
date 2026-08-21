#!/usr/bin/env Rscript

# Serial, cap-enforcing MV8-O orchestrator. It runs exactly the three authorized
# sources plus required repeats, writes private caches, and publishes only
# aggregate source/view/resource evidence. It does not run PH or landscapes.

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% 6:7)) stop(paste(
  "usage: run_mv08o_residual_source_sentinel.R <mv08n-audit-dir> <primary-raw-dir>",
  "<added-raw-dir> <hca-bm002.h5> <private-output-dir> <public-output-dir> [--single-unit=<unit-id>|--single-run=<unit-id>:<role>]"), call. = FALSE)
audit_dir <- normalizePath(args[[1L]], mustWork = TRUE)
primary_raw <- normalizePath(args[[2L]], mustWork = TRUE)
added_raw <- normalizePath(args[[3L]], mustWork = TRUE)
hca_path <- normalizePath(args[[4L]], mustWork = TRUE)
private_dir <- normalizePath(args[[5L]], mustWork = FALSE)
public_dir <- normalizePath(args[[6L]], mustWork = FALSE)
single_unit <- NULL; single_role <- NULL
if (length(args) == 7L && startsWith(args[[7L]], "--single-unit=")) {
  single_unit <- sub("^--single-unit=", "", args[[7L]]); single_role <- "primary"
} else if (length(args) == 7L && startsWith(args[[7L]], "--single-run=")) {
  pieces <- strsplit(sub("^--single-run=", "", args[[7L]]), ":", fixed = TRUE)[[1L]]
  if (length(pieces) != 2L || !nzchar(pieces[[1L]]) || !(pieces[[2L]] %in% c("primary", "repeat"))) {
    stop("--single-run must be --single-run=<unit-id>:<primary|repeat>", call. = FALSE)
  }
  single_unit <- pieces[[1L]]; single_role <- pieces[[2L]]
} else if (length(args) == 7L) {
  stop("optional MV8-O argument must select one authorized source run", call. = FALSE)
}
if (dir.exists(private_dir) || dir.exists(public_dir)) stop("refusing to overwrite MV8-O output roots", call. = FALSE)
dir.create(private_dir, recursive = TRUE); dir.create(public_dir, recursive = TRUE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(warn = 2)
for (package in c("processx", "ps", "digest")) if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required")
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
classify_stderr <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!nzchar(text)) return("empty")
  known <- grepl("vst.flavor.*glmGamPoi", text) &&
    grepl("Falling back to native \\(slower\\) implementation", text) &&
    !grepl("(^|\\n)(Error|Execution halted)", text)
  if (known) "known_glmGamPoi_native_fallback" else "unexpected_nonempty"
}
atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial"); utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE), error = function(e) list()))
  sum(vapply(handles, function(h) tryCatch(as.numeric(ps::ps_memory_info(h)[["rss"]]), error = function(e) 0), numeric(1L)))
}
queue <- utils::read.csv(file.path(audit_dir, "mv08n-residual-source-queue.csv"), check.names = FALSE,
                         stringsAsFactors = FALSE)
jobs <- queue[queue$authorization_state == "source_view_sentinel_authorized", , drop = FALSE]
jobs <- jobs[match(c("internal_minimum_cell_sentinel", "internal_maximum_cell_sentinel", "external_mv08m_bridge_sentinel"),
                    jobs$sentinel_role), , drop = FALSE]
if (nrow(jobs) != 3L || anyNA(jobs$unit_id) || !all(jobs$elapsed_cap_seconds == 1800L) ||
    !all(jobs$rss_cap_bytes == 12 * 1024^3) || !all(jobs$workers == 1L) || !all(jobs$retries == 0L)) {
  stop("MV8-N source sentinel queue drift", call. = FALSE)
}
run_jobs <- rbind(jobs, jobs[jobs$repeat_required, , drop = FALSE])
run_jobs$run_role <- c(rep("primary", nrow(jobs)), rep("repeat", sum(jobs$repeat_required)))
if (nrow(run_jobs) != 5L) stop("expected three primary plus two repeat jobs", call. = FALSE)
if (!is.null(single_unit)) {
  run_jobs <- run_jobs[run_jobs$unit_id == single_unit & run_jobs$run_role == single_role, , drop = FALSE]
  if (nrow(run_jobs) != 1L) stop("requested MV8-O source run is not authorized", call. = FALSE)
}
worker <- normalizePath("scripts/run_mv08o_residual_source_worker.R", mustWork = TRUE)
rows <- vector("list", nrow(run_jobs))
for (i in seq_len(nrow(run_jobs))) {
  job <- run_jobs[i, , drop = FALSE]
  is_hca <- job$dataset_scope == "external8"
  source_kind <- if (is_hca) "h5" else "raw"
  source_path <- if (is_hca) hca_path else file.path(if (job$source_tier == "primary90") primary_raw else added_raw,
                                                       paste0(job$unit_id, "__raw.rds"))
  if (!file.exists(source_path)) stop("required source is absent for ", job$unit_id, call. = FALSE)
  stem <- paste0(job$unit_id, "__", job$run_role)
  cache_path <- file.path(private_dir, "cache", paste0(stem, ".rds"))
  worker_audit <- file.path(private_dir, "worker-audit", paste0(stem, ".csv"))
  stdout_path <- file.path(private_dir, "logs", paste0(stem, "-stdout.txt"))
  stderr_path <- file.path(private_dir, "logs", paste0(stem, "-stderr.txt"))
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(worker_audit), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(stdout_path), recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time(); peak <- 0
  proc <- processx::process$new(Sys.which("Rscript"), c("--vanilla", worker, audit_dir, source_kind,
    source_path, job$unit_id, cache_path, worker_audit, job$run_role), stdout = stdout_path,
    stderr = stderr_path, cleanup_tree = TRUE)
  disposition <- "running"
  while (proc$is_alive()) {
    Sys.sleep(0.25); elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs")); peak <- max(peak, tree_rss(proc$get_pid()))
    if (elapsed > job$elapsed_cap_seconds || peak > job$rss_cap_bytes) {
      disposition <- if (elapsed > job$elapsed_cap_seconds) "elapsed_cap_exceeded" else "rss_cap_exceeded"
      proc$kill_tree(); break
    }
  }
  proc$wait(timeout = 10000); elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs")); status <- proc$get_exit_status()
  if (disposition == "running") disposition <- if (identical(status, 0L)) "completed" else "child_failed"
  stderr_bytes <- if (file.exists(stderr_path)) file.info(stderr_path)$size else NA_real_
  stderr_class <- classify_stderr(stderr_path)
  audit <- if (disposition == "completed" && file.exists(worker_audit) && file.exists(cache_path))
    utils::read.csv(worker_audit, check.names = FALSE, stringsAsFactors = FALSE) else NULL
  rows[[i]] <- data.frame(contract_id = "mv08o_source_sentinel_resource_v1", job_order = i,
    unit_id = job$unit_id, dataset_scope = job$dataset_scope, sentinel_role = job$sentinel_role,
    run_role = job$run_role, disposition = disposition, exit_status = status,
    elapsed_seconds = elapsed, peak_rss_bytes = peak, elapsed_cap_seconds = job$elapsed_cap_seconds,
    rss_cap_bytes = job$rss_cap_bytes, cache_bytes = if (file.exists(cache_path)) file.info(cache_path)$size else NA_real_,
    cache_sha256 = if (file.exists(cache_path)) sha_file(cache_path) else NA_character_,
    worker_audit_sha256 = if (file.exists(worker_audit)) sha_file(worker_audit) else NA_character_,
    stderr_bytes = stderr_bytes, stderr_class = stderr_class, worker_rows = if (is.null(audit)) 0L else nrow(audit),
    all_geometry_valid = if (is.null(audit)) FALSE else {
      required <- audit$dataset_scope == "internal124" |
        audit$panel_id == "common475" |
        (audit$panel_id == "exact500" &
          audit$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
      all(audit$values_finite[required]) &&
        all(audit$zero_variance_gene_count[required] == 0L) &&
        all(audit$correlation_chord_valid[required])
    },
    persistence_computed = FALSE, landscapes_computed = FALSE, clustering_computed = FALSE,
    fusion_computed = FALSE, outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE)
  message("MV8-O ", i, "/", nrow(run_jobs), " ", job$unit_id, " ", job$run_role, ": ", disposition)
  if (disposition != "completed") break
}
resource <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
atomic_csv(resource, file.path(public_dir, "mv08o-source-sentinel-resource.csv"))
if (any(resource$disposition != "completed") ||
    any(resource$elapsed_seconds > resource$elapsed_cap_seconds) || any(resource$peak_rss_bytes > resource$rss_cap_bytes) ||
    any(!resource$all_geometry_valid) || any(resource$worker_rows != ifelse(resource$dataset_scope == "internal124", 10L, 4L)) ||
    (is.null(single_unit) && (nrow(resource) != 5L ||
      any(!resource$stderr_class %in% c("empty", "known_glmGamPoi_native_fallback"))))) {
  stop("MV8-O source sentinel did not pass; partial artifacts preserved", call. = FALSE)
}
if (!is.null(single_unit)) {
  cat("MV8-O single-source resource check passed for ", single_unit, " ", single_role,
      "; global closure remains pending\n", sep = "")
  quit(status = 0L)
}
audits <- lapply(seq_len(nrow(resource)), function(i) utils::read.csv(file.path(private_dir, "worker-audit",
  paste0(resource$unit_id[[i]], "__", resource$run_role[[i]], ".csv")), check.names = FALSE, stringsAsFactors = FALSE))
views <- do.call(rbind, audits)
atomic_csv(views, file.path(public_dir, "mv08o-source-sentinel-view-summary.csv"))
primary <- views[views$run_role == "primary", , drop = FALSE]
repeat_views <- views[views$run_role == "repeat", , drop = FALSE]
repeat_rows <- merge(primary[primary$unit_id %in% unique(resource$unit_id[resource$run_role == "repeat"]),
                                      c("unit_id", "seed", "representation_id", "panel_id", "matrix_sha256", "distance_sha256")],
                     repeat_views[, c("unit_id", "seed", "representation_id", "panel_id", "matrix_sha256", "distance_sha256")],
                     by = c("unit_id", "seed", "representation_id", "panel_id"), suffixes = c("_primary", "_repeat"))
determinism <- data.frame(contract_id = "mv08o_source_sentinel_determinism_v1", repeat_rows,
  passed = repeat_rows$matrix_sha256_primary == repeat_rows$matrix_sha256_repeat &
    ((is.na(repeat_rows$distance_sha256_primary) & is.na(repeat_rows$distance_sha256_repeat)) |
      repeat_rows$distance_sha256_primary == repeat_rows$distance_sha256_repeat),
  stringsAsFactors = FALSE)
atomic_csv(determinism, file.path(public_dir, "mv08o-source-sentinel-determinism.csv"))
if (nrow(primary) != 24L || nrow(repeat_views) != 14L || nrow(determinism) != 14L || !all(determinism$passed)) {
  stop("MV8-O deterministic repeat gate failed", call. = FALSE)
}
validation <- data.frame(check_id = c("authorized_units", "five_runs_complete", "resource_caps", "stderr_policy",
  "exact500_geometry", "internal_five_seed_axes", "hca_bridge_axis", "repeat_determinism", "downstream_firewall"),
  passed = c(setequal(unique(resource$unit_id), jobs$unit_id), nrow(resource) == 5L,
    all(resource$elapsed_seconds <= resource$elapsed_cap_seconds & resource$peak_rss_bytes <= resource$rss_cap_bytes),
    all(resource$stderr_class %in% c("empty", "known_glmGamPoi_native_fallback")), {
      required <- views$dataset_scope == "internal124" | views$panel_id == "common475" |
        (views$panel_id == "exact500" & views$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
      all(views$values_finite[required] & views$zero_variance_gene_count[required] == 0L & views$correlation_chord_valid[required])
    },
    nrow(primary[primary$dataset_scope == "internal124", ]) == 20L,
    nrow(primary[primary$unit_id == "HCA_BM_002", ]) == 4L, all(determinism$passed),
    !any(resource$persistence_computed | resource$landscapes_computed | resource$clustering_computed | resource$fusion_computed | resource$biological_outcomes_computed)),
  evidence = c("exactly MV8-N's three authorized units", "three primary plus two required repeats",
    "serial 1800-second/12-GiB policy", "only documented glmGamPoi-native fallback diagnostics or empty stderr", "internal exact500 plus external common475 data/residual and exact500 residual geometry valid; external SCT-data exact500 remains diagnostic-only",
    "two internal samples by five seeds by two representations", "HCA_BM_002 one seed by two representations and two panels",
    "maximum internal and HCA bridge hashes match repeats", "no PH, landscapes, clustering, fusion, or outcomes"), stringsAsFactors = FALSE)
atomic_csv(validation, file.path(public_dir, "mv08o-source-sentinel-validation.csv"))
if (!all(validation$passed)) stop("MV8-O validation failed", call. = FALSE)
cat("MV8-O source sentinel checks=", sum(validation$passed), "/", nrow(validation), "; all source/view gates pass; PH remains closed\n", sep = "")
