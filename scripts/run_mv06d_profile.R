#!/usr/bin/env Rscript

for (package in c("processx", "ps", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " is required for monitored MV6-D execution.",
         call. = FALSE)
  }
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
args <- getOption("mv06d.args", commandArgs(trailingOnly = TRUE))
if (length(args) != 7L) {
  stop(
    "usage: run_mv06d_profile.R SENTINEL_CSV PANEL_CSV CANDIDATE_CSV ",
    "RESOURCE_CSV CACHE_DIR PRIVATE_ROOT PUBLIC_DIR", call. = FALSE
  )
}
sentinel_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
panel_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
candidate_path <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
resource_path <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)
cache_dir <- normalizePath(args[[5L]], winslash = "/", mustWork = TRUE)
private_root <- args[[6L]]
public_dir <- args[[7L]]
dir.create(private_root, recursive = TRUE, showWarnings = FALSE)
dir.create(public_dir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("source", "ph", "landscape", "repeat", "logs")) {
  dir.create(file.path(private_root, subdir), recursive = TRUE,
             showWarnings = FALSE)
}

source("R/toy_baseline.R")
source("R/dual_view_topology.R")
source("R/landscape_contract.R")
source("R/landscape_reference.R")
source("R/landscape_public_api.R")
source("R/mv06d_matched_profile.R")

sentinels <- utils::read.csv(sentinel_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
if (nrow(sentinels) != 10L || !setequal(sentinels$stage, c(1L, 2L)) ||
    any(sentinels$outcome_label_state != "closed") ||
    any(as.logical(sentinels$biological_outcomes_computed))) {
  stop("MV6-D sentinel manifest violates the frozen staged contract.",
       call. = FALSE)
}

process_tree_rss <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(e) NULL)
  if (is.null(root)) return(0)
  handles <- c(list(root), tryCatch(ps::ps_children(root, recursive = TRUE),
                                   error = function(e) list()))
  sum(vapply(handles, function(handle) tryCatch(
    as.numeric(ps::ps_memory_info(handle)[["rss"]]),
    error = function(e) 0
  ), numeric(1L)))
}

run_unit <- function(unit_id, script, script_args, output_path,
                     elapsed_cap, rss_cap) {
  log_stem <- gsub("[^A-Za-z0-9_.-]", "_", unit_id)
  stdout_path <- file.path(private_root, "logs", paste0(log_stem, "__stdout.txt"))
  stderr_path <- file.path(private_root, "logs", paste0(log_stem, "__stderr.txt"))
  started <- Sys.time()
  process <- processx::process$new(
    command = Sys.which("Rscript"),
    args = c("--vanilla", script, script_args), stdout = stdout_path,
    stderr = stderr_path, cleanup_tree = TRUE
  )
  peak <- 0
  cap <- ""
  while (process$is_alive()) {
    Sys.sleep(0.25)
    peak <- max(peak, process_tree_rss(process$get_pid()))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (elapsed > elapsed_cap) {
      cap <- "elapsed_cap_exceeded"
      process$kill_tree()
    } else if (peak > rss_cap) {
      cap <- "rss_cap_exceeded"
      process$kill_tree()
    }
  }
  process$wait(timeout = 5000)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  status <- process$get_exit_status()
  disposition <- if (nzchar(cap)) cap else if (
    identical(status, 0L) && file.exists(output_path)
  ) "completed" else "failed"
  data.frame(
    unit_id = unit_id, disposition = disposition,
    exit_status = if (is.null(status)) NA_integer_ else status,
    elapsed_seconds = elapsed, peak_process_tree_rss_bytes = peak,
    output_bytes = if (file.exists(output_path)) file.info(output_path)$size else NA_real_,
    output_sha256 = if (file.exists(output_path))
      mv06d_file_sha256_v1(output_path) else NA_character_,
    elapsed_cap_seconds = elapsed_cap, rss_cap_bytes = rss_cap,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
}

file_stem <- function(fold_id, seed) paste0(
  gsub("[^A-Za-z0-9_.-]", "_", fold_id), "__", seed
)
source_rows <- list()
ph_rows <- list()
landscape_rows <- list()
repeat_rows <- list()

bind_rows_fill <- function(rows) {
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (name in missing) row[[name]] <- NA
    row[columns]
  })
  do.call(rbind, rows)
}

private_bytes <- function() {
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (!length(files)) 0 else sum(file.info(files)$size, na.rm = TRUE)
}

measured_seconds <- function() {
  rows <- c(source_rows, ph_rows, landscape_rows, repeat_rows)
  if (!length(rows)) 0 else sum(vapply(rows, function(row) {
    as.numeric(row$elapsed_seconds[[1L]])
  }, numeric(1L)), na.rm = TRUE)
}

execute_fold <- function(fold_rows, repeat_mode = FALSE) {
  fold_id <- fold_rows$fold_id[[1L]]
  seed <- fold_rows$seed[[1L]]
  stage <- fold_rows$stage[[1L]]
  stem <- file_stem(fold_id, seed)
  base_dir <- if (repeat_mode) file.path(private_root, "repeat") else private_root
  source_path <- file.path(base_dir, if (repeat_mode) "" else "source",
                           paste0(stem, "__source.rds"))
  unit <- run_unit(
    paste0(stem, if (repeat_mode) "__source_repeat" else "__source"),
    "scripts/run_mv06d_source_entry.R",
    c(sentinel_path, panel_path, candidate_path, resource_path, cache_dir,
      fold_id, as.character(seed), source_path), source_path,
    elapsed_cap = 1800, rss_cap = 8 * 1024^3
  )
  unit$contract_id <- "mv06d_source_resource_metric_v1"
  unit$stage <- stage; unit$fold_id <- fold_id; unit$seed <- seed
  if (unit$disposition != "completed") return(list(source = unit))
  source_record <- readRDS(source_path)
  mv06d_validate_source_record_v1(source_record)
  unit$shared_and_two_pair_bytes <- unit$output_bytes
  unit$cell_view_object_bytes <- sum(vapply(
    source_record$payload$views, function(pair) as.numeric(object.size(
      pair$cell_topology_v1)), numeric(1L)
  ))
  unit$gene_view_object_bytes <- sum(vapply(
    source_record$payload$views, function(pair) as.numeric(object.size(
      pair$gene_topology_v1)), numeric(1L)
  ))
  ph <- list()
  for (role in c("held_out", "training")) {
    for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
      ph_stem <- paste(stem, role, view_id, sep = "__")
      ph_path <- file.path(base_dir, if (repeat_mode) "" else "ph",
                           paste0(ph_stem, "__ph.rds"))
      cap_time <- if (view_id == "cell_topology_v1") 600 else 1800
      cap_rss <- if (view_id == "cell_topology_v1") 4 * 1024^3 else 8 * 1024^3
      metric <- run_unit(
        paste0(ph_stem, if (repeat_mode) "__repeat" else ""),
        "scripts/run_mv06d_ph_entry.R",
        c(source_path, role, view_id, ph_path), ph_path, cap_time, cap_rss
      )
      metric$contract_id <- "mv06d_ph_resource_metric_v1"
      metric$stage <- stage; metric$fold_id <- fold_id; metric$seed <- seed
      metric$role <- role; metric$view_id <- view_id
      if (metric$disposition == "completed") {
        record <- readRDS(ph_path)
        mv06d_validate_ph_record_v1(record)
        metric$sample_id <- record$identity$sample_id
        metric$point_count <- record$h0_mst_oracle$point_count
        metric$h0_intervals <- nrow(record$topology_result$diagram[
          record$topology_result$diagram[, "dimension"] == 0, , drop = FALSE])
        metric$h1_intervals <- record$h0_mst_oracle$finite_h1_intervals
        metric$h0_mst_error <- record$h0_mst_oracle$maximum_absolute_error
        metric$h0_mst_passed <- record$h0_mst_oracle$passed
        metric$ph_cache_key <- record$cache_key
      }
      ph[[paste(role, view_id, sep = "\r")]] <- metric
      if (metric$disposition != "completed") {
        return(list(source = unit, ph = ph))
      }
    }
  }
  landscapes <- list()
  for (view_id in c("cell_topology_v1", "gene_topology_v1")) {
    held_key <- paste("held_out", view_id, sep = "\r")
    train_key <- paste("training", view_id, sep = "\r")
    held_path <- file.path(base_dir, if (repeat_mode) "" else "ph",
      paste0(paste(stem, "held_out", view_id, sep = "__"), "__ph.rds"))
    train_path <- file.path(base_dir, if (repeat_mode) "" else "ph",
      paste0(paste(stem, "training", view_id, sep = "__"), "__ph.rds"))
    landscape_stem <- paste(stem, view_id, sep = "__")
    landscape_path <- file.path(base_dir,
      if (repeat_mode) "" else "landscape",
      paste0(landscape_stem, "__landscape.rds"))
    metric <- run_unit(
      paste0(landscape_stem, if (repeat_mode) "__repeat" else ""),
      "scripts/run_mv06d_landscape_entry.R",
      c(held_path, train_path, landscape_path), landscape_path,
      elapsed_cap = 1800, rss_cap = 8 * 1024^3
    )
    metric$contract_id <- "mv06d_landscape_resource_metric_v1"
    metric$stage <- stage; metric$fold_id <- fold_id; metric$seed <- seed
    metric$view_id <- view_id
    if (metric$disposition == "completed") {
      record <- readRDS(landscape_path)
      metric$h0_distance <- record$result$distances[["H0"]]
      metric$h1_distance <- record$result$distances[["H1"]]
      metric$combined_distance <- record$result$distances[["combined"]]
      metric$h1_squared_fraction <- record$result$h1_squared_distance_fraction
      metric$h0_exact <- record$result$dimensions$H0$exact
      metric$h1_exact <- record$result$dimensions$H1$exact
      metric$h0_first_intervals <-
        record$result$dimensions$H0$first_finite_intervals
      metric$h0_second_intervals <-
        record$result$dimensions$H0$second_finite_intervals
      metric$h1_first_intervals <-
        record$result$dimensions$H1$first_finite_intervals
      metric$h1_second_intervals <-
        record$result$dimensions$H1$second_finite_intervals
      metric$landscape_cache_key <- record$cache_key
    }
    landscapes[[view_id]] <- metric
    if (metric$disposition != "completed") break
  }
  list(source = unit, ph = ph, landscape = landscapes,
       source_path = source_path)
}

append_result <- function(result, repeat_mode = FALSE) {
  target <- if (repeat_mode) "repeat_rows" else "source_rows"
  if (repeat_mode) {
    repeat_rows[[length(repeat_rows) + 1L]] <<- result$source
  } else {
    source_rows[[length(source_rows) + 1L]] <<- result$source
  }
  if (length(result$ph)) {
    for (row in result$ph) {
      if (repeat_mode) repeat_rows[[length(repeat_rows) + 1L]] <<- row else
        ph_rows[[length(ph_rows) + 1L]] <<- row
    }
  }
  if (length(result$landscape)) {
    for (row in result$landscape) {
      if (repeat_mode) repeat_rows[[length(repeat_rows) + 1L]] <<- row else
        landscape_rows[[length(landscape_rows) + 1L]] <<- row
    }
  }
}

stage1 <- sentinels[sentinels$stage == 1L, , drop = FALSE]
primary1 <- execute_fold(stage1)
append_result(primary1)
stage1_complete <- primary1$source$disposition == "completed" &&
  length(primary1$ph) == 4L &&
  all(vapply(primary1$ph, function(x) x$disposition == "completed", logical(1L))) &&
  length(primary1$landscape) == 2L &&
  all(vapply(primary1$landscape, function(x) x$disposition == "completed",
             logical(1L)))
repeat_complete <- FALSE
repeat_exact <- FALSE
if (stage1_complete) {
  repeated1 <- execute_fold(stage1, repeat_mode = TRUE)
  append_result(repeated1, repeat_mode = TRUE)
  repeat_complete <- repeated1$source$disposition == "completed" &&
    length(repeated1$ph) == 4L && length(repeated1$landscape) == 2L &&
    all(vapply(c(repeated1$ph, repeated1$landscape), function(x)
      x$disposition == "completed", logical(1L)))
  if (repeat_complete) {
    primary_hashes <- c(primary1$source$output_sha256,
      vapply(primary1$ph, `[[`, character(1L), "output_sha256"),
      vapply(primary1$landscape, `[[`, character(1L), "output_sha256"))
    repeat_hashes <- c(repeated1$source$output_sha256,
      vapply(repeated1$ph, `[[`, character(1L), "output_sha256"),
      vapply(repeated1$landscape, `[[`, character(1L), "output_sha256"))
    primary_bytes <- c(primary1$source$output_bytes,
      vapply(primary1$ph, `[[`, numeric(1L), "output_bytes"),
      vapply(primary1$landscape, `[[`, numeric(1L), "output_bytes"))
    repeat_bytes <- c(repeated1$source$output_bytes,
      vapply(repeated1$ph, `[[`, numeric(1L), "output_bytes"),
      vapply(repeated1$landscape, `[[`, numeric(1L), "output_bytes"))
    repeat_exact <- identical(primary_hashes, repeat_hashes) &&
      identical(primary_bytes, repeat_bytes)
  }
}
gene_stage1 <- if (length(primary1$ph)) do.call(rbind, primary1$ph) else NULL
continue_stage2 <- stage1_complete && repeat_complete && repeat_exact &&
  all(gene_stage1$elapsed_seconds[gene_stage1$view_id ==
                                  "gene_topology_v1"] <= 900) &&
  all(gene_stage1$peak_process_tree_rss_bytes[gene_stage1$view_id ==
                                  "gene_topology_v1"] <= 6 * 1024^3) &&
  measured_seconds() <= 14400 && private_bytes() <= 10 * 1024^3

if (continue_stage2) {
  stage2 <- sentinels[sentinels$stage == 2L, , drop = FALSE]
  groups <- split(stage2, paste(stage2$fold_id, stage2$seed, sep = "\r"))
  groups <- groups[order(vapply(groups, function(x) x$seed[[1L]], integer(1L)))]
  for (group in groups) {
    result <- execute_fold(group)
    append_result(result)
    complete <- result$source$disposition == "completed" &&
      length(result$ph) == 4L && length(result$landscape) == 2L &&
      all(vapply(c(result$ph, result$landscape), function(x)
        x$disposition == "completed", logical(1L)))
    if (!complete || measured_seconds() > 14400 ||
        private_bytes() > 10 * 1024^3) break
  }
}

source_metrics <- bind_rows_fill(source_rows)
ph_metrics <- bind_rows_fill(ph_rows)
landscape_metrics <- bind_rows_fill(landscape_rows)
repeat_metrics <- bind_rows_fill(repeat_rows)
utils::write.csv(source_metrics, file.path(public_dir, "mv06d-source-metrics.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(ph_metrics, file.path(public_dir, "mv06d-ph-metrics.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(landscape_metrics,
                 file.path(public_dir, "mv06d-landscape-metrics.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(repeat_metrics, file.path(public_dir, "mv06d-repeat-metrics.csv"),
                 row.names = FALSE, na = "")

full_complete <- continue_stage2 && nrow(source_metrics) == 5L &&
  nrow(ph_metrics) == 20L && nrow(landscape_metrics) == 10L &&
  all(source_metrics$disposition == "completed") &&
  all(ph_metrics$disposition == "completed") &&
  all(landscape_metrics$disposition == "completed")
projection <- if (full_complete) mv06d_project_workload_v1(
  source_metrics, ph_metrics, landscape_metrics
) else data.frame()
if (nrow(projection)) utils::write.csv(
  projection, file.path(public_dir, "mv06d-worker-projection.csv"),
  row.names = FALSE, na = ""
)

storage_projection <- data.frame()
projected_storage <- NA_real_
if (full_complete) {
  mean_cell_view <- mean(source_metrics$cell_view_object_bytes / 2)
  mean_gene_view <- mean(source_metrics$gene_view_object_bytes / 2)
  mean_shared <- mean(pmax(
    0, source_metrics$shared_and_two_pair_bytes -
      source_metrics$cell_view_object_bytes -
      source_metrics$gene_view_object_bytes
  ))
  components <- c(
    fold_shared = 75 * mean_shared,
    cell_views = 6750 * mean_cell_view,
    gene_views = 6750 * mean_gene_view,
    cell_ph = 6750 * mean(ph_metrics$output_bytes[
      ph_metrics$view_id == "cell_topology_v1"]),
    gene_ph = 6750 * mean(ph_metrics$output_bytes[
      ph_metrics$view_id == "gene_topology_v1"]),
    cell_landscape_pairs = 35350 * mean(landscape_metrics$output_bytes[
      landscape_metrics$view_id == "cell_topology_v1"]),
    gene_landscape_pairs = 35350 * mean(landscape_metrics$output_bytes[
      landscape_metrics$view_id == "gene_topology_v1"])
  )
  storage_projection <- data.frame(
    contract_id = "mv06d_private_storage_projection_v1",
    component = names(components), projected_bytes = unname(components),
    projection_basis = "bounded_mean_serialized_or_object_size_v1",
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  projected_storage <- sum(components)
  utils::write.csv(storage_projection,
    file.path(public_dir, "mv06d-storage-projection.csv"),
    row.names = FALSE, na = "")
}

maximum_hours <- if (nrow(projection)) sum(projection$projected_worker_hours[
  projection$scenario == "observed_maximum"]) else NA_real_
correctness_pass <- stage1_complete && repeat_complete && repeat_exact &&
  all(ph_metrics$h0_mst_passed)
decision_id <- if (!stage1_complete || !repeat_complete || !repeat_exact) {
  "stop_matched_sct_scaleup"
} else if (!full_complete) {
  "insufficient_profile"
} else if (correctness_pass && maximum_hours <= 72 &&
           projected_storage <= 10 * 1024^3) {
  "go_prefreeze_full_matched_sct"
} else {
  "revise_for_targeted_acceleration"
}
decision <- data.frame(
  contract_id = "mv06d_matched_sct_profile_decision_v1",
  decision = decision_id, stage1_complete = stage1_complete,
  stage1_repeat_complete = repeat_complete,
  stage1_seven_files_byte_identical = repeat_exact,
  stage2_authorized = continue_stage2, full_profile_complete = full_complete,
  bounded_worker_seconds = measured_seconds(),
  bounded_private_storage_bytes = private_bytes(),
  source_units = nrow(source_metrics), ph_units = nrow(ph_metrics),
  landscape_pair_units = nrow(landscape_metrics),
  projected_maximum_worker_hours = maximum_hours,
  projected_private_storage_bytes = projected_storage,
  full_production_authorized = FALSE, fusion_jobs_executed = 0L,
  clustering_jobs_executed = 0L, outcome_jobs_executed = 0L,
  outcome_label_state = "closed", biological_outcomes_computed = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(decision, file.path(public_dir, "mv06d-decision.csv"),
                 row.names = FALSE, na = "")
message("MV6-D disposition: ", decision_id)
