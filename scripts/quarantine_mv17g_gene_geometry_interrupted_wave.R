#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "usage: quarantine_mv17g_gene_geometry_interrupted_wave.R",
      "<original-private-prefreeze> <parallel-public-prefreeze>",
      "<production-private-root> <quarantine-root>"
    ),
    call. = FALSE
  )
}

original_private <- normalizePath(args[[1L]], mustWork = TRUE)
parallel_public <- normalizePath(args[[2L]], mustWork = TRUE)
production <- normalizePath(args[[3L]], mustWork = TRUE)
quarantine <- args[[4L]]
if (dir.exists(quarantine)) {
  stop("MV17-G gene-geometry interruption quarantine already exists", call. = FALSE)
}

source("R/mv08z_landscape_production.R")
source("R/mv17_null_models.R")
source("R/mv17_calibration.R")
source("R/mv17_full_calibration.R")
source("R/mv17g_parallel_recovery.R")
read_csv <- .mv08z_read_csv
write_csv <- .mv08z_atomic_csv
sha256 <- .mv08z_sha256_file

plan_text <- paste(readLines("PROJECT_PLAN.md", warn = FALSE), collapse = "\n")
if (!grepl(
  "mv17g_gene_geometry_interrupted_wave_quarantine_authorized_v1",
  plan_text, fixed = TRUE
)) {
  stop("MV17-G interrupted-wave quarantine lacks authorization", call. = FALSE)
}
processes <- system2("ps", c("-eo", "args="), stdout = TRUE)
if (any(grepl("run_mv17g_calibration_parallel_recovery[.]R", processes)) ||
    any(grepl("run_mv17g_calibration_group_worker[.]R", processes))) {
  stop("MV17-G controller/workers must be absent before quarantine", call. = FALSE)
}

queue <- read_csv(file.path(original_private, "mv17g-primary-grouped-queue.csv"))
contract <- read_csv(file.path(parallel_public, "mv17g-parallel-contract.csv"))
progress <- read_csv(file.path(production, "mv17g-private-progress.csv"))
scan <- mv17g_checkpoint_scan_v1(queue, production)
partials <- list.files(
  production, pattern = "[.]partial$", recursive = TRUE, full.names = TRUE
)
admitted_prefix <- nrow(progress)
artifact_prefix <- mv17g_complete_prefix_v1(scan, require_incomplete = TRUE)
wave_orders <- seq.int(admitted_prefix + 1L, artifact_prefix)
if (length(partials) != 0L || admitted_prefix != 1018L ||
    artifact_prefix != 1026L || !identical(wave_orders, 1019:1026) ||
    !identical(as.integer(progress$job_order), seq_len(admitted_prefix)) ||
    any(queue$view[wave_orders] != "gene")) {
  stop("MV17-G interrupted-wave boundary drift", call. = FALSE)
}

inspect_one <- function(i) {
  q <- queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, production)
  if (!all(file.exists(paths)) ||
      file.info(paths[["stdout"]])$size != 0L ||
      file.info(paths[["stderr"]])$size != 0L) {
    stop("MV17-G interrupted-wave artifact drift at order ", i, call. = FALSE)
  }
  result <- readRDS(paths[["result"]])
  matrix_path <- file.path(
    production, "matrices", sprintf("gene__%03d.rds", q$unit_order)
  )
  expected_seeds <- if (q$null_family == "observed") {
    0L
  } else {
    seq.int(q$seed_first, length.out = q$replicate_count)
  }
  if (!identical(result$contract_id, "mv17g_group_result_v1") ||
      !identical(result$null_family, q$null_family) ||
      result$matrix_sha256 != sha256(matrix_path) ||
      result$replicate_count != length(expected_seeds) ||
      !setequal(unique(result$metrics$seed), expected_seeds) ||
      nrow(result$metrics) != length(expected_seeds) * 8L ||
      any(result$metrics$h0_mst_maximum_absolute_error > 1e-8) ||
      !isTRUE(result$finite) || result$labels_opened || result$outcomes_opened) {
    stop("MV17-G interrupted-wave payload drift at order ", i, call. = FALSE)
  }
  resource <- mv17c_parse_gnu_time_v1(paths[["time"]])
  if (resource$exit_status != 0L ||
      resource$wall_seconds > contract$child_timeout_seconds ||
      resource$maximum_RSS_bytes > contract$child_RSS_cap_bytes) {
    stop("MV17-G interrupted-wave resource drift at order ", i, call. = FALSE)
  }
  do.call(rbind, lapply(names(paths), function(role) {
    data.frame(
      contract_id = "mv17g_gene_geometry_interrupted_wave_quarantine_v1",
      admitted_prefix = admitted_prefix,
      artifact_prefix = artifact_prefix,
      job_order = q$job_order,
      view = q$view,
      unit_order = q$unit_order,
      null_family = q$null_family,
      role = role,
      original_artifact = normalizePath(paths[[role]]),
      bytes = as.numeric(file.info(paths[[role]])$size),
      sha256 = sha256(paths[[role]]),
      time_receipt_valid = TRUE,
      disposition = "preserved_unadmitted_rejected_raw_euclidean_gene_evidence_v1",
      stringsAsFactors = FALSE
    )
  }))
}

inventory <- do.call(rbind, lapply(wave_orders, inspect_one))
dir.create(file.path(quarantine, "jobs"), recursive = TRUE)
dir.create(file.path(quarantine, "logs"), recursive = TRUE)
destinations <- file.path(
  quarantine,
  ifelse(inventory$role == "result", "jobs", "logs"),
  basename(inventory$original_artifact)
)
moved <- logical(nrow(inventory))
rollback <- TRUE
on.exit({
  if (rollback && any(moved)) {
    for (i in rev(which(moved))) {
      file.rename(destinations[[i]], inventory$original_artifact[[i]])
    }
  }
}, add = TRUE)
for (i in seq_len(nrow(inventory))) {
  if (!file.rename(inventory$original_artifact[[i]], destinations[[i]])) {
    stop("MV17-G interrupted-wave quarantine move failed", call. = FALSE)
  }
  moved[[i]] <- TRUE
}
inventory$quarantine_artifact <- normalizePath(destinations)
if (any(file.exists(inventory$original_artifact)) ||
    !all(file.exists(inventory$quarantine_artifact)) ||
    !identical(
      unname(as.numeric(file.info(inventory$quarantine_artifact)$size)),
      unname(as.numeric(inventory$bytes))
    ) ||
    !identical(
      unname(vapply(inventory$quarantine_artifact, sha256, character(1L))),
      unname(tolower(inventory$sha256))
    )) {
  stop("MV17-G interrupted-wave quarantine verification failed", call. = FALSE)
}
post_scan <- mv17g_checkpoint_scan_v1(queue, production)
if (mv17g_complete_prefix_v1(post_scan, require_incomplete = TRUE) != admitted_prefix) {
  stop("MV17-G production root did not return to admitted prefix", call. = FALSE)
}
write_csv(
  inventory,
  file.path(quarantine, "mv17g-gene-geometry-interrupted-wave-quarantine.csv")
)
rollback <- FALSE
message(
  "Quarantined MV17-G unadmitted orders 1019--1026; artifacts=",
  nrow(inventory), "; admitted_prefix=", admitted_prefix
)
