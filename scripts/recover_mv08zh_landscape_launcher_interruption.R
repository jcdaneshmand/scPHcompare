#!/usr/bin/env Rscript

# Adopt one fully written but unreceipted MV8-ZF child after its parent launcher
# was interrupted. This script performs no landscape computation and no retry.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: recover_mv08zh_landscape_launcher_interruption.R",
  "<mv08zh-prefreeze> <mv08zf-prefreeze> <private-root> <public-root>",
  "<recovery-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zh_root <- normalizePath(args[[1L]], mustWork = TRUE)
zf_root <- normalizePath(args[[2L]], mustWork = TRUE)
private_root <- normalizePath(args[[3L]], mustWork = TRUE)
public_root <- normalizePath(args[[4L]], mustWork = TRUE)
recovery_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) stop("invalid MV8-ZH recovery head", call. = FALSE)

source("R/mv08z_landscape_production.R")
read_csv <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
truth <- .mv08z_truth
atomic_csv <- .mv08z_atomic_csv

verify_manifest <- function(root, name) {
  manifest <- read_csv(file.path(root, name))
  paths <- file.path(root, manifest$artifact)
  if (!all(file.exists(paths)) ||
      !all(as.numeric(file.info(paths)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(paths, sha_file, character(1L)) == manifest$sha256)) {
    stop("MV8-ZH manifest drift: ", name, call. = FALSE)
  }
  manifest
}
verify_manifest(zh_root, "mv08zh-artifact-manifest.csv")
verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")

decision <- read_csv(file.path(zh_root, "mv08zh-decision.csv"))
snapshot <- read_csv(file.path(zh_root, "mv08zh-stopped-snapshot.csv"))
orphan <- read_csv(file.path(zh_root, "mv08zh-orphan-binding.csv"))
implementation <- read_csv(file.path(zh_root, "mv08zh-implementation-bindings.csv"))
queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
contract <- read_csv(file.path(zf_root, "mv08zf-contract.csv"))
ledger_path <- file.path(public_root, "mv08zf-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zf-chunk-completions.csv")
progress_path <- file.path(public_root, "mv08zf-progress.csv")
ledger <- read_csv(ledger_path)
completed <- read_csv(completion_path)
progress <- read_csv(progress_path)

current_head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
recovery_file <- "scripts/recover_mv08zh_landscape_launcher_interruption.R"
bound_recovery <- implementation[implementation$role == "recovery_executor", , drop = FALSE]
if (nrow(decision) != 1L || !truth(decision$orphan_adoption_authorized) ||
    decision$orphan_production_order != 164L || decision$resume_at_production_order != 165L ||
    decision$automatic_retries != 0L || truth(decision$landscape_recomputation_authorized) ||
    current_head != recovery_head || nrow(bound_recovery) != 1L ||
    bound_recovery$file != recovery_file || sha_file(recovery_file) != bound_recovery$sha256) {
  stop("MV8-ZH committed recovery authorization drift", call. = FALSE)
}

prefix_n <- as.integer(snapshot$completed_prefix)
if (prefix_n != 163L || nrow(queue) != 628L || nrow(ledger) != prefix_n ||
    nrow(completed) != prefix_n || nrow(progress) != 1L ||
    progress$state != "running" || progress$completed_chunks != prefix_n ||
    sha_file(ledger_path) != snapshot$ledger_sha256 ||
    sha_file(completion_path) != snapshot$completion_sha256 ||
    sha_file(progress_path) != snapshot$progress_sha256 ||
    !identical(as.integer(ledger$production_order), seq_len(prefix_n)) ||
    !identical(as.integer(completed$production_order), seq_len(prefix_n)) ||
    any(ledger$disposition != "completed") || any(ledger$exit_status != 0L) ||
    any(ledger$retries != 0L) || any(completed$retries != 0L)) {
  stop("MV8-ZH accepted public prefix drift", call. = FALSE)
}

chunk_paths <- function(row) {
  root <- file.path(private_root, "production",
                    .mv08z_safe_group(row$group_order),
                    .mv08z_safe_chunk(row$chunk_order))
  c(distance = file.path(root, "distances.csv"), status = file.path(root, "status.csv"))
}
for (index in seq_len(prefix_n)) {
  paths <- chunk_paths(queue[index, , drop = FALSE])
  if (!all(file.exists(paths)) ||
      sha_file(paths[["distance"]]) != completed$distances_sha256[[index]] ||
      sha_file(paths[["status"]]) != completed$status_sha256[[index]]) {
    stop("MV8-ZH accepted-prefix private evidence drift at ", index, call. = FALSE)
  }
}

index <- 164L
row <- queue[index, , drop = FALSE]
paths <- chunk_paths(row)
stdout <- file.path(private_root, "logs", sprintf("chunk_%04d.stdout", index))
stderr <- file.path(private_root, "logs", sprintf("chunk_%04d.stderr", index))
bound_paths <- c(paths, stdout = stdout, stderr = stderr)
if (!all(file.exists(bound_paths)) ||
    !all(as.numeric(file.info(bound_paths)$size) ==
         as.numeric(unlist(orphan[paste0(names(bound_paths), "_bytes")], use.names = FALSE))) ||
    !all(vapply(bound_paths, sha_file, character(1L)) ==
         as.character(unlist(orphan[paste0(names(bound_paths), "_sha256")], use.names = FALSE)))) {
  stop("MV8-ZH orphan artifact or log drift", call. = FALSE)
}
status <- read_csv(paths[["status"]])
distances <- read_csv(paths[["distance"]])
stderr_text <- trimws(paste(readLines(stderr, warn = FALSE), collapse = "\n"))
expected_stderr <- paste0("Completed MV8-Z ", .mv08z_safe_group(row$group_order), "/",
                          .mv08z_safe_chunk(row$chunk_order), "; pairs=", row$pair_count)
if (nrow(status) != 1L || nrow(distances) != row$pair_count ||
    status$completion_state != "complete" || status$mode != "production" ||
    status$execution_head != snapshot$execution_head ||
    status$group_order != row$group_order || status$chunk_order != row$chunk_order ||
    status$pair_subset_sha256 != row$pair_subset_sha256 ||
    status$distances_sha256 != sha_file(paths[["distance"]]) ||
    !truth(status$workers == 1L) || status$retries != 0L || truth(status$fallback_used) ||
    !identical(distances$pair_identity_sha256,
               read_csv(file.path(
                 "docs/audits/mv08z-landscape-execution-prefreeze-v1",
                 "mv08z-pair-queue.csv"
               ))$pair_identity_sha256[row$pair_start:row$pair_end]) ||
    stderr_text != expected_stderr) {
  stop("MV8-ZH orphan scientific validation failed", call. = FALSE)
}

logs <- list.files(file.path(private_root, "logs"), full.names = TRUE)
expected_logs <- unlist(lapply(seq_len(index), function(i) {
  file.path(private_root, "logs", sprintf("chunk_%04d.%s", i, c("stdout", "stderr")))
}))
partials <- list.files(private_root, pattern = "partial", recursive = TRUE,
                       full.names = TRUE, all.files = TRUE)
if (!setequal(normalizePath(logs, mustWork = TRUE),
              normalizePath(expected_logs, mustWork = TRUE)) || length(partials)) {
  stop("MV8-ZH ambiguous logs or partial artifacts remain", call. = FALSE)
}

# The child completed under the live parent monitor, but the launcher was lost
# before publishing its outer telemetry. Use the already-frozen caps as explicit
# conservative upper bounds; do not represent them as recovered measurements.
elapsed_upper <- as.numeric(contract$child_elapsed_cap_seconds)
rss_upper <- as.numeric(contract$child_rss_cap_bytes)
ledger_row <- ledger[1L, , drop = FALSE]
ledger_row[1L, ] <- NA
ledger_row$contract_id <- "mv08zf_resource_ledger_v1"
ledger_row$execution_head <- snapshot$execution_head
ledger_row$production_order <- index
ledger_row$global_chunk_order <- row$global_chunk_order
ledger_row$group_order <- row$group_order
ledger_row$chunk_order <- row$chunk_order
ledger_row$pair_count <- row$pair_count
ledger_row$disposition <- "completed"
ledger_row$exit_status <- 0L
ledger_row$elapsed_seconds <- elapsed_upper
ledger_row$peak_process_tree_rss_bytes <- rss_upper
ledger_row$elapsed_cap_seconds <- elapsed_upper
ledger_row$rss_cap_bytes <- rss_upper
ledger_row$distances_bytes <- as.numeric(file.info(paths[["distance"]])$size)
ledger_row$distances_sha256 <- sha_file(paths[["distance"]])
ledger_row$status_bytes <- as.numeric(file.info(paths[["status"]])$size)
ledger_row$status_sha256 <- sha_file(paths[["status"]])
ledger_row$stdout_bytes <- as.numeric(file.info(stdout)$size)
ledger_row$stderr_bytes <- as.numeric(file.info(stderr)$size)
ledger_row$stderr_class <- "expected_completion"
ledger_row$workers <- 1L
ledger_row$retries <- 0L
ledger_row$outcome_label_state <- "closed"
ledger_row$biological_outcomes_computed <- FALSE

completion_row <- completed[1L, , drop = FALSE]
completion_row[1L, ] <- NA
completion_row$contract_id <- "mv08zf_chunk_completion_v1"
completion_row$execution_head <- snapshot$execution_head
completion_row$production_order <- index
completion_row$global_chunk_order <- row$global_chunk_order
completion_row$group_order <- row$group_order
completion_row$chunk_order <- row$chunk_order
completion_row$pair_count <- row$pair_count
completion_row$pair_subset_sha256 <- row$pair_subset_sha256
completion_row$distances_bytes <- ledger_row$distances_bytes
completion_row$distances_sha256 <- ledger_row$distances_sha256
completion_row$status_bytes <- ledger_row$status_bytes
completion_row$status_sha256 <- ledger_row$status_sha256
completion_row$exact <- TRUE
completion_row$all_active_levels <- TRUE
completion_row$grid_points <- 0L
completion_row$level_cap_applied <- FALSE
completion_row$workers <- 1L
completion_row$retries <- 0L
completion_row$outcome_label_state <- "closed"
completion_row$biological_outcomes_computed <- FALSE

ledger <- rbind(ledger, ledger_row)
completed <- rbind(completed, completion_row)
progress$completed_chunks <- index
progress$completed_pairs <- sum(completed$pair_count)
progress$aggregate_child_seconds <- sum(ledger$elapsed_seconds)
progress$private_bytes <- sum(as.numeric(file.info(list.files(
  private_root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
  no.. = TRUE
))$size), na.rm = TRUE)
progress$state <- "running"

atomic_csv(ledger, ledger_path)
atomic_csv(completed, completion_path)
atomic_csv(progress, progress_path)
cat("MV8-ZH adopted completed unreceipted chunk 164 without retry; resume_at=165\n")
