#!/usr/bin/env Rscript

# Adopt the exact completed-but-unreceipted MV8-ZT engine-v2 child at order
# 326. This executor performs no landscape computation and no retry. It is
# re-entrant across ledger, completion, and progress receipt promotion.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(paste(
  "usage: recover_mv08zt_landscape_interruption.R",
  "<recovery-prefreeze> <mv08zt-prefreeze> <private-bindings>",
  "<private-root> <public-root> <recovery-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

recovery_root <- normalizePath(args[[1L]], mustWork = TRUE)
zt_root <- normalizePath(args[[2L]], mustWork = TRUE)
bindings_path <- normalizePath(args[[3L]], mustWork = TRUE)
private_root <- normalizePath(args[[4L]], mustWork = TRUE)
public_root <- normalizePath(args[[5L]], mustWork = TRUE)
recovery_head <- tolower(trimws(args[[6L]]))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) stop("invalid recovery head", call. = FALSE)

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
    stop("MV8-ZT recovery manifest drift: ", name, call. = FALSE)
  }
}
verify_manifest(recovery_root, "mv08zt-r1-artifact-manifest.csv")
verify_manifest(zt_root, "mv08zt-artifact-manifest.csv")

snapshot <- read_csv(file.path(recovery_root, "mv08zt-r1-stopped-snapshot.csv"))
orphan <- read_csv(file.path(recovery_root, "mv08zt-r1-orphan-binding.csv"))
decision <- read_csv(file.path(recovery_root, "mv08zt-r1-decision.csv"))
implementation <- read_csv(file.path(recovery_root, "mv08zt-r1-implementation-bindings.csv"))
frozen_ledger <- read_csv(file.path(recovery_root, "mv08zt-r1-prefix-ledger.csv"))
frozen_completed <- read_csv(file.path(recovery_root, "mv08zt-r1-prefix-completions.csv"))
frozen_progress <- read_csv(file.path(recovery_root, "mv08zt-r1-prefix-progress.csv"))
contract <- read_csv(file.path(zt_root, "mv08zt-contract.csv"))
queue <- read_csv(file.path(zt_root, "mv08zt-production-queue.csv"))
inputs <- read_csv(file.path(zt_root, "mv08zt-input-bindings.csv"))
bindings <- read_csv(bindings_path)
ledger_path <- file.path(public_root, "mv08zt-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zt-chunk-completions.csv")
progress_path <- file.path(public_root, "mv08zt-progress.csv")
ledger <- read_csv(ledger_path)
completed <- read_csv(completion_path)
progress <- read_csv(progress_path)

environment_head <- tolower(trimws(Sys.getenv(
  "MV08ZT_RECOVERY_GIT_HEAD", unset = ""
)))
executor_file <- "scripts/recover_mv08zt_landscape_interruption.R"
bound_executor <- implementation[implementation$role == "recovery_executor", , drop = FALSE]
binding_hash <- inputs$sha256[inputs$role == "private_unit_bindings"]
if (nrow(decision) != 1L || !truth(decision$orphan_adoption_authorized) ||
    decision$orphan_production_order != 326L ||
    decision$resume_at_production_order != 327L ||
    decision$automatic_retries != 0L ||
    truth(decision$landscape_recomputation_authorized) ||
    truth(decision$scientific_contract_changed) ||
    !grepl("^[0-9a-f]{40}$", environment_head) ||
    environment_head != recovery_head ||
    nrow(bound_executor) != 1L || bound_executor$file != executor_file ||
    sha_file(executor_file) != bound_executor$sha256 ||
    length(binding_hash) != 1L || sha_file(bindings_path) != binding_hash) {
  stop("MV8-ZT committed recovery authorization drift", call. = FALSE)
}
runner_lines <- suppressWarnings(system2(
  "pgrep", c("-f", "[r]un_mv08zt_full_landscape_production[.]R"),
  stdout = TRUE, stderr = FALSE
))
if (length(runner_lines)) stop("MV8-ZT runner is active", call. = FALSE)

prefix_n <- 325L
index <- 326L
if (nrow(queue) != 628L || nrow(frozen_ledger) != prefix_n ||
    nrow(frozen_completed) != prefix_n || nrow(frozen_progress) != 1L ||
    snapshot$completed_prefix != prefix_n ||
    snapshot$unreceipted_complete_order != index ||
    snapshot$execution_head != frozen_progress$execution_head ||
    !identical(as.integer(frozen_ledger$production_order), seq_len(prefix_n)) ||
    !identical(as.integer(frozen_completed$production_order), seq_len(prefix_n))) {
  stop("MV8-ZT frozen recovery prefix is invalid", call. = FALSE)
}
if (!(nrow(ledger) %in% c(prefix_n, index)) ||
    !(nrow(completed) %in% c(prefix_n, index)) ||
    !identical(ledger[seq_len(prefix_n), , drop = FALSE], frozen_ledger) ||
    !identical(completed[seq_len(prefix_n), , drop = FALSE], frozen_completed) ||
    nrow(progress) != 1L || !(progress$completed_chunks %in% c(prefix_n, index))) {
  stop("MV8-ZT public recovery state is not an accepted boundary", call. = FALSE)
}
state_key <- paste(nrow(ledger), nrow(completed), progress$completed_chunks, sep = "/")
if (!(state_key %in% c("325/325/325", "326/325/325", "326/326/325", "326/326/326"))) {
  stop("MV8-ZT public recovery receipts are non-monotone", call. = FALSE)
}

chunk_paths <- function(value) {
  root <- file.path(private_root, "production",
                    .mv08z_safe_group(value$group_order),
                    .mv08z_safe_chunk(value$chunk_order))
  c(distance = file.path(root, "distances.csv"), status = file.path(root, "status.csv"))
}
for (prefix_index in seq_len(prefix_n)) {
  value <- chunk_paths(queue[prefix_index, , drop = FALSE])
  if (!all(file.exists(value)) ||
      sha_file(value[["distance"]]) != frozen_completed$distances_sha256[[prefix_index]] ||
      sha_file(value[["status"]]) != frozen_completed$status_sha256[[prefix_index]]) {
    stop("MV8-ZT accepted-prefix private evidence drift at ", prefix_index, call. = FALSE)
  }
}

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
  stop("MV8-ZT orphan artifact or log drift", call. = FALSE)
}
status <- read_csv(paths[["status"]])
distances <- read_csv(paths[["distance"]])
group_bindings <- bindings[as.integer(bindings$group_order) == row$group_order, , drop = FALSE]
expected_pairs <- .mv08z_add_pair_identities(.mv08z_group_pairs(group_bindings), row$group_id)
expected_pairs <- expected_pairs[
  expected_pairs$pair_ordinal >= row$pair_start &
    expected_pairs$pair_ordinal <= row$pair_end, , drop = FALSE
]
stderr_text <- trimws(paste(readLines(stderr, warn = FALSE), collapse = "\n"))
expected_stderr <- paste0("Completed MV8-Z ", .mv08z_safe_group(row$group_order), "/",
                          .mv08z_safe_chunk(row$chunk_order), "; pairs=", row$pair_count)
if (nrow(status) != 1L || nrow(distances) != row$pair_count ||
    status$completion_state != "complete" || status$mode != "production" ||
    status$execution_head != snapshot$execution_head ||
    status$scientific_engine_version != 2L ||
    status$group_order != row$group_order || status$chunk_order != row$chunk_order ||
    status$pair_subset_sha256 != row$pair_subset_sha256 ||
    status$distances_sha256 != sha_file(paths[["distance"]]) ||
    !truth(status$workers == 1L) || status$retries != 0L || truth(status$fallback_used) ||
    .mv08z_sha256_text(expected_pairs$pair_identity_sha256) != row$pair_subset_sha256 ||
    !identical(distances$pair_identity_sha256, expected_pairs$pair_identity_sha256) ||
    stderr_text != expected_stderr) {
  stop("MV8-ZT orphan scientific validation failed", call. = FALSE)
}
logs <- list.files(file.path(private_root, "logs"), full.names = TRUE)
expected_logs <- unlist(lapply(seq_len(index), function(value) {
  file.path(private_root, "logs", sprintf("chunk_%04d.%s", value, c("stdout", "stderr")))
}))
partials <- list.files(private_root, pattern = "[.]partial$", recursive = TRUE,
                       full.names = TRUE, all.files = TRUE)
if (!setequal(normalizePath(logs, mustWork = TRUE),
              normalizePath(expected_logs, mustWork = TRUE)) || length(partials)) {
  stop("MV8-ZT ambiguous logs or partial artifacts remain", call. = FALSE)
}

# Parent telemetry was lost after successful child completion. Record the
# frozen caps as conservative upper bounds, never as recovered measurements.
elapsed_upper <- as.numeric(contract$child_elapsed_cap_seconds)
rss_upper <- as.numeric(contract$child_rss_cap_bytes)
ledger_row <- frozen_ledger[1L, , drop = FALSE]
ledger_row[1L, ] <- NA
ledger_row$contract_id <- "mv08zt_resource_ledger_v1"
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
ledger_row$scientific_engine_version <- 2L
ledger_row$workers <- 1L
ledger_row$retries <- 0L
ledger_row$outcome_label_state <- "closed"
ledger_row$biological_outcomes_computed <- FALSE

completion_row <- frozen_completed[1L, , drop = FALSE]
completion_row[1L, ] <- NA
completion_row$contract_id <- "mv08zt_chunk_completion_v1"
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
completion_row$scientific_engine_version <- 2L
completion_row$workers <- 1L
completion_row$retries <- 0L
completion_row$outcome_label_state <- "closed"
completion_row$biological_outcomes_computed <- FALSE

if (nrow(ledger) == prefix_n) {
  ledger <- rbind(frozen_ledger, ledger_row)
  atomic_csv(ledger, ledger_path)
} else if (!identical(ledger[index, , drop = FALSE], ledger_row)) {
  stop("MV8-ZT existing recovered ledger row drift", call. = FALSE)
}
if (nrow(completed) == prefix_n) {
  completed <- rbind(frozen_completed, completion_row)
  atomic_csv(completed, completion_path)
} else if (!identical(completed[index, , drop = FALSE], completion_row)) {
  stop("MV8-ZT existing recovered completion row drift", call. = FALSE)
}
if (progress$completed_chunks == prefix_n) {
  progress$completed_chunks <- index
  progress$completed_pairs <- sum(completed$pair_count)
  progress$aggregate_child_seconds <- sum(ledger$elapsed_seconds)
  files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[!file.info(files)$isdir]
  progress$private_bytes <- sum(as.numeric(file.info(files)$size))
  progress$state <- "running"
  atomic_csv(progress, progress_path)
}
cat("MV8-ZT adopted completed unreceipted engine-v2 chunk 326 without retry; resume_at=327\n")
