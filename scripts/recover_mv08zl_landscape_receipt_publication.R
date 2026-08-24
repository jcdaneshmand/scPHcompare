#!/usr/bin/env Rscript

# Promote the exact completed order-280 receipt already present in the frozen
# public partial, then refresh progress. No landscape, retry, or ledger row is
# computed. The executor is re-entrant across its two publication operations.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "usage: recover_mv08zl_landscape_receipt_publication.R",
  "<mv08zl-prefreeze> <mv08zf-prefreeze> <private-root> <public-root>",
  "<recovery-head>"
), call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required", call. = FALSE)

zl_root <- normalizePath(args[[1L]], mustWork = TRUE)
zf_root <- normalizePath(args[[2L]], mustWork = TRUE)
private_root <- normalizePath(args[[3L]], mustWork = TRUE)
public_root <- normalizePath(args[[4L]], mustWork = TRUE)
recovery_head <- tolower(trimws(args[[5L]]))
if (!grepl("^[0-9a-f]{40}$", recovery_head)) stop("invalid MV8-ZL head", call. = FALSE)

source("R/mv08z_landscape_production.R")
read_csv <- .mv08z_read_csv
sha_file <- .mv08z_sha256_file
truth <- .mv08z_truth
atomic_csv <- .mv08z_atomic_csv
verify_manifest <- function(root, name) {
  manifest <- read_csv(file.path(root, name))
  files <- file.path(root, manifest$artifact)
  if (!all(file.exists(files)) ||
      !all(as.numeric(file.info(files)$size) == as.numeric(manifest$bytes)) ||
      !all(vapply(files, sha_file, character(1L)) == manifest$sha256)) {
    stop("MV8-ZL manifest drift: ", name, call. = FALSE)
  }
}
verify_manifest(zl_root, "mv08zl-artifact-manifest.csv")
verify_manifest(zf_root, "mv08zf-artifact-manifest.csv")

snapshot <- read_csv(file.path(zl_root, "mv08zl-stopped-snapshot.csv"))
binding <- read_csv(file.path(zl_root, "mv08zl-order-280-binding.csv"))
policy <- read_csv(file.path(zl_root, "mv08zl-recovery-policy.csv"))
decision <- read_csv(file.path(zl_root, "mv08zl-decision.csv"))
implementation <- read_csv(file.path(zl_root, "mv08zl-implementation-bindings.csv"))
queue <- read_csv(file.path(zf_root, "mv08zf-production-queue.csv"))
ledger_path <- file.path(public_root, "mv08zf-resource-ledger.csv")
completion_path <- file.path(public_root, "mv08zf-chunk-completions.csv")
completion_partial_path <- paste0(completion_path, ".partial")
progress_path <- file.path(public_root, "mv08zf-progress.csv")
progress_partial_path <- paste0(progress_path, ".partial")
executor_file <- "scripts/recover_mv08zl_landscape_receipt_publication.R"
bound_executor <- implementation[
  implementation$role == "receipt_recovery_executor", , drop = FALSE
]
current_head <- tolower(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
if (nrow(decision) != 1L || !truth(decision$receipt_promotion_authorized) ||
    decision$adopted_production_order != 280L ||
    decision$resume_at_production_order != 281L ||
    truth(decision$landscape_recomputation_authorized) ||
    decision$automatic_retries != 0L || current_head != recovery_head ||
    nrow(bound_executor) != 1L || bound_executor$file != executor_file ||
    sha_file(executor_file) != bound_executor$sha256) {
  stop("MV8-ZL committed authorization drift", call. = FALSE)
}
runner_lines <- suppressWarnings(system2(
  "pgrep", c("-f", "[r]un_mv08zf_full_landscape_production[.]R"),
  stdout = TRUE, stderr = FALSE
))
if (length(runner_lines)) stop("MV8-ZL runner is active", call. = FALSE)

ledger <- read_csv(ledger_path)
progress <- read_csv(progress_path)
if (nrow(ledger) != 280L ||
    !identical(as.integer(ledger$production_order), 1:280) ||
    sha_file(ledger_path) != snapshot$ledger_sha256 ||
    any(ledger$disposition != "completed") || any(ledger$exit_status != 0L) ||
    any(ledger$workers != 1L) || any(ledger$retries != 0L)) {
  stop("MV8-ZL durable ledger drift", call. = FALSE)
}

validate_completion <- function(value, rows) {
  nrow(value) == rows &&
    identical(as.integer(value$production_order), seq_len(rows)) &&
    all(value$workers == 1L) && all(value$retries == 0L) &&
    all(value$outcome_label_state == "closed") &&
    !any(truth(value$biological_outcomes_computed))
}
final <- read_csv(completion_path)
partial_exists <- file.exists(completion_partial_path)
partial <- if (partial_exists) read_csv(completion_partial_path) else data.frame()
state <- if (validate_completion(final, 279L) && partial_exists &&
             validate_completion(partial, 280L)) {
  "pre_promotion"
} else if (validate_completion(final, 280L) && !partial_exists &&
           progress$completed_chunks == 279L) {
  "completion_promoted_progress_pending"
} else if (validate_completion(final, 280L) && !partial_exists &&
           progress$completed_chunks == 280L) {
  "already_complete"
} else {
  stop("MV8-ZL publication state is not an accepted recovery boundary", call. = FALSE)
}

if (state == "pre_promotion") {
  if (sha_file(completion_path) != snapshot$completion_sha256 ||
      sha_file(completion_partial_path) != snapshot$completion_partial_sha256 ||
      sha_file(progress_path) != snapshot$progress_sha256 ||
      !identical(final, partial[seq_len(279L), , drop = FALSE])) {
    stop("MV8-ZL frozen pre-promotion state drift", call. = FALSE)
  }
  recovery_root <- file.path(private_root, "recovery", "mv08zl")
  dir.create(recovery_root, recursive = TRUE, showWarnings = FALSE)
  preserved_prefix <- file.path(recovery_root, "mv08zf-chunk-completions-prefix-279.csv")
  if (file.exists(preserved_prefix)) {
    if (sha_file(preserved_prefix) != snapshot$completion_sha256) {
      stop("MV8-ZL preserved prefix drift", call. = FALSE)
    }
  } else if (!file.copy(completion_path, preserved_prefix, overwrite = FALSE) ||
             sha_file(preserved_prefix) != snapshot$completion_sha256) {
    stop("MV8-ZL could not preserve the 279-row prefix", call. = FALSE)
  }
  if (!file.rename(completion_partial_path, completion_path)) {
    stop("MV8-ZL exact completion promotion failed", call. = FALSE)
  }
  final <- read_csv(completion_path)
  if (!validate_completion(final, 280L) ||
      final$distances_sha256[[280L]] != binding$distances_sha256 ||
      final$status_sha256[[280L]] != binding$status_sha256) {
    stop("MV8-ZL promoted completion failed validation", call. = FALSE)
  }
}

if (progress$completed_chunks == 279L) {
  progress$completed_chunks <- 280L
  progress$completed_pairs <- sum(as.integer(final$pair_count))
  progress$aggregate_child_seconds <- sum(as.numeric(ledger$elapsed_seconds))
  private_files <- list.files(private_root, recursive = TRUE, full.names = TRUE,
                              all.files = TRUE, no.. = TRUE)
  private_files <- private_files[!file.info(private_files)$isdir]
  progress$private_bytes <- sum(as.numeric(file.info(private_files)$size))
  progress$state <- "running"
  atomic_csv(progress, progress_path)
}

final <- read_csv(completion_path)
progress <- read_csv(progress_path)
partials <- c(
  list.files(public_root, pattern = "[.]partial$", full.names = TRUE),
  list.files(private_root, pattern = "(__partial__|[.]partial$)", recursive = TRUE,
             full.names = TRUE, all.files = TRUE)
)
if (!validate_completion(final, 280L) || nrow(ledger) != 280L ||
    progress$completed_chunks != 280L ||
    progress$completed_pairs != sum(final$pair_count) ||
    progress$aggregate_child_seconds != sum(ledger$elapsed_seconds) ||
    length(partials) != 0L) {
  stop("MV8-ZL final receipt recovery validation failed", call. = FALSE)
}
cat("MV8-ZL adopted exact completed receipt 280 without retry; resume_at=281\n")
