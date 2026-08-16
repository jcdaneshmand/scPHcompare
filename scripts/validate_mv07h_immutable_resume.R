#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "usage: validate_mv07h_immutable_resume.R MODE PRIVATE_ROOT",
    "PUBLIC_DIR STATE_RDS OUTPUT"))
}
mode <- match.arg(args[[1L]], c("snapshot", "compare"))
private_root <- args[[2L]]
public_dir <- args[[3L]]
state_path <- args[[4L]]
output <- args[[5L]]
public_relative <- c(
  "mv07h-source-metrics.csv", "mv07h-ph-metrics.csv",
  "mv07h-ph-engine-attempts.csv", "mv07h-sentinel-equivalence.csv",
  "mv07h-repeat-validation.csv", "mv07h-full-ph-decision.csv",
  "mv07h-independent-validation.csv", "mv07h-cross-engine-validation.csv",
  "mv07h-validation-decision.csv"
)
inventory <- function() {
  private_paths <- sort(list.files(
    private_root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
    no.. = TRUE), method = "radix")
  private_paths <- private_paths[file.info(private_paths)$isdir %in% FALSE]
  normalized_private_root <- normalizePath(private_root, winslash = "/",
                                           mustWork = TRUE)
  private_paths <- normalizePath(private_paths, winslash = "/",
                                 mustWork = TRUE)
  private_relative <- substring(private_paths,
                                nchar(normalized_private_root) + 2L)
  private_relative <- chartr("\\", "/", private_relative)
  public_paths <- file.path(public_dir, public_relative)
  relative <- c(private_relative, file.path("public", public_relative))
  paths <- c(private_paths, public_paths)
  if (!length(private_paths) || anyDuplicated(relative) ||
      !all(file.exists(paths))) {
    stop("MV7-H immutable-resume axis is incomplete.")
  }
  info <- file.info(paths)
  data.frame(
    contract_id = "mv07h_immutable_resume_inventory_v1",
    scope = c(rep("private_artifact_ledger_or_log", length(private_paths)),
              rep("public_production_or_validation_evidence",
                  length(public_paths))),
    relative_file = relative, bytes = as.numeric(info$size),
    sha256 = vapply(paths, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE), character(1L)),
    mtime_unix_seconds = as.numeric(info$mtime), stringsAsFactors = FALSE
  )
}
if (mode == "snapshot") {
  state <- inventory()
  if (file.exists(state_path)) stop("Refusing overwrite: ", state_path)
  dir.create(dirname(state_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(state, state_path, compress = FALSE, version = 3)
  message("MV7-H immutable-resume snapshot: ", nrow(state), " files")
} else {
  if (!file.exists(state_path)) stop("Missing MV7-H resume snapshot.")
  if (file.exists(output)) stop("Refusing overwrite: ", output)
  before <- readRDS(state_path)
  after <- inventory()
  axis_ok <- identical(before$scope, after$scope) &&
    identical(before$relative_file, after$relative_file)
  hash_ok <- axis_ok && identical(before$sha256, after$sha256)
  bytes_ok <- axis_ok && identical(before$bytes, after$bytes)
  mtime_ok <- axis_ok &&
    identical(before$mtime_unix_seconds, after$mtime_unix_seconds)
  private <- before$scope == "private_artifact_ledger_or_log"
  public <- before$scope == "public_production_or_validation_evidence"
  maximum_mtime_delta <- if (axis_ok) {
    max(abs(before$mtime_unix_seconds - after$mtime_unix_seconds))
  } else Inf
  checks <- data.frame(
    contract_id = "mv07h_immutable_resume_validation_v1",
    category = c("artifact_axis", "private_hashes", "private_bytes",
                 "private_mtimes", "public_hashes", "public_bytes",
                 "public_mtimes"),
    passed = c(
      axis_ok,
      hash_ok && all(before$sha256[private] == after$sha256[private]),
      bytes_ok && all(before$bytes[private] == after$bytes[private]),
      mtime_ok && all(before$mtime_unix_seconds[private] ==
                        after$mtime_unix_seconds[private]),
      hash_ok && all(before$sha256[public] == after$sha256[public]),
      bytes_ok && all(before$bytes[public] == after$bytes[public]),
      mtime_ok && all(before$mtime_unix_seconds[public] ==
                        after$mtime_unix_seconds[public])),
    private_files = sum(private), public_files = sum(public),
    maximum_mtime_delta_seconds = maximum_mtime_delta,
    outcome_label_state = "closed", biological_outcomes_computed = FALSE,
    stringsAsFactors = FALSE
  )
  if (!all(checks$passed)) stop(
    "MV7-H immutable-resume validation failed: ",
    paste(checks$category[!checks$passed], collapse = ", "))
  write.csv(checks, output, row.names = FALSE, na = "")
  message("MV7-H immutable-resume validation: 7/7 pass across ",
          nrow(before), " files")
}
