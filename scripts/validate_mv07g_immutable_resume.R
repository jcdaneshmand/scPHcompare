#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("usage: validate_mv07g_immutable_resume.R MODE PREFREEZE PRIVATE_ROOT PUBLIC_DIR STATE_RDS OUTPUT")
}
mode <- match.arg(args[[1]], c("snapshot", "compare"))
prefreeze <- args[[2]]
private_root <- args[[3]]
public_dir <- args[[4]]
state_path <- args[[5]]
output <- args[[6]]
queue <- read.csv(file.path(prefreeze, "mv07g-queue.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
repeat_queue <- queue[queue$seed == min(queue$seed), , drop = FALSE]
private_relative <- c(
  queue$output_file,
  file.path("repeat", repeat_queue$output_file)
)
public_relative <- c(
  "mv07g-source-metrics.csv", "mv07g-ph-metrics.csv",
  "mv07g-repeat-validation.csv", "mv07g-full-ph-projection.csv",
  "mv07g-decision.csv"
)
inventory <- function() {
  relative <- c(private_relative, file.path("public", public_relative))
  paths <- c(file.path(private_root, private_relative),
             file.path(public_dir, public_relative))
  if (length(paths) != 83L || anyDuplicated(relative) ||
      !all(file.exists(paths))) {
    stop("MV7-G immutable-resume axis is incomplete.")
  }
  info <- file.info(paths)
  data.frame(
    contract_id = "mv07g_immutable_resume_inventory_v1",
    scope = c(rep("private_scientific_artifact", length(private_relative)),
              rep("public_production_evidence", length(public_relative))),
    relative_file = relative,
    bytes = as.numeric(info$size),
    sha256 = vapply(paths, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE), character(1L)),
    mtime_unix_seconds = as.numeric(info$mtime),
    stringsAsFactors = FALSE
  )
}
if (mode == "snapshot") {
  state <- inventory()
  dir.create(dirname(state_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(state, state_path, compress = FALSE, version = 3)
  message("MV7-G immutable-resume snapshot: 83 artifacts")
} else {
  if (!file.exists(state_path)) stop("Missing MV7-G resume snapshot.")
  before <- readRDS(state_path)
  after <- inventory()
  axis_ok <- identical(before$scope, after$scope) &&
    identical(before$relative_file, after$relative_file)
  hash_ok <- identical(before$sha256, after$sha256)
  bytes_ok <- identical(before$bytes, after$bytes)
  mtime_delta <- max(abs(before$mtime_unix_seconds - after$mtime_unix_seconds))
  mtime_ok <- identical(before$mtime_unix_seconds, after$mtime_unix_seconds)
  private <- before$scope == "private_scientific_artifact"
  public <- before$scope == "public_production_evidence"
  checks <- data.frame(
    contract_id = "mv07g_immutable_resume_validation_v1",
    category = c("artifact_axis", "private_hashes", "private_bytes",
                 "private_mtimes", "public_hashes", "public_bytes",
                 "public_mtimes"),
    passed = c(axis_ok, hash_ok && all(before$sha256[private] ==
      after$sha256[private]), bytes_ok && all(before$bytes[private] ==
      after$bytes[private]), mtime_ok && all(before$mtime_unix_seconds[private] ==
      after$mtime_unix_seconds[private]), hash_ok && all(before$sha256[public] ==
      after$sha256[public]), bytes_ok && all(before$bytes[public] ==
      after$bytes[public]), mtime_ok && all(before$mtime_unix_seconds[public] ==
      after$mtime_unix_seconds[public])),
    detail = c("78 private plus 5 public", "78 unchanged", "78 unchanged",
               "78 unchanged", "5 unchanged", "5 unchanged", "5 unchanged"),
    maximum_mtime_delta_seconds = mtime_delta,
    stringsAsFactors = FALSE
  )
  if (!all(checks$passed)) {
    stop("MV7-G immutable-resume validation failed: ",
         paste(checks$category[!checks$passed], collapse = ", "))
  }
  write.csv(checks, output, row.names = FALSE, na = "")
  message("MV7-G immutable-resume validation: 7/7 pass")
}
