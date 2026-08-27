#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    paste(
      "usage: quarantine_mv17g_interrupted_wave.R",
      "<original-private-prefreeze> <recovery-public-prefreeze>",
      "<recovery-private-prefreeze> <production-private-root> <quarantine-root>"
    ),
    call. = FALSE
  )
}
original_private <- normalizePath(args[[1]], mustWork = TRUE)
recovery_public <- normalizePath(args[[2]], mustWork = TRUE)
recovery_private <- normalizePath(args[[3]], mustWork = TRUE)
production <- normalizePath(args[[4]], mustWork = TRUE)
quarantine <- args[[5]]
if (dir.exists(quarantine)) stop("MV17-G interrupted-wave quarantine already exists", call. = FALSE)

source("R/mv08z_landscape_production.R")
source("R/mv17g_parallel_recovery.R")
r <- .mv08z_read_csv
w <- .mv08z_atomic_csv
h <- .mv08z_sha256_file

manifest <- r(file.path(recovery_public, "mv17g-controller-recovery-artifact-manifest.csv"))
manifest_paths <- file.path(recovery_public, manifest$artifact)
if (!all(file.exists(manifest_paths)) ||
    !identical(unname(as.numeric(file.info(manifest_paths)$size)), unname(as.numeric(manifest$bytes))) ||
    !identical(unname(vapply(manifest_paths, h, character(1L))), unname(tolower(manifest$sha256)))) {
  stop("MV17-G recovery public manifest drift", call. = FALSE)
}
private_binding <- r(file.path(recovery_public, "mv17g-controller-recovery-private-binding.csv"))
private_paths <- file.path(recovery_private, private_binding$artifact)
if (!all(file.exists(private_paths)) ||
    !identical(unname(as.numeric(file.info(private_paths)$size)), unname(as.numeric(private_binding$bytes))) ||
    !identical(unname(vapply(private_paths, h, character(1L))), unname(tolower(private_binding$sha256)))) {
  stop("MV17-G recovery private binding drift", call. = FALSE)
}
contract <- r(file.path(recovery_public, "mv17g-controller-recovery-contract.csv"))
decision <- r(file.path(recovery_public, "mv17g-controller-recovery-decision.csv"))
if (!contract$owner_authorization_recorded || !contract$launch_authorized_after_prefreeze_commit ||
    contract$manual_replay_exception_children != 8L ||
    decision$disposition != "preserve_and_replay_interrupted_wave_787_794_then_continue") {
  stop("MV17-G interrupted-wave replay is not authorized", call. = FALSE)
}

queue <- r(file.path(original_private, "mv17g-primary-grouped-queue.csv"))
wave <- r(file.path(recovery_private, "mv17g-controller-interrupted-wave-artifacts.csv"))
if (nrow(wave) != 32L || !identical(sort(unique(as.integer(wave$job_order))), 787:794) ||
    !all(file.exists(wave$artifact)) ||
    !identical(unname(as.numeric(file.info(wave$artifact)$size)), unname(as.numeric(wave$bytes))) ||
    !identical(unname(vapply(wave$artifact, h, character(1L))), unname(tolower(wave$sha256)))) {
  stop("MV17-G interrupted-wave source drift", call. = FALSE)
}
expected <- do.call(rbind, lapply(787:794, function(i) {
  q <- queue[i, , drop = FALSE]
  paths <- mv17g_job_artifacts_v1(q, production)
  data.frame(job_order = i, role = names(paths), artifact = unname(paths), stringsAsFactors = FALSE)
}))
expected <- expected[order(expected$job_order, match(expected$role, c("result", "time", "stdout", "stderr"))), ]
wave <- wave[order(wave$job_order, match(wave$role, c("result", "time", "stdout", "stderr"))), ]
if (!identical(normalizePath(expected$artifact), normalizePath(wave$artifact))) {
  stop("MV17-G interrupted-wave production path drift", call. = FALSE)
}

dir.create(file.path(quarantine, "jobs"), recursive = TRUE)
dir.create(file.path(quarantine, "logs"), recursive = TRUE)
destinations <- file.path(
  quarantine,
  ifelse(wave$role == "result", "jobs", "logs"),
  basename(wave$artifact)
)
moved <- logical(nrow(wave))
rollback <- TRUE
on.exit({
  if (rollback && any(moved)) {
    for (i in rev(which(moved))) file.rename(destinations[i], wave$artifact[i])
  }
}, add = TRUE)
for (i in seq_len(nrow(wave))) {
  if (!file.rename(wave$artifact[i], destinations[i])) {
    stop("MV17-G interrupted-wave quarantine move failed", call. = FALSE)
  }
  moved[i] <- TRUE
}
if (!all(file.exists(destinations)) || any(file.exists(wave$artifact)) ||
    !identical(unname(as.numeric(file.info(destinations)$size)), unname(as.numeric(wave$bytes))) ||
    !identical(unname(vapply(destinations, h, character(1L))), unname(tolower(wave$sha256)))) {
  stop("MV17-G interrupted-wave quarantine verification failed", call. = FALSE)
}
quarantine_manifest <- data.frame(
  contract_id = "mv17g_interrupted_wave_quarantine_v1",
  job_order = wave$job_order,
  role = wave$role,
  original_artifact = wave$artifact,
  quarantine_artifact = normalizePath(destinations),
  bytes = wave$bytes,
  sha256 = wave$sha256,
  time_receipt_valid = wave$time_receipt_valid,
  stringsAsFactors = FALSE
)
w(quarantine_manifest, file.path(quarantine, "mv17g-interrupted-wave-quarantine.csv"))
rollback <- FALSE
message("Quarantined 32 MV17-G artifacts for replay of orders 787--794")
