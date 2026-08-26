#!/usr/bin/env Rscript
options(warn = 2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste("usage: build_mv17b_null_qualification_closure.R",
  "<prefreeze> <recovery-prefreeze> <primary> <repeat> <output>"), call. = FALSE)
prefreeze <- normalizePath(args[[1L]], mustWork = TRUE)
recovery <- normalizePath(args[[2L]], mustWork = TRUE)
primary <- normalizePath(args[[3L]], mustWork = TRUE)
repeat_root <- normalizePath(args[[4L]], mustWork = TRUE); output <- args[[5L]]
if (dir.exists(output)) stop("MV17-B closure output exists", call. = FALSE)
dir.create(output, recursive = TRUE)
source("R/mv08z_landscape_production.R")
readc <- .mv08z_read_csv; atomic <- .mv08z_atomic_csv; sha <- .mv08z_sha256_file
verify_manifest <- function(root, name) {
  manifest <- readc(file.path(root, name)); paths <- file.path(root, manifest$artifact)
  if (any(!file.exists(paths)) ||
      !identical(unname(as.numeric(file.info(paths)$size)), unname(as.numeric(manifest$bytes))) ||
      !identical(unname(vapply(paths, sha, character(1L))), unname(manifest$sha256)))
    stop("MV17-B bound manifest drift", call. = FALSE)
  manifest
}
pre_manifest <- verify_manifest(prefreeze, "mv17b-artifact-manifest.csv")
recovery_manifest <- verify_manifest(recovery, "mv17b-recovery-artifact-manifest.csv")
contract <- readc(file.path(prefreeze, "mv17b-contract.csv"))
recovery_contract <- readc(file.path(recovery, "mv17b-recovery-contract.csv"))
primary_paths <- file.path(primary, c("mv17b-qualification.csv", "mv17b-status.csv"))
repeat_paths <- file.path(repeat_root, c("mv17b-qualification.csv", "mv17b-status.csv"))
bindings <- readc(file.path(recovery, "mv17b-recovery-production-binding.csv"))
actual_paths <- c(primary_paths, repeat_paths)
if (!identical(unname(vapply(actual_paths, sha, character(1L))), unname(bindings$sha256)) ||
    !identical(unname(as.numeric(file.info(actual_paths)$size)), unname(as.numeric(bindings$bytes))))
  stop("MV17-B preserved production binding drift", call. = FALSE)
x <- readc(primary_paths[[1L]]); y <- readc(repeat_paths[[1L]])
s1 <- readc(primary_paths[[2L]]); s2 <- readc(repeat_paths[[2L]])
required <- c("fixture", "seed", "null_family", "invariant_error", "deterministic",
              "pair_distance_spearman", "neighbor_jaccard", "passed_invariant")
grid <- expand.grid(fixture = c("circle_in_3d", "correlated_gaussian", "independent_gaussian"),
  seed = c(17001L, 17002L, 17003L), null_family = c("coordinate_permutation",
  "covariance_gaussian", "radial_density_cloud", "within_row_axis_shuffle"),
  stringsAsFactors = FALSE)
key <- function(z) paste(z$fixture, z$seed, z$null_family, sep = "\r")
cardinality_ok <- nrow(x) == 36L && all(required %in% names(x)) &&
  setequal(key(x), key(grid)) && !anyDuplicated(key(x))
numeric_ok <- all(is.finite(x$invariant_error)) && all(is.finite(x$pair_distance_spearman)) &&
  all(is.finite(x$neighbor_jaccard)) && all(x$neighbor_jaccard >= 0 & x$neighbor_jaccard <= 1)
repeat_ok <- identical(x, y) && identical(s1, s2) &&
  sha(primary_paths[[1L]]) == sha(repeat_paths[[1L]]) &&
  sha(primary_paths[[2L]]) == sha(repeat_paths[[2L]])
circle_ok <- all(x$neighbor_jaccard[x$fixture == "circle_in_3d"] <= contract$maximum_circle_neighbor_jaccard)
status_ok <- s1$jobs == 36L && s1$all_invariants && s1$all_deterministic &&
  s1$real_corpus_jobs == 0L && s1$PH_jobs == 0L && !s1$labels_opened && !s1$outcomes_opened
families <- split(seq_len(nrow(x)), x$null_family)
summary <- do.call(rbind, lapply(names(families), function(family) {
  z <- x[families[[family]], , drop = FALSE]
  data.frame(contract_id = "mv17b_family_summary_v2", null_family = family,
    jobs = nrow(z), maximum_invariant_error = max(z$invariant_error),
    minimum_pair_distance_spearman = min(z$pair_distance_spearman),
    maximum_pair_distance_spearman = max(z$pair_distance_spearman),
    minimum_neighbor_jaccard = min(z$neighbor_jaccard),
    maximum_neighbor_jaccard = max(z$neighbor_jaccard),
    maximum_circle_neighbor_jaccard = max(z$neighbor_jaccard[z$fixture == "circle_in_3d"]),
    all_invariants = all(z$passed_invariant), all_deterministic = all(z$deterministic))
})); rownames(summary) <- NULL
validation <- data.frame(contract_id = "mv17b_closure_validation_v2",
  check_id = c("prefreeze_manifest", "recovery_manifest", "production_binding",
    "36_exact_jobs", "complete_grid", "numeric_metrics", "all_invariants",
    "all_deterministic", "byte_exact_independent_repeat", "circle_separation",
    "real_corpus_PH_closed", "labels_outcomes_closed", "downstream_firewall"),
  passed = c(nrow(pre_manifest)>0L, nrow(recovery_manifest)>0L, nrow(bindings)==4L,
    cardinality_ok, cardinality_ok, numeric_ok, all(x$passed_invariant),
    all(x$deterministic), repeat_ok, circle_ok, status_ok, status_ok,
    !recovery_contract$clustering_authorized && !recovery_contract$fusion_authorized &&
      !recovery_contract$biology_authorized && !recovery_contract$manuscript_claims_authorized))
if (!all(validation$passed)) stop("MV17-B v2 closure failed", call. = FALSE)
decision <- data.frame(contract_id = "mv17b_closure_decision_v2", MV17B_closed = TRUE,
  admitted_null_families = paste(sort(unique(x$null_family)), collapse = ";"),
  MV17C_implementation_eligible = TRUE, real_calibration_authorized = FALSE,
  real_localization_authorized = FALSE, real_H2_authorized = FALSE,
  labels_outcomes_clustering_fusion_claims = "closed")
artifacts <- list("mv17b-complete-qualification.csv"=x, "mv17b-family-summary.csv"=summary,
  "mv17b-production-binding.csv"=bindings, "mv17b-validation.csv"=validation,
  "mv17b-decision.csv"=decision)
for (name in names(artifacts)) atomic(artifacts[[name]], file.path(output, name))
writeLines(c("# MV17-B null-model qualification closure v2", "",
  "All 36 synthetic jobs pass their fixed invariants, exact repeat, and circle-separation gate.",
  "Complete metrics and all primary/repeat hashes are bound.", "",
  "All four families qualify for their declared compatible view only. Real calibration remains closed."),
  file.path(output, "MV17B_NULL_QUALIFICATION_CLOSURE_2026-08-26.md"))
files <- sort(list.files(output)); manifest <- data.frame(
  contract_id="mv17b_closure_artifact_manifest_v2", artifact=files,
  bytes=as.numeric(file.info(file.path(output,files))$size),
  sha256=vapply(file.path(output,files),sha,character(1L)))
atomic(manifest,file.path(output,"mv17b-artifact-manifest.csv"))
message("Built MV17-B v2 closure; checks=",nrow(validation))
