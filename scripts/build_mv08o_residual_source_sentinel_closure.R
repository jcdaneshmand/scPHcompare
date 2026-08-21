#!/usr/bin/env Rscript

# Bind independently completed, resource-capped MV8-O source attempts into one
# public closure. Input paths stay in a private manifest; this script publishes
# only hashes, aggregate resource records, and geometry/determinism evidence.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop(
  "usage: build_mv08o_residual_source_sentinel_closure.R <private-attempt-manifest.csv> <public-output-dir>",
  call. = FALSE
)
manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
public_dir <- normalizePath(args[[2L]], mustWork = FALSE)
if (dir.exists(public_dir)) stop("refusing to overwrite MV8-O closure output", call. = FALSE)
for (package in c("digest")) if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required", call. = FALSE)
dir.create(public_dir, recursive = TRUE)

sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  utils::write.csv(x, partial, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
atomic_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("failed to publish ", basename(path), call. = FALSE)
}
classify_stderr <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!nzchar(text)) return("empty")
  known <- grepl("vst.flavor.*glmGamPoi", text) &&
    grepl("Falling back to native \\(slower\\) implementation", text) &&
    !grepl("(^|\\n)(Error|Execution halted)", text)
  if (known) "known_glmGamPoi_native_fallback" else "unexpected_nonempty"
}

attempts <- utils::read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("attempt_id", "unit_id", "run_role", "private_root", "public_root")
if (!identical(names(attempts), needed) || nrow(attempts) != 5L || anyDuplicated(attempts[c("unit_id", "run_role")])) {
  stop("invalid MV8-O private attempt manifest", call. = FALSE)
}
expected <- data.frame(
  unit_id = c("SRA701877_SRS3279688", "SRA742961_SRS3565197", "HCA_BM_002", "SRA742961_SRS3565197", "HCA_BM_002"),
  run_role = c("primary", "primary", "primary", "repeat", "repeat"), stringsAsFactors = FALSE
)
if (!identical(attempts[c("unit_id", "run_role")], expected)) stop("MV8-O closure attempt order drift", call. = FALSE)

resource_rows <- vector("list", nrow(attempts)); audit_rows <- vector("list", nrow(attempts))
for (i in seq_len(nrow(attempts))) {
  a <- attempts[i, , drop = FALSE]
  private_root <- normalizePath(a$private_root[[1L]], mustWork = TRUE)
  public_root <- normalizePath(a$public_root[[1L]], mustWork = TRUE)
  resource_path <- file.path(public_root, "mv08o-source-sentinel-resource.csv")
  resource <- utils::read.csv(resource_path, check.names = FALSE, stringsAsFactors = FALSE)
  resource <- resource[resource$unit_id == a$unit_id & resource$run_role == a$run_role, , drop = FALSE]
  if (nrow(resource) != 1L) stop("attempt resource identity drift for ", a$attempt_id, call. = FALSE)
  stem <- paste0(a$unit_id, "__", a$run_role)
  cache_path <- file.path(private_root, "cache", paste0(stem, ".rds"))
  audit_path <- file.path(private_root, "worker-audit", paste0(stem, ".csv"))
  stderr_path <- file.path(private_root, "logs", paste0(stem, "-stderr.txt"))
  if (!all(file.exists(c(cache_path, audit_path, stderr_path)))) stop("private evidence absent for ", a$attempt_id, call. = FALSE)
  if (!identical(tolower(sha_file(cache_path)), tolower(resource$cache_sha256[[1L]])) ||
      !identical(tolower(sha_file(audit_path)), tolower(resource$worker_audit_sha256[[1L]]))) {
    stop("private evidence hash drift for ", a$attempt_id, call. = FALSE)
  }
  audit <- utils::read.csv(audit_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(audit) != resource$worker_rows[[1L]]) stop("worker row-count drift for ", a$attempt_id, call. = FALSE)
  resource$attempt_id <- a$attempt_id
  resource$stderr_class <- classify_stderr(stderr_path)
  resource_rows[[i]] <- resource
  audit_rows[[i]] <- audit
}
resource <- do.call(rbind, resource_rows)
views <- do.call(rbind, audit_rows)
required <- views$dataset_scope == "internal124" | views$panel_id == "common475" |
  (views$panel_id == "exact500" & views$representation_id == "sct_pearson_residual_all_qc_fit_selected384")
diagnostic <- views$dataset_scope == "external8" & views$panel_id == "exact500" &
  views$representation_id == "sct_data_all_qc_fit_selected384"
primary <- views[views$run_role == "primary", , drop = FALSE]
repeat_views <- views[views$run_role == "repeat", , drop = FALSE]
repeat_rows <- merge(
  primary[primary$unit_id %in% c("SRA742961_SRS3565197", "HCA_BM_002"),
          c("unit_id", "seed", "representation_id", "panel_id", "matrix_sha256", "distance_sha256")],
  repeat_views[, c("unit_id", "seed", "representation_id", "panel_id", "matrix_sha256", "distance_sha256")],
  by = c("unit_id", "seed", "representation_id", "panel_id"), suffixes = c("_primary", "_repeat")
)
determinism <- data.frame(
  contract_id = "mv08o_source_sentinel_determinism_v1", repeat_rows,
  passed = repeat_rows$matrix_sha256_primary == repeat_rows$matrix_sha256_repeat &
    ((is.na(repeat_rows$distance_sha256_primary) & is.na(repeat_rows$distance_sha256_repeat)) |
      repeat_rows$distance_sha256_primary == repeat_rows$distance_sha256_repeat),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  check_id = c("authorized_runs", "five_runs_complete", "resource_caps", "stderr_policy",
    "private_evidence_rehashed", "required_geometry", "external_data_exact500_diagnostic_only",
    "internal_five_seed_axes", "hca_bridge_axis", "repeat_determinism", "downstream_firewall"),
  passed = c(
    identical(as.character(resource$unit_id), as.character(expected$unit_id)) &&
      identical(as.character(resource$run_role), as.character(expected$run_role)),
    nrow(resource) == 5L && all(resource$disposition == "completed"),
    all(resource$elapsed_seconds <= resource$elapsed_cap_seconds & resource$peak_rss_bytes <= resource$rss_cap_bytes),
    all(resource$stderr_class %in% c("empty", "known_glmGamPoi_native_fallback")),
    TRUE,
    all(views$values_finite[required] & views$zero_variance_gene_count[required] == 0L & views$correlation_chord_valid[required]),
    sum(diagnostic) == 2L && all(!views$correlation_chord_valid[diagnostic]),
    nrow(primary[primary$dataset_scope == "internal124", ]) == 20L,
    nrow(primary[primary$unit_id == "HCA_BM_002", ]) == 4L,
    nrow(determinism) == 14L && all(determinism$passed),
    !any(resource$persistence_computed | resource$landscapes_computed | resource$clustering_computed |
      resource$fusion_computed | resource$biological_outcomes_computed)
  ),
  evidence = c(
    "three MV8-N authorized units, including two required repeats", "all five source attempts completed",
    "serial 1,800-second and 12-GiB caps", "only documented glmGamPoi-native fallback diagnostics; no errors",
    "private cache and worker-audit hashes independently rechecked", "all required internal exact500 and external common475/residual views pass",
    "external SCT-data exact500 is recorded as invalid diagnostic-only, never PH eligible",
    "two internal samples x five seeds x two representations", "HCA_BM_002 one seed x two representations x two panels",
    "maximum internal and HCA matrix/distance hashes match repeats", "no PH, landscapes, clustering, fusion, labels, or outcomes"
  ), stringsAsFactors = FALSE
)
if (!all(validation$passed)) stop("MV8-O closure validation failed", call. = FALSE)
atomic_csv(resource, file.path(public_dir, "mv08o-source-sentinel-resource.csv"))
atomic_csv(views, file.path(public_dir, "mv08o-source-sentinel-view-summary.csv"))
atomic_csv(determinism, file.path(public_dir, "mv08o-source-sentinel-determinism.csv"))
atomic_csv(validation, file.path(public_dir, "mv08o-source-sentinel-validation.csv"))
report <- c(
  "# MV8-O Pearson-residual source/view sentinel closure", "",
  "## Result", "",
  "The bounded source/view sentinel passed. It validates source-scale feasibility and deterministic representation construction only; persistence, landscapes, clustering, fusion, labels, outcomes, and default adoption remain closed.", "",
  "## Contract evidence", "",
  "- Five serial runs: internal minimum primary; internal maximum primary/repeat; HCA_BM_002 primary/repeat.",
  "- Every run stayed inside the 1,800-second and 12-GiB limits.",
  "- Standard SCTransform used the documented native fallback because glmGamPoi is unavailable; this is logged as a known performance diagnostic, not an error or a changed transform.",
  "- Required geometry passed for internal exact500 and external common475 plus exact500 Pearson residuals.",
  "- External SCT-data exact500 remained the expected diagnostic-only invalid view; it was not eligible for PH.",
  "- Repeats matched matrix and eligible distance hashes exactly.", "",
  "## Scope firewall", "",
  "This closure does not compute persistence diagrams, landscapes, pairwise comparisons, clustering, fusion, labels, biological outcomes, or claims. It does not make the Pearson-residual representation the project default.", "",
  "## Next gate", "",
  "A separately prefrozen full-source production decision is required before any remaining source fits or topology execution. The maximum-source margin should be considered explicitly in that decision."
)
atomic_text(report, file.path(public_dir, "MV08O_RESIDUAL_SOURCE_SENTINEL_CLOSURE_2026-08-21.md"))
artifacts <- list.files(public_dir, full.names = TRUE)
artifacts <- artifacts[basename(artifacts) != "mv08o-artifact-manifest.csv"]
manifest <- data.frame(artifact = basename(artifacts), bytes = file.info(artifacts)$size,
  sha256 = vapply(artifacts, sha_file, character(1L)), stringsAsFactors = FALSE)
atomic_csv(manifest, file.path(public_dir, "mv08o-artifact-manifest.csv"))
cat("MV8-O closure checks=", sum(validation$passed), "/", nrow(validation),
  "; all source/view gates pass; PH remains closed\n", sep = "")
