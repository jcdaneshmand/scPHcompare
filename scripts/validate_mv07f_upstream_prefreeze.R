#!/usr/bin/env Rscript

options(warn = 2)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest required.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: validate_mv07f_upstream_prefreeze.R AUDIT_DIR SOURCE_ROOT ",
       "EXPECTED_HEAD OUTPUT", call. = FALSE)
}
audit_dir <- args[[1L]]
source_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
expected_head <- tolower(trimws(args[[3L]])); output <- args[[4L]]
readc <- function(name) read.csv(file.path(audit_dir, name),
  stringsAsFactors = FALSE, check.names = FALSE)
truth <- function(x) if (is.logical(x)) !is.na(x) & x else
  tolower(trimws(as.character(x))) == "true"
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
required <- c("mv07f-source-freeze.csv", "mv07f-source-coverage.csv",
  "mv07f-upstream-queue.csv", "mv07f-resource-caps.csv",
  "mv07f-repeat-target.csv", "mv07f-acceptance-criteria.csv",
  "mv07f-decision.csv")
if (any(!file.exists(file.path(audit_dir, required)))) stop("MV7-F prefreeze incomplete.")
s <- readc(required[1]); coverage <- readc(required[2]); queue <- readc(required[3])
caps <- readc(required[4]); target <- readc(required[5]); criteria <- readc(required[6])
decision <- readc(required[7])
source_ok <- nrow(s) == 19L && all(s$accepted_head == expected_head) &&
  !any(truth(s$private_source)) && !any(grepl("^/|^[A-Za-z]:", s$artifact_locator)) &&
  all(file.exists(s$artifact_locator)) &&
  all(vapply(seq_len(nrow(s)), function(i)
    sha(s$artifact_locator[[i]]) == s$sha256[[i]], logical(1L)))
paths <- list.files(source_root, pattern = "[.]rds$", recursive = TRUE,
  full.names = TRUE, ignore.case = TRUE)
ids <- tools::file_path_sans_ext(basename(paths))
samples <- sort(unique(queue$sample_id), method = "radix")
coverage_ok <- nrow(coverage) == 1L && coverage$expected_samples == 34L &&
  coverage$uniquely_located_samples == 34L && length(samples) == 34L &&
  all(vapply(samples, function(id) sum(ids == id), integer(1L)) == 1L) &&
  !truth(coverage$source_root_published) && !truth(coverage$source_paths_published)
queue_ok <- nrow(queue) == 204L && identical(queue$queue_order, 1:204) &&
  sum(queue$stage == "raw") == 34L && sum(queue$stage == "sct") == 170L &&
  !anyDuplicated(queue$sample_id[queue$stage == "raw"]) &&
  !anyDuplicated(paste(queue$sample_id[queue$stage == "sct"],
                       queue$seed[queue$stage == "sct"])) &&
  identical(sort(unique(queue$seed[queue$stage == "sct"])),
            20260805:20260809) &&
  all(queue$selected_cells[queue$stage == "sct"] == 384L) &&
  !any(truth(queue$biological_outcomes_computed)) &&
  !any(truth(queue$panel_fit) | truth(queue$pca) | truth(queue$ph) |
       truth(queue$landscape))
caps_ok <- nrow(caps) == 4L &&
  caps$elapsed_cap_seconds[caps$scope == "raw_child"] == 1800 &&
  caps$elapsed_cap_seconds[caps$scope == "sct_child"] == 1800 &&
  caps$elapsed_cap_seconds[caps$scope == "aggregate_worker"] == 14400 &&
  caps$storage_cap_bytes[caps$scope == "aggregate_storage"] == 4 * 1024^3
target_ok <- nrow(target) == 1L && target$seed == 20260809L &&
  target$sample_id == sort(samples, method = "radix")[[1L]]
criteria_ok <- nrow(criteria) == 10L && all(truth(criteria$passed))
decision_ok <- nrow(decision) == 1L &&
  decision$decision == "authorize_serial_atomic_upstream_production_only" &&
  decision$raw_jobs == 34L && decision$sct_jobs == 170L && decision$workers == 1L &&
  decision$retries == 0L && !truth(decision$panel_fit_authorized) &&
  !truth(decision$pca_authorized) && !truth(decision$ph_authorized) &&
  !truth(decision$landscape_authorized) && !truth(decision$labels_authorized) &&
  !truth(decision$outcomes_authorized)
checks <- data.frame(
  contract_id = "mv07f_prefreeze_independent_validation_v1",
  category = c("source_freeze", "private_source_coverage", "exact_queue",
    "resource_caps", "repeat_target", "acceptance_criteria", "decision_gate"),
  passed = c(source_ok, coverage_ok, queue_ok, caps_ok, target_ok, criteria_ok,
             decision_ok),
  detail = c("19 public-safe source identities", "34 unique private sources",
    "34 raw and 170 SCT jobs", "1800s 8GiB child and aggregate gates",
    "lexicographic first sample at final seed", "10 of 10 criteria",
    "serial upstream only"), stringsAsFactors = FALSE)
if (!all(checks$passed)) stop("MV7-F independent prefreeze validation failed: ",
  paste(checks$category[!checks$passed], collapse = ", "), call. = FALSE)
if (file.exists(output)) stop("Refusing overwrite: ", output)
write.table(checks, output, sep = ",", row.names = FALSE, col.names = TRUE,
            quote = TRUE, na = "")
message("MV7-F prefreeze independent validation: 7/7 categories pass")
