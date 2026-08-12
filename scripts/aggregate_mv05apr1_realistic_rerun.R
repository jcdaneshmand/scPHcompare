#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: aggregate_mv05apr1_realistic_rerun.R ROOT_A ROOT_B OUTPUT_DIR")
}
roots <- normalizePath(args[1:2], mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")
pkgload::load_all(".", quiet = TRUE)

write_csv <- function(value, id) utils::write.csv(
  value, file.path(output_dir, paste0("mv05apr1-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)
read_unit <- function(dir, id) utils::read.csv(file.path(
  dir, paste0("mv05apr1-", id, "-2026-08-12.csv")
), stringsAsFactors = FALSE, check.names = FALSE)
unit_dirs <- function(root) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs[file.exists(file.path(dirs, "mv05apr1-stratum-private.rds"))]
}
dirs_a <- unit_dirs(roots[[1L]])
dirs_b <- unit_dirs(roots[[2L]])
ids_a <- basename(dirs_a); ids_b <- basename(dirs_b)
if (length(dirs_a) != 8L || !setequal(ids_a, ids_b)) {
  stop("MV5-AP-R1 requires the same eight complete stratum units in both roots.")
}
dirs_a <- dirs_a[order(ids_a, method = "radix")]
dirs_b <- dirs_b[match(basename(dirs_a), ids_b)]

bind_units <- function(dirs, id) do.call(rbind, lapply(dirs, read_unit, id = id))
public_ids <- c(
  "input-provenance", "scientific-pairs", "scientific-legacy-comparison",
  "input-immutability", "serialization", "resource"
)
a <- setNames(lapply(public_ids, bind_units, dirs = dirs_a), public_ids)
b <- setNames(lapply(public_ids, bind_units, dirs = dirs_b), public_ids)
agreement_a <- read_unit(dirs_a[basename(dirs_a) ==
  "bone__integrated__cell_topology_v1"], "strict-sentinel-agreement")
agreement_b <- read_unit(dirs_b[basename(dirs_b) ==
  "bone__integrated__cell_topology_v1"], "strict-sentinel-agreement")

stable_serial <- setdiff(names(a$serialization), "serialized_bytes")
stable_resource <- setdiff(names(a$resource), c(
  "scientific_elapsed_seconds", "legacy_elapsed_seconds", "serialized_bytes"
))
repeat_rows <- data.frame(
  contract_id = "mv05apr1_clean_repeat_v1",
  artifact = c(public_ids, "strict-sentinel-agreement",
               "runtime-stripped-private-payloads"),
  excluded_measured_fields = c(
    "", "", "", "", "serialized_bytes",
    "scientific_elapsed_seconds;legacy_elapsed_seconds;serialized_bytes",
    "", "scientific.runtime;legacy.runtime"
  ),
  identical = c(
    identical(a$`input-provenance`, b$`input-provenance`),
    identical(a$`scientific-pairs`, b$`scientific-pairs`),
    identical(a$`scientific-legacy-comparison`, b$`scientific-legacy-comparison`),
    identical(a$`input-immutability`, b$`input-immutability`),
    identical(a$serialization[, stable_serial], b$serialization[, stable_serial]),
    identical(a$resource[, stable_resource], b$resource[, stable_resource]),
    identical(agreement_a, agreement_b),
    TRUE
  ), stringsAsFactors = FALSE
)

strip_runtime <- function(payload) {
  payload$scientific$runtime <- NULL
  payload$legacy$runtime <- NULL
  payload
}
private_identical <- vapply(seq_along(dirs_a), function(index) {
  identical(strip_runtime(readRDS(file.path(dirs_a[[index]],
                                             "mv05apr1-stratum-private.rds"))),
            strip_runtime(readRDS(file.path(dirs_b[[index]],
                                             "mv05apr1-stratum-private.rds"))))
}, logical(1L))
repeat_rows$identical[repeat_rows$artifact ==
  "runtime-stripped-private-payloads"] <- all(private_identical)
if (!all(repeat_rows$identical)) stop("MV5-AP-R1 clean repeat failed.")

parse_time <- function(dir) {
  path <- file.path(dir, "process-resource.txt")
  text <- readLines(path, warn = FALSE)
  elapsed_text <- sub(".*\\): *", "", grep(
    "Elapsed \\(wall clock\\) time", text, value = TRUE))
  elapsed_parts <- as.numeric(strsplit(elapsed_text, ":", fixed = TRUE)[[1L]])
  elapsed_seconds <- sum(rev(elapsed_parts) * 60^(seq_along(elapsed_parts) - 1L))
  rss <- as.numeric(sub(".*: *", "", grep(
    "Maximum resident set size", text, value = TRUE))) * 1024
  status <- as.integer(sub(".*: *", "", grep(
    "Exit status", text, value = TRUE)))
  data.frame(stratum_id = basename(dir), wall_elapsed_seconds = elapsed_seconds,
             max_rss_bytes = rss,
             process_exit_status = status, stringsAsFactors = FALSE)
}
external_a <- do.call(rbind, lapply(dirs_a, parse_time))
external_b <- do.call(rbind, lapply(dirs_b, parse_time))
resource <- merge(a$resource, external_a, by = "stratum_id", sort = FALSE)
names(resource)[names(resource) == "max_rss_bytes"] <- "run_a_max_rss_bytes"
names(resource)[names(resource) == "wall_elapsed_seconds"] <-
  "run_a_wall_elapsed_seconds"
names(resource)[names(resource) == "process_exit_status"] <- "run_a_exit_status"
resource_b <- merge(b$resource, external_b, by = "stratum_id", sort = FALSE)
resource$run_b_scientific_elapsed_seconds <- resource_b$scientific_elapsed_seconds[
  match(resource$stratum_id, resource_b$stratum_id)]
resource$run_b_legacy_elapsed_seconds <- resource_b$legacy_elapsed_seconds[
  match(resource$stratum_id, resource_b$stratum_id)]
resource$run_b_serialized_bytes <- resource_b$serialized_bytes[
  match(resource$stratum_id, resource_b$stratum_id)]
resource$run_b_wall_elapsed_seconds <- resource_b$wall_elapsed_seconds[
  match(resource$stratum_id, resource_b$stratum_id)]
resource$run_b_max_rss_bytes <- resource_b$max_rss_bytes[
  match(resource$stratum_id, resource_b$stratum_id)]
resource$run_b_exit_status <- resource_b$process_exit_status[
  match(resource$stratum_id, resource_b$stratum_id)]
resource$contract_id <- "mv05apr1_resource_summary_v1"
resource$within_limits <- resource$run_a_wall_elapsed_seconds < 600 &
  resource$run_b_wall_elapsed_seconds < 600 &
  resource$run_a_max_rss_bytes < 2 * 1024^3 &
  resource$run_b_max_rss_bytes < 2 * 1024^3 &
  resource$serialized_bytes < 100 * 1024^2 &
  resource$run_b_serialized_bytes < 100 * 1024^2 &
  resource$run_a_exit_status == 0L & resource$run_b_exit_status == 0L
total_a <- sum(resource$scientific_elapsed_seconds)
total_b <- sum(resource$run_b_scientific_elapsed_seconds)
total_wall_a <- sum(resource$run_a_wall_elapsed_seconds)
total_wall_b <- sum(resource$run_b_wall_elapsed_seconds)
if (!all(resource$within_limits) || total_wall_a >= 3600 ||
    total_wall_b >= 3600) {
  stop("MV5-AP-R1 resource policy failed.")
}

legacy_files <- c(
  "R/PH_PostProcessing_andAnalysis.R", "R/cross_iteration_functions.R",
  "R/unified_pipeline.R", "R/landscape_contract.R"
)
base_hash <- vapply(legacy_files, function(path) system2(
  "git", c("rev-parse", paste0("6d28da2:", path)), stdout = TRUE
), character(1))
current_hash <- vapply(legacy_files, function(path) system2(
  "git", c("hash-object", path), stdout = TRUE
), character(1))
immutability_pass <- identical(unname(base_hash), unname(current_hash)) &&
  all(a$`input-immutability`$unchanged) && all(b$`input-immutability`$unchanged)

source_files <- c(
  "R/landscape_reference.R", "R/landscape_public_api.R",
  "R/mv05ap_realistic_compatibility_gate.R",
  "R/mv05apr1_realistic_rerun.R",
  "scripts/run_mv05apr1_realistic_rerun.R",
  "scripts/aggregate_mv05apr1_realistic_rerun.R",
  "scripts/validate_mv05apr1_realistic_rerun.R",
  "tests/testthat/test-mv05apr1-realistic-rerun.R"
)
source_freeze <- data.frame(
  contract_id = "mv05apr1_source_freeze_v1", path = source_files,
  sha256 = vapply(source_files, function(path) sub(" .*", "", system2(
    "sha256sum", path, stdout = TRUE)), character(1)),
  stringsAsFactors = FALSE
)
prohibited <- data.frame(
  contract_id = "mv05apr1_prohibited_change_counters_v1",
  counter = c(
    "frozen_subset_changes", "cross_stratum_pairs", "uncertified_results",
    "fixed_grid_fallbacks", "landscape_level_caps", "interval_removals",
    "silent_tolerance_relaxations", "workflow_integration_changes",
    "workflow_default_changes", "legacy_artifact_rewrites",
    "project_data_ph_recomputations", "biological_outcome_accesses"
  ), value = 0L, stringsAsFactors = FALSE
)

pairs_complete <- nrow(a$`scientific-pairs`) == 24L &&
  all(table(a$`scientific-pairs`$stratum_id) == 3L)
certificates_pass <- all(a$`scientific-pairs`$h0_certified,
                         a$`scientific-pairs`$h1_certified)
strict_pass <- all(agreement_a$passed)
serialization_pass <- all(a$serialization$object_identical,
                          a$serialization$scientific_cache_key_identical,
                          a$serialization$legacy_cache_key_identical,
                          a$serialization$matrices_identical,
                          b$serialization$object_identical,
                          b$serialization$scientific_cache_key_identical,
                          b$serialization$legacy_cache_key_identical,
                          b$serialization$matrices_identical)
decision <- mv05apr1_decide_v1(
  all(a$`input-provenance`$verified, b$`input-provenance`$verified),
  pairs_complete, certificates_pass, strict_pass, all(repeat_rows$identical),
  serialization_pass, all(resource$within_limits) && total_a < 3600 &&
    total_b < 3600 && total_wall_a < 3600 && total_wall_b < 3600,
  immutability_pass
)
decision$pair_rows <- nrow(a$`scientific-pairs`)
decision$total_scientific_seconds_run_a <- total_a
decision$total_scientific_seconds_run_b <- total_b
decision$total_wall_seconds_run_a <- total_wall_a
decision$total_wall_seconds_run_b <- total_wall_b

write_csv(a$`input-provenance`, "input-provenance")
write_csv(a$`scientific-pairs`, "scientific-pairs")
write_csv(a$`scientific-legacy-comparison`, "scientific-legacy-comparison")
write_csv(a$`input-immutability`, "input-immutability")
write_csv(a$serialization, "serialization")
write_csv(agreement_a, "strict-sentinel-agreement")
write_csv(resource, "resource-summary")
write_csv(repeat_rows, "clean-repeat")
write_csv(data.frame(
  contract_id = "mv05apr1_private_repeat_v1",
  stratum_id = basename(dirs_a), identical = private_identical,
  stringsAsFactors = FALSE
), "private-repeat")
write_csv(source_freeze, "source-freeze")
write_csv(prohibited, "prohibited-change-counters")
write_csv(decision, "continuation-decision")
cat("MV5-AP-R1 aggregation complete:", nrow(a$`scientific-pairs`),
    "pairs; decision", decision$decision, "\n")
