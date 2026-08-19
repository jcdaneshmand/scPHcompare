#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: validate_mv08i_hca_validation.R <audit-dir> <output.csv>", call. = FALSE)
}
audit_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_path <- args[[2L]]
read_csv <- function(name) utils::read.csv(file.path(audit_dir, name), check.names = FALSE,
                                           stringsAsFactors = FALSE)
hash_file <- function(path) toupper(digest::digest(file = path, algo = "sha256", serialize = FALSE))
topology <- read_csv("mv08i-topology-summary.csv")
validation <- read_csv("mv08i-topology-validation.csv")
identity <- read_csv("mv08i-input-identity.csv")
execution <- read_csv("mv08i-landscape-execution.csv")
levels <- read_csv("mv08i-landscape-level-summary.csv")
ledger <- read_csv("mv08i-repeat-ledger.csv")
manifest <- read_csv("mv08i-public-artifact-manifest.csv")
expected_panel <- "b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba"
units <- sprintf("HCA_BM_%03d", 1:8)
checks <- list()
add <- function(id, passed, evidence) checks[[length(checks) + 1L]] <<- data.frame(
  check_id = id, passed = isTRUE(passed), evidence = evidence, stringsAsFactors = FALSE)
add("topology_unit_count", identical(sort(unique(topology$unit_id)), units) && nrow(topology) == 32L,
    "8 units x 2 views x H0/H1")
add("topology_point_axes", all(topology$point_count[topology$view_id == "cell_topology_v1"] == 384L) &&
      all(topology$point_count[topology$view_id == "gene_topology_v1"] == 475L),
    "384-cell and 475-gene views")
add("topology_h0_connectivity", all(topology$positive_intervals[topology$homology_dimension == "H0"] ==
      topology$point_count[topology$homology_dimension == "H0"] - 1L), "finite positive H0 = points - 1")
add("topology_panel_identity", all(tolower(topology$panel_sha256) == expected_panel) &&
      all(tolower(identity$panel_sha256) == expected_panel), "common475 SHA-256")
add("topology_selection", nrow(identity) == 8L && all(identity$selected_cells == 384L) &&
      all(identity$outcome_label_state == "closed"), "8/8 deterministic 384-cell selections")
add("topology_firewall", all(!topology$labels_outcomes_opened) && all(!topology$landscapes_computed) &&
      all(!topology$fusion_computed) && all(!identity$labels_outcomes_opened) &&
      all(!identity$landscapes_computed) && all(!identity$fusion_computed), "labels/outcomes, landscapes, and fusion closed during PH")
add("landscape_profiles", nrow(execution) == 32L && all(execution$active_level_count > 0L) &&
      all(execution$integration_method == "exact_critical_breakpoint_v1") &&
      all(execution$level_policy == "all_active_consecutive_levels") && all(execution$grid_policy == "none"),
    "32 exact all-active-level profiles; no grid or universal cap")
add("landscape_levels", nrow(levels) == sum(execution$active_level_count) &&
      all(vapply(split(levels, interaction(levels$unit_id, levels$view_id, levels$homology_dimension, drop = TRUE)),
        function(x) identical(as.integer(x$level), seq_len(unique(x$active_level_count))), logical(1))),
    "every retained landscape level is consecutive and accounted for")
add("landscape_firewall", all(!execution$labels_outcomes_opened) && all(!execution$fusion_computed) &&
      all(!levels$labels_outcomes_opened) && all(!levels$fusion_computed), "landscape downstream firewall closed")
add("repeat_ledger", nrow(ledger) >= 5L && all(ledger$byte_identical) &&
      all(ledger$primary_sha256 == ledger$repeat_sha256), "primary/repeat scientific artifacts agree")
manifest_pass <- all(vapply(seq_len(nrow(manifest)), function(i) {
  path <- file.path(audit_dir, manifest$artifact[[i]])
  file.exists(path) && identical(hash_file(path), toupper(manifest$sha256[[i]]))
}, logical(1)))
add("public_manifest", manifest_pass, "all public audit artifacts rehash to the bound manifest")
result <- do.call(rbind, checks)
partial <- paste0(output_path, ".partial")
utils::write.csv(result, partial, row.names = FALSE, quote = TRUE)
file.rename(partial, output_path)
if (!all(result$passed)) stop("MV8-I validation failed", call. = FALSE)
cat("MV8-I independent validation passed ", sum(result$passed), "/", nrow(result), " checks\n", sep = "")
