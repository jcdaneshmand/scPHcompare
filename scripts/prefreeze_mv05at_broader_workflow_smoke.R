#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: prefreeze_mv05at_broader_workflow_smoke.R OUTPUT_DIR")
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve the repository root.")
setwd(normalizePath(file.path(
  dirname(gsub("~+~", " ", sub("^--file=", "", script_arg[[1L]]), fixed = TRUE)), ".."
), mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)
write_csv <- function(x, id) utils::write.csv(
  x, file.path(output_dir, paste0("mv05at-", id, "-2026-08-12.csv")),
  row.names = FALSE, na = "", quote = TRUE
)

subset <- utils::read.csv("docs/audits/mv05ap-frozen-subset-2026-08-12.csv",
                          stringsAsFactors = FALSE, check.names = FALSE)
manifest <- utils::read.csv("docs/audits/mv04-input-diagram-manifest-2026-08-05.csv",
                            stringsAsFactors = FALSE, check.names = FALSE)
binding_columns <- c("diagram_id", "stratum_id", "sample_id", "view_id",
  "h0_finite_intervals", "h1_finite_intervals", "diagram_sha256",
  "result_file_sha256", "result_file", "selection_role", "selection_rule")
frozen <- subset[order(subset$diagram_id, method = "radix"), binding_columns]
selected <- mv05ap_select_depth_triplets_v1(manifest)
selected <- selected[order(selected$diagram_id, method = "radix"), binding_columns]
rownames(frozen) <- rownames(selected) <- NULL
if (!identical(frozen, selected)) stop("MV5-AT frozen input binding drifted.")

strata <- sort(unique(subset$stratum_id), method = "radix")
scope <- do.call(rbind, lapply(strata, function(id) {
  rows <- subset[subset$stratum_id == id, ]
  data.frame(contract_id = "mv05at_scope_v1", stratum_id = id,
    diagrams = nrow(rows), pairs = choose(nrow(rows), 2L),
    h0_min = min(rows$h0_finite_intervals), h0_max = max(rows$h0_finite_intervals),
    h1_min = min(rows$h1_finite_intervals), h1_max = max(rows$h1_finite_intervals),
    expected_h0_method = "exact_breakpoint_stream_v1",
    expected_h1_method = if (max(rows$h1_finite_intervals) > 500L)
      "adaptive_global_error_v1" else "exact_breakpoint_stream_v1",
    stringsAsFactors = FALSE)
}))
stopifnot(nrow(scope) == 8L, sum(scope$diagrams) == 24L,
          sum(scope$pairs) == 24L)
write_csv(scope, "scope")

prior <- utils::read.csv("docs/audits/mv05apr1-resource-summary-2026-08-12.csv",
                         stringsAsFactors = FALSE, check.names = FALSE)
resources <- merge(scope[, c("stratum_id", "expected_h1_method")],
  prior[, c("stratum_id", "run_a_wall_elapsed_seconds", "run_a_max_rss_bytes",
            "run_b_wall_elapsed_seconds", "run_b_max_rss_bytes")],
  by = "stratum_id", sort = TRUE)
resources$contract_id <- "mv05at_resource_admission_v1"
resources$planned_wall_seconds <- ifelse(
  resources$expected_h1_method == "adaptive_global_error_v1", 750, 120)
resources$max_wall_seconds <- 750
resources$max_pairs <- 3L
resources$max_rss_bytes <- 2 * 1024^3
resources$workers <- 1L
resources$admitted <- pmax(resources$run_a_wall_elapsed_seconds,
                           resources$run_b_wall_elapsed_seconds) < 750 &
  pmax(resources$run_a_max_rss_bytes, resources$run_b_max_rss_bytes) < 2 * 1024^3
resources <- resources[, c("contract_id", "stratum_id", "expected_h1_method",
  "planned_wall_seconds", "max_wall_seconds", "max_pairs", "max_rss_bytes",
  "workers", "run_a_wall_elapsed_seconds", "run_b_wall_elapsed_seconds",
  "run_a_max_rss_bytes", "run_b_max_rss_bytes", "admitted")]
if (!all(resources$admitted)) stop("MV5-AT prior resource evidence refuses a unit.")
write_csv(resources, "resource-admission")

abort_rules <- data.frame(contract_id = "mv05at_abort_rules_v1",
  rule_id = sprintf("AT-ABORT-%02d", 1:12),
  condition = c("frozen binding drift", "input hash or eligibility mismatch",
    "interval envelope drift", "pair-count drift", "resource-plan refusal",
    "process wall/RSS/exit breach", "uncertified dimension", "artifact collision",
    "completion verification failure", "resume mutation", "legacy/default/downstream drift",
    "test, check, or scope failure"), action = "stop_without_authorization",
  stringsAsFactors = FALSE)
write_csv(abort_rules, "abort-rules")

sources <- c("R/corrected_landscape_workflow.R", "R/PH_PostProcessing_andAnalysis.R",
  "R/unified_pipeline.R", "R/landscape_public_api.R", "R/landscape_reference.R",
  "docs/audits/mv05ap-frozen-subset-2026-08-12.csv",
  "docs/audits/mv05apr1-resource-summary-2026-08-12.csv",
  "docs/specifications/MV05AT_BROADER_REALISTIC_WORKFLOW_SMOKE_SPECIFICATION_V1.md",
  "scripts/prefreeze_mv05at_broader_workflow_smoke.R",
  "scripts/run_mv05at_broader_workflow_smoke.R")
freeze <- data.frame(contract_id = "mv05at_source_freeze_v1", path = sources,
  sha256 = vapply(sources, function(path) sub(" .*", "", system2(
    "sha256sum", path, stdout = TRUE)), character(1)), stringsAsFactors = FALSE)
write_csv(freeze, "source-freeze")

decision <- data.frame(contract_id = "mv05at_prefreeze_decision_v1",
  decision = "execute_eight_bounded_existing_data_units",
  strata = 8L, diagrams = 24L, pairs = 24L, execute_authorized = TRUE,
  downstream_consumption_authorized = FALSE, default_change_authorized = FALSE,
  legacy_rewrite_authorized = FALSE, new_data_authorized = FALSE,
  optimization_authorized = FALSE, stringsAsFactors = FALSE)
write_csv(decision, "prefreeze-decision")
cat("MV5-AT prefreeze complete:", nrow(scope), "admitted strata\n")
