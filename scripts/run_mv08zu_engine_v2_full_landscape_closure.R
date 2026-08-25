#!/usr/bin/env Rscript

# Launch only the independent MV8-ZU rehashing closure after MV8-ZT has
# published its complete terminal receipt. This launcher cannot run landscape
# workers and binds every argument to its prospectively reviewed path.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0L) {
  stop("run_mv08zu_engine_v2_full_landscape_closure.R takes no arguments",
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest required", call. = FALSE)
}

sha_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

builder <- "scripts/build_mv08zu_engine_v2_full_landscape_closure.R"
expected_builder_sha256 <- "933bc0660cb8992f56db071dd729a0bde9db07ec221e35ed9f48ede6d23fb3b5"
closure_args <- c(
  "docs/audits/mv08zt-engine-v2-full-landscape-prefreeze-v3",
  "docs/audits/mv08z-landscape-execution-prefreeze-v1",
  "tmp/mv08z-landscape-execution-private-v1/mv08z-private-unit-bindings.csv",
  "tmp/mv08s-ph-sentinel-v3",
  "tmp/mv08va-full-ph-production-v2",
  paste0(
    "tmp/mv08zp-landscape-kernel-repair-private-v1/",
    "scph-landscape-kernel-v2-final-candidate_x86_64-unknown-linux-gnu.so"
  ),
  "tmp/mv08zt-full-landscape-private-v2",
  "tmp/mv08zt-full-landscape-public-v2",
  "docs/audits/mv08zt-engine-v2-interruption-recovery-prefreeze-v2"
)
output_dir <- "docs/audits/mv08zu-engine-v2-full-landscape-closure-v1"
progress_path <- file.path(
  closure_args[[8L]], "mv08zt-progress.csv"
)
ledger_path <- file.path(
  closure_args[[8L]], "mv08zt-resource-ledger.csv"
)
completion_path <- file.path(
  closure_args[[8L]], "mv08zt-chunk-completions.csv"
)

if (!file.exists(builder) || sha_file(builder) != expected_builder_sha256) {
  stop("MV8-ZU closure builder hash drift", call. = FALSE)
}
if (!all(file.exists(closure_args)) || !all(file.exists(c(
  progress_path, ledger_path, completion_path
)))) {
  stop("MV8-ZU closure input missing", call. = FALSE)
}
if (dir.exists(output_dir)) {
  stop("refusing to overwrite MV8-ZU output", call. = FALSE)
}

progress <- read_csv(progress_path)
ledger <- read_csv(ledger_path)
completed <- read_csv(completion_path)
terminal_ok <- nrow(progress) == 1L &&
  progress$state == "landscape_production_complete_closure_pending" &&
  progress$completed_chunks == 628L && progress$total_chunks == 628L &&
  progress$completed_pairs == 152744L && progress$total_pairs == 152744L &&
  progress$workers == 1L && progress$retries == 0L &&
  nrow(ledger) == 628L && nrow(completed) == 628L &&
  identical(as.integer(ledger$production_order), 1:628) &&
  identical(as.integer(completed$production_order), 1:628)
if (!terminal_ok) {
  stop("MV8-ZT is not terminal; MV8-ZU closure remains closed", call. = FALSE)
}

rscript <- file.path(R.home("bin"), "Rscript")
status <- system2(
  rscript,
  args = vapply(c(builder, closure_args, output_dir), shQuote, character(1L))
)
if (!identical(as.integer(status), 0L)) {
  stop("MV8-ZU independent closure failed", call. = FALSE)
}
cat("MV8-ZU terminal launcher completed the independent closure\n")
