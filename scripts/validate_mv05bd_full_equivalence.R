args <- commandArgs(trailingOnly = TRUE)
private_run <- if (length(args)) args[[1L]] else "tmp/mv05bd/run-b"

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required artifact: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

equivalence <- read_required(file.path(private_run, "private-equivalence.csv"))
self <- read_required(file.path(private_run, "private-self.csv"))
environment <- read_required(file.path(private_run, "environment.csv"))
summary <- read_required(
  "docs/audits/mv05bd-full-equivalence-summary-2026-08-13.csv"
)
resources <- read_required(
  "docs/audits/mv05bd-invariant-resource-summary-2026-08-13.csv"
)
decision <- read_required(
  "docs/audits/mv05bd-continuation-decision-2026-08-13.csv"
)

checks <- c(
  nrow(equivalence) == 408L,
  sum(equivalence$tier == "D") == 318L,
  sum(equivalence$tier == "E") == 90L,
  all(as.logical(equivalence$equivalent)),
  all(equivalence$status == 0L),
  all(equivalence$engine_version == 1L),
  all(as.logical(equivalence$finite_counts_match)),
  all(as.logical(equivalence$reverse_bit_identical)),
  all(as.logical(equivalence$reverse_counts_swap)),
  all(as.logical(equivalence$reverse_diagnostics_match)),
  nrow(self) == 112L,
  all(as.logical(self$exactly_zero)),
  all(as.logical(self$counts_match)),
  all(self$status == 0L),
  identical(summary$tier, c("D", "E")),
  identical(summary$passed_results, summary$required_results),
  all(as.logical(summary$gate_passed)),
  isTRUE(resources$normalized_runs_identical),
  resources$reverse_bit_identical == 408L,
  resources$self_exact_zero == 112L,
  resources$peak_rss_bytes <= resources$rss_limit_bytes,
  isTRUE(resources$resource_gate_passed),
  environment$labels_opened == FALSE,
  environment$outcomes_computed == FALSE,
  environment$defaults_changed == FALSE,
  isTRUE(decision$full_numerical_equivalence_accepted),
  isTRUE(decision$r_engine_canonical),
  !isTRUE(decision$rust_production_adoption_authorized),
  isTRUE(decision$linux_build_runtime_certified),
  !isTRUE(decision$windows_build_runtime_certified),
  !isTRUE(decision$macos_build_runtime_certified),
  !isTRUE(decision$artifact_distribution_certified),
  !isTRUE(decision$public_default_changed),
  identical(decision$next_sprint,
            "MV5-BE_cross_platform_build_and_distribution_prefreeze")
)

if (!all(checks)) {
  stop("MV5-BD validation failed at checks: ",
       paste(which(!checks), collapse = ", "), call. = FALSE)
}
cat("MV5-BD independent validation: ", sum(checks), "/",
    length(checks), " passed\n", sep = "")
