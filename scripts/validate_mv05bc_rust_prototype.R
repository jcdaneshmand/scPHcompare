args <- commandArgs(trailingOnly = TRUE)
private_run <- if (length(args)) args[[1L]] else "tmp/mv05bc/run-c"

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required artifact: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

equivalence <- read_required(
  "docs/audits/mv05bc-equivalence-summary-2026-08-13.csv"
)
performance <- read_required(
  "docs/audits/mv05bc-performance-summary-2026-08-13.csv"
)
validation <- read_required(
  "docs/audits/mv05bc-independent-validation-2026-08-13.csv"
)
decision <- read_required(
  "docs/audits/mv05bc-continuation-decision-2026-08-13.csv"
)
fixtures <- read_required(file.path(private_run, "fixture-validation.csv"))
tier_b <- read_required(file.path(private_run, "tier-b-equivalence.csv"))
tier_c <- read_required(file.path(private_run, "tier-c-equivalence.csv"))
speed <- read_required(file.path(private_run, "speed.csv"))

checks <- c(
  identical(equivalence$tier, c("A", "B", "C")),
  identical(equivalence$required_results, c(3L, 20L, 12L)),
  identical(equivalence$passed_results, equivalence$required_results),
  all(equivalence$gate_passed),
  nrow(fixtures) == 3L && all(as.logical(fixtures$passed)),
  nrow(tier_b) == 20L && all(as.logical(tier_b$equivalent)),
  nrow(tier_c) == 12L && all(as.logical(tier_c$equivalent)),
  max(as.double(tier_b$absolute_error)) <=
    max(as.double(tier_b$acceptance_threshold)),
  max(as.double(tier_c$absolute_error)) <=
    max(as.double(tier_c$acceptance_threshold)),
  nrow(speed) == 6L && all(as.logical(speed$rust_no_slower)),
  stats::median(as.double(speed$speedup_vs_r)) >= 3,
  isTRUE(performance$speed_gate_passed) &&
    isTRUE(performance$memory_gate_passed),
  performance$peak_rss_bytes <= performance$maximum_allowed_rss_bytes,
  nrow(validation) == 18L && all(as.logical(validation$result)),
  isTRUE(decision$prototype_accepted),
  isTRUE(decision$r_engine_canonical),
  !isTRUE(decision$rust_production_adoption_authorized),
  !isTRUE(decision$tier_d_e_complete),
  !isTRUE(decision$public_default_changed),
  !isTRUE(decision$additional_seed_production_authorized),
  !isTRUE(decision$partitions_authorized),
  identical(decision$next_sprint,
            "MV5-BD_full_Rust_equivalence_and_adoption_gate")
)

if (!all(checks)) {
  stop("MV5-BC independent validation failed at checks: ",
       paste(which(!checks), collapse = ", "), call. = FALSE)
}
cat("MV5-BC independent validation: ", sum(checks), "/",
    length(checks), " passed\n", sep = "")
