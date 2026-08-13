args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "Usage: run_mv05bc_rust_prototype.R <intervals.csv> <references.csv> ",
    "<speed-references.csv> <r-speed.csv> <rust-library> <output-dir>",
    call. = FALSE
  )
}

interval_path <- args[[1L]]
reference_path <- args[[2L]]
speed_reference_path <- args[[3L]]
r_speed_path <- args[[4L]]
rust_library <- args[[5L]]
output_dir <- args[[6L]]

source("R/landscape_rust_prototype.R")

atomic_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  utils::write.csv(value, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not atomically publish ", path, ".", call. = FALSE)
  }
}

stable_identity <- function(tier, row, candidate) {
  payload <- paste(
    "mv05bc_rust_result_v1", tier, row$stratum_id, row$pair_order,
    row$pair_cache_key, row$dimension,
    sprintf("%.17g", candidate$squared_distance), candidate$active_levels,
    candidate$event_segments, sep = "|"
  )
  paste0(
    "mv05bc_rust_result_v1:",
    digest::digest(payload, algo = "sha256", serialize = FALSE)
  )
}

intervals <- utils::read.csv(interval_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
references <- utils::read.csv(reference_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
speed_references <- utils::read.csv(
  speed_reference_path, stringsAsFactors = FALSE, check.names = FALSE
)
r_speed <- utils::read.csv(r_speed_path, stringsAsFactors = FALSE,
                           check.names = FALSE)

intervals <- intervals[order(intervals$diagram_id, intervals$dimension,
                             intervals$interval_order), , drop = FALSE]
interval_keys <- paste(intervals$diagram_id, intervals$dimension, sep = "|")
interval_map <- split(intervals[, c("birth", "death"), drop = FALSE],
                      interval_keys, drop = TRUE)

interval_for <- function(diagram_id, dimension) {
  key <- paste(diagram_id, as.integer(sub("^H", "", dimension)), sep = "|")
  value <- interval_map[[key]]
  if (is.null(value)) {
    return(matrix(numeric(), nrow = 0L, ncol = 2L,
                  dimnames = list(NULL, c("birth", "death"))))
  }
  as.matrix(value)
}

run_row <- function(row) {
  dimension <- as.integer(sub("^H", "", row$dimension))
  landscape_rust_prototype_dimension(
    interval_for(row$first_diagram_id, row$dimension),
    interval_for(row$second_diagram_id, row$dimension),
    dimension = dimension, library = rust_library
  )
}

threshold_for <- function(row) {
  reference <- as.double(row$reference_squared_distance)
  if (isTRUE(as.logical(row$reference_exact))) {
    1e-10 + 1e-10 * abs(reference)
  } else {
    as.double(row$achieved_absolute_error_estimate) +
      100 * .Machine$double.eps * max(1, abs(reference))
  }
}

equivalence_row <- function(row, tier) {
  candidate <- run_row(row)
  reference <- as.double(row$reference_squared_distance)
  threshold <- threshold_for(row)
  error <- abs(candidate$squared_distance - reference)
  data.frame(
    contract_id = "mv05bc_private_equivalence_v1", tier = tier,
    stratum_id = row$stratum_id, pair_order = row$pair_order,
    pair_cache_key = row$pair_cache_key, dimension = row$dimension,
    reference_exact = as.logical(row$reference_exact),
    reference_squared_distance = reference,
    candidate_squared_distance = candidate$squared_distance,
    absolute_error = error, acceptance_threshold = threshold,
    equivalent = isTRUE(candidate$rust_used) && error <= threshold,
    active_levels = candidate$active_levels,
    event_segments = candidate$event_segments,
    first_finite_intervals = candidate$first_finite_intervals,
    second_finite_intervals = candidate$second_finite_intervals,
    engine_version = candidate$engine_version, status = candidate$status,
    result_identity = stable_identity(tier, row, candidate),
    stringsAsFactors = FALSE
  )
}

empty <- matrix(numeric(), nrow = 0L, ncol = 2L)
fixture_definitions <- list(
  single_tent = list(matrix(c(0, 2), nrow = 1L), empty, 2 / 3),
  sign_changing_tents = list(
    matrix(c(0, 2), nrow = 1L), matrix(c(0.25, 2.25), nrow = 1L),
    7 / 64
  ),
  narrow_feature = list(
    matrix(c(0.499, 0.501), nrow = 1L), empty, 0.002 ^ 3 / 12
  )
)
fixture_rows <- lapply(names(fixture_definitions), function(case) {
  definition <- fixture_definitions[[case]]
  candidate <- landscape_rust_prototype_dimension(
    definition[[1L]], definition[[2L]], 0L, library = rust_library
  )
  error <- abs(candidate$squared_distance - definition[[3L]])
  data.frame(
    contract_id = "mv05bc_fixture_v1", case = case,
    expected_squared_distance = definition[[3L]],
    observed_squared_distance = candidate$squared_distance,
    absolute_error = error, acceptance_threshold = 1e-12,
    passed = isTRUE(candidate$rust_used) && error <= 1e-12,
    active_levels = candidate$active_levels,
    event_segments = candidate$event_segments,
    engine_version = candidate$engine_version, status = candidate$status,
    stringsAsFactors = FALSE
  )
})
fixture_rows <- do.call(rbind, fixture_rows)

tractable <- references[
  as.logical(references$reference_exact) &
    pmax(references$first_finite_intervals,
         references$second_finite_intervals) <= 500L,
  , drop = FALSE
]
tractable <- tractable[order(
  tractable$stratum_id, tractable$pair_order, tractable$dimension,
  tractable$pair_cache_key
), , drop = FALSE]
if (nrow(tractable) < 20L) {
  stop("Fewer than 20 deterministic tractable exact references exist.",
       call. = FALSE)
}
tier_b_selection <- tractable[seq_len(20L), , drop = FALSE]
tier_b <- do.call(rbind, lapply(seq_len(nrow(tier_b_selection)), function(i) {
  equivalence_row(tier_b_selection[i, , drop = FALSE], "B")
}))

speed_references <- speed_references[order(
  speed_references$stratum_id, speed_references$pair_order,
  speed_references$dimension
), , drop = FALSE]
if (nrow(speed_references) != 12L) {
  stop("Frozen Tier C must contain exactly 12 H0/H1 references.",
       call. = FALSE)
}
tier_c <- do.call(rbind, lapply(seq_len(nrow(speed_references)), function(i) {
  equivalence_row(speed_references[i, , drop = FALSE], "C")
}))

pair_keys <- unique(speed_references[, c("stratum_id", "pair_order"),
                                     drop = FALSE])
speed_rows <- lapply(seq_len(nrow(pair_keys)), function(i) {
  key <- pair_keys[i, , drop = FALSE]
  pair_rows <- speed_references[
    speed_references$stratum_id == key$stratum_id &
      speed_references$pair_order == key$pair_order,
    , drop = FALSE
  ]
  started <- proc.time()[["elapsed"]]
  values <- lapply(seq_len(nrow(pair_rows)), function(j) {
    run_row(pair_rows[j, , drop = FALSE])
  })
  elapsed <- proc.time()[["elapsed"]] - started
  accepted <- r_speed[
    r_speed$stratum_id == key$stratum_id &
      r_speed$pair_order == key$pair_order,
    , drop = FALSE
  ]
  if (nrow(accepted) != 1L) stop("Missing accepted R speed row.", call. = FALSE)
  h0 <- values[[match("H0", pair_rows$dimension)]]$squared_distance
  h1 <- values[[match("H1", pair_rows$dimension)]]$squared_distance
  data.frame(
    contract_id = "mv05bc_private_speed_v1",
    panel_order = accepted$panel_order, stratum_id = key$stratum_id,
    pair_order = key$pair_order, rust_h0_squared_distance = h0,
    rust_h1_squared_distance = h1, rust_elapsed_seconds = elapsed,
    accepted_r_elapsed_seconds = accepted$elapsed_seconds,
    speedup_vs_r = accepted$elapsed_seconds / elapsed,
    rust_no_slower = elapsed <= accepted$elapsed_seconds,
    stringsAsFactors = FALSE
  )
})
speed_rows <- do.call(rbind, speed_rows)
speed_rows <- speed_rows[order(speed_rows$panel_order), , drop = FALSE]

environment <- data.frame(
  contract_id = "mv05bc_environment_v1",
  engine_id = "rust_exact_critical_pairs_segment_tree_v1",
  rust_toolchain = "1.97.1-x86_64-unknown-linux-gnu",
  engine_version = 1L,
  external_rust_crates = 0L,
  tier_a_results = nrow(fixture_rows), tier_b_results = nrow(tier_b),
  tier_c_results = nrow(tier_c), labels_opened = FALSE,
  outcomes_computed = FALSE, production_adoption_authorized = FALSE,
  stringsAsFactors = FALSE
)

atomic_csv(fixture_rows, file.path(output_dir, "fixture-validation.csv"))
atomic_csv(tier_b, file.path(output_dir, "tier-b-equivalence.csv"))
atomic_csv(tier_c, file.path(output_dir, "tier-c-equivalence.csv"))
atomic_csv(speed_rows, file.path(output_dir, "speed.csv"))
atomic_csv(environment, file.path(output_dir, "environment.csv"))

if (!all(fixture_rows$passed) || !all(tier_b$equivalent) ||
    !all(tier_c$equivalent)) {
  stop("MV5-BC equivalence gate failed.", call. = FALSE)
}

cat(
  "MV5-BC prototype: Tier A ", sum(fixture_rows$passed), "/",
  nrow(fixture_rows), "; Tier B ", sum(tier_b$equivalent), "/",
  nrow(tier_b), "; Tier C ", sum(tier_c$equivalent), "/",
  nrow(tier_c), "; median speedup ",
  format(stats::median(speed_rows$speedup_vs_r), digits = 6), "x\n",
  sep = ""
)
