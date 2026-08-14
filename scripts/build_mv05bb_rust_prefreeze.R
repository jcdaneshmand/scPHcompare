#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: build_mv05bb_rust_prefreeze.R OUTPUT_DIR")
out <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_public <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)
write_public <- function(value, name) utils::write.csv(
  value, file.path(out, name), row.names = FALSE, quote = TRUE
)

decision <- read_public("docs/audits/mv05ba-continuation-decision-2026-08-13.csv")
summary <- read_public("docs/audits/mv05ba-benchmark-summary-2026-08-13.csv")
equivalence <- read_public("docs/audits/mv05ba-equivalence-summary-2026-08-13.csv")
if (nrow(decision) != 1L || !decision$rust_trigger_satisfied ||
    !decision$rust_prefreeze_authorized || decision$rust_implementation_authorized ||
    decision$corrected_persim_adopted || !decision$accepted_r_engine_retained ||
    nrow(summary) != 1L || summary$median_candidate_speedup >= 1 ||
    !equivalence$corrected_persim_mathematically_supported) {
  stop("MV5-BA does not authorize the MV5-BB prefreeze.")
}

kernel <- data.frame(
  contract_id = "mv05bb_rust_kernel_contract_v1",
  boundary_id = c("input", "canonicalization", "dimensions", "intervals",
                  "levels", "integration", "error", "output", "ordering",
                  "memory", "parallelism", "fallback", "forbidden"),
  frozen_requirement = c(
    "two canonical numeric persistence diagrams plus requested homology dimension; no paths, metadata, labels, or outcomes",
    "R validates dimensions, finite births, positive finite persistence, stable row ordering, and source hashes before FFI",
    "one H0 or H1 call at a time; R retains separate H0/H1 result and certification objects",
    "finite deaths only; essential interval counts remain R provenance and are never silently passed as finite",
    "all consecutive active landscape levels with zero padding at missing depth; no universal cap",
    "squared L2 over exact piecewise-linear differences; no uniform grid; compensated or pairwise summation required",
    "prototype is exact critical-event arithmetic; no adaptive tolerance claim until separately specified and certified",
    "squared distance, active-level count, segment/event count, finite input counts, engine version, and status code",
    "deterministic lexical pair order and stable event tie policy; normalized outputs byte-identical across repeats",
    "pair-bounded allocation; 1 GiB prototype RSS cap and no unbounded landscape tensor or all-corpus retention",
    "one internal thread; process-level concurrency remains at most two independent pairs",
    "any nonzero status, panic boundary, nonfinite value, negative value beyond tolerance, or failed equivalence returns control to canonical R; no partial adoption",
    "no PH generation, normalization, view construction, cache naming, file IO, clustering, k selection, labels, outcomes, combined-primary distance, or workflow-default change"
  ), stringsAsFactors = FALSE
)

ffi <- data.frame(
  contract_id = "mv05bb_ffi_contract_v1",
  field_order = 1:10,
  field = c("abi_version", "dimension", "first_interval_count",
            "first_birth_death", "second_interval_count", "second_birth_death",
            "squared_distance", "active_levels", "event_segments", "status"),
  type = c("u32", "u8", "usize", "const_f64_ptr_interleaved", "usize",
           "const_f64_ptr_interleaved", "mut_f64_ptr", "mut_u64_ptr",
           "mut_u64_ptr", "i32"),
  ownership = c("value", "value", "value", "R_borrowed_read_only", "value",
                "R_borrowed_read_only", "R_owned_write", "R_owned_write",
                "R_owned_write", "value"),
  requirement = c(
    "must equal 1", "0 or 1", "bounded by usize and preflight resource gate",
    "length exactly 2*n; finite; immutable during call",
    "bounded by usize and preflight resource gate",
    "length exactly 2*n; finite; immutable during call",
    "written only on status zero", "written only on status zero",
    "written only on status zero", "zero success; stable named nonzero failures"
  ), stringsAsFactors = FALSE
)

corpus <- data.frame(
  contract_id = "mv05bb_equivalence_corpus_v1",
  tier = c("A", "B", "C", "D", "E"),
  scope = c("analytical_and_sign_crossing", "tractable_R_exact",
            "MV5BA_worst_depth", "MV5AY_exact_only", "MV5AY_adaptive_H1"),
  required_results = c(3L, 20L, 12L, 318L, 90L),
  acceptance = c(
    "squared distance absolute error <= 1e-12",
    "squared distance absolute error <= 1e-10 plus 1e-10 relative",
    "within each accepted exact/adaptive MV5-AY certificate",
    "squared distance absolute error <= 1e-10 plus 1e-10 relative",
    "within accepted R achieved absolute error plus 100 machine eps times scale"
  ),
  required_before = c("prototype", "prototype", "speed_decision",
                      "production_adoption", "production_adoption"),
  stringsAsFactors = FALSE
)

gate <- data.frame(
  contract_id = "mv05bb_prototype_gate_v1",
  gate_id = c("toolchain", "source", "panic", "sanitizers", "fixtures",
              "worst_depth", "determinism", "speed", "memory", "fallback",
              "package", "adoption"),
  requirement = c(
    "Rust toolchain version and installer/source hashes frozen before any build; no toolchain is currently present",
    "minimal crate and lockfile committed; dependency count justified; unsafe code denied except reviewed FFI shim",
    "panic cannot unwind across FFI; stable status code and R fallback",
    "Miri or equivalent memory-safety test where applicable plus address/undefined sanitizer or valgrind boundary test",
    "all Tier A and B equivalence checks pass",
    "all six frozen MV5-BA worst-depth pairs pass 12/12 certificates",
    "two clean builds/runs yield identical normalized outputs and result identities",
    "at least 3x median speedup versus accepted R on same six-pair panel; no pair slower than R",
    "peak RSS <= 1 GiB per pair-bounded process",
    "forced status/error fixtures prove canonical R fallback and no partial artifact",
    "exact-staged R package check remains Status OK with Rust unavailable and with prototype available",
    "production adoption requires separate gate plus complete Tier D/E corpus; prototype success alone cannot alter defaults"
  ),
  abort_on_failure = TRUE,
  stringsAsFactors = FALSE
)

decision_out <- data.frame(
  contract_id = "mv05bb_decision_v1",
  prefreeze_accepted = TRUE,
  rust_prototype_authorized = TRUE,
  rust_toolchain_install_authorized = FALSE,
  rust_production_adoption_authorized = FALSE,
  r_engine_retained = TRUE,
  corrected_persim_role = "independent_oracle_only",
  additional_seed_production_authorized = FALSE,
  partitions_authorized = FALSE,
  next_sprint = "MV5-BC_bounded_rust_kernel_prototype_requires_toolchain_approval",
  stringsAsFactors = FALSE
)

write_public(kernel, "mv05bb-rust-kernel-contract-2026-08-13.csv")
write_public(ffi, "mv05bb-ffi-contract-2026-08-13.csv")
write_public(corpus, "mv05bb-equivalence-corpus-2026-08-13.csv")
write_public(gate, "mv05bb-prototype-gate-2026-08-13.csv")
write_public(decision_out, "mv05bb-continuation-decision-2026-08-13.csv")
cat("MV5-BB prefreeze: bounded Rust prototype specified; toolchain install closed\n")
