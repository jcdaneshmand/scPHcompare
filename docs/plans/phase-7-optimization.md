# Phase 7 Subplan — Profiling and Optimization

## Objective

Improve end-to-end runtime and memory without altering approved scientific behavior. Rust is a possible implementation choice, not a predetermined outcome.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P7-01 | Freeze correctness corpus | Approved fixtures, reference outputs, tolerances, failure tests | Optimization branch must pass corpus before benchmark comparison | `in_progress` — exact/adaptive H0/H1 corpus now independently cross-checked and includes sign-changing landscape differences; eligible corrected-diagram corpus remains |
| P7-02 | Define benchmark protocol | Hardware, software, warmup, repetitions, inputs, timing/memory/statistics | Protocol produces stable estimates and comparable runs | `in_progress` — candidate scaling records kernel/end-to-end time, process-tree RSS, determinism, timeouts, versions, and SciPy error; eligible-diagram pairwise protocol remains |
| P7-03 | Profile end to end | Stage time, CPU, peak memory, copies/dense conversions, I/O, parallel overhead | Dominant hotspots and scaling curves are measured | `in_progress` — R exact/adaptive, Persim, corrected critical pairs, SciPy, and GUDHI profiled through 5,000 controlled intervals; corrected PH-to-LDM profiling remains |
| P7-04 | Remove/reuse work | Invalid/redundant computation removal and dependency-aware caching | End-to-end improvement measured; outputs remain equivalent | `in_progress` — removed the unconditional 60-second PH polling floor with process-aware waiting; reference hashes unchanged |
| P7-05 | Improve memory/parallelism | Sparse preservation, copy reduction, bounded worker strategy | Peak memory and/or time improves without instability | `in_progress` — corrected critical-pair prototype held near 155 MB through 5,000 controlled intervals; batching and worst-case eligible-diagram behavior remain |
| P7-06 | Evaluate mature libraries | Ripser/sparse/GPU/alternative complex candidates against bottleneck | Keep/replace decision has benchmark and maintenance rationale | `in_progress` — isolated Persim/GUDHI/SciPy trial complete; built-in Persim norm rejected, corrected critical-pair path leads, Python GUDHI grid limited to sensitivity/display |
| P7-07 | Evaluate approximation | Accuracy/stability versus runtime across subsampling/landmarks/thresholds | Trade-off is quantified and acceptable domain tolerance is approved | `not_started` |
| P7-08 | Apply Rust gate | Architecture decision record addressing all Rust criteria | Rust work begins only if every gate criterion passes | `in_progress` — ADR-001/002 retain a negative decision; reconsider only if corrected critical-pair batching remains an eligible-diagram bottleneck |
| P7-09 | Implement narrow Rust component if approved | Versioned API, R bindings, error handling, equivalence tests, cross-platform CI | Material end-to-end benefit and maintainable installation demonstrated | `not_started` |
| P7-10 | Evaluate G7 | Before/after report | Benefit is reproducible and behavior preserved | `not_started` |

## Rust decision criteria

- Stable dominant hotspot remains after algorithmic and caching fixes.
- Existing optimized libraries are insufficient.
- Narrow typed input/output contract exists.
- Numerical equivalence/error tolerances are approved.
- Material end-to-end—not only microbenchmark—benefit is demonstrated.
- Windows, macOS, and Linux installation/maintenance burden is acceptable.

## Gate G7 checklist

- [ ] Benchmark protocol is reproducible.
- [ ] Correctness corpus passes unchanged.
- [ ] Runtime and peak-memory changes are reported with uncertainty.
- [ ] Approximation trade-offs are explicit.
- [ ] Any Rust decision has an architecture record and maintenance owner.
