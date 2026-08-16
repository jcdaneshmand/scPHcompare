# MV7-H landscape-stress reference-oracle rebind

Date: 2026-08-16

Status: accepted; clean independent-validation rerun authorized

## Failed validation attempt

The first independent stress validation stopped after approximately 121 seconds
while evaluating one of its three prospectively selected R reference pairs.
QUADPACK reported `roundoff error is detected in the extrapolation table`.
The validator published no decision or validation output, so the other 19
landscape groups and every clustering, label, outcome, combination, and claim
stage remained closed.

The production stress result itself did not fail. The frozen seed-20260807
gene-H1 group completed all 7,626 unordered pairs twice, and its primary and
repeat distance files are byte-identical. This rebind changes only the
independent R reference engine used to validate those Rust results.

## Correction and unchanged estimand

Commit `f633d98` versions the reference engine as `landscape_reference_v3` and
extends its existing error-budget-preserving recursive partition bisection to
the two narrowly recognized QUADPACK diagnostics:

- `extremely bad integrand behaviour`; and
- `roundoff error is detected`.

Each recursive split divides the parent absolute-error allowance equally
between its children. Child error estimates are summed, the independent
coarse/fine refinement delta remains part of the achieved error estimate, and
the original global `1e-8` absolute and relative tolerances are unchanged.
Unrelated errors and exhaustion of the maximum split depth remain fatal. The
engine now records the coarse- and fine-pass fallback split counts and the
partition-failure policy in every adaptive result.

No persistence diagram, finite-interval policy, landscape level, integration
estimand, H0/H1 separation, Rust kernel, pair axis, resource cap, or acceptance
tolerance changed. The definition remains all active consecutive levels,
finite intervals only, H0 and H1 evaluated separately, and exact or
error-controlled squared-L2 integration without a fixed grid or level cap.

## Verification before validator rerun

The package-aware focused reference/API/workflow/MV7-H suite passes. The
canonical source suite passes with only its two established optional MV5-BC
skips.

The repaired engine was then exercised directly against the frozen production
stress artifacts:

| Stratum | Finite intervals | R squared distance | Absolute Rust/R error | Acceptance tolerance | Splits (coarse/fine) | Result |
|---|---:|---:|---:|---:|---:|---|
| median depth | 773 + 1,459 | 0.0008206492 | 8.101708e-13 | 8.495744e-11 | 0 / 0 | pass |
| maximum depth | 3,172 + 2,635 | 0.001074764 | 8.969279e-13 | 6.188235e-12 | 1 / 1 | pass |

The maximum-depth pair therefore demonstrates the repaired path on the actual
full-corpus H1 workload while retaining a certified error bound. A clean
independent-validation rerun is authorized in new private and public output
roots. The remaining 19 landscape groups may be authorized only if that rerun
passes every frozen check.
