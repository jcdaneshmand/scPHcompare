# MV7-H complete-landscape oracle v4 rebind

Date: 2026-08-16

Status: accepted; clean complete-corpus validator rerun authorized

## Failed-closed validator attempt

The first complete-corpus validation attempt stopped after 167.2 seconds while
evaluating the prospectively frozen cell-H1 R oracle. QUADPACK reported
`roundoff error was detected`. Reference v3 recognized the equivalent
`roundoff error is detected ...` wording but not this tense variant. No public
completion-validation directory or authorization decision was written. All 20
production groups and four byte-identical component repeats remained unchanged,
and MV7-I stayed closed.

## Narrow correction

Commit `aca6716` versions the reference engine as `landscape_reference_v4` and
matches QUADPACK's common `roundoff error` message stem. The recoverable set is
still restricted to QUADPACK bad-integrand and roundoff diagnostics. Unrelated
errors, invalid midpoints, and exhaustion of the maximum split depth remain
fatal. Recursive children still divide the parent absolute-error allowance in
half, and their error estimates are summed; the global `1e-8` absolute and
relative tolerances are unchanged.

The exact frozen cell-H1 pair that triggered the failure contains 533 and 492
finite intervals. Under v4 it passes with:

- Rust squared distance: 0.1876789;
- absolute Rust/R discrepancy: `1.053326e-10`;
- certified acceptance tolerance: `1.062933e-09`;
- coarse/fine fallback splits: 6/27; and
- reference `within_requested_tolerance`: true.

The focused landscape/API/workflow/MV7-H suite passes, and the canonical source
suite passes with only the two established optional MV5-BC skips. No persistence
diagram, landscape level, pair, distance, Rust output, error tolerance, or
scientific estimand changed. A clean rerun of the same frozen complete-corpus
validator is authorized; MV7-I remains closed pending an explicit full pass.
