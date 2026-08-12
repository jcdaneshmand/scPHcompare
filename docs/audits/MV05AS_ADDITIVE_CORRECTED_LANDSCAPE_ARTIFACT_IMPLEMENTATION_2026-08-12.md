# MV5-AS additive corrected-landscape artifact implementation

Date: 2026-08-12

Authorization: MV5-AR completion `a16e8c2`

Implementation: `2aa502e`

Realistic-smoke runner and validator: `5db4247`

Aggregate-validation hardening: `4f06337`

## Outcome

MV5-AS passes. The postprocessing and unified workflows now expose an explicit
`corrected_landscape_control = NULL` boundary. When enabled, it produces a
versioned corrected-landscape sidecar without populating, rewriting, or
redirecting the historical landscape fields. When omitted, the established
workflow path and defaults remain unchanged.

This implementation preserves the dissertation-aligned definition: every
finite persistence interval contributes; every active consecutive landscape
level is used; H0 and H1 are calculated and stored separately; the combined
matrix is descriptive; and numerical integration is exact or independently
error-controlled. There is no fixed grid fallback, level cap, interval
removal, or tolerance relaxation.

## Implemented boundary

The additive producer includes:

1. strict normalization of the MV5-AR control contract;
2. admission-only resource planning before calculation;
3. a hash-bound input manifest;
4. one create-only atomic pair shard per canonical sample pair;
5. a pair index reconstructed from verified shards;
6. separate H0, H1, and descriptive combined matrices;
7. provenance and a hash-bound completion marker written last; and
8. verified resume that refuses mismatched or corrupted artifacts.

The corrected-only postprocessing route deliberately returns a sidecar and
does not create the legacy landscape list or legacy combined-matrix fields.
`run_modular_analysis()` and cross-iteration analysis do not accept or consume
the corrected artifacts in this sprint.

## Atomic interruption and resume evidence

The focused test contract interrupts after the first pair shard, resumes from
that verified shard, reconstructs the complete matrix, and checks it directly
against the public pairwise API. It also verifies immutable completed resume,
corruption refusal, cache/signature boundaries, serialization, resource
refusal, strict-control rejection, and the absence of legacy landscape fields.

Forty focused expectations pass. The complete repository suite also passes,
with only the two established CRAN-guard skips. A package tarball built from
the exact Git index passes `R CMD check --no-manual --no-build-vignettes` with
`Status: OK`.

## Realistic smoke

Two clean processes used the frozen `bone__integrated__cell_topology_v1`
stratum: three existing persistence diagrams, three canonical pairs, H0 depth
383 for all diagrams, and H1 depth 130--146. The smoke entered through the
actual corrected-only `run_postprocessing_pipeline()` path.

All six H0/H1 pair calculations used `exact_breakpoint_stream_v1` and passed
strict certification. Scientific calculation time was 14.757 and 14.671
seconds. Total process wall time was 46.05 and 31.78 seconds; peak resident
memory was 942,567,424 and 943,030,272 bytes. Both runs stayed below the frozen
120-second and 1.5-GiB bounds.

The public pair, matrix, resume, and artifact-manifest ledgers are identical
across clean runs. The execution ledger is identical after excluding elapsed
time. Completed resume leaves every artifact path, byte count, modification
time, and SHA-256 unchanged.

## Independent validation

The independent validator passes twice with byte-identical output across 15
categories:

- frozen realistic scope and strict exact certification;
- artifacts-only and public-repeat boundaries;
- private versioned-contract verification;
- pair-shard-to-matrix reconstruction;
- immutable resume and completion-bound artifact schema;
- default-NULL and legacy-source immutability boundaries;
- clean repeat and wall/RSS limits;
- implementation scope and 11-file source freeze;
- 15 prohibited-change counters at zero; and
- a bounded continuation decision.

Private persistence diagrams, pair shards, matrix RDS files, process logs, and
cache directories remain under ignored `tmp/` storage. Only public-safe CSV
evidence and this report are tracked. The dissertation and paper PDFs remain
untracked as previously directed, and `example_run.r` remains untouched and
untracked.

## Decision and next sprint

MV5-AS accepts the additive producer and authorizes MV5-AT only: a broader
realistic workflow smoke across the already frozen existing-data strata. That
sprint should test orchestration, resource admission, atomic resume, and
artifact identity at broader but still bounded scope before any consumer is
allowed to use corrected matrices.

Corrected downstream clustering, visualization, modular/cross-iteration
consumption, workflow-default changes, legacy artifact rewrites, new data,
biological claims, and Rust or other optimization remain unauthorized.
