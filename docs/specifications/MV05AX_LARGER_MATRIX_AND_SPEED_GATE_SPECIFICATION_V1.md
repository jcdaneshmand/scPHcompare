# MV5-AX larger corrected-matrix and speed gate v1

Date: 2026-08-12

Authorization: owner direction after MV5-AW `5b1dae2`.

## Frozen scope

Use every row of the accepted MV-04 manifest: eight strata, 56 diagrams, and
204 within-stratum pairs. Bone strata contain four samples (six pairs; candidate
`k=2:3`); large strata contain ten samples (45 pairs; candidate `k=2:9`). Cell
and gene views and H0/H1 remain separate. No cross-stratum pair is allowed.

## Execution and speed decision

Six strata are exact-only. The two large gene strata contain 90 adaptive H1
pairs. Accepted three-pair measurements project approximately 4--5 hours of
adaptive work sequentially. Execute the two gene strata as separate concurrent
processes, each internally one-worker, with independent atomic sidecars and a
2-GiB per-process cap. Execute exact-only strata in bounded processes.

This process-level concurrency changes scheduling only, not arithmetic,
ordering within a sidecar, tolerance, landscapes, or cache identity. Rust is a
high-value later candidate for the repeated sort/evaluation kernel, but it must
pass exact/adaptive reference equivalence and is not required before this
bounded production run.

## Stability boundary

Complete matrices make multiple `k` values possible but do not themselves
create stability replicates. Partitioning remains closed after production until
a separate label-closed perturbation/resampling contract is accepted.

Abort on input/hash drift, resource breach, uncertified output, collision,
source mutation, or scientific/default/legacy drift.
