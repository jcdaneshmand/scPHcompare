# MV5-AP-R1 realistic landscape rerun specification v1

Date: 2026-08-12

Authorization: MV5-AQ completion `49bce6b`

## Purpose

Complete the read-only compatibility and resource evaluation that MV5-AP
stopped. This sprint tests whether the repaired public landscape engine is
ready for a later, separately prefrozen opt-in integration design. It does not
integrate the engine into a workflow and does not change any historical
artifact.

## Frozen inputs

The exact tracked MV5-AP depth-triplet table is authoritative. Its 24 rows are
the deterministic minimum, middle-order, and maximum H1-depth diagrams in each
of eight accepted MV-04 strata. Every diagram ID, path, diagram hash, result
file hash, stored provenance hash, and eligibility flag must reproduce.

The subset may not be reselected, enlarged, reduced, or replaced based on
numerical or biological outcomes.

## Comparison units

The complete evaluation contains the three unordered pairs among the three
frozen diagrams inside each stratum: 8 strata x 3 pairs = 24 pairs. No
cross-stratum pair is part of this gate because cohort, representation, view,
and filtration scale define the compatibility strata and cross-stratum pairs
were not frozen by MV5-AP.

## Numerical contract

- scientific mode only for the primary result;
- `method = "auto"`, `exact_max_intervals = 500`;
- absolute and relative tolerances `1e-8`;
- all finite intervals and all active consecutive levels;
- H0 and H1 separate; combined is descriptive;
- squared-L2 integration;
- exact critical-breakpoint integration when routed;
- adaptive output accepted only with its global error certificate;
- no grid fallback, level cap, interval removal, or tolerance relaxation;
- any failed certificate aborts without a partial accepted matrix.

Explicit `legacy_k1_unit_grid_v0` is run on the same pairs only as descriptive
compatibility evidence. It cannot select a winner or alter a default.

## Execution and resource policy

Each stratum is an atomic unit with one scientific matrix, one explicit legacy
matrix, public pair diagnostics, private versioned RDS, and serialization
check. Run strata sequentially with a 600-second hard wall-time ceiling per
unit. The complete eight-stratum pass must remain below 3,600 seconds, each
process below 2 GiB maximum RSS, and each serialized unit below 100 MiB.

Execute two clean roots. Scientific fields, matrices, diagnostics, cache keys,
certificates, input identities, and runtime-stripped private payloads must be
identical. Wall time, maximum RSS, and serialized byte size are measured fields
and are compared against bounds rather than byte equality.

## Independent validation

The validator must reconstruct the 24-pair plan without using a decision
helper; verify all input hashes; compare every public row to both private matrix
objects; prove matrix symmetry, zero diagonals, canonical order, error
certificates, exact/adaptive routing, explicit legacy isolation, serialization,
clean repeat, resource limits, source immutability, and prohibited-change
counters. The repaired strict sentinel must also be rerun exact and adaptive
and agree within `1e-8` in every reported dimension.

## Decision

Only a complete pass may authorize a later opt-in workflow-integration
prefreeze. Even a pass does not authorize integration, a workflow default
change, artifact rewriting, project-data recomputation, clustering, fusion,
new data, Rust optimization, or manuscript claims.

Any provenance drift, missing pair, uncertified result, exact/adaptive
disagreement, repeat/serialization failure, resource-limit breach,
workflow/legacy mutation, test/check failure, or scope drift stops the sprint.
