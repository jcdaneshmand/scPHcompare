# MV5-AP-R1 realistic landscape compatibility rerun

Date: 2026-08-12

Authorization: MV5-AQ completion `49bce6b`

Prospective contract: `a74e523`

Validator hardening: `8a6be9f`

## Outcome

MV5-AP-R1 passes. The repaired dissertation-aligned persistence-landscape
engine completed the full frozen realistic gate twice with certified results,
identical scientific evidence, and bounded resources. The outcome authorizes
only a later opt-in workflow-integration **prefreeze**. It does not authorize
integration itself, a workflow default change, a legacy artifact rewrite,
project-data recomputation, or a manuscript claim.

## Frozen scope

The exact MV5-AP subset blob at completion `6d28da2` remains unchanged. It
contains the minimum, middle-order, and maximum-H1-depth diagrams in each of
eight accepted MV-04 strata: 24 diagrams spanning cell and gene views, bone
and large cohorts, and the accepted representations.

Before execution, MV5-AP-R1 froze all three unordered pairs inside each
three-diagram stratum: 24 total pairs. No cross-stratum comparison was added.
All 24 file hashes, diagram hashes, stored hashes, classes, and scientific
eligibility flags reproduced.

## Landscape definition and routing

The calculation retains every finite interval and every active consecutive
landscape level. H0 and H1 remain separate primary outputs; combined distance
is descriptive. Scientific distance is the squared-L2 landscape integral.
There is no fixed grid, level cap, interval deletion, or tolerance relaxation.

`method = "auto"` with an exact guard of 500 chooses only the numerical engine:

- 18 pairs used exact critical-breakpoint integration for both H0 and H1;
- six high-depth gene-view pairs used exact H0 and certified adaptive H1;
- zero pair used adaptive H0 because all frozen H0 counts are at most 499.

All 48 dimension-specific results passed their strict absolute/relative
`1e-8` global certificate. The largest adaptive H1 achieved error estimate was
`7.39880071420669e-11`, below its pair-specific threshold. A failed certificate
would have aborted the atomic stratum without an accepted matrix.

## Strict sentinel

The original MV5-AP failure pair was rerun with exact and adaptive methods:

| Dimension | Exact | Adaptive | Absolute difference |
|---|---:|---:|---:|
| H0 | 54.1735947138941 | 54.1735947138941 | 7.1054e-15 |
| H1 | 6.23387432164545 | 6.23387432146666 | 1.7879e-10 |
| Combined | 54.5310879525003 | 54.5310879524799 | 2.0435e-11 |

All differences pass the frozen `1e-8` comparison limit.

## Compatibility evidence

Explicit `legacy_k1_unit_grid_v0` was evaluated for the same 24 pairs. Twelve
legacy combined distances were zero, consistent with the historical unit grid
missing informative filtration ranges in part of this realistic corpus. These
comparisons are descriptive only: all 24 rows record `winner_selected = FALSE`,
and legacy mode never entered a scientific calculation or routing decision.

## Resources

Two complete clean runs produced:

| Measure | Run A | Run B |
|---|---:|---:|
| Total scientific time | 1,095.556 s | 1,107.943 s |
| Total external wall time | 1,285.96 s | 1,280.65 s |
| Deepest unit scientific time | 542.960 s | 546.429 s |
| Deepest unit wall time | 564.58 s | 567.94 s |

Peak measured RSS across all units and both runs was 990,363,648 bytes. Every
unit stayed below the prospectively frozen 600-second wall, 2-GiB RSS, and
100-MiB serialized-size limits; each complete run stayed below 3,600 wall
seconds.

The deepest gene-view unit's roughly 32-second margin is acceptable under the
frozen gate but operationally tight. This is evidence for later exact-contract
optimization—potentially including Rust—not permission to simplify the
landscape definition.

## Repeat, serialization, and independent validation

- All stable public scientific/provenance fields repeated exactly.
- All eight runtime-stripped private stratum payloads repeated exactly.
- All 16 within-run serialization checks preserved objects, matrices, and
  scientific/legacy cache keys exactly.
- Every public CSV distance agreed with both private matrix values within the
  explicit `1e-12` decimal round-trip limit; the measured maximum discrepancy
  was `4.97e-14`.
- All input sizes, timestamps, and hashes remained unchanged.
- Seventeen independent validation categories passed twice, and the two
  validator outputs were byte-identical.
- Twelve prohibited-change counters are zero.

Private matrix payloads and process-level raw evidence remain under ignored
`tmp/`. Public aggregate tables contain no private file paths.

## Decision and next sprint

MV5-AP-R1 authorizes MV5-AR, a prospective opt-in workflow-integration
prefreeze only. MV5-AR must map versioned inputs/outputs, cache and resume
behavior, coexistence with explicit legacy mode, resource controls, migration
and rollback, validation, and abort rules before any workflow source is
changed. It should explicitly account for the narrow maximum-depth runtime
margin and keep optimization separate from scientific semantics.
