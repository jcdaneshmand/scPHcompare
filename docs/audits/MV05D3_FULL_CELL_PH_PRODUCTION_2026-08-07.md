# MV5-D3 full cell-PH production audit

| Field | Result |
|---|---|
| Date | 2026-08-07 |
| Scope | All 6,750 frozen SCT cell views |
| Coverage | 75 fold-seed groups; 15 folds; 5 seeds |
| Roles | 450 held-out and 6,300 training views |
| Outcome-label state | Closed |
| Result | MV5-D3 complete; landscape stage now eligible for separate specification |

## Outcome

MV5-D3 completed the corrected SCT cell persistence-diagram cache without
opening labels or starting landscapes. All 6,750 384-cell by 30-PC views
produced complete typed H0/H1 records. Independent validation passed every
record identity, hash, and diagram invariant. No landscape, distance,
clustering, integration, gene-view, or biological-outcome job ran.

## Production scope and execution

The public manifest contains exactly 75 groups of 90 views. Grouped execution
loads each MV5-D1 fold cache once, while every view retains an independent
immutable result and provenance identity. Two admission groups completed before
the remaining 73 groups were allowed to run. The first group's 90 records and
checkpoint remained unchanged through an explicit resume attempt.

| Measurement | Result |
|---|---:|
| Groups completed / failed | 75 / 0 |
| Views completed / failed | 6,750 / 0 |
| H0 intervals | 2,592,000 |
| H1 intervals | 2,096,835 |
| Group worker time | 3,767.759 s (1.047 h) |
| Typed PH plus MST operation time | 980.385 s |
| Median / maximum group elapsed | 50.072 / 52.673 s |
| Maximum process-tree RSS | 286,371,840 B (273.1 MiB) |
| Private diagram storage | 196,293,225 B (187.2 MiB) |
| Implementation SHA-256 | `c1f4de704278d4f6fe971b9a434e1f98733cf5290aeb1509aea602de4875764a` |
| Manifest SHA-256 | `af1385b4fa11e3b747103882e9ba2c74ce680b44ba06a29f980d1b737cc281f5` |

The 900-second, 4-GiB per-group and 43,200-second stage guards were not
approached. The measured grouped PH time is substantially below the earlier
one-view-per-process pilot projection because fold-cache deserialization and
process startup are amortized safely across each group.

## Correctness and immutability

| Validation | Result |
|---|---:|
| Result-file hashes and static record invariants | 6,750 / 6,750 pass |
| Stored full-view H0 MST checks | 6,750 / 6,750 pass |
| Fresh full-view H0 MST recomputations | 75 / 75 pass |
| Maximum recorded fresh MST absolute error | 0 |
| Independent repeated production group | 90 / 90 pass |
| Object-identical repeated records | 90 / 90 |
| Byte-identical repeated files | 90 / 90 |
| Resume-preserved first-group records/checkpoint | 91 / 91 files |

The fresh MST oracle selected one deterministic held-out-priority view from
each group, reloaded its original MV5-D1 coordinates, and compared all 383
finite H0 deaths against the full 384-cell Euclidean MST. This is independent
of the oracle evidence stored during production.

## Measured pre-landscape projection

| Component | Worker-hours |
|---|---:|
| Measured SCT normalization | 2.562 |
| Measured training-only cell coordinates | 2.376 |
| Measured full cell PH | 1.047 |
| Projected exact all-level landscape distances | 3.572 |
| Measured plus projected total | 9.556 |
| Planning cap with reserve | 21.600 |
| Remaining margin | 12.044 |

The PH stage is complete, but the 3.572-hour landscape term remains a
projection. It must be replaced by measured evidence in a separate monitored
sprint. Passing this resource gate does not authorize outcomes or claims.

## Landscape contract carried forward

The next stage must use the project-owner-approved dissertation-aligned
definition: separate H0 and H1; exclude the one infinite H0 class and any
non-positive-persistence interval; retain all active consecutive levels; and
use exact or error-controlled L2 integration on each dimension's natural
support. No fixed level cap or fixed uniform grid is permitted. A combined
H0/H1 distance is secondary and must preserve both component distances and the
H1 squared-distance contribution.

## Privacy and scope

- Private records, checkpoints, repeats, and logs remain under ignored
  `tmp/mv05d3` paths.
- PDFs, reviewer correspondence, and `example_run.r` remain untracked.
- No dependency or lockfile was changed and nothing was pushed.
- Outcome labels remained closed and biological outcomes computed remained
  false in every public completion artifact.
- The complete current-source test suite passed with only its two intentional
  CRAN skips.
- The isolated source-package check under R 4.4.1 reported `Status: OK` with
  zero errors, warnings, or notes.

## Decision table

| Question | Disposition |
|---|---|
| Complete corrected cell diagram cache? | Yes: 6,750/6,750 independently valid |
| Correctness and determinism demonstrated? | Yes: all-record invariants, MST oracles, exact group repeat, and resume |
| Resource gate passed? | Yes: 9.556 h measured-plus-projected, 12.044 h below cap |
| Biological interpretation permitted? | Prohibited |
| Next action | Specify and admit a bounded label-closed exact all-level SCT cell landscape-distance stage |
