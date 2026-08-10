# MV5-I integrated cell landscape-distance audit

Date: 2026-08-09
Source revision: `7486c8d` (accepted MV5-H)
Representation: corrected integrated cells as observations
Outcome state: closed

## Outcome

MV5-I is complete. The accepted MV5-D4 exact persistence-landscape definition
has been transferred without approximation to all 6,750 validated MV5-H
integrated cell H0/H1 records. All 35,350 frozen held-out-query-to-training
biological pairs produced separate H0 and H1 rows, for 70,700 exact distances
in 360 immutable chunks.

This sprint establishes validated topology-distance inputs only. It did not run
retrieval or baselines, open biological outcomes, rank methods, cluster samples,
construct gene topology, fuse views, add data, or make manuscript claims.

## Landscape definition and interval boundary

The production definition is:

- finite strictly positive-persistence intervals only;
- one essential H0 interval excluded per source diagram;
- H0 and H1 constructed and reported separately;
- every consecutive active landscape level retained, with no universal cap;
- exact piecewise-linear subtraction on critical pairs; and
- exact L2 integration across every linear segment, with no fixed grid or
  approximate quadrature.

Staging independently reopened every MV5-H record and matched its file hash,
record cache key, and diagram hash to public evidence.

| Staged item | Count |
|---|---:|
| Groups | 75 |
| MV5-H records | 6,750 |
| Finite positive H0 intervals | 2,585,250 |
| Finite positive H1 intervals | 1,545,943 |
| Essential H0 intervals excluded | 6,750 |
| Total staged finite intervals | 4,131,193 |
| Private staged request/interval bytes | 1,188,929,209 |

MV5-H had already certified zero invalid and zero-persistence diagram
intervals. MV5-I nevertheless applies the finite/positive filter explicitly
before construction. Every result records query and training active-level
counts, their total processed count, difference-level/segment/critical-point
counts, `all_active_levels = TRUE`, and `level_cap_applied = FALSE`.
Across production, 29,559,529 source active levels were processed.

## Frozen pair and chunk population

The deterministic manifest contains:

| Axis | Result |
|---|---:|
| LOSO folds | 15 |
| Seeds | 5 |
| Fold-seed groups | 75 |
| Biological query-to-training pairs | 35,350 |
| H0 rows | 35,350 |
| H1 rows | 35,350 |
| Total dimension rows | 70,700 |
| Chunks | 360 |
| Maximum rows per chunk | 250 |

Each request binds the group/fold/seed/dimension, query and training job IDs,
record cache keys, diagram hashes, result-file hashes, and landscape-definition
identity. Chunk IDs bind ordered request IDs. The compressed public manifest
reproduces the exact uncompressed CSV content used for execution.

## Admission and independent exact oracles

The smallest and largest real production groups were admitted prospectively:

| Admission group | Rows | Chunks | Chunk seconds | Maximum chunk | Peak RSS |
|---|---:|---:|---:|---:|---:|
| SRA713577 / 20260805 | 178 | 2 | 8.622 | 7.451 s | 234,946,560 B |
| SRA779509 / 20260805 | 3,250 | 14 | 48.780 | 9.851 s | 245,092,352 B |

One eligible H0 and one eligible H1 pair in each group were recomputed with the
independent R exact breakpoint-stream oracle. All four passed. Absolute
differences were 0 to (7.11	imes10^{-14}), far inside the prospective
tolerances. A representative admission resume rebuilt zero chunks and left all
four files unchanged.

## Complete production resources

Two heavy workers executed all 75 groups under 1,800-second group, 8-GiB RSS,
43,200-worker-second stage, and 10-GB projected-storage guards.

| Measure | Result | Guard |
|---|---:|---:|
| Groups | 75 / 75 | 75 |
| Distance rows | 70,700 / 70,700 | 70,700 |
| Chunks | 360 / 360 | 360 |
| Group-worker time | 2,072.893 s (0.576 h) | 43,200 s |
| Exact pair-operation time | 1,120.978 s | descriptive |
| Median group time | 25.486 s | 1,800 s |
| Maximum group time | 61.261 s | 1,800 s |
| Peak process-tree RSS | 255,774,720 B (243.9 MiB) | 8 GiB |
| Private output/status bytes | 66,554,027 | recorded |
| Private staging plus output bytes | 1,255,483,236 | included in projection |

There were zero group failures and no guard violations.

## Global correctness, component, and matrix-boundary checks

Independent validation matched every request ID, source record identity,
diagram hash, chunk identity, implementation hash, manifest hash, status hash,
and output hash. Every distance and squared distance is finite, nonnegative,
and algebraically consistent. The 35,350 H0 and 35,350 H1 rows pair one-to-one.
A secondary component artifact reports H0, H1, and
(sqrt{d_{H0}^2+d_{H1}^2}), while preserving H0/H1 as primary components.

The frozen workload is a directed rectangular held-out-query by
training-reference relation. It contains zero self rows and zero reciprocal
rows and assembles zero square matrices. Symmetry and diagonal checks are
therefore explicitly not applicable here rather than silently asserted. Any
later method that assembles a square matrix must validate symmetry and a zero
diagonal in its own contract.

Retrieval, clustering, gene-view, fusion, new-data, and biological-outcome
counters are zero in every scientific result and completion record.

## Determinism and resume

A fresh complete repeat of the maximum 3,250-row group rebuilt all 14 chunks.
All 14 scientific distance files are byte-identical to production, and all
static status identity fields match. The complete production tree was then
snapshotted, the 75-group monitor rerun, and snapshotted again. It rebuilt zero
groups; all 720 distance/status files retained identical sizes and SHA-256
values.

## Post-distance projection and decision

The projection uses complete measured storage, including staged intervals.

| Component | Worker hours | State |
|---|---:|---|
| Integrated coordinates | 9.107 | measured |
| Integrated cell PH | 1.098 | measured |
| Integrated landscape distances | 0.576 | measured |
| Integrated retrieval inputs | 0.507 | projected with 25% reserve |
| **Total** | **11.289** | **passes 21.6-hour cap** |

Measured coordinate/PH/landscape storage is 2,241,808,051 bytes. With the
reserved retrieval-input projection, total storage is 2,443,617,957 bytes,
below the 10-GB cap. The existing upstream peak of 3,303,829,504 bytes and
790.95-second maximum group remain below their caps.

Decision: authorize a separate label-closed integrated retrieval-input sprint.
That later sprint must stop before biological evaluation and cannot reuse the
negative SCT retrieval outcome to tune or select methods.

## Verification and public boundary

The focused MV5-I suite passes 8/8 expectations. The complete source-loaded
suite passes 644/644 expectations with zero failures, errors, warnings, or
skipped blocks. The isolated source-package build and installed-package
`R CMD check --no-manual` pass with `Status: OK`, zero errors, zero warnings,
and zero notes. The tarball has 141 entries and zero paths matching docs,
temporary data, PDFs, reviewer material, or `example_run.r`. No dependency was
installed and no lockfile or dependency declaration changed.

Private staged intervals, production chunks/status/logs, repeats, and MV5-H
records remain under ignored `tmp/`. PDFs, reviewer correspondence, private
caches, and `example_run.r` remain untracked. Nothing was pushed.

MV5-I proves that the dissertation-aligned all-level landscape calculation is
well-defined, exact, reproducible, and feasible on the corrected integrated
cell view. It does not establish biological usefulness or superiority.
