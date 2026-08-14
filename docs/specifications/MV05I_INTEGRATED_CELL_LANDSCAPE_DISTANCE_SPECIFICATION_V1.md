# MV5-I integrated cell landscape-distance specification v1

Date: 2026-08-09
Status: frozen before admission or production execution
Predecessor: accepted MV5-H integrated cell persistence diagrams at local revision `7486c8d`

## Purpose

MV5-I computes the label-closed persistence-landscape distances needed for the
integrated cells-as-observations view. It transfers the accepted MV5-D4 exact
landscape definition to the corrected MV5-H representation without changing
the frozen LOSO folds, seeds, samples, normalization maps, coordinates, or
persistence diagrams.

This is a precomputation and method-integrity sprint. It does not run retrieval,
rank samples, inspect biological outcomes, cluster samples, build gene
topology, fuse views, introduce new data, or make manuscript claims.

## Frozen scientific population

- 75 fold-seed groups: 15 leave-one-study-out folds crossed with seeds
  20260805 through 20260809.
- 6,750 immutable MV5-H views: 90 samples per group.
- Each view is the accepted 384-cell by 30-coordinate
  `inductive_integrated` cells-as-observations representation.
- Eligible biological pairs are directed held-out queries against training
  references within the same fold-seed group.
- The manifest must contain exactly 35,350 biological pairs and two separate
  dimension rows per pair, for exactly 70,700 requests.
- H0 and H1 remain separate primary components. Their Euclidean component
  combination may be emitted only as a secondary, explicitly identified
  convenience value.

## Frozen persistence-landscape definition

For each MV5-H diagram and homology dimension:

1. Select only finite intervals with strictly positive persistence.
2. Exclude the single essential H0 interval before landscape construction.
3. Retain every consecutive active landscape level. No universal level cap is
   permitted.
4. Construct exact piecewise-linear landscapes from critical pairs.
5. For a pair, subtract landscapes level by level, retaining all active
   difference levels.
6. Integrate the squared difference exactly on every consecutive linear
   segment:
   `width * (y0^2 + y0*y1 + y1^2) / 3`.
7. Report the nonnegative square root as the L2 landscape distance.

No fixed uniform grid, interpolation grid, truncation level, or approximate
quadrature is part of the primary definition. Source active-level counts,
difference-level counts, segment counts, and critical-point counts are recorded
for auditable accounting. `level_cap_applied` must be false for every result.

## Identity and immutability

Every request binds:

- group, fold, seed, H0/H1 dimension, query job, and training job;
- query and training MV5-H record cache keys;
- query and training diagram SHA-256 values;
- query and training result-file SHA-256 values; and
- the landscape definition identifier.

The request identifier is a SHA-256 digest of that identity. Requests are
ordered deterministically and partitioned into immutable chunks of at most 250
rows. Each chunk identifier is a digest of its ordered request identifiers.
Existing input, output, and status artifacts are reusable only after their
identities and hashes validate; partial or stale artifact sets are errors and
must never be overwritten.

## Execution and resource guards

Production is authorized by MV5-H evidence and must use:

- at most two simultaneous heavy Python workers;
- at most 250 requests per chunk;
- 900 seconds per chunk;
- 1,800 seconds per fold-seed group;
- 8 GiB process-tree RSS per group;
- 43,200 aggregate group-worker seconds for the distance stage; and
- 10 GB total projected project storage.

The monitor records group elapsed time, exact pair-operation time, peak
process-tree RSS, private bytes, chunk count, and completed H0/H1 rows. A guard
breach stops new launches but lets active workers finish naturally so that
their state remains inspectable. It never kills WSL or worker processes.

## Admission gates

Before full production:

- build and reproduce the 70,700-row pair and 360-chunk manifests;
- stage all 75 input bundles while independently validating all 6,750 MV5-H
  record/file/diagram identities;
- run one representative group and the maximum-request group;
- validate every admission request/result/status and all zero-scope counters;
- compare one eligible H0 and one eligible H1 result per admission group with
  the independent exact R breakpoint-stream oracle;
- prove immutable resume for an admission group; and
- confirm elapsed, RSS, and storage feasibility under the production guards.

## Completion gates

MV5-I is complete only if all of the following hold:

- 75 groups, 360 chunks, and 70,700 dimension rows complete;
- exactly 35,350 H0 and 35,350 H1 rows exist and pair one-to-one;
- every request/source/chunk/file identity validates independently;
- every distance and squared distance is finite, nonnegative, and mutually
  consistent;
- every result reports exact integration, all active levels, no level cap, and
  complete source-level accounting;
- essential H0 and nonpositive intervals are excluded before construction;
- retrieval, clustering, gene-view, fusion, new-data, and biological-outcome
  counters remain zero;
- no self or reciprocal rows occur in the frozen directed rectangular request
  scope; therefore no square matrix is assembled and symmetry/diagonal checks
  are explicitly not applicable rather than silently claimed;
- one complete group repeat is byte-identical for every distance artifact and
  status identity;
- a complete production resume rebuilds zero chunks and changes zero immutable
  artifact bytes;
- focused tests, the complete test suite, and a staged public-only package
  check pass; and
- a measured post-distance projection explicitly authorizes or rejects the
  later label-closed integrated retrieval-input sprint.

## Stop boundary

Stop after producing validated distance and component artifacts plus the
post-distance authorization. Do not compute retrieval rankings, baselines,
performance metrics, biological outcomes, clustering, gene-side topology,
fusion, new-data analyses, or manuscript conclusions in MV5-I.
