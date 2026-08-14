# MV5-V streamed full-robustness production prefreeze specification v1

Date frozen: 2026-08-10

Accepted admission base: `2c844f6`

Labels opened: prohibited

Robustness or biological outcomes: prohibited

## Purpose

MV5-V converts the passing MV5-U admission into an exact, resource-bounded,
streamed execution contract. It freezes the full label-closed calculation
scope, atomic publication and resume behavior, independent validation, and
configuration-stratified resource envelope before any full robustness work.

This sprint is a prefreeze. It cannot execute the 600-group calculation,
evaluate retrieval, open labels, compare robustness effects, rank a
configuration or representation, or make a manuscript claim.

## Scientific scope

The later robustness estimand remains the change in the already frozen
held-out-query retrieval effects under each perturbation. It reuses the
MV5-E/MV5-K endpoint, aggregation, study-blocked uncertainty, missingness, and
complete-reporting rules. It does not reopen clustering or add a clustering
robustness search.

Only the four MV5-T one-factor-at-a-time configurations are eligible:

1. nested 192 cells, 30 PCs, Euclidean geometry;
2. nested 256 cells, 30 PCs, Euclidean geometry;
3. 384 cells, first 20 accepted PCs, Euclidean geometry; and
4. 384 cells, 30 PCs, cosine geometry represented as Euclidean chord distance
   after row-unit normalization.

The accepted 384-cell/30-PC Euclidean results are the unchanged reference.
The five original seeds, SCT/integrated representations, H0/H1 separation,
raw unweighted H0/H1 combination, and pseudobulk context are already covered.
No new seed, feature panel, PC count, distance, clustering algorithm,
integration method, or factorial interaction enters this contract.

## Exact complete scope

Cross 15 frozen LOSO folds, five seeds (`20260805` through `20260809`), two
representations, and four configurations:

| Object | Exact count |
|---|---:|
| Atomic group units | 600 |
| Transformed sample views | 54,000 |
| Heldout-training biological pairs | 282,800 |
| H0/H1 landscape request rows | 565,600 |
| Deterministic landscape subchunks (at most 250 rows) | 2,880 |
| Matched energy rows | 282,800 |
| Label-closed assembled method rows | 1,131,200 |
| Maximum-fold clean-repeat groups | 8 |

Each atomic group contains all 90 sample views for one fold, seed,
representation, and configuration. Pair identities are the exact accepted
MV5-D4/MV5-I heldout-query-to-training axes, rebound to the perturbation group
identity. Each biological pair has one H0 and one H1 request. The four assembled
methods are H0 landscape, H1 landscape, descriptive raw H0/H1 Euclidean
combination, and matched cell-distribution energy. The unchanged pseudobulk
comparator is reused from accepted evidence and is not recalculated.

## Coordinate and PH contract

- Read only the 150 accepted coordinate files whose hashes were frozen by
  MV5-T and revalidated by MV5-U.
- Reproduce the representation-neutral SHA-256 nested cell order for 192/256.
- Truncate the accepted 30-coordinate matrix for PC20; do not refit PCA.
- Reject zero or nonfinite norms before cosine-chord normalization.
- Use the typed `cell_topology_v1` view, `ripserr` Euclidean Vietoris-Rips PH,
  dimensions H0/H1, threshold `-1`, and field 2.
- Require 90 transformed views and one essential H0 interval per view.
- Independently validate every finite H0 death multiset against an Euclidean
  minimum-spanning-tree oracle.

## Landscape contract

The revised dissertation-aligned definition is unchanged:

- H0 and H1 stay separate;
- exclude the essential H0 interval;
- retain every consecutive active landscape level;
- apply no universal level cap and no uniform integration grid;
- integrate squared L2 differences exactly over linear critical-pair segments.

One group holds at most 90 landscapes in memory. Pair requests are processed in
deterministic subchunks of at most 250 rows and appended to a private temporary
group artifact. No dense landscape matrix is retained. The completed group is
published only after all subchunks and method rows validate.

## Atomic artifacts and resume

The 600 group directories are the only resumable production units.

- Missing group: eligible to build.
- Complete group with matching status, manifest, queue/source/implementation
  hashes, artifact hashes, row counts, and sizes: reusable without mutation.
- Partial temporary group: never published; clean process-owned temporary data
  may be discarded by the creating process only.
- Published partial, stale, hash-invalid, or identity-mismatched group: hard
  abort; never overwrite or auto-repair.

Each group binds source identity, transformed-view metrics, finite H0/H1
intervals, landscape summaries, complete H0/H1 pair distances, matched energy,
assembled four-method label-closed rows, internal resources, artifact manifest,
and terminal status. Publication is a same-filesystem atomic rename.

A validation-only resume must reuse all completed groups and preserve every
path, hash, byte size, and modification timestamp. The eight frozen maximum-
fold representation/configuration groups are rebuilt in a separate private
root and must match all deterministic scientific artifacts byte-for-byte.

## Resource envelope

Accepted historical measurements cover the same 75 groups per representation:

| Stage | SCT seconds/config | Integrated seconds/config |
|---|---:|---:|
| PH | 3,767.759 | 3,952.673 |
| Exact landscapes | 4,194.152 | 2,072.893 |
| Retrieval-input assembly | 1,461.448 | 1,896.369 |

The four configurations project to 69,381.174 worker-seconds, or 19.273
worker-hours, before the explicit validation/repeat reserve. MV5-V freezes:

- one heavy worker;
- configuration-stratified launch and stop decisions;
- 600 seconds and 4 GiB process-tree RSS per group;
- 8 worker-hours per configuration;
- 30 worker-hours across all four configurations, validation, and repeats;
- 4 GiB new private storage per configuration;
- 16 GiB total new private storage;
- at most 250 exact-landscape rows per deterministic subchunk.

The older conservative four-setting storage projection is 10.18 GB. MV5-U's
measured group artifacts extrapolate to about 7.22 GB before complete pair-row
expansion. The 16-GiB cap covers the conservative projection, complete pair
outputs, eight repeats, manifests, logs, and temporary group staging without
authorizing uncontrolled growth.

## Independent validation

The execution sprint must independently verify:

1. every committed and private source hash and the exact engine/runtime hashes;
2. the 600-group/54,000-view axes and all configuration-isolation invariants;
3. exact accepted heldout-training pair axes, two dimensions per biological
   pair, 2,880 deterministic subchunks, and complete group row counts;
4. all transformed point/coordinate identities and shapes;
5. all-view finite PH invariants and H0 MST agreement;
6. analytic square H1 and exact-landscape oracles;
7. direct-loop energy checks sampled prospectively across every configuration,
   representation, fold-size stratum, and seed;
8. exact-distance square identities and all-active/no-cap flags;
9. atomic manifests, complete failure reporting, resource caps, eight-group
   byte-identical clean repeat, and immutable full resume;
10. zero label access, zero outcomes, and public artifact safety.

## Abort rules

Abort the current configuration on source, queue, implementation, runtime,
shape, cell, coordinate, pair-axis, dimension, subchunk, PH, MST, landscape,
energy, artifact, repeat, resume, label-firewall, or resource mismatch. Preserve
completed immutable groups. Do not proceed to a later configuration until the
current configuration has complete resources and independent validation.

Abort the whole program if any failure indicates a shared implementation or
source defect, if the projected/observed total crosses 30 worker-hours or 16
GiB, or if labels/outcomes are accessed.

## Decision boundary

MV5-V can only freeze a prospective engine and exact queue. Even a passing
prefreeze records `full_execution_authorized=FALSE` until the engine and source
identities are committed and rebound in a separate launch-readiness commit.
The next eligible sprint may perform that final binding and a one-group real-
runner smoke check; it must not silently combine prefreeze and full execution.

Robustness outcomes, method ranking, clustering reopening, spectral promotion,
gene topology, cell/gene fusion, new data, optimization or Rust work, package-
default changes, manuscript claims, PDFs/reviewer material/private-cache
tracking, `example_run.r`, pushing, and modification of accepted MV5-D through
MV5-U artifacts remain prohibited.
