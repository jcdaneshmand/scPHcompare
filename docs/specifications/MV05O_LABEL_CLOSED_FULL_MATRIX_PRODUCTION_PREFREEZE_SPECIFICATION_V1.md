# MV5-O label-closed complete-training-matrix production prefreeze specification v1

## Document control

| Field | Value |
|---|---|
| Status | Frozen prospective execution contract; production not executed |
| Date | 2026-08-10 |
| Base revision | `36305383a706b1e3699cb90f086e9b4bfe0b6e8a` |
| Upstream gate | MV5-N complete and canonically verified |
| Outcome-label state | Closed |
| Required stop | Before complete-matrix execution, production clustering, or label opening |

## 1. Purpose and boundary

MV5-O converts the cap-passing MV5-N resource decision into an executable,
immutable production plan. It freezes every source, implementation, group,
chunk, baseline unit, cap, output/status rule, validation sample, repeat,
resume, and abort condition before full work begins.

This prefreeze does not calculate a complete training matrix, fit a production
cluster, assign a held-out sample, open tissue/study/approach labels, or
compute ARI/NMI. It changes no accepted MV5-D through MV5-N artifact.

## 2. Scientific object remains unchanged

The only authorized topological distances remain:

- corrected cells-as-observations `cell_topology_v1` diagrams;
- finite positive-persistence intervals, with the one essential H0 interval
  excluded;
- H0 and H1 calculated and stored separately;
- every consecutive active persistence-landscape level;
- exact piecewise-linear critical-pair L2 integration;
- no universal level cap and no uniform evaluation grid.

The raw `sqrt(H0^2 + H1^2)` combination remains descriptive and is not a
production calculation in this stage. Energy is representation-specific;
training-standardized pseudobulk is computed once per fold/seed and reused
across representations.

## 3. Frozen sources and implementations

The source-freeze root is
`541e7d3aa8acce5d512bbb4819c034735eef47387e91a63abccfa259f53d6de1`.
It binds 18 accepted manifests, resource records, identity generators,
formula engines, and production implementations by path, size, and SHA-256.

| Implementation | SHA-256 |
|---|---|
| Prefreeze/queue contract | `1b7b7d9553ca3028dea886e2de1ec0265d38be60eefd3102b521c6d597bf5a7f` |
| Per-group request/interval stager | `608a50a94fe21bac89319f79afb7369a91f85a3b7fc0f226e9cd02ef88c059b0` |
| Cache-preserving exact landscape runner | `ffee3d0884f1bbf84194320a3263a943945425e46ef45799160fb90201e4da1d` |
| Baseline-group runner | `af10c597c9d49051f91f55d90139f369d7cc55db99816a3ec7ce440be8f19e1f` |

Production must refuse any mismatched base revision, source-freeze root,
implementation hash, group root, chunk root, diagram hash, result-file hash,
or input-cache identity. A changed file requires a new prefreeze; it cannot be
accepted by editing a status record.

## 4. Exact production scope

| Work family | Groups/units | Pair or component rows |
|---|---:|---:|
| Fold-seed-representation groups | 150 | — |
| Exact landscape chunks | 4,340 | 1,050,700 H0/H1 rows |
| Representation-specific energy groups | 150 | 525,350 pairs |
| Shared pseudobulk groups | 75 | 262,675 pairs |
| **All resumable units** | **4,565** | **1,838,725 values** |

There are 15 LOSO folds, five seeds, and two representations. Each group
regenerates the exact MV5-N canonical unordered training-pair manifest and
must reproduce its published group identity root before staging intervals.
Each landscape chunk must reproduce its published chunk root. No held-out to
held-out or query-to-query work is permitted.

Large pair-row manifests remain private execution inputs. Git stores the
complete group/chunk roots, first/last pair identities, counts, queues, and
source hashes needed to regenerate and verify them without publishing
million-row duplicates.

## 5. Execution order and concurrency

The 150 groups have one deterministic order: representation, fold ID, then
seed, all under radix ordering. Within a group, landscape chunks are ordered
by H0/H1 and first pair-request ID. Baseline units follow their group.

At most two group workers may run concurrently. Each landscape worker stages
one group, constructs each sample/dimension landscape once, and processes that
group's chunks sequentially so MV5-N's measured caching assumptions remain
valid. Starting a fresh process per chunk is prohibited because it would
invalidate the resource projection.

Pseudobulk is calculated only for the SCT-anchored fold/seed group and reused
for the corresponding integrated stratum after exact pair/sample-axis
identity validation.

## 6. Resource and storage guards

| Guard | Frozen value |
|---|---:|
| Aggregate worker-hour cap | 21.6 h |
| Projected work including 10% reserve | 16.117047 h |
| Per-group/per-unit elapsed cap | 900 s |
| Per-process-tree RSS cap | 4 GiB |
| Additional private-storage cap | 10 GiB |
| Landscape output projection | 618,845,884 bytes |
| Conservative baseline output projection | 403,468,800 bytes (512 bytes/value) |
| Status/validation/slow-storage reserve | 255,578,671 bytes (25%) |
| Total private-storage projection | 1,277,893,355 bytes |

Resource accounting must be updated after every completed unit. The queue may
launch a new unit only while observed plus projected remaining work stays
inside every cap. No cap may be met by dropping H1, landscape levels, samples,
folds, seeds, or validation.

## 7. Atomic artifacts and resume

Every unit has one deterministic output path and one deterministic status
path. Output is written to a same-directory partial path, flushed and closed,
then atomically renamed. Status is written only after output SHA-256 and byte
size are known.

Each distance row retains contract/engine, group, chunk or baseline unit,
fold, seed, representation, method or homology dimension, canonical pair ID,
sample IDs, distance, source-freeze root, implementation hash, closed-label
state, and zero clustering/outcome counters. Landscape rows additionally
retain squared distance, active-level/segment/critical-point counts, exact
status, all-active-level status, and no-cap status.

Status binds request file and request-subset hashes, implementation and source
roots, output hash/size, row count, elapsed/operation time, peak RSS, and
completion state. A rerun may reuse a unit only after all identities and hashes
validate. One-sided output/status, stale hashes, or a partial artifact is an
abort—not an invitation to overwrite.

## 8. Validation and deterministic evidence

All 4,565 units receive schema, count, identity, numerical-finiteness,
nonnegativity, source-root, label-firewall, and completion validation.

The frozen sampled plan additionally requires:

- 12 independent exact R landscape oracles: minimum, representative, and
  maximum training sizes × two representations × H0/H1;
- a clean byte-for-byte repeat of every landscape and energy output in the
  maximum group for each representation (33 units per representation);
- pseudobulk formula and cross-representation identity checks;
- energy V-statistic formula checks;
- a snapshot/rerun/snapshot comparison across all 4,565 units showing
  unchanged hashes, sizes, and timestamps with zero rebuilds; and
- complete accounting reconciliation against 150 groups, 1,838,725 values,
  16.117047 projected worker-hours, and the storage projection.

The frozen landscape runner already passes a 32-row corrected real-diagram
fixture: exact agreement with accepted MV5-N output, one independent R oracle,
byte-identical clean repeat, and immutable resume.

## 9. Abort and recovery

The ten public abort rules cover source/implementation drift, identity-root
drift, missing/duplicate/query-query requests, label leakage or nonzero
downstream counters, time/RSS/work/storage caps, partial or stale artifacts,
and any oracle/formula/repeat/resume failure.

No automatic retry is authorized. Cap or correctness failure stops launches,
preserves already validated immutable artifacts, quarantines partial state,
and requires a separately audited corrective sprint. Recovery may resume only
after the cause and resulting source/implementation hashes are documented.

## 10. Label firewall and future estimands

Production inputs and outputs accept no tissue, study outcome, approach,
class, label, biological-label, technical-label, or downstream-result column.
Study IDs remain split identifiers only. Production clustering is not part of
this prefreeze or the matrix run.

After matrix and prediction locks, a separate label-open sprint may evaluate:

1. descriptive alignment of frozen training partitions; and
2. inductive generalization of immutable held-out assignments.

Those estimands must remain separate and neither may retune `k`, method,
representation, component, fold, or tissue. Oracle class-count `k` remains
prohibited except as an explicitly historical sensitivity in a separately
approved analysis.

## 11. Prefreeze acceptance and stop

The prefreeze passes only when all source hashes, queue counts, unit identities,
caps, validation samples, abort rules, runner fixtures, focused/full tests,
package check, public artifact hashes, and deterministic repeats validate and
the contract is committed locally.

It must then stop. The subsequent production goal may execute only these
label-closed matrices. It may not fit clusters, open labels, report ARI/NMI,
select winners, run robustness, add gene topology/fusion/new data, optimize or
rewrite in Rust, change package defaults, publish private artifacts, or push.
