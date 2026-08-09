# MV5-D3 full cell-PH production specification v1

| Field | Frozen value |
|---|---|
| Date | 2026-08-07 |
| Stage | MV5-D3 |
| Input | 75 independently validated MV5-D1 fold-seed coordinate caches |
| Jobs | 6,750 typed SCT `cell_topology_v1` views |
| PH | Complete Vietoris-Rips H0/H1, Euclidean, field 2 |
| Execution unit | One 90-view fold-seed group |
| Maximum heavy workers | 2 |
| Outcome-label state | Closed |
| Required stop | Before landscapes, distances, clustering, integration, gene views, or outcomes |

## Purpose

MV5-D3 creates the complete immutable persistence-diagram cache needed by the
corrected cells-as-observations SCT primary analysis. It promotes the validated
MV5-D2 diagram contract to all frozen cell views without constructing or
approximating a persistence landscape.

## Frozen scope and identity

The public manifest contains exactly 75 fold-seed groups and 90 views per
group: 15 leave-one-study-out folds by five frozen seeds. It binds each job to
its MV5-D1 fold-cache key and SHA-256, typed-view key and payload SHA-256,
sample role, mapping stratum, representation, orientation, shape, metric, and
PH parameters. The manifest contains no tissue, assay-approach, or outcome
field.

Every immutable record identity binds:

- the full production-manifest SHA-256;
- source fold-cache and typed-view identities;
- complete Euclidean Vietoris-Rips H0/H1 parameters;
- implementation-file SHA-256;
- runtime and package versions; and
- the single-thread numerical environment.

An existing result is reusable only when its matching checkpoint, record
identity, and hashes validate. A stale, unreadable, ambiguous, or mismatched
result may not be silently overwritten.

## Diagram contract

Each 384-cell by 30-training-PC view produces one typed record containing
separate H0 and H1 diagrams. The contract requires:

1. exactly 384 H0 intervals, comprising 383 finite component merges and one
   explicit essential interval with infinite death;
2. finite H0 deaths equal to the full-view Euclidean minimum-spanning-tree edge
   multiset within `max(1e-7, 1e-7 * maximum_edge)`;
3. every H1 interval finite and of positive persistence; and
4. no invalid or zero-persistence interval.

## Execution and resource guards

Work is grouped by fold cache so each source cache is loaded once for its 90
independently saved results. Each result is serialized atomically. A completed
group receives a private per-view checkpoint and one public-safe group metric.

| Guard | Frozen value |
|---|---:|
| Per-group elapsed cap | 900 seconds |
| Per-group process-tree RSS cap | 4 GiB |
| Whole-stage elapsed cap | 43,200 seconds |
| Maximum concurrent heavy groups | 2 |

The first two groups form the admission pilot. Production may continue only if
both complete all 90 views without a time, RSS, storage, identity, or diagram
failure. Resume validation must prove that all 90 records and the checkpoint of
the first group remain byte-unchanged.

## Independent validation

Completion requires a fresh validator to:

- validate all 6,750 manifest rows and public/private execution metrics;
- verify every result-file hash, cache identity, and static diagram invariant;
- confirm stored full-view H0 MST evidence for all 6,750 records;
- reload the corresponding MV5-D1 coordinates and independently recompute one
  full 384-cell H0 MST oracle per group;
- compare all 90 records in one independently executed production-group repeat
  for diagram hash, R-object identity, serialized-file hash, and bytes;
- confirm 75 groups, 6,750 completed views, zero failures, and zero downstream
  counters; and
- replace the prior PH projection with measured production time and storage.

## Landscape boundary

MV5-D3 stores diagrams only. The project-owner-approved dissertation-aligned
landscape definition is unchanged:

- calculate and report H0 and H1 separately;
- exclude the one essential H0 interval and every non-positive-persistence
  interval from landscape input;
- retain all active consecutive landscape levels, with no universal level cap;
- use exact or error-controlled L2 integration on dimension-specific support,
  with no universal fixed uniform-grid count; and
- treat an unweighted H0/H1 combination as secondary while retaining both
  component distances and the H1 squared-distance contribution.

## Completion gate

MV5-D3 is complete only after the full production cache, independent all-record
validation, 75 recomputed MST checks, 90-view exact repeat, resume proof,
measured resource reconciliation, full source-loaded test suite, and isolated
source-package check pass. Completion authorizes specification of a separate
landscape-distance stage; it does not itself authorize landscapes, clustering,
labels, gene topology, integration, or biological claims.
