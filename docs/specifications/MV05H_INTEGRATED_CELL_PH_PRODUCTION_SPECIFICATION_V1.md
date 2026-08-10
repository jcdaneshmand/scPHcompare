# MV5-H label-closed integrated cell-PH production specification v1

## Document control

| Field | Frozen value |
|---|---|
| Contract ID | `mv05h_integrated_cell_ph_production_v1` |
| Date | 2026-08-09 |
| Source revision | `dcd975e` |
| Input | 75 independently validated MV5-G coordinate groups |
| Jobs | 6,750 typed integrated `cell_topology_v1` views |
| PH | Complete Vietoris–Rips H0/H1, Euclidean, field 2 |
| Execution unit | One 90-view fold-seed group |
| Maximum heavy workers | 2 |
| Outcome-label state | Closed |
| Required stop | Before landscapes, distances, retrieval, clustering, gene views, fusion, or outcomes |

## 1. Purpose and stop boundary

MV5-H creates the complete immutable persistence-diagram cache for the accepted
MV5-G inductive integrated cell coordinates. It changes only the coordinate
representation relative to MV5-D3; corrected cells-as-observations orientation,
complete Euclidean H0/H1, and all diagram correctness rules remain unchanged.

This sprint permits source validation, typed coordinate-view construction,
H0/H1 PH, immutable private caching, resource measurement, independent
validation, full-group deterministic repeat, resume proof, and post-PH
projection. It prohibits landscapes, landscape distances, baselines, retrieval,
clustering, gene topology, fusion, new data, biological outcomes, and claims.

## 2. Frozen input and identity

The public manifest contains exactly 75 fold-seed groups and 90 views per group:
15 leave-one-study-out folds crossed with seeds `20260805:20260809`. Every view
binds to the accepted MV5-G private file SHA-256, group cache key, payload hash,
coordinate-set hash, sample role, sample ID, coordinate payload hash, and typed
view cache key.

Each coordinate matrix must retain its exact ordered 384 cell IDs and 30 columns
`PC1` through `PC30`. It becomes a scientific-eligible integrated
`cell_topology_v1` view with rows as points and Euclidean geometry. It must not
be labeled `sct_whole`, and no expression matrix is reconstructed or invented.

Each PH identity additionally binds the complete manifest hash,
implementation-file hashes, runtime/package versions, single-thread numerical
environment, and complete VR parameters. Existing records are reusable only
when checkpoint, record identity, and file hash all validate.

## 3. Diagram contract

Each integrated view produces separate H0 and H1 diagrams and must contain:

1. exactly 384 H0 intervals: 383 finite merges and one explicit essential
   interval with infinite death;
2. finite H0 deaths equal to the full-view Euclidean minimum-spanning-tree edge
   multiset within `max(1e-7, 1e-7 * maximum_edge)`;
3. only finite, positive-persistence H1 intervals; and
4. zero invalid or zero-persistence intervals.

The essential H0 class remains stored for provenance and is excluded only in a
later landscape sprint.

## 4. Execution and resource guards

Groups load one MV5-G coordinate bundle and atomically serialize 90 independent
PH records. The first two groups form an admission run in a separate private
root. Full production may begin only if both pass every source, identity,
diagram, time, memory, storage, and scope check.

| Guard | Frozen value |
|---|---:|
| Per-group elapsed cap | 900 seconds |
| Per-group process-tree RSS cap | 4 GiB |
| Whole-stage worker-time cap | 43,200 seconds |
| PH result storage cap | 1 GB |
| Maximum concurrent heavy groups | 2 |

## 5. Required validation

An independent validator must verify all 6,750 manifest rows, source identities,
typed coordinate identities, result-file hashes, static diagram invariants, and
stored full-view MST evidence. It must reload MV5-G coordinates and independently
recompute one complete 384-cell H0 MST oracle per group.

One separately executed 90-view group must be object- and byte-identical to its
production counterpart. A completed-queue rerun must rebuild zero records and
preserve all record/checkpoint hashes. Every landscape, distance, retrieval,
clustering, gene-view, fusion, new-data, and outcome counter must remain zero.

## 6. Post-PH decision

After full validation, measured coordinate-plus-PH resources replace the prior
PH projection. A later exact all-active-level landscape-distance sprint is
authorized only if the reserved total through retrieval remains below 21.6
worker-hours, 8 GiB RSS, 1,800 seconds per upstream mapping group, and 10 GB
combined private storage. Authorization does not execute landscapes or open
outcomes.
