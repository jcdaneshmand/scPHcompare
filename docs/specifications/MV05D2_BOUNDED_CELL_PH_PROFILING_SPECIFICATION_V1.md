# MV5-D2 bounded cell-PH profiling specification v1

| Field | Frozen value |
|---|---|
| Date | 2026-08-07 |
| Stage | MV5-D2 |
| Input | 75 independently validated MV5-D1 fold-seed coordinate caches |
| View | SCT `cell_topology_v1` only |
| PH | Complete Vietoris-Rips H0/H1, Euclidean, field 2 |
| Primary pilot | 30 diagrams |
| Exact repeats | 5 diagrams, one per seed |
| Outcome-label state | Closed |
| Required stop | Before all 6,750 diagrams, landscapes, distances, clustering, integration, gene views, or outcomes |

## Purpose

MV5-D2 measures and validates persistence-diagram construction on the corrected
cells-as-observations coordinates before a full cell-primary PH stage is
authorized. It resolves the unmeasured PH term left by MV5-D1. It does not
construct a persistence landscape or make a biological comparison.

## Frozen input and pilot selection

Only the 6,750 typed 384-cell by 30-PC views inside the 75 validated MV5-D1
caches are eligible. The public pilot manifest binds every selected view to the
fold-cache key and hash, typed-view key and payload hash, fit scope, seed,
sample, and execution role.

The deterministic 30-job pilot contains:

- one held-out view and one training control in every one of the 15 LOSO folds;
- exactly six jobs for each of the five frozen seeds;
- mapped held-out views from all nine folds in which any seed has a missing
  training-selected feature;
- unmapped held-out controls from the other six folds; and
- five repeat jobs, one per seed, spanning held-out and training roles.

Mapped held-out selection uses the greatest observed training-schema mapping
burden within the seed assigned by a deterministic capacity-constrained
matching. Training controls use a deterministic fold-index rule. Neither
tissue, assay approach, nor any outcome may influence selection.

## PH and diagram contract

Every job must use the typed corrected entry point with:

- observations: 384 cells;
- coordinates: 30 shared, training-fitted PCs;
- metric: Euclidean distance on the frozen coordinates;
- filtration: complete Vietoris-Rips with `threshold = -1`;
- maximum homology dimension: 1; and
- coefficient field: 2.

Each cached result contains H0 and H1 in one typed diagram plus immutable
provenance. Exactly one essential H0 interval is retained with infinite death,
so H0 has 384 records: 383 finite merge intervals and one essential class. All
H1 records must be finite and have positive persistence. Invalid and
zero-persistence intervals are prohibited.

## Correctness requirements

1. The finite H0 death multiset must equal the Euclidean minimum-spanning-tree
   edge multiset within `max(1e-7, 1e-7 * maximum_edge)`.
2. The five independent repeats must have identical diagram hashes, record
   keys, R objects, serialized-file hashes, and bytes.
3. For each repeat view, the first 32 ordered cells must be checked with both
   Ripserr and `TDA::ripsDiag`/GUDHI. Finite H0 and H1 intervals must agree
   within `1e-6` after explicitly removing GUDHI's maximum-scale-capped
   representation of the essential H0 class.
4. Existing stale, unreadable, or identity-mismatched results may not be
   overwritten.

The MST oracle checks the full 384-cell H0 result. The 32-cell cross-engine
check is analytical corroboration only and is not a scientific result.

## Cache and execution rules

The identity binds the pilot-manifest hash, source fold-cache hash and key,
typed-view hash and key, PH parameters, implementation hash, runtime/package
versions, and thread environment. The immutable record excludes elapsed time
and wall-clock timestamps so exact repeats can be compared byte-for-byte.
Resource metrics are stored separately.

| Guard | Frozen value |
|---|---:|
| Maximum heavy workers used by the pilot | 1 |
| Per-job elapsed cap | 600 seconds |
| Per-job process-tree RSS cap | 4 GiB |
| Pilot admission cap | 3,600 seconds |
| Primary pilot storage cap | 100 MB |

The queue writes a completed public metric after each atomic private cache.
Safe resume requires a completed matching metric; a cache without that metric
is treated as ambiguous and is not silently reused for resource projection.

## Landscape boundary

MV5-D2 produces diagrams, not landscapes. The dissertation-aligned definition
remains unchanged:

- H0 and H1 are constructed and compared separately;
- only finite positive-persistence intervals enter a landscape, so the one
  essential H0 interval is excluded from the landscape input;
- every active consecutive landscape level is retained;
- primary L2 integration is exact or error-controlled on dimension-specific
  support; and
- any unweighted H0/H1 combination is secondary and must retain its component
  distances and H1 energy contribution.

No fixed landscape level cap or fixed uniform-grid count is introduced here.

## Completion and decision gate

MV5-D2 completes only if all 30 jobs pass, all full-view H0 MST checks pass,
all five repeats are byte-identical, all ten reduced H0/H1 cross-engine checks
pass, no resource guard is breached, and an honest 6,750-view projection is
recorded. A projection below the planning cap permits owner review of a later
full cell-PH cache stage; it does not launch that stage automatically.

