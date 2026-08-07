# ADR-005: MV5-A/MV5-B execution scaffold

## Status

Accepted on 2026-08-06 for label-closed staged execution. This decision does
not authorize biological endpoint calculation, method ranking, fusion,
full-cohort execution, or a public-default change.

## Context

ADR-004 froze a sample-level leave-one-study-out benchmark but left two
implementation questions open: whether folds and matched baselines could be
made immutable before labels were opened, and whether the locked Seurat stack
could map held-out cells into a reference-trained coordinate system without
jointly refitting on the held-out study.

The legacy `perform_integration()` function remains transductive. The locked
environment contains Seurat 5.3.0, whose reference-mapping API supports PCA
projection from a reference to a query and correction of the projected query
embedding.

## Decision

1. Use `mv05_loso_manifest_v1` as the executable split artifact. It contains
   study identifiers needed for fold assignment but contains no tissue or
   approach columns and declares `outcome_label_state=closed`.
2. Bind every fold to a unique training `fit_scope_id`, all 124 large-cohort
   sample IDs, one held-out study, disjoint reference/query roles, a source
   identity digest, and an immutable cache key.
3. Implement three matched sample-distance baselines:
   - square-root empirical V-statistic energy divergence over the same
     shared-PCA cell coordinates;
   - Euclidean distance between training-scaled gene means; and
   - RMS Frobenius distance between gene-correlation matrices.
4. Require common fit scope, representation, subsample seed, axes, identifiers,
   input cache keys, exact zero diagonals, symmetry, matrix hashes, and bundle
   hashes for every baseline result.
5. Use a separate inductive route based on reference-only feature selection and
   PCA, `FindTransferAnchors(reduction="pcaproject")`, and
   `IntegrateEmbeddings()`. Do not call `TransferData()` or transfer any label.
6. A mapping failure returns `held_out_mapping_unavailable` with its reason; it
   never falls back to the transductive integration object.

## Evidence

The generated LOSO artifact contains 2,232 rows: 124 samples in each of 18
folds. Every sample is held out exactly once, training/query IDs are disjoint,
and outcome labels are absent.

Analytical fixtures pass for all three baseline formulas. A synthetic
scientific-shape fixture with 500 genes, 384 cells per sample, 30 PCs, and three
samples completed all baseline matrices with exact zero diagonals and symmetry.

The label-closed Seurat fixture used 80 reference and 40 held-out query cells,
200 reference-selected features, and 10 reference PCs. Two independent mapping
runs each found 111 anchors and produced the identical query-embedding SHA-256
hash while preserving the reference identity hash.

## Consequences

Inductive mapping is technically available rather than categorically blocked.
This is synthetic feasibility, not evidence that every existing-data fold will
map or that integration improves biology. The next gate is MV5-C: process one
eligible existing-data tissue across all five frozen seeds and methods, retain
all failures, hash matrices/predictions before outcome labels are opened, and
project full-run cost. G-MV5 remains open.
