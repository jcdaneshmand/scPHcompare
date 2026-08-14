# MV5-AJ label-closed nested-256 execution specification v1

Date frozen: 2026-08-11

Accepted authorization base: `23640a48b31b569dd52878c6f7358e6082b297ba`

## Scope

MV5-AJ may bind and calculate only the 150 frozen
`nested_cells_256_pc30_euclidean_v1` groups: 15 LOSO folds, five accepted seeds,
and two separate representations. The exact totals are 13,500 views, 70,700
directed held-out-query/training-reference pairs, 141,400 separate H0/H1 exact-
landscape rows in 720 deterministic subchunks, 70,700 matched-energy rows, and
282,800 four-method rows.

The calculation consumes the exact fourth-position queue authorized by MV5-AI.
Rankings, tissue labels, outcomes, clustering, within-training pairs, other
configurations, and method or subgroup selection are outside this sprint.

## Scientific implementation

The accepted MV5-U/W primitives and MV5-AF streamed execution architecture are
reused. Each accepted 384-by-30 training-fit coordinate view is reduced by the
frozen `sha256_sample_seed_cell_nested_v1` order. Cells are ordered by a SHA-256
key over sample ID, accepted seed, and cell ID; the first 256 are retained. This
is not source-matrix row order and no redraw is allowed.

For every view, the first 192 selected IDs must exactly equal the accepted
nested-192 selection, and the 256 IDs must be a strict subset of the same
384-cell realization. Production and independent validators both reconstruct
these identities. All 30 coordinates and Euclidean geometry remain unchanged;
there is no PCA refit, normalization change, coordinate truncation, or outcome-
informed selection.

Typed Euclidean Vietoris-Rips PH uses H0 and H1, threshold -1, and field 2.
Every view must contain one essential H0 interval; finite H0 deaths must match
an independently reconstructed Euclidean MST. H0 and H1 remain separate, the
essential H0 interval is excluded, every consecutive active landscape level is
retained, and squared L2 distance is integrated exactly over critical-pair
segments. No universal landscape-level cap or uniform grid is introduced.

Energy is recomputed on the same transformed views and directed pair axis. Four
label-closed rows are emitted per pair: H0 landscape, H1 landscape, descriptive
raw H0/H1 Euclidean combination, and matched energy.

## Execution and reproducibility

Before production, the prospective engine must bind its Git identity, the
complete MV5-AI authorization and repeat evidence, all 150 coordinate hashes,
the private coordinate inventory, source/runtime/implementation hashes, and
exact queue counts. One heavy worker is allowed. Hard caps are 600 seconds and
4 GiB process-tree RSS per group, six worker-hours and 4 GiB new private storage
for the configuration.

Groups publish atomically. Existing groups are reusable only when every source,
queue, implementation, runtime, artifact hash, size, and completed status
matches. Partial, stale, or invalid publications hard-fail without repair.

Independent validation must reconstruct configuration isolation, source and
nested point identities, all-view H0 MST invariants, selected H1 diagrams,
exact-landscape semantics, stratified direct-loop energy values, complete pair
and method rows, manifests, resources, and prohibited-operation counters without
production scientific helpers. Prespecified minimum/maximum groups must repeat
exactly. A complete clean validation and full immutable resume, including paths,
hashes, sizes, and timestamps, are required.

## Interpretation and stop boundary

Successful calculation establishes only that the complete nested-256 distance
panel was generated reproducibly under the frozen definition. It does not say
whether retrieval changed, and it does not authorize robustness, equivalence,
invariance, superiority, or default claims.

MV5-AJ does not rank samples, read `Tissue` or `Approach`, calculate outcomes,
compare nested 256 with another setting, cluster, or open the next stage. Gene
topology, cell/gene fusion, new data, Rust optimization, defaults, manuscript
claims, private/PDF/reviewer tracking, `example_run.r`, pushing, and any result-
driven repair remain prohibited.
