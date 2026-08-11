# MV5-AF label-closed nested192 execution specification v1

Date frozen: 2026-08-11

Accepted authorization base: `0dc460d`

## Scope

MV5-AF may bind and calculate only the 150 frozen
`nested_cells_192_pc30_euclidean_v1` groups: 15 LOSO folds, five accepted seeds and
two representations. The exact totals are 13,500 views, 70,700 directed
heldout-training pairs, 141,400 H0/H1 exact-landscape rows in 720 deterministic
subchunks, 70,700 matched-energy rows and 282,800 four-method rows.

## Scientific implementation

Reuse the accepted MV5-U/W calculation primitives. Each accepted 384-by-30
training-fit coordinate view is reduced to 192 cells by the frozen
`sha256_sample_seed_cell_nested_v1` order: cells are ordered by a SHA-256 key
over sample ID, accepted seed and cell ID, and the first 192 cells in that
deterministic order are retained. This is the exact code-level meaning of
"nested 192"; it is not the first 192 rows of the source matrix. The same order
makes the 192-cell set a strict subset of the still-closed 256-cell set.
All 30 unchanged coordinates and Euclidean geometry are retained. No PCA refit,
normalization, coordinate truncation, or outcome-informed cell selection occurs.

Typed Euclidean Vietoris-Rips PH uses dimensions H0/H1, threshold -1 and field
2. Every view must contain one essential H0 interval, and finite H0 deaths must
agree with an independent Euclidean MST oracle. H0/H1 remain separate; the
essential H0 interval is excluded; every consecutive active landscape level is
retained; squared L2 is integrated exactly over critical-pair segments. Energy
uses the same transformed views and directed pair axis. Four label-closed rows
are emitted per biological pair: H0, H1, descriptive raw H0/H1 combination,
and matched energy.

## Execution and reproducibility

Bind the committed queue, 150 coordinate hashes, source/runtime/implementation
hashes and exact counts before production. One heavy worker is allowed. Hard
caps are 600 seconds and 4 GiB process-tree RSS per group, 6 worker-hours and
4 GiB new private storage for the configuration. Groups publish atomically and
existing groups are reusable only when every identity, hash, size and status
matches. Partial, stale or invalid published groups hard-fail without repair.

Independent validation must reconstruct queue/configuration isolation,
source and deterministic nested point identities, all-view PH/MST invariants, exact-landscape
semantics, stratified direct-loop energy values, complete pair/method rows,
artifacts/resources and prohibited-operation counters without production
scientific helpers. Frozen clean repeats and a full immutable resume are
required.

## Stop boundary

MV5-AF does not rank training samples, access labels, calculate outcomes,
compare nested 192 with another setting, cluster, authorize nested 256,
select methods/subgroups, or alter accepted evidence. Gene/fusion/new-data,
optimization/Rust/default/claim work, private/PDF/reviewer tracking,
`example_run.r`, and pushing remain prohibited.

## Execution record

The prospective engine was committed at `08f8332` and its runtime/source queue
at `e9ba7d2` before production. All 150 groups completed in 4,439.486 seconds
(1.233 worker-hours), with a 53.302-second maximum group, 454,893,568-byte peak
RSS, and 1,254,700,479-byte private result tree.

Independent checks reconstructed all 13,500 deterministic nested selections,
proved every 192-cell set is contained in its frozen 256-cell superset,
recomputed all 2,578,500 H0 MST deaths, 60 H1 diagrams, 60 exact landscape
distances, and 30 energy distances, and validated all manifests and method
rows. Sixteen deterministic group artifacts and 11/11 clean validator ledgers
match byte-for-byte. All 1,650 production files retain exact paths, hashes,
sizes, and timestamps through a 150-group validation-only resume.

See `docs/audits/MV05AF_LABEL_CLOSED_NESTED192_CONFIGURATION_EXECUTION_2026-08-11.md`.
