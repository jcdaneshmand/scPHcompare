# MV5-AB label-closed cosine execution specification v1

Date frozen: 2026-08-11

Accepted authorization base: `89be63d`

## Scope

MV5-AB may bind and calculate only the 150 frozen
`cells384_pc30_cosine_chord_v1` groups: 15 LOSO folds, five accepted seeds and
two representations. The exact totals are 13,500 views, 70,700 directed
heldout-training pairs, 141,400 H0/H1 exact-landscape rows in 720 deterministic
subchunks, 70,700 matched-energy rows and 282,800 four-method rows.

## Scientific implementation

Reuse the accepted MV5-U/W calculation primitives. Each accepted 384-by-30
training-fit coordinate view is row-L2-normalized after rejecting every zero or
nonfinite norm. Euclidean distance on those unit vectors is chord distance.
No PCA refit, cell reselection or coordinate truncation occurs.

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
caps are 600 seconds and 4 GiB process-tree RSS per group, 8 worker-hours and
4 GiB new private storage for the configuration. Groups publish atomically and
existing groups are reusable only when every identity, hash, size and status
matches. Partial, stale or invalid published groups hard-fail without repair.

Independent validation must reconstruct queue/configuration isolation,
normalization and point identities, all-view PH/MST invariants, exact-landscape
semantics, stratified direct-loop energy values, complete pair/method rows,
artifacts/resources and prohibited-operation counters without production
scientific helpers. Frozen clean repeats and a full immutable resume are
required.

## Stop boundary

MV5-AB does not rank training samples, access labels, calculate outcomes,
compare cosine with another setting, cluster, authorize nested-cell settings,
select methods/subgroups, or alter accepted evidence. Gene/fusion/new-data,
optimization/Rust/default/claim work, private/PDF/reviewer tracking,
`example_run.r`, and pushing remain prohibited.

## Execution record

The engine was frozen at `20ec50e` and its exact queue/runtime/source binding at
`6b37eac`. All 150 groups completed in 2.608 worker-hours within every cap.
Independent validation passed 15 artifact categories, all 13,500 raw-coordinate
normalization/MST oracles, and 30 stratified direct energy oracles. Sixteen
repeat artifacts match byte-for-byte; all 1,650 production files retain exact
paths, hashes, sizes, and timestamps through a 150-group validation-only resume.
See
`docs/audits/MV05AB_LABEL_CLOSED_COSINE_CONFIGURATION_EXECUTION_2026-08-11.md`.
