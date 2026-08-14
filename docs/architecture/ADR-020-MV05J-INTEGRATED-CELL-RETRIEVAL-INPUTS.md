# ADR-020: Freeze label-closed integrated cell retrieval inputs

- Status: Accepted and executed
- Date: 2026-08-09
- Decision scope: MV5-J only

## Context

MV5-G produced 75 validated integration groups, MV5-H converted every accepted
integrated coordinate matrix into the corrected typed cell-topology view, and
MV5-I completed exact all-active-level H0 and H1 landscape distances for all
35,350 frozen held-out-query-to-training-reference pairs. Those components are
not yet an immutable prediction input: the matched integrated-coordinate energy
baseline, shared pseudobulk context, deterministic rankings, completion ledger,
and full replay evidence remain to be constructed.

The accepted MV5-D5 SCT retrieval-input contract supplies the appropriate
statistical shape, but its SCT energy distances and rankings cannot be reused
for an integrated view. Conversely, its seed-specific mean-profile bundles and
training-derived pseudobulk definition are representation-independent and can
be reused read-only.

## Decision

1. Create a separate integrated retrieval-input artifact family; never mix SCT
   and integrated rankings or identities.
2. Keep raw integrated H0 and H1 as separate confirmatory topology methods.
3. Retain `sqrt(H0^2 + H1^2)` only as descriptive secondary output.
4. Do not fit a topology-component scale from held-out-query rows and do not
   expand to within-training topology in this sprint.
5. Recompute energy distance from the exact MV5-G 384-cell by 30-coordinate
   matrices used by integrated topology.
6. Reuse the accepted MV5-D5 mean-profile bundles and MV5-D1 training-selected,
   training-standardized 500-gene pseudobulk definition as context.
7. Bind every group to MV5-D1, MV5-D5 mean-profile, MV5-G coordinate, MV5-I
   component, and implementation hashes; rank by distance then canonical sample
   ID; publish atomically; and refuse stale reuse.
8. Keep labels and all retrieval evaluation, clustering, gene topology, fusion,
   new-data, and outcome work closed.
9. Admit representative and maximum groups before the monitored 75-group run,
   then require independent all-row validation, repeatability, immutable resume,
   resource evidence, and a clean package check.

## Consequences

MV5-J can produce a fair integrated counterpart to the accepted SCT prediction
inputs while preserving the same frozen folds, cells, seeds, query/reference
scope, and pseudobulk context. Integrated energy is recalculated rather than
borrowed, so topology and its matched distributional baseline share the exact
coordinate geometry.

The sprint cannot reveal whether integrated topology performs better, select a
method, or revise biological conclusions. A later sprint may open labels only
after MV5-J explicitly authorizes it and the later specification locks endpoints
and verifies immutable public hashes.

## Evidence

All 75 groups completed 35,350 biological pairs across five methods, producing
176,750 unique ranked rows and 375 successful method-group records. Independent
validation reproduced every H0/H1 value, checked 450 direct baseline pairs, and
found a maximum numerical difference of `9.805490e-13`. Representative and
maximum complete-group repeats, all 150 resume-tracked private files, and all
six public artifacts are byte-stable. Production used 1,896.369 aggregate
worker-seconds, 301,789,184 bytes peak process-tree RSS, and 209,877,213 private
result bytes. A later separately prediction-locked integrated retrieval
evaluation is authorized; every other stop boundary remains closed.
