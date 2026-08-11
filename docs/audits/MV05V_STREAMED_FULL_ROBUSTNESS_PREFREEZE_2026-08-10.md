# MV5-V streamed full-robustness production prefreeze

Date: 2026-08-10

Status: complete

Accepted MV5-U base: `2c844f6`

Prospective gate commit: `9d14c16`

Calculation-primitives SHA-256: `a62050b36ee6c618cf171280043df47c9504f342b151c2f9066a8ba94530bb16`

Labels opened: no

Outcomes computed: no

Full execution authorized: no

## Decision

MV5-V passes as an exact scope, streaming, resource, validation, and stop-rule
prefreeze. It does not pass launch readiness: the calculation primitives are
bound, but a dedicated full-group runner/monitor and its runtime identity have
not yet been implemented or smoke-tested. The generated queue therefore keeps
all 600 `execution_authorized` and `execution_completed` fields false.

The next action is a separate MV5-W launch-readiness sprint: implement and bind
the atomic group orchestrator, independently test its full-pair semantics, and
run one real label-closed group smoke. No configuration-wide execution is
permitted until that gate passes.

## Scientific scope frozen

The eventual estimand is sensitivity of the already frozen held-out-query
retrieval effects. MV5-V does not reopen clustering. The accepted 384-cell,
30-PC Euclidean results remain the reference, and only these four
one-factor-at-a-time configurations enter the prospective queue:

- nested 192 cells/30 PCs/Euclidean;
- nested 256 cells/30 PCs/Euclidean;
- 384 cells/first 20 accepted PCs/Euclidean;
- 384 cells/30 PCs/cosine chord after row-unit normalization.

No new seed, panel, coordinate fit, PC count, metric, representation,
integration method, clustering method, or factorial interaction is added.

## Exact scope evidence

The builder reads and validates both accepted 70,700-row H0/H1 directed pair
manifests, summarizes their exact label-closed biological pair axes, joins them
one-to-one to all 150 frozen coordinate sources, and crosses only the four
accepted configurations.

| Object | Frozen count |
|---|---:|
| Source identities | 176 |
| Private coordinate identities | 150 |
| Base fold/seed/representation scopes | 150 |
| Robustness group units | 600 |
| Sample views | 54,000 |
| Heldout-training biological pairs | 282,800 |
| H0/H1 landscape rows | 565,600 |
| Landscape subchunks, at most 250 rows | 2,880 |
| Matched energy rows | 282,800 |
| Four-method assembled rows | 1,131,200 |
| Maximum-fold clean-repeat groups | 8 |

The group queue spans exactly 15 folds, five seeds, two representations, and
four configurations. Every group contains 90 views. It carries no tissue,
approach, cluster, or outcome field.

## Landscape and streaming definition

The corrected dissertation-aligned definition is unchanged: H0 and H1 are
separate, essential H0 is excluded, all consecutive active levels are retained,
there is no fixed level cap or uniform grid, and L2 integration is exact over
linear critical-pair segments.

The planned runner holds at most one 90-view group in memory. It processes
pair requests in deterministic subchunks of at most 250 rows and never saves a
dense landscape matrix. A group is published atomically only after every
subchunk, energy row, four-method row, file hash, and row count validates.

## Resource model

The projection is based on the actual accepted 75-group production stages for
each representation, not on admission startup timing alone.

| Stage | SCT seconds/config | Integrated seconds/config |
|---|---:|---:|
| PH | 3,767.759 | 3,952.673 |
| Exact landscapes | 4,194.152 | 2,072.893 |
| Retrieval-input assembly | 1,461.448 | 1,896.369 |

All four configurations project to 69,381.174 worker-seconds or 19.273 worker-
hours before the explicit reserve. The frozen envelope is one heavy worker,
600 seconds/4 GiB per group, 8 worker-hours and 4 GiB new private storage per
configuration, and 30 worker-hours/16 GiB for the full program including
validation and eight repeat groups.

The 16-GiB storage cap covers both MV5-T's conservative 10.18-GB projection
and MV5-U's measured approximately 7.22-GB group-artifact extrapolation before
complete pair expansion. Crossing a configuration or program cap stops future
launches while preserving completed immutable groups.

## Atomicity, repeat, and resume

One fold/seed/representation/configuration group is the publication and resume
unit. Missing groups are buildable; complete hash-valid groups are reusable;
published partial, stale, or hash-invalid groups are hard failures and are
never overwritten automatically. Process-owned temporary groups are never
published before validation.

The eight maximum-training groups selected by configuration and representation
must reproduce all deterministic scientific artifacts byte-for-byte in a
separate private root. Full validation-only resume must reuse every completed
group without changing any path, hash, size, or modification timestamp.

## Independent validation frozen

Fifteen required categories cover source/runtime hashes, the complete group and
view axes, configuration isolation, accepted pair-axis rebinding, subchunk
coverage, all-view H0 MST checks, analytic square H1 and exact-landscape
oracles, stratified direct energy checks, atomic manifests/failures, resources,
eight-group repeat, immutable full resume, and public label safety.

Fourteen hard abort families cover identity drift, pair/configuration leakage,
shape/norm failures, PH/MST/landscape/energy disagreement, subchunk gaps or
duplicates, invalid publication, unit/configuration/program resource breaches,
repeat/resume failure, and any label or outcome access.

## Artifact index

- `mv05v-source-freeze-2026-08-10.csv`
- `mv05v-base-pair-scope-2026-08-10.csv`
- `mv05v-full-group-queue-2026-08-10.csv`
- `mv05v-resource-projection-2026-08-10.csv`
- `mv05v-prefreeze-decision-2026-08-10.csv`
- `mv05v-prefreeze-summary-2026-08-10.csv`
- `mv05v-artifact-schema-2026-08-10.csv`
- `mv05v-validation-plan-2026-08-10.csv`
- `mv05v-abort-rules-2026-08-10.csv`

Private coordinate caches, labels, outcomes, PDFs, reviewer material,
`example_run.r`, and all later production artifacts remain untracked.

## Stop boundary

MV5-V stops before orchestration implementation, real-runner smoke, full
calculation, labels, robustness outcomes, uncertainty, comparison, ranking, or
claims. Spectral promotion, gene topology, cell/gene fusion, new data,
optimization/Rust, package-default changes, and pushing also remain outside
the next sprint.
