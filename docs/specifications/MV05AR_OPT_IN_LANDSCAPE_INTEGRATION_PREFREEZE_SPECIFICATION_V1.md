# MV5-AR opt-in corrected-landscape integration prefreeze v1

Date: 2026-08-12

Authorization: MV5-AP-R1 completion `f6df036`

## Decision boundary

MV5-AR freezes an additive artifact-production integration only. A later
implementation may add an explicit corrected-landscape control to
`run_postprocessing_pipeline()` and a pass-through control to
`run_unified_pipeline()`. The default remains `NULL`, meaning no corrected
artifact work and byte-for-byte historical behavior.

Opting in creates a new versioned corrected-landscape artifact tree beside the
legacy landscape list and combined matrix. It does not redirect
`process_iteration_calculate_matrices()`, replace an iteration's
`landscape_l2_distance_matrix`, feed corrected H0/H1/combined matrices to
clustering or visualization, or pass corrected objects to legacy Betti or
cross-iteration landscape-curve analyses.

Downstream corrected-matrix consumption requires a later separate prefreeze.

## Proposed caller control

The later additive implementation may accept
`corrected_landscape_control = NULL` on the two entrypoints above. A non-NULL
named list must contain:

- `contract_id = "scph_corrected_landscape_workflow_control_v1"`;
- `enabled = TRUE`;
- `max_wall_seconds`, a positive explicit caller budget;
- `max_pairs`, a positive explicit pair-count ceiling;
- `max_rss_bytes`, at least 1.5 GiB and defaulting to 2 GiB;
- `workers = 1` in v1;
- `existing = "resume_or_fail"`;
- `downstream_use = "artifacts_only"`.

Scientific settings are fixed in v1: scientific mode, `auto` routing, exact
guard 500, absolute/relative tolerance `1e-8`, and no fixed grid, level cap,
interval removal, weighting, normalization, or tolerance relaxation. A caller
may not override these through `...`.

## Input and identity

Each iteration consumes its named `pd_list` only after validating that every
entry is a typed finite persistence diagram with a unique non-empty sample ID.
The canonical input-set identity binds the control contract, ordered sample
IDs, per-diagram SHA-256 hashes, scientific specification, engine version,
method, exact guard, tolerances, and subdivisions. Runtime, path, hostname, and
timestamps are excluded from cache identity.

The implementation must not accept legacy landscape-list or combined-matrix
objects as scientific inputs. Custom corrected input, if added later, requires
a distinct `corrected_landscape_v1_path` key and exact versioned schema; all
existing custom landscape keys remain explicitly legacy.

## Artifact tree

For each iteration and input-set hash, create a non-colliding directory under
`results/corrected_landscape_v1/`:

`<iteration-token>--<input-set-sha256>/`

It contains:

- `resource-plan-v1.csv` written before calculation;
- `input-manifest-v1.csv` with ordered IDs and diagram hashes;
- `pairs/<pair-id>--<pair-cache-key>.rds`, one versioned
  `scph_landscape_distance_v1` shard per canonical pair;
- `pair-index-v1.csv` with status, hashes, methods, certificates, and sizes;
- `distance-matrix-v1.rds`, a versioned
  `scph_landscape_distance_matrix_v1` object with H0, H1, and descriptive
  combined matrices;
- `provenance-v1.csv` with contract and engine identity;
- `completion-v1.csv`, written last and binding every prior artifact hash.

No legacy filename may be reused. No corrected artifact may overwrite any
existing file. Runtime/process evidence remains separate from cache identity.

## Atomic execution and resume

Pair IDs use lexical sample order. Each pair is computed through
`persistence_landscape_distance()`, written to a sibling temporary file,
reloaded and validated, then atomically renamed. An existing shard is reusable
only when class, cache key, pair IDs, input hashes, contract, method provenance,
certificate, file hash, and size all verify. Any mismatch hard-fails.

The matrix is assembled only after every canonical pair verifies. It is
cross-checked against every shard, serialized/reloaded, and atomically renamed.
`completion-v1.csv` is the sole completeness marker and is written last.
Interruption leaves resumable verified shards but no apparently complete
matrix. Resume may never rewrite a verified artifact or silently delete a
conflict.

## Resource admission

MV5-AP-R1 measured a maximum 567.94-second wall time and 990,363,648-byte RSS
for a three-pair maximum-depth gene stratum. V1 therefore runs one worker and
plans pair shards before execution.

The audited planning model is deliberately conservative:

- exact/exact pair: 30 wall seconds;
- any adaptive dimension: 240 wall seconds;
- iteration startup/finalization: 30 wall seconds;
- required RSS budget: at least 1.5 GiB;
- H0/H1 interval counts outside the observed 383--499 / 79--2802 envelope are
  `profiling_required`, not silently capped;
- execution is admitted only when planned pairs do not exceed `max_pairs` and
  planned wall time does not exceed `max_wall_seconds`.

Rejecting an unprofiled or over-budget run is computational protection, not a
scientific change. The plan must still report the full requested interval and
pair counts.

## Legacy coexistence and rollback

All existing legacy functions, files, custom keys, defaults, clustering paths,
Betti summaries, cross-iteration curves, and plots remain unchanged. Schema
detection is read-only and silent conversion is forbidden.

Because integration is additive, rollback is disabling the explicit control.
Incomplete corrected directories may be retained for audited resume or moved
aside by an explicit owner action; implementation may not delete them
automatically. Legacy outputs remain the workflow's downstream inputs until a
separate consumption prefreeze passes.

## Validation and stop rules

The implementation sprint must prove default-off equivalence, analytic
oracles, pair-shard/matrix equivalence, H0/H1 separation, cache identity,
failure without uncertified output, atomic interruption/resume, immutable
resume, serialization, legacy coexistence, resource admission/rejection,
clean repeat, and source/artifact immutability. It must include a realistic
bounded smoke but may not rerun the complete MV5-AP-R1 corpus without a new
execution authorization.

Stop for any ambiguous input or downstream consumer, legacy overwrite,
default-on behavior, uncertified result, combined-as-primary use, non-atomic
write, cache collision, unprofiled resource admission, grid/cap/tolerance
shortcut, workflow behavior drift, or test/check failure.
