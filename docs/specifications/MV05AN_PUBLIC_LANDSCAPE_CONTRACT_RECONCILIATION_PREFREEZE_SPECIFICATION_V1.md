# MV5-AN public landscape contract reconciliation prefreeze specification v1

Date frozen: 2026-08-12

Accepted base: `1b29ed6`

Behavior/default/export/calculation changes: zero

## Purpose

MV5-AN reconciles the corrected dissertation-aligned persistence-landscape
definition with the package's actual functions, exported workflows, legacy
artifacts, diagnostics, production engine, reference implementation, tests,
and documentation. It classifies every landscape-named top-level R function
and every public workflow that can create, accept, or analyze landscape data.

## Existing state

The accepted robustness analyses use exact all-active critical-pair distances.
The R reference engine independently supports exact integration and an
error-controlled adaptive fallback. In contrast, the exported postprocessing
workflow reaches an internal historical implementation that constructs only
level one on 100 points over [0,1], combines H0/H1, and writes unversioned
landscape-list and distance-matrix RDS/CSV artifacts. Cross-iteration summaries
also use fixed 500-point unit grids and pointwise curve reductions.

No individual landscape function is currently exported. This makes in-place
redirection unnecessary and unsafe: existing workflow artifacts may be relied
upon despite their legacy semantics.

## Frozen scientific contract

The later public scientific API must use every finite positive-persistence
interval; exclude essential H0; retain all consecutive active levels; compute
H0 and H1 separately; integrate squared level differences over filtration;
use exact linear critical-pair segments within an explicit guard; use only a
partitioned, refinement-checked, error-controlled fallback outside that guard;
and fail rather than silently return partial or uncontrolled results. The
unweighted H0/H1 Euclidean norm is secondary and descriptive. No universal
uniform grid, level cap, post-hoc normalization, or weighting is allowed.

## API and compatibility decision

MV5-AO may add two new exports only:

- `persistence_landscape_distance()` for one diagram pair; and
- `persistence_landscape_distance_matrix()` for a complete named diagram set.

Both require versioned result classes and full provenance. MV5-AO may add
legacy-schema detection/read-only compatibility and an explicit legacy mode.
It may not redirect an existing function, overwrite an existing artifact,
change an exported workflow default, or recompute project results. A later
gate must evaluate realistic compatibility, resource, and migration evidence
before any workflow default changes.

## Artifact policy

Legacy landscape-list RDS, combined matrix RDS/CSV, and custom override inputs
must be identified as legacy grid artifacts. Exact pair ledgers, summary
ledgers, and versioned in-memory reference objects remain distinct. New matrix
artifacts require separate H0/H1 matrices, a descriptive combined matrix, and
a versioned provenance sidecar; they may not reuse ambiguous legacy filenames.

## Validation and stop boundary

Independent validation must rescan R definitions and NAMESPACE, confirm every
landscape pathway is classified, verify source hashes, reproduce all prefreeze
ledgers byte-for-byte, and prove all behavior/default/export/calculation
counters remain zero.

MV5-AN excludes code behavior changes, exports, project calculations,
clustering, gene/fusion work, new data, Rust/optimization, manuscript claims,
PDF/private/reviewer tracking, pushing, and `example_run.r`.
