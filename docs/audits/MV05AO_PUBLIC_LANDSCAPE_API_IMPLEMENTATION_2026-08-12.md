# MV5-AO public persistence-landscape API implementation

Date: 2026-08-12  
Accepted base: `c452932`  
Public API engine: `1deeec6`  
Independent validator: `0b85c20`

## Result

MV5-AO implements the dissertation-aligned, all-active persistence-landscape
definition as two additive public APIs:

- `persistence_landscape_distance()` for one diagram pair; and
- `persistence_landscape_distance_matrix()` for a complete named collection.

The implementation adds versioned pair and matrix result classes, immutable
diagram hashes and cache keys, separate H0/H1 results, a secondary descriptive
combined result, per-dimension integration diagnostics, and complete method and
source provenance. No existing API was removed or redirected.

## Scientific behavior

Scientific mode accepts typed H0/H1 persistence diagrams, rejects malformed or
nonpositive finite intervals, records and excludes essential intervals, uses
every finite interval and every active landscape level, and zero-pads unequal
landscape depths. Exact critical-breakpoint squared-L2 integration is the
default. Adaptive integration is available only explicitly (or through
explicit `auto`) and hard-fails unless its independent refinement estimate is
within the requested error bound. There is no public grid parameter and no
landscape-level cap.

H0 and H1 remain the primary outputs. `combined` is explicitly the secondary,
unweighted Euclidean norm of their distances; it does not erase the components.

## Compatibility boundary

Historical level-1, 100-point `[0, 1]` behavior is available only through
`mode = "legacy_k1_unit_grid_v0"`. The new read-only detector classifies
unversioned legacy shapes as candidates, never as losslessly convertible
scientific results, and never silently converts them. Existing workflow source,
defaults, filenames, and artifact writers are unchanged.

## Validation

The focused suite passes 54 expectations, including analytic single-tent and
shifted-tent oracles, a 12-level overlap oracle, essential-interval accounting,
identity/symmetry/triangle behavior, exact/adaptive agreement, exact-guard
failure, explicit legacy reproduction, canonical hashing/order, matrix/pair
consistency, deterministic cache identity, schema detection, and bounded
resource use. The complete repository suite passes with only its two expected
CRAN skips.

Two independent validator builds pass all 16 categories. Five deterministic
ledgers reproduce byte-for-byte. The resource smoke completes 10 diagrams and
45 pairs in 0.335 and 0.331 seconds respectively; both produce an identical
34,720-byte result and remain below the 30-second and 5-MB limits. Runtime is
recorded rather than misrepresented as byte-deterministic.

## Boundary and continuation

MV5-AO performs no project-data calculation, clustering analysis, workflow
default migration, legacy artifact rewrite, gene/cell fusion, new-data work,
Rust optimization, or biological inference.

The immediate next sprint should be MV5-AP: a prospective, read-only realistic
compatibility and resource prefreeze. It should exercise the new APIs against a
frozen representative subset of existing persistence diagrams, quantify exact
guard feasibility and adaptive requirements, compare scientific versus legacy
results without selecting a winner, verify versioned serialization/reload, and
decide whether a later opt-in workflow integration sprint is safe. MV5-AP must
not change any workflow default or rewrite any artifact.
