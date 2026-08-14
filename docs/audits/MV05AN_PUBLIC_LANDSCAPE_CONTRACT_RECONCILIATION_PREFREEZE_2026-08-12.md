# MV5-AN public landscape contract reconciliation prefreeze

Date: 2026-08-12  
Accepted base: `1b29ed6`  
Prospective engine: `f2aab37`

## Result

MV5-AN classified every landscape-named top-level R function plus the accepted
Python production engine: 46 pathways in total. It also classified six exported
workflow entrypoints and six artifact schemas. No pathway remained ambiguous,
and the sprint changed no behavior, export, default, artifact, or calculation.

The audit found two coexisting systems. Accepted robustness production uses
exact, all-active persistence landscapes. Exported historical workflows can
still use a level-1, 100-point grid on `[0, 1]`, combined H0/H1 distances, and
unversioned artifacts; cross-iteration summaries also use fixed 500-point grids.

## Frozen scientific contract

The public contract is frozen to:

1. accept typed persistence-diagram inputs;
2. include every finite interval;
3. exclude the essential H0 interval;
4. retain every active landscape level with zero-padding across diagrams;
5. calculate H0 and H1 separately;
6. integrate squared landscape differences over filtration scale;
7. use exact critical-pair integration when exact mode is requested;
8. use checked adaptive refinement only as an explicit error-controlled fallback;
9. expose any H0/H1 composite as a secondary descriptive result;
10. impose no universal fixed level cap or uniform-grid count;
11. return method, dimension, tolerance, interval, and version provenance; and
12. fail clearly when the requested guarantees cannot be met.

## Compatibility decision

MV5-AO may add two public APIs:

- `persistence_landscape_distance()`
- `persistence_landscape_distance_matrix()`

They must return versioned, non-colliding result schemas. Historical artifacts
remain read-only, may be identified by legacy-schema detection, and may be
reproduced only through an explicit legacy mode. Existing workflow defaults are
not redirected and existing artifacts are not overwritten. Level-1, 100-point
legacy artifacts cannot be converted losslessly without original persistence
diagrams or recomputation.

## Validation

The audit used 20 source artifacts and rescanned 45 R landscape functions plus
the accepted Python engine, `NAMESPACE`, workflow callers, documentation, tests,
and artifact writers. All 12 independent validation categories passed. All 15
production and validation ledgers reproduced byte-for-byte across clean builds.

## Boundary and continuation

MV5-AN made no behavior, default, export, artifact, or project-calculation
change. Its sole authorized continuation is MV5-AO: implement the two additive,
versioned APIs and result/provenance schemas; analytic oracles; exact/adaptive
agreement checks; legacy-schema detection and explicit legacy mode;
documentation; and bounded resource smoke evidence. Workflow-default changes,
legacy artifact rewrites, clustering expansion, gene/cell fusion, new data,
Rust optimization, and new biological claims remain outside MV5-AO.
