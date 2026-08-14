# MV5-AT broader realistic corrected-landscape workflow smoke specification v1

Date frozen: 2026-08-12

Authorization: MV5-AS completion `4d427ce`

## Purpose

MV5-AT tests the additive corrected-landscape producer through the actual
postprocessing orchestration across every already accepted MV5-AP-R1 realistic
stratum. It is a workflow and artifact-contract smoke, not a new scientific
analysis and not a downstream-consumption sprint.

## Frozen inputs

The input population is exactly `mv05ap-frozen-subset-2026-08-12.csv`: eight
strata, three deterministic min/middle/max-H1 diagrams per stratum, 24 diagrams
and 24 within-stratum pairs. The binding must reproduce from the accepted MV-04
manifest. Every result-file SHA-256, diagram SHA-256, stored provenance hash,
and scientific-eligibility flag must verify before calculation.

No alternate diagram, new data, cross-stratum pair, biological label, outcome,
or result-driven selection is permitted.

## Scientific and workflow contract

Each stratum is an independent corrected-only invocation of
`run_postprocessing_pipeline()` with all legacy clustering, visualization,
Betti, and cross-iteration flags false. The explicit v1 control remains
default-off outside the invocation and fixes:

- scientific `auto` routing with exact guard 500;
- all finite intervals and all active consecutive landscape levels;
- separate H0 and H1 squared-L2 distances plus descriptive combined distance;
- absolute and relative tolerance `1e-8`, 200 subdivisions, and one worker;
- artifacts-only downstream use and resume-or-fail behavior; and
- create-only atomic shards with completion written last.

The smoke must not populate legacy landscape fields or alter any existing
diagram file.

## Resource admission

Each unit contains exactly three pairs, receives a 750-second planned-wall
ceiling and a 2-GiB RSS capacity, and runs sequentially. The frozen planner
therefore admits both exact-only units (120-second conservative estimate) and
units with adaptive H1 (750-second conservative estimate).

This is supported by two accepted MV5-AP-R1 runs: the largest observed unit was
567.94 seconds and the largest observed RSS was 990,363,648 bytes. A unit is
rejected before calculation if its interval counts leave the profiled H0
383--499 / H1 79--2802 envelope, its pair count changes, or its planner exceeds
the frozen limits.

## Resume and validation

After first completion, each unit is invoked again against the same directory.
Every artifact path, size, modification time, and SHA-256 must remain unchanged,
and the returned sidecar must report resume. An independent validator must:

1. reproduce all 24 input bindings and hashes;
2. verify eight completion markers and all 24 pair shards;
3. reconstruct every H0/H1/combined matrix entry from its shard;
4. verify strict certificates and expected exact/adaptive routing;
5. verify all corrected-only legacy-field absences;
6. verify immutable resume and input-file immutability;
7. enforce process wall/RSS/exit bounds;
8. reproduce the public aggregate ledgers; and
9. confirm default, legacy, downstream, optimization, outcome, and new-data
   change counters remain zero.

The accepted MV5-AS clean repeat already covers independent three-diagram
recalculation. MV5-AT uses exact completed resume for all eight strata to avoid
duplicating more than 20 minutes of unchanged high-depth numerical work.

## Abort and decision boundary

Abort on source or manifest drift, input mismatch, resource refusal or breach,
uncertified output, artifact collision/corruption, non-atomic completion,
resume mutation, legacy/default/downstream drift, or test/check failure.

Passing MV5-AT may authorize only a separate prospective corrected-matrix
consumer prefreeze. It does not authorize clustering, visualization, Betti or
cross-iteration consumption; workflow defaults; legacy rewrites; new data;
biological claims; or Rust/optimization.
