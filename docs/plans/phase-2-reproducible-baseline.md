# Phase 2 Subplan — Reproducible Baseline and Repository Health

## Objective

Turn the approved Phase 1 method into a deterministic, testable package and workflow that runs from a clean clone.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P2-01 | Define pipeline contracts | Typed schemas for stage inputs/outputs, identifiers, units, and errors | Each stage can be tested independently without hidden global state | `in_progress` — sample-flow, attempt, status, and reconciliation contracts implemented |
| P2-02 | Build analytical fixtures | Tiny point clouds with known H0/H1 and summary distances | Expected outputs are derived and asserted with tolerances | `in_progress` — maintained square/line baseline asserts H0/H1 birth/death values at `1e-10`; summary-distance fixtures remain |
| P2-03 | Build realistic fixture | Small redistributable/de-identified subset plus expected outputs | Completes quickly on CI and exercises the full path | `complete` — deterministic Harmony route covers H0 and H1 with dimension-specific reference hashes |
| P2-04 | Lock environments | R/system dependency inventory, lockfile verification, platform notes | Clean environment resolves declared dependencies reproducibly | `in_progress` — empty-cache Ubuntu restore verified 265/265 records; CI and other platforms remain |
| P2-05 | Establish determinism | Seed hierarchy, parallel RNG rules, provenance metadata | Repeated runs match declared deterministic/tolerance policy | `in_progress` — analytical baseline seed/RNG/tolerance contract verified; full-pipeline hierarchy remains |
| P2-06 | Add workflow orchestration | Evaluation and prototype using `targets` or justified alternative | Invalidated inputs rerun only dependent stages; hashes prevent stale checkpoints | `not_started` |
| P2-07 | Add CI | Parse, lint if adopted, unit tests, docs, package check, privacy check, toy smoke run | Required checks pass on supported platforms from clean checkout | `in_progress` — local required workflow covers privacy, check, and installed H0/H1 fixtures; hosted observation awaits push |
| P2-08 | Improve observability | Structured bounded logs, timing/memory fields, configuration and revision IDs | A failed run identifies stage, input, configuration, and restart point | `in_progress` — structured PH attempts, stage timings, and sampled monitor/child/process-tree RSS implemented; revision fields remain |
| P2-09 | Audit public repository | README/install/example/license/citation/data-policy/size/ignore remediation ledger | A new user can install and run the toy example without private data | `not_started` |
| P2-10 | Evaluate G2 | Clean-clone reproduction report | Independent clean run passes or exceptions are approved | `not_started` |

## Gate G2 checklist

- [x] Toy analytical and realistic fixtures pass.
- [x] Corrected reference outputs are versioned with tolerances.
- [x] CI exercises the public route and privacy guard locally; hosted evidence remains.
- [ ] Dependency and seed policies are documented.
- [ ] Clean clone completes the toy workflow.
- [ ] No confidential material is included.
