# scPHcompare Auditable Execution Subplans

## Purpose

This directory converts `PROJECT_PLAN.md` into executable, reviewable work. The main plan defines strategy; these files define the evidence required to complete each phase.

## Control rules

1. Use the task IDs in filenames, commits, analysis outputs, and decisions.
2. Do not mark a task complete until its acceptance test is satisfied and its evidence path is recorded.
3. Record diagnoses before changing scientific behavior. A correction must have a regression test or a documented reason why one is impossible.
4. Preserve raw evidence. Derived files may be regenerated; source data, historical outputs, and confidential materials must not be overwritten.
5. Record failed and null results. They are evidence, not disposable attempts.
6. Gates are author-team decisions. Completing tasks does not automatically approve a gate.
7. Reviewer correspondence remains under the Git-ignored `docs/private/` tree and must never be copied into these public plans.
8. Credit for Dr. Eric Rouchka and Dr. Akshitkumar Mistry must be preserved; final author order and CRediT roles require author-team approval.
9. Use existing project data first. A new dataset enters scope only through the Phase 4 evidence-gap decision, not merely because additional data are available.

## Status vocabulary

| Status | Meaning |
|---|---|
| `not_started` | No task evidence has been collected. |
| `in_progress` | Work has begun and the current evidence path is logged. |
| `blocked` | A named dependency or decision prevents progress. |
| `needs_review` | Acceptance evidence exists and awaits review. |
| `complete` | Acceptance test passed and reviewer/date are recorded. |
| `superseded` | Replaced by a documented decision; retained for audit history. |

## Phase index

| Phase | Subplan | Depends on | Current state | Gate |
|---|---|---|---|---|
| 0 | [Preservation and provenance](phase-0-preservation.md) | None | `in_progress` | G0 |
| 1 | [Scientific and implementation audit](phase-1-scientific-audit.md) | Preservation controls from Phase 0 | `in_progress` | G1 |
| 2 | [Reproducible baseline and repository health](phase-2-reproducible-baseline.md) | Approved Phase 1 specifications | `in_progress` | G2 |
| 3 | [Literature, references, and figures](phase-3-literature-and-figures.md) | Can overlap Phases 0–2 | `in_progress` | G3 |
| 4 | [Primary benchmark redesign](phase-4-primary-benchmark.md) | G1; statistical plan | `in_progress` | G4 |
| 5 | [Method expansion](phase-5-method-expansion.md) | G4 | `in_progress` | G5 |
| 6 | [Biological and practical validation](phase-6-biological-validation.md) | G4; may overlap Phase 5 | `not_started` | G6 |
| 7 | [Profiling and optimization](phase-7-optimization.md) | Stable corrected pipeline; normally G5/G6 | `in_progress` | G7 |
| 8 | [Manuscript and response package](phase-8-manuscript.md) | Stable confirmatory results | `not_started` | G8 |
| 9 | [Release and archive](phase-9-release.md) | G8 | `not_started` | G9 |

## Audit record format

For every task, append a row to [WORK_LOG.md](WORK_LOG.md) containing:

- date and task ID;
- person or agent performing the work;
- exact input revision/configuration;
- evidence paths;
- outcome and unresolved issues;
- reviewer and status transition.

Decisions that alter scope, mathematical definitions, primary endpoints, authorship, or release behavior must also be added to the decision log in `PROJECT_PLAN.md`.

## Cross-phase roadmap

The project-owner-directed [dual-view and multiview topology sprint roadmap](DUAL_VIEW_MULTIVIEW_SPRINT_ROADMAP.md) coordinates work across Phases 1, 4, 5, 7, and 8. It treats cell topology as the corrected confirmatory view, gene topology as a separately specified secondary view, and fusion as exploratory until both component views pass independent correctness and feasibility gates.

## Recommended starting sequence

1. **P0-01:** freeze and record the current canonical repository state.
2. **P0-02:** review the remote/local divergence without pulling or overwriting anything.
3. **P0-03/P0-04:** classify and checksum irreplaceable artifacts.
4. **P1-01:** construct the observation-unit trace.
5. **P1-02/P1-03:** reproduce current metrics and recompute them once per sample.
6. **P1-05:** lock the persistence-landscape mathematical specification.

The first high-value scientific decision is G1. Phase 0 prevents loss; Phase 1 determines which prior results can still be trusted.

## New-data trigger

The default is to complete Phases 0–3 and the initial Phase 4 design using the existing data. A new dataset is considered only when all of the following are recorded:

1. the exact claim or estimand the existing data cannot evaluate;
2. why reanalysis, blocking, subsampling, or simulation cannot close the gap;
3. the required biological groups, study structure, annotations, and sample count;
4. candidate datasets and their access, license, quality, comparability, and confounding risks;
5. estimated preprocessing, compute, storage, and validation cost;
6. the decision to acquire, defer, or reject, with owner/author-team approval.
