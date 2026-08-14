# Phase 4 Subplan — Primary Benchmark Redesign

## Objective

Prespecify and execute a fair sample-level benchmark testing whether PH representations add stable information beyond matched non-topological representations.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P4-01 | Freeze estimands/endpoints | Statistical analysis plan defining primary/secondary estimands and units | Approved before confirmatory computation | `needs_review` — prospective MV-05 endpoint and unit registry frozen; no outcome computation |
| P4-02 | Specify matched baselines | Pseudobulk, centroid, composition, distributional/OT, and graph-summary contracts | Same samples, splits, clustering/task, and endpoint used across representations | `in_progress` — MV5-C baselines remain exact; MV5-D1 cell coordinates are complete, but baseline production and endpoint comparison remain closed |
| P4-03 | Specify PH primary method | Approved cell-view primary contract plus separately labeled gene-view secondary contract; distance, space, homology dimensions, filtration, and vectorization | Exactly one confirmatory cell-topology configuration, one approved secondary gene-topology configuration, and labeled sensitivities | `in_progress` — MV-02 typed implementation/analytical gate complete; MV-03 biological-scale eligibility/feasibility and later author-team approval remain |
| P4-04 | Define validation splits | Study-blocked/leave-one-study-out rules and leakage audit | No sample/study information crosses a prohibited boundary | `in_progress` — 75/75 label-closed LOSO fold-seed transforms independently validate; feature panel/scaling/PCA are training-only; outcomes remain closed pending later immutable prediction gate |
| P4-05 | Define balance/subsampling | Cell-count-matched and tissue-balanced repeated sampling plan | Repetitions, seeds, stopping, and aggregation are prespecified | `needs_review` — five fixed 384-cell seeds, tissue macro-averaging, and study-block aggregation frozen |
| P4-06 | Evaluate the new-data trigger | Existing-data gap analysis; alternatives; candidate/quality/access/cost assessment only if triggered | Continue with existing data unless a named estimand cannot be evaluated and owner/author team approves expansion | `needs_review` — 90 samples across five multi-study tissues support the current attempt; no new-data trigger opened |
| P4-07 | Build controlled simulations | Generative scenarios separating biology, batch, composition, and geometry | Expected direction/failure modes are stated before results | `not_started` |
| P4-08 | Build robustness controls | Outliers, shuffles, composition controls, subsampling, dimensions/features/scales | Each control targets a named alternative explanation | `not_started` |
| P4-09 | Execute benchmark | Immutable result bundle with manifests, logs, failures, and uncertainty | All prespecified methods/results are reported, including nulls | `not_started` |
| P4-10 | Evaluate G4 | Benchmark report and gate decision | Fairness, blocking, reproducibility, and practical question approved | `not_started` |

## Analysis guardrails

- Separate biological-conservation and batch-removal endpoints.
- Estimate direct paired differences with uncertainty.
- Label oracle-`k` analyses separately; primary results cannot depend on known class count unless justified by task definition.
- Evaluate cell and gene topology independently before any fused result; fusion must be compared with both component views.
- Freeze confirmatory choices before inspecting full benchmark outcomes.
- Do not add a dataset simply to increase volume. New data must resolve a specific limitation that materially affects the paper's defensibility or utility.

## Gate G4 checklist

- [ ] Primary estimand and endpoints are approved.
- [ ] PH and non-PH methods operate on matched observation units.
- [ ] Blocking, leakage, balance, and subsampling checks pass.
- [ ] Existing-data sufficiency is documented; any new-data decision satisfies the recorded trigger.
- [ ] Confirmatory results are complete and reproducible.
