# Phase 4 Subplan — Primary Benchmark Redesign

## Objective

Prespecify and execute a fair sample-level benchmark testing whether PH representations add stable information beyond matched non-topological representations.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P4-01 | Freeze estimands/endpoints | Statistical analysis plan defining primary/secondary estimands and units | Approved before confirmatory computation | `in_progress` — MV5-E/MV5-K retrieval and the MV5-L paired representation DID were frozen and executed; technical-mixing, clustering, robustness, and external-validation estimands remain open |
| P4-02 | Specify matched baselines | Pseudobulk, centroid, composition, distributional/OT, and graph-summary contracts | Same samples, splits, clustering/task, and endpoint used across representations | `in_progress` — both SCT and integrated cell views include matched same-coordinate energy and shared pseudobulk context; no SCT-versus-integrated comparison has been run |
| P4-03 | Specify PH primary method | Approved cell-view primary contract plus separately labeled gene-view secondary contract; distance, space, homology dimensions, filtration, and vectorization | Exactly one confirmatory cell-topology configuration, one approved secondary gene-topology configuration, and labeled sensitivities | `in_progress` — MV-02 typed implementation/analytical gate complete; MV-03 biological-scale eligibility/feasibility and later author-team approval remain |
| P4-04 | Define validation splits | Study-blocked/leave-one-study-out rules and leakage audit | No sample/study information crosses a prohibited boundary | `in_progress` — both 75-group SCT and integrated LOSO pipelines and their prediction locks independently validate; any paired representation comparison remains separately gated |
| P4-05 | Define balance/subsampling | Cell-count-matched and tissue-balanced repeated sampling plan | Repetitions, seeds, stopping, and aggregation are prespecified | `needs_review` — five fixed 384-cell seeds, tissue macro-averaging, and study-block aggregation frozen |
| P4-06 | Evaluate the new-data trigger | Existing-data gap analysis; alternatives; candidate/quality/access/cost assessment only if triggered | Continue with existing data unless a named estimand cannot be evaluated and owner/author team approves expansion | `in_progress` — MV5-M ranks external validation third but defers it under existing-data-first; current clustering/robustness work remains answerable without new data |
| P4-07 | Build controlled simulations | Generative scenarios separating biology, batch, composition, and geometry | Expected direction/failure modes are stated before results | `not_started` |
| P4-08 | Build robustness controls | Outliers, shuffles, composition controls, subsampling, dimensions/features/scales | Each control targets a named alternative explanation | `in_progress` - MV5-X completes and independently validates all 150 exact PC20 calculation groups; PC20 outcome estimands/reporting and the other three configurations remain separately gated |
| P4-09 | Execute benchmark | Immutable result bundle with manifests, logs, failures, and uncertainty | All prespecified methods/results are reported, including nulls | `in_progress` - retrieval/DID and prediction-locked clustering are complete; the label-closed PC20 robustness calculation is complete, but its outcome-evaluation contract and execution remain pending |
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
