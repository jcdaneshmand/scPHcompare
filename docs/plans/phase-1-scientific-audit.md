# Phase 1 Subplan — Scientific and Implementation Audit

## Objective and boundary

Determine what the current code actually computes, whether it matches the paper/dissertation, and which claims survive. Diagnose first; corrections occur only after specifications and regression tests are approved.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P1-01 | Observation-unit trace | Table mapping every input row, cluster label, distance, and metric to cell/sample/study/group | No metric has an ambiguous unit or weighting rule | `in_progress` — PH conflict traced and historical diagrams invalidated; MV-02 typed cell/gene objects now enforce point/coordinate axes, IDs, metric, fit scope, and result provenance; broader downstream metric/weighting trace remains |
| P1-02 | Reproduce legacy metrics | Script/configuration that reconstructs current reported metrics and identifies irreproducible values | Values match within declared tolerance or discrepancy is localized | `not_started` |
| P1-03 | Same-unit recalculation | Once-per-sample PH evaluation plus comparison with cell-replicated evaluation | Effect of cell-count weighting is quantified for every headline metric | `not_started` |
| P1-04 | PCA/integration-space trace | Fit/transform provenance, feature sets, dimensions, centering/scaling, shared versus per-sample basis | Every comparison states whether coordinates are commensurate and why | `in_progress` — MV5-B reference-only feature/PCA plus held-out `pcaproject`/`IntegrateEmbeddings` route passes a deterministic label-closed fixture; full existing-data fold trace remains |
| P1-05 | Landscape specification audit | Paper/dissertation/code comparison; approved mathematical definition; toy fixtures | One definition covers levels, grids, aggregation, H0/H1 combination, and distance | `in_progress` — definition, feasibility, owner approval, R oracle, independent cross-check, and historical-diagram eligibility complete; eligible-diagram production validation and broader author-team confirmation remain |
| P1-06 | BDM/SDM/LDM verification | Hand-computable fixtures and automated tests for formulas, orientation, scaling, infinite deaths, censoring | Tests agree with analytical results and fail under known wrong implementations | `in_progress` — LDM oracle independently verified and a Persim sign-crossing regression added; historical BDM/SDM/LDM invalidated upstream; corrected BDM/SDM and censoring verification remain |
| P1-07 | Statistical reconstruction | Null-generation trace, seeds, blocks, cluster-size rules, p-value correction, multiplicity families, functional-test assumptions | Every p-value/interval has a documented estimand and valid procedure or is withdrawn | `in_progress` — prospective MV-05 blocks, seeds, nonzero randomization p-values, bootstrap, Holm families, cluster rules, and failure policy frozen; legacy-result reconstruction and benchmark execution remain |
| P1-08 | Dataset integrity audit | Regenerated sample/cell table, accession/metadata checks, pre/post-filter reconciliation, imbalance summaries | Counts reconcile to source objects; unresolved discrepancies are isolated | `in_progress` — historical 127-to-124 and validation 25-to-25 flows reconciled; future rerun instrumentation implemented |
| P1-09 | Execution/package audit | Locked-environment test/check logs; missing-function trace; regression tests for discovered failures | Package checks and toy execution have explicit pass/fail records | `not_started` |
| P1-10 | Claim disposition | Claim-to-evidence table: supported, narrowed, untested, contradicted | Every abstract/result/conclusion claim has one disposition | `not_started` |
| P1-11 | Evaluate G1 | Approved observation unit, landscape definition, valid-result list, invalidation list | Author team makes and records a gate decision | `not_started` |

## Mandatory comparisons

- Cell-replicated versus one-row-per-sample metrics.
- Unweighted sample estimates versus cell-count-weighted estimates.
- First-landscape-only, all-level L2, and legacy `rowMeans()` behavior.
- Oracle `k` versus a non-oracle selection/stability rule.
- Bone-marrow `k = n` as a technical distinctness test, not biological clustering validation.
- Common filtration scale versus sample-adaptive scale.

## Gate G1 checklist

- [ ] Primary observation unit and weighting are approved.
- [ ] Mathematical landscape specification is approved.
- [ ] Dataset and metadata discrepancies are resolved or bounded.
- [ ] Current results are partitioned into valid, regenerable, and withdrawn.
- [ ] Statistical analysis plan requirements are listed.
- [ ] Optimization remains frozen until corrected reference behavior exists.

## Stop condition

If the headline advantage disappears under matched sample-level evaluation, do not conceal or tune around it. Record the result and reframe the project around integration diagnostics, failure characterization, or another evidence-supported contribution.
