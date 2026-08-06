# Phase 5 Subplan — Method Expansion

## Objective

Test generalizability with a prespecified representative method set, not an unconstrained parameter search.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P5-01 | Select integration panel | Paradigm coverage, rationale, native output contract, feasibility/licensing | Panel is approved before full runs; omissions are justified | `not_started` |
| P5-02 | Select topology panel | Cell/gene views; within-sample metrics; H0/H1; landscape, diagram, and curve distances; limited H2 feasibility | Confirmatory cell, secondary gene, sensitivity, and exploratory fusion configurations are labeled | `in_progress` — MV-01 pilot panel frozen; implementation/feasibility and later exploratory admission remain |
| P5-03 | Select task/clustering panel | Hierarchical, spectral, PAM, consensus/fusion methods, graph communities, and continuous tasks where valid | Same algorithms/selection rules apply to comparable distance representations; fusion is benchmarked against both components | `in_progress` — PAM/stability primary and average/spectral sensitivities frozen prospectively; MV-05 execution remains |
| P5-04 | Define tuning budget | Parameter grids, defaults, compute budget, failure handling | No method receives outcome-informed preferential tuning | `not_started` |
| P5-05 | Pilot feasibility | Runtime/memory/failure and output-validity report on fixture data | Invalid combinations are removed by documented assumptions, not performance results | `not_started` |
| P5-06 | Execute integration matrix | Versioned outputs and logs for approved methods | Native documented representation is used and all failures retained | `not_started` |
| P5-07 | Execute topology/task matrix | Prespecified comparisons with multiplicity families | Complete result matrix or documented technical exclusions | `not_started` |
| P5-08 | Run ablations | Components removed/varied one at a time | Conclusions identify which components drive behavior | `not_started` |
| P5-09 | Synthesize boundaries | Works/fails/uncertain table by data and method conditions | Generalization claims match observed scope | `not_started` |
| P5-10 | Evaluate G5 | Approved generalizability report | Seurat-specific versus general findings are explicit | `not_started` |

## Gate G5 checklist

- [ ] Method panel represents justified paradigms.
- [ ] Confirmatory and exploratory analyses are separated.
- [ ] Tuning and multiplicity rules were followed.
- [ ] Failures and null results are reported.
- [ ] Claims state where evidence works, fails, or is uncertain.
