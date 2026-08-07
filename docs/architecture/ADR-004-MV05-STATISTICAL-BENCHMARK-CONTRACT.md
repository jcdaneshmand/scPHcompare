# ADR-004: MV-05 statistical benchmark contract

## Status

Accepted on 2026-08-06 as a prospective implementation contract. Biological
execution, method ranking, fusion, manuscript claims, and public-default
changes remain unauthorized.

## Context

The dissertation intended sample-level comparison of persistence landscapes,
but the legacy evaluation code selects cluster counts from known tissue,
study, or approach labels and then copies each sample assignment to every cell
before calculating agreement. That changes the inferential unit and weights
samples by their cell counts. The immutable MV-04 pilot is valid technical
evidence, but its 14 samples were selected for feasibility, no pilot tissue
spans studies, and its transformations were fit over all pilot samples.

The existing large-cohort metadata contains 124 samples. Five tissues occur in
at least two studies, yielding a prospective 90-sample, 15-study candidate set.
The separate 25-sample bone-marrow cohort contains one study and can support
only within-study technical approach diagnostics. Current Seurat integration
is transductive because it jointly calls `FindIntegrationAnchors()` and
`IntegrateData()` and has no audited held-out query mapping path.

## Decision

1. Treat the biological sample as the observation, study as the outer fold and
   dependence block, cells as within-sample subsamples, and subsampling seeds
   as repeated technical realizations.
2. Use leave-one-study-out validation for the large cohort. Restrict the
   confirmatory tissue endpoint to tissues represented in at least two studies
   and present in the training portion of a fold.
3. Fit feature selection, scaling, PCA, integration/reference mapping, and
   distance normalization on training samples only. Hash matrices and
   predictions before evaluation labels are opened.
4. Make cross-study tissue retrieval macro mean reciprocal rank the primary
   continuous endpoint. Report biological conservation and residual technical
   structure separately.
5. Compare cell H0 and H1 separately against a matched energy-distance baseline
   on the same cells and shared training-fit coordinates. Use pseudobulk as a
   context baseline and gene-correlation distance as the matched secondary
   gene-view baseline.
6. Make PAM with label-free repeated-subsample stability the primary clustering
   analysis. Keep average linkage as a sensitivity; admit spectral clustering
   only after its graph, eigensolver, seed, eigengap, and provenance contract is
   implemented. Keep known-label `k` historical-only.
7. Use paired tissue-stratified study-block bootstrap intervals and optional
   study-level sign flips with nonzero Monte Carlo p-values. Apply Holm
   adjustment within three prespecified families and make no inferential claim
   with fewer than four contributing study blocks.
8. Retain every failure and non-estimable case. Do not replace methods, samples,
   folds, or seeds after outcomes are opened.
9. Use existing data first. Do not open a new-data search unless an executed
   existing-data analysis leaves a named, decision-relevant estimand gap.

## Activation boundary

This decision authorizes contract validators, fold-fit scaffolding, matched
baseline implementations, deterministic clustering helpers, provenance, and
technical smoke tests. Confirmatory integrated execution remains blocked until
a training-reference/held-out-query mapping contract passes. G-MV5 remains
open until the frozen single-view benchmark is executed and cell and gene views
receive independent works/fails/uncertain dispositions.

## Consequences

The MV-04 pilot cannot be repurposed as biological validation, and legacy
cell-replicated cluster scores cannot enter the new primary analysis. The next
implementation stage is MV5-A/MV5-B: build fold manifests and baseline fixtures,
then prove inductive integration on a bounded two-study technical case without
opening biological outcomes.

The complete contract and evidence are in
`docs/specifications/MV05_STATISTICAL_BENCHMARK_PLAN_V1.md` and
`docs/audits/MV05_STATISTICAL_PLAN_FREEZE_2026-08-06.md`.
