# MV7-A robustness/confounding evidence-map audit

## Disposition

`authorize_no_new_ph_confounding_prefreeze_only`

MV7-A reconciles the completed MV5 cell-view robustness, representation, and
clustering evidence with the completed MV6 cell/gene/fusion evidence. The map
contains 14 robustness axes and 10 confounding axes. It authorizes only a
separately committed MV7-B specification for diagnostics that reuse accepted
sample-level outcomes and authoritative metadata. It does not authorize new
PH, new data, method selection, fusion expansion, default changes, or claim
promotion.

## Landscape definition

The landscapes remain appropriate and unchanged:

- all finite positive-persistence intervals are retained;
- the essential H0 interval is excluded;
- every consecutive active landscape level is used;
- H0 and H1 are calculated and reported separately;
- squared-L2 integration is exact or error-controlled on each dimension's
  support;
- no universal grid or level cap is introduced; and
- the unweighted H0/H1 Euclidean composite is descriptive only.

The completed sensitivity panels change cell count, PCA dimension, or point
geometry while holding this landscape definition fixed. Consequently, their
differences are not artifacts of changing grid resolution or truncating
landscape levels.

## Coverage result

Complete or adequately bounded evidence already exists for five fixed seeds,
separate H0/H1, the no-grid/no-cap landscape contract, complete
Vietoris–Rips filtration, study-blocked validation, equal-tissue aggregation,
cell-view depth (192/256/384), cell PC20/PC30, cell Euclidean/cosine-chord,
cell SCT/integrated comparison, and secondary PAM/average-linkage clustering.

The final gene view is narrower: it uses 384 cells, the fixed technically
transductive 500-gene global core, the 30-PC cell transform, SCT only, and
Pearson-correlation chord geometry. Gene-panel-size, gene-metric, and
integrated-gene sensitivities have not been run and would require new PH. They
are real gaps, but they are not yet justified as the next computation.

Library size and harmonized cell-type composition are unavailable in the
accepted metadata. Retained cell count is not used as a library-size proxy.
Technology/approach is recorded but partly nested within study and tissue, so
only descriptive association—not a causal technology effect—is identifiable.

## Selection-resistant next action

The fixed prerequisite order admits three related no-new-PH questions:

1. leave-one-study and leave-one-tissue influence on the final separate-view
   summaries and fixed contrasts;
2. association with retained cell count; and
3. sequencing-approach stratification with explicit nesting diagnostics.

These analyses can reuse the locked 4,050 MV6-H query-method outcome rows and
the authoritative metadata. The MV7-B contract must freeze exact methods,
transforms, missingness, influence thresholds, multiplicity, and reporting
before execution. All studies, tissues, views, dimensions, and fixed methods
must be reported; no favorable subgroup or fusion weight can be selected.

External data remain deferred. The trigger is a named manuscript-critical
estimand that remains unresolved after existing-data diagnostics—not a general
desire for a larger dataset.

## Reproducibility

- prospective implementation commit: to be bound by the final source freeze;
- 23/23 source files exist and were independently rehashed;
- independent validation: 10/10 categories pass;
- repeat: 10/10 production artifacts are byte-identical;
- source-freeze, decision, validation, and repeat hashes are bound after the
  exact-commit build.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | approve; corrected landscape contract preserved |
| Correctness demonstrated? | pass; independent 10/10 and byte repeat 10/10 |
| Computation feasible? | yes; next diagnostics reuse existing CSV artifacts and require no PH |
| Biological interpretation permitted? | secondary/descriptive with explicit confounding limits |
| Next action | prefreeze MV7-B no-new-PH influence/confounding diagnostics |
