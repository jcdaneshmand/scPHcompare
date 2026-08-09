# MV5-D5 label-closed SCT cell retrieval-input specification v1

## Document control

| Field | Value |
|---|---|
| Contract | `mv05d5_retrieval_group_bundle_v1` |
| Date frozen | 2026-08-08 |
| Status | Executed and independently validated |
| Input scope | Accepted MV5-D1 fold coordinates, MV5-D0 SCT caches, and MV5-D4 component distances |
| Output scope | Immutable query-to-training distances and rankings only |
| Outcome state | Closed |

## 1. Purpose and stop boundary

MV5-D5 converts the corrected cell-topology precomputation into the immutable
prediction-side inputs required by the frozen MV-05 retrieval benchmark. It
does not evaluate tissue, study, or approach labels. It does not choose a
winning method, cluster samples, run integration, construct a gene view, fuse
views, or acquire new data.

Every group is one frozen LOSO fold and cell-subsample seed. Its query axis is
the held-out study and its reference axis is every training sample in that
fold. The output is sufficient for continuous held-out retrieval, not for a
full sample-distance matrix or clustering.

## 2. Immutable inputs

Each group binds all of the following:

- one accepted `mv05d1_sct_cell_fold_v1` cache and payload hash;
- one seed-specific `mv05d5_mean_profiles_v1` bundle derived from the exact 90
  accepted MV5-D0 SCT caches;
- the accepted MV5-D4 component rows for the same fold and seed;
- the training/query partition and training-fit panel, centering, scaling, PCA,
  normalization-cache keys, and view cache keys;
- the implementation SHA-256.

The input tables contain no `tissue` or `approach` columns. `study` is used only
through the already frozen LOSO partition. All artifacts declare
`outcome_label_state = closed` and
`biological_outcomes_computed = FALSE`.

## 3. Frozen method panel

| Method | Role | Distance |
|---|---|---|
| `cell_landscape_h0_v1` | Confirmatory topology component | Raw exact all-active-level H0 landscape L2 |
| `cell_landscape_h1_v1` | Confirmatory topology component | Raw exact all-active-level H1 landscape L2 |
| `cell_landscape_h0_h1_raw_euclidean_v1` | Descriptive secondary only | `sqrt(H0^2 + H1^2)` |
| `cell_distribution_energy_shared_pca_v1` | Matched cell baseline | Square root of empirical V-statistic energy divergence in the same 384-cell, training-fit 30-PC space |
| `pseudobulk_shared_panel_euclidean_v1` | Context baseline | Euclidean distance between sample means on the same training-selected and training-standardized 500-gene panel |

The matched comparison for the cell-topology methods is the energy baseline.
Pseudobulk is retained as context and cannot substitute for the matched
distributional baseline.

## 4. Energy-distance definition

For point clouds `X` and `Y`, each containing the same 384-cell realization in
the common 30-PC coordinate system,

```text
sqrt(max(0,
  2 mean(||X_i - Y_j||)
  - mean(||X_i - X_j||)
  - mean(||Y_i - Y_j||)))
```

The within-sample means include the zero diagonal and both ordered directions,
matching `sqrt_v_statistic_energy_divergence_v1`. Within-sample terms are
cached once per group; this is an algebraic optimization, not a changed
formula.

## 5. Pseudobulk definition and missing features

For each accepted MV5-D0 sample/seed cache, MV5-D5 stores the named SCT
row-mean vector. For each fold it then:

1. selects the exact MV5-D1 training-derived 500-feature panel;
2. subtracts the MV5-D1 training-only feature means;
3. divides by the MV5-D1 training-only feature scales;
4. maps a held-out absent feature to zero after standardization, exactly as in
   the accepted cell-fold contract; and
5. computes Euclidean query-to-training distance.

No held-out value refits the panel, center, scale, or PCA.

## 6. Component-scale disposition

H0 and H1 remain separate primary methods. Multiplying all distances within one
method/fold/seed by a positive scalar cannot change its retrieval ordering, so
no scale is needed for either primary ranking.

The accepted MV5-D4 scope contains held-out-query-to-training-reference pairs
only. It contains zero within-training topology pairs. Fitting a median
off-diagonal component scale from the available rows would therefore use held-
out queries and violate the label/firewall fit boundary. Computing the missing
scope would require 262,675 additional biological topology pairs, or 525,350
new H0/H1 distance rows—7.431 times the completed MV5-D4 scope.

MV5-D5 consequently does not invent a fitted scale and does not launch that
expansion. The raw H0/H1 Euclidean combination already declared in MV5-D4 is
retained only as descriptive secondary output. It cannot override separate H0
and H1 results or the matched energy-baseline comparison.

## 7. Ranking and tie contract

Within each fold, seed, method, and held-out query:

1. sort training references by ascending distance;
2. break exact numerical ties by canonical training sample ID;
3. assign consecutive ranks beginning at one; and
4. record exact tie size, tie status, and the tie policy.

The prediction payload contains ranked training sample IDs and distances only.
No label-derived prediction or endpoint is stored.

## 8. Failure, identity, and resume contract

Every group contains five method-completion records. A method is either
`completed` with its exact expected pair count or has a structured failure
code and zero silent replacement. Each group identity binds its fold cache,
mean-profile bundle, MV5-D4 component subset, and implementation.

Publication is atomic. An existing bundle is reopened, fully validated, and
reused only when its identity is exact. A stale artifact causes refusal rather
than overwrite. Runtime and RSS measurements remain outside the scientific
payload so replay can be byte-identical.

## 9. Accepted production evidence

- 75 fold-seed groups and five fixed seeds;
- 35,350 biological query-to-training pairs;
- five methods and 176,750 unique ranked rows;
- 375 completed method groups and zero failures;
- all topology rows reproduced from MV5-D4;
- 225 independent energy and 225 independent pseudobulk formula checks;
- maximum independent numerical difference `1.136868e-13`;
- two admission bundles and all six public artifacts byte-identical;
- all 150 private production bundle/audit files unchanged after resume;
- zero exact-distance ties in the accepted data; and
- zero clustering, integration, gene-view, or biological-outcome jobs.

## 10. Next permitted stage

MV5-E may open the already frozen labels only after verifying the public
ranking hash and all MV5-D5 completion records. It may compute prespecified
held-out retrieval endpoints and blocked uncertainty. It may not refit,
rescale, reorder, select, or replace any MV5-D5 method or prediction artifact.
