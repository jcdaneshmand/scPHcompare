# MV5-J label-closed integrated cell retrieval-input specification v1

## Document control

| Field | Value |
|---|---|
| Contract | `mv05j_integrated_retrieval_group_bundle_v1` |
| Date frozen | 2026-08-09 |
| Status | Executed and independently validated |
| Input scope | Accepted MV5-D1 folds, MV5-D5 mean-profile bundles, MV5-G integrated coordinates, and MV5-I integrated H0/H1 components |
| Output scope | Immutable integrated query-to-training distances and rankings only |
| Outcome state | Closed |

## 1. Purpose and stop boundary

MV5-J constructs the prediction-side inputs needed to evaluate the integrated
cell-topology view in a later, separately prediction-locked sprint. Every group
is one of the already frozen 15 leave-one-study-out folds and five cell-subsample
seeds. Its query axis is the held-out study and its reference axis is every
training sample in that fold.

This sprint may calculate distances and deterministic neighbor rankings. It
must not open tissue or other biological outcomes, calculate retrieval
performance endpoints, compare the integrated view with the accepted SCT
result, tune or select methods, cluster samples, construct gene topology, fuse
views, acquire data, or make manuscript claims.

## 2. Immutable group inputs

Each group binds:

- one accepted `mv05d1_sct_cell_fold_v1` record, including its frozen query and
  training partitions, normalization keys, training-selected panel, center,
  and scale;
- one accepted seed-specific `mv05d5_mean_profile_bundle_v1`, reused read-only
  for the same-panel pseudobulk context baseline;
- one accepted MV5-G integrated-coordinate group, including its private file
  hash, group cache key, payload hash, coordinate-set hash, and all 90 sample
  coordinate-view identities;
- the accepted MV5-I component rows for the identical fold, seed, query axis,
  and training axis; and
- the implementation SHA-256.

The group identity and cache key change if any bound input changes. Existing
outputs are reopened and reused only when the complete identity validates.
Stale artifacts are refused, never silently replaced.

## 3. Frozen method panel

| Method | Role | Distance |
|---|---|---|
| `integrated_cell_landscape_h0_v1` | Confirmatory integrated topology component | Raw exact all-active-level H0 landscape L2 from MV5-I |
| `integrated_cell_landscape_h1_v1` | Confirmatory integrated topology component | Raw exact all-active-level H1 landscape L2 from MV5-I |
| `integrated_cell_landscape_h0_h1_raw_euclidean_v1` | Descriptive secondary only | `sqrt(H0^2 + H1^2)` |
| `integrated_cell_distribution_energy_v1` | Matched integrated-cell baseline | Square root of empirical V-statistic energy divergence in the exact same 384-cell, 30-dimensional MV5-G coordinate space |
| `pseudobulk_training_standardized_panel_v1` | Context baseline | Euclidean distance between sample means on the frozen training-selected and training-standardized 500-gene panel |

H0 and H1 remain separate primary topology methods. The raw composite is
descriptive and cannot override either component. Energy is the matched
same-cell/same-coordinate distributional baseline. Pseudobulk is context only.

## 4. Integrated energy definition

For MV5-G coordinate clouds `X` and `Y`:

```text
sqrt(max(0,
  2 mean(||X_i - Y_j||)
  - mean(||X_i - X_j||)
  - mean(||Y_i - Y_j||)))
```

The within-sample means include diagonal zeros and both ordered directions.
All operands come from the exact integrated coordinate matrices used to produce
the corresponding MV5-H/MV5-I topology. Within-sample terms may be cached as an
algebraic optimization. No coordinate, cell, or dimension may be substituted.

## 5. Pseudobulk context definition

MV5-J reuses the accepted MV5-D5 mean-profile bundles without modification.
For each fold it selects the exact MV5-D1 500-feature panel, applies MV5-D1
training-only centers and scales, maps held-out absent features to zero after
standardization as frozen in MV5-D1, and calculates Euclidean query-to-training
distance. This baseline is intentionally shared with the SCT retrieval input
contract so later view comparisons do not confound topology representation
with a changed pseudobulk comparator.

## 6. Component-scale disposition

The available MV5-I scope contains directed held-out-query-to-training-reference
pairs and no within-training topology pairs. A fitted H0/H1 component scale
therefore cannot be estimated from the available topology rows without using
held-out queries. MV5-J forbids held-out-derived component-scale fitting and
does not expand topology scope.

H0 and H1 need no scalar for their own within-method rankings. The raw
`sqrt(H0^2 + H1^2)` value is retained only to make the already published MV5-I
combined field inspectable; it is not a calibrated multiview fusion method.

## 7. Ranking, identities, and completion records

Within each fold, seed, method, and query, references are ordered by ascending
distance and exact ties are broken by canonical training sample ID. Consecutive
ranks begin at one, and every row records tie size, tie status, and policy.

Every ranking row has an immutable ID binding group, method, query, reference,
distance policy, and its source identity. Every group has exactly five
structured method-completion rows. A failed method records a code and produces
no silent fallback.

## 8. Prospective admission and resource guards

Before full production, MV5-J must admit:

1. the representative `SRA550660` / `20260805` group; and
2. the prospectively identified maximum-pair group from the frozen fold table.

Each admission is computed twice into fresh paths and must be byte-identical.
The independent validator must reproduce all topology rows, directly recompute
spread-out integrated-energy and pseudobulk pairs, and verify partitions,
rankings, ties, hashes, and stop counters.

Full production is limited to two workers, 600 elapsed seconds and 4 GiB
process-tree RSS per group, 7,200 aggregate worker-seconds, and 2 GiB for new
private MV5-J results. The monitor records a guard breach and allows an active
worker to finish; it does not kill R or WSL processes. Any guard breach blocks
authorization of the next stage pending audit.

## 9. Required completion evidence

Completion requires:

- 75 validated groups over 15 folds and five seeds;
- 35,350 biological query-to-training pairs and 176,750 ranking rows;
- 375 completed method groups and zero silent failures;
- exact agreement of H0, H1, and raw-composite rows with MV5-I;
- eligible direct formula oracles for integrated energy and pseudobulk;
- canonical ranks, ties, partitions, method identities, and unique row IDs;
- zero held-out-derived scale fitting;
- zero retrieval-evaluation, clustering, gene-topology, fusion, new-data, and
  biological-outcome jobs;
- byte-identical complete-group and public-assembly repeats;
- a completed-queue resume with zero rebuilt or changed private files;
- resource and storage evidence within the frozen guards; and
- passing focused and full tests plus a clean staged package check.

## 10. Authorization decision

The required evidence passed: 75 groups, 35,350 biological pairs, 176,750
rankings, 375/375 completed method groups, exact H0/H1 agreement, 450 direct
baseline oracles, byte-identical admission and public repeats, and a 150-file
zero-rebuild resume. The authorization disposition is
`approve_separate_prediction_locked_integrated_retrieval_evaluation`.

MV5-J must end with an explicit approve/hold decision for a later integrated
retrieval evaluation. Approval authorizes only a separately specified,
prediction-locked label-opening sprint that verifies the public ranking hashes
before calculating prespecified endpoints. It does not authorize method
selection, comparison with SCT results, clustering, gene topology, fusion, new
data, or manuscript interpretation.
