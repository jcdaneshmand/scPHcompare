# MV5-N label-closed sample-clustering and complete-matrix resource gate specification v1

## Document control

| Field | Value |
|---|---|
| Contract ID | `mv05n_label_closed_clustering_resource_gate_v1` |
| Date | 2026-08-10 |
| Status | Frozen before real admission execution |
| Source revision | `e5db294ee185326e182b19aca76701df260038ca` |
| Primary unit | biological sample |
| Fit scope | training samples within one LOSO fold and one cell-subsample seed |
| Required stop | before full matrix production, clustering outcomes, robustness, gene, fusion, new-data, or optimization work |

## 1. Purpose

MV5-N converts the dissertation-era sample-clustering idea into a leakage-safe,
auditable inductive design. It freezes the clustering and held-out prediction
rules, constructs immutable identities for every missing training-training
distance request, and profiles a small outcome-closed subset. It does not run
the full matrices and it does not ask whether clusters recover tissue, study,
or approach.

The sprint uses the corrected cell-as-observation persistence diagrams already
accepted in MV5-D3 and MV5-H. It does not rehabilitate legacy gene-as-point
results or the original transductive/oracle-k analyses.

## 2. Landscape definition

Every topological request uses the revised dissertation-aligned definition:

- finite positive-persistence intervals only;
- the one essential H0 interval excluded at landscape construction;
- H0 and H1 constructed and compared separately;
- every consecutive active persistence-landscape level retained;
- exact piecewise-linear integration over critical-pair segments;
- no fixed level cap and no uniform evaluation grid; and
- zero reported numerical integration error for the exact engine.

The unscaled `sqrt(H0^2 + H1^2)` distance may be reconstructed only as a
descriptive secondary matrix. It cannot replace the separate H0 and H1
components because no outcome-closed component scaling rule was fitted.

## 3. Inductive split and required matrices

For each of 15 leave-one-study-out folds and seeds `20260805:20260809`:

1. fit all transforms using training samples only;
2. calculate a complete symmetric distance matrix among training samples;
3. choose `k` from partitions generated across the five frozen seeds, without
   labels;
4. fit the seed-specific training partition at the selected `k`; and
5. assign each held-out sample using its already accepted query-to-training
   distances.

Held-out samples never participate in transform fitting, distance scaling,
`k` selection, medoid selection, dendrogram construction, or cluster naming.
Distances among held-out samples are neither required nor authorized.

The exact missing scope per representation is:

| Quantity | Count |
|---|---:|
| unordered training-training sample pairs | 262,675 |
| H0 request rows | 262,675 |
| H1 request rows | 262,675 |
| H0/H1 request rows | 525,350 |
| already accepted query-to-training pairs | 35,350 |

Each pair is canonicalized as `first_sample_id < second_sample_id` and bound to
fold, seed, representation, source record keys, diagram hashes, result-file
hashes, homology dimension, and landscape definition. The gate instantiates,
validates, and digests every identity group, then discards non-admitted rows.
Public global/group/chunk inventories retain counts, first/last request IDs,
source hashes, and SHA-256 identity-set digests. A later authorized production
run can stream the exact rows again from the same generator without adding very
large generated manifests to Git or paying pointless maximum-compression cost.

## 4. Distance panel

The same clustering rules apply independently to both representations and the
following distance methods:

| Distance | Role |
|---|---|
| exact all-active-level cell H0 landscape L2 | confirmatory topology component |
| exact all-active-level cell H1 landscape L2 | confirmatory topology component |
| raw `sqrt(H0^2 + H1^2)` | descriptive secondary |
| cell-distribution energy in the same 384-cell/30-coordinate space | matched cell baseline |
| training-standardized 500-gene pseudobulk Euclidean | context baseline |

SCT energy uses the MV5-D1 shared training-fit PCA coordinates. Integrated
energy uses the MV5-G inductive integrated coordinates. Pseudobulk uses the
same training-selected and training-standardized panel in both representation
comparisons and remains an identity/context control.

## 5. Primary PAM contract

PAM/k-medoids is primary because it consumes a general dissimilarity without
pretending the matrix is a Euclidean feature table.

For each fold, representation, and distance method:

- candidate `k` is every integer from 2 through `min(10, n_training - 1)`;
- PAM is fit separately to each of the five seed-specific training matrices;
- stability at each `k` is the mean of the ten pairwise adjusted Rand indices
  between seed partitions over the identical training-sample axis;
- Monte Carlo uncertainty is the leave-one-seed-out jackknife SE;
- the best row is the first maximum in ascending `k`;
- the threshold is `best_mean - best_SE`; and
- selected `k` is the smallest candidate whose mean stability reaches the
  threshold.

This is the already implemented `mv05_select_stable_k_v1` one-SE rule. The
MV5-N wrapper rejects missing seeds, mismatched sample axes, or an incomplete
candidate grid.

Cluster labels are canonical rather than inherited from an algorithm: clusters
are ordered by their sorted member-ID signatures and numbered `1..k`. This
makes assignments stable to arbitrary source label permutations.

## 6. Held-out assignment

For PAM, each held-out sample is assigned to the nearest frozen training
medoid. Exact distance ties are resolved by lexicographically smallest medoid
sample ID. No held-out sample can become a medoid.

Predictions must retain fold, seed, representation, distance method, selected
`k`, training-matrix identity, partition identity, held-out sample ID, assigned
canonical cluster, assignment distance, winning medoid, and tie disposition.
They must be made immutable before any biological or technical label is joined.

## 7. Average-linkage sensitivity

Average-linkage hierarchical clustering is the sole currently eligible
clustering sensitivity. It uses the same training matrix and the `k` selected
by the primary PAM stability rule; it does not select its own `k` from labels
or outcomes.

A held-out sample is assigned to the training cluster with minimum mean
dissimilarity from that sample to all members of the frozen cluster. Exact ties
are resolved by the lexicographically smallest sorted training-member
signature. This is the out-of-sample analogue of average linkage and avoids
rebuilding a transductive dendrogram containing held-out samples.

## 8. Ineligible and excluded methods

- Spectral clustering remains **ineligible**. The legacy Gaussian-kernel code
  lacks a separately frozen and validated distance-to-affinity, bandwidth,
  eigengap, disconnected-graph, sign, and held-out extension contract.
- Ward linkage remains excluded for an arbitrary topological dissimilarity
  unless a valid Euclidean embedding is independently demonstrated.
- K-means on a distance matrix is excluded because a distance matrix is not a
  feature representation.
- Known tissue/study/approach class counts are prohibited for primary `k`.
  Oracle `k` may only return later as an explicitly historical sensitivity.

## 9. Label firewall

The pair generator, matrix builder, `k` selector, training partition, and
held-out assignment accept no tissue, approach, class, label, outcome, or
biological/technical-label columns. Study IDs are used only to define the
already frozen outer folds; they are not a clustering outcome.

During MV5-N:

- outcome-label state is `closed`;
- biological outcomes computed is `FALSE`;
- clustering outcome jobs executed is zero;
- ARI/NMI against tissue/study/approach is prohibited; and
- neither a method, representation, component, fold, nor tissue may be selected
  from already known retrieval results.

Any later label-open evaluation must keep two estimands separate:

1. **Descriptive training-partition alignment** compares a frozen training
   partition with labels for those same training samples. It describes the
   partition already fit and does not estimate out-of-study generalization.
2. **Inductive held-out generalization** compares the immutable held-out
   assignments with labels opened only after prediction lock. This is the
   relevant out-of-study generalization estimand.

The two may not be pooled, substituted for one another, or used to retune
`k`, the clustering method, representation, distance component, or tissue.
Neither estimand is computed in MV5-N.

## 10. Bounded real admission

Three folds are selected solely from training-matrix size, with canonical fold
ID resolving ties:

- minimum training size;
- median training size (representative); and
- maximum training size.

Only seed `20260805`, the canonical first frozen seed, is admitted. For each
profile, representation, and H0/H1 dimension, 32 requests are selected at
equally spaced ranks across the canonical request-ID ordering. The admission
therefore contains exactly 384 landscape rows in 12 resource groups. This is a
computational profile, not a biological sample and not a partial matrix on
which clustering may be run.

Each group has a 900-second elapsed guard. Existing execution wrappers must
enforce a 4-GiB process-tree RSS cap. The run records elapsed and pair-operation
seconds, peak RSS, output bytes, exact status, active levels, segments, and
critical points.

Correctness requires:

1. all 384 requested rows complete exactly once;
2. H0/H1 diagram and file hashes match accepted sources;
3. one eligible H0 and H1 result per profile/representation agrees with an
   independent R exact-linear-segment oracle within `1e-10` absolute error;
4. a second clean output directory is byte-identical;
5. an in-place rerun validates resume without rewriting outputs;
6. synthetic PAM, average-linkage, stable-k, canonical-label, tie, matrix, and
   baseline formula fixtures pass; and
7. no outcome column or outcome job appears.

## 11. Resource decision

Admission measurements are projected to all 1,050,700 topological component
rows across both representations, then combined with measured energy,
pseudobulk, validation, I/O, and deterministic-repeat costs. Full production
is not authorized by this specification.

A later production authorization may pass only if:

- all correctness and repeat checks pass;
- projected aggregate work, including reserve, is at most 21.6 worker-hours;
- projected additional private storage is at most 10 GB;
- no profile exceeds 900 seconds or 4 GiB;
- no large profile shows an unexplained nonlinear slow tail; and
- execution can be chunked, resumed, and validated without dense landscape
  matrices in memory.

If the projection fails, the next action is a numerical-equivalence-preserving
engine optimization/resource sprint (possibly Rust) or a narrower approved
execution design—not silent truncation of landscape levels, H1, samples, or
folds.

## 12. MV5-N acceptance and stop

MV5-N passes when the contract, exact pair identities, public inventories,
bounded real admission, independent oracle, deterministic repeat/resume,
synthetic clustering/baseline fixtures, full test suite, package check, and
public artifact hashes all validate.

It must then stop. A separate prospective authorization is required before
full training matrices, cluster predictions, label opening, ARI/NMI, robustness
execution, additional integration, gene topology, fusion, new data,
optimization, or manuscript claim promotion.
