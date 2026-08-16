# MV7-I descriptive topology and clustering prefreeze v1

Date: 2026-08-16

## Purpose

MV7-I evaluates the corrected 124-sample dual-view topology without allowing
tissue, study, approach, or prior outcome values to choose a view, homology
dimension, seed, composite, clustering algorithm, or cluster count. This stage
first freezes and produces label-closed matrices and partitions. Biological
metadata may be joined only after those artifacts are independently validated
and a separate outcome-execution prefreeze is accepted.

## Representations

The primary scientific representations are cell-H0, cell-H1, gene-H0, and
gene-H1. H0 and H1 remain individually reportable. Within each view, the
unweighted composite

`d_combined = sqrt(d_H0^2 + d_H1^2)`

is a secondary descriptive representation only. It cannot replace or suppress
either component. Cross-view fusion remains outside MV7-I.

For every view and seed, also report the pairwise H1 squared-distance fraction
`d_H1^2 / (d_H0^2 + d_H1^2)`, defined as zero only when both components are
zero. This is a contribution measure, not evidence that an H1 feature is a
biological cycle.

## Seed aggregation and uncertainty

The five seed-specific 124 × 124 matrices are the fundamental inputs. For each
representation and sample pair, report all five values plus median, minimum,
maximum, IQR, and raw median absolute deviation (`constant = 1`). Combined
distances and H1 fractions are computed within seed before aggregation.

The median distance matrix is a descriptive heatmap/ordination summary. It does
not replace seed-specific matrices for cluster-count selection or uncertainty.
No favorable seed may be selected.

## Label-closed clustering

For each of the six representations independently:

1. fit PAM to every seed-specific matrix for `k = 2:10`;
2. canonicalize clusters from sorted member-ID signatures;
3. compute all ten pairwise seed adjusted Rand indices for each `k`;
4. use mean pairwise ARI as stability;
5. estimate uncertainty with the five-seed delete-one-seed jackknife;
6. select the smallest `k` within one jackknife standard error of the maximum;
7. if the full grid or five matched axes are unavailable, report
   `no_stable_k`; and
8. at the PAM-selected `k`, fit average-linkage clustering as the sole
   algorithm sensitivity.

PAM is primary because it consumes a general dissimilarity matrix. Average
linkage does not select its own `k`. Ward linkage and direct k-means remain
prohibited. Spectral clustering remains deferred until a distance-to-affinity
and eigengap/stability contract is prospectively validated.

## Populations and interpretation

The 124-sample set is a descriptive population spanning eight tissues. The
existing 90-sample, five-tissue, multi-study population remains the only
cross-study primary-claim population. The 34 added samples from pancreatic
islets, prostate, and substantia nigra describe corpus coverage and topology
location but cannot establish cross-study tissue generalization because each
tissue comes from one study. The three below-250-cell samples remain excluded
and may enter only a separately authorized threshold sensitivity.

After label-closed artifacts are immutable, a separate outcome prefreeze may
authorize complete, fixed-order descriptive association summaries for tissue,
study, and canonical approach. It must use the MV7-E canonical approach field,
report all representations/seeds/algorithms, preserve study and tissue
confounding warnings, and prohibit causal technology, external-generalization,
or single-study cross-study claims.

## Resource, validation, and abort rules

Execution uses one worker, a 3,600-second wall cap, a 4-GiB RSS cap, atomic
private artifacts, public manifests only, deterministic repeat, and immutable
resume. Stop on axis/hash drift, missing matrix entries, nonfinite/asymmetric
distance, label access before authorization, incomplete candidate `k`, unstable
tie handling, resource breach, or any attempt to use outcomes for selection.

Passing this prefreeze authorizes label-closed matrix and clustering artifact
production only. It does not authorize metadata outcomes, manuscript claims,
new distances, cross-view fusion, or external data.
