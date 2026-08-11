# MV5-R prediction-locked clustering-outcome prefreeze specification v1

Date frozen: 2026-08-10

Accepted clustering-artifact base: `f16321c`

Outcome state: not computed

## Purpose and boundary

MV5-R defines the complete evaluation contract for the already immutable MV5-Q
sample clusters and held-out assignments before any clustering outcome is
calculated. It binds the external metadata by SHA-256 but does not copy it into
Git and does not publish sample labels.

This prefreeze authorizes contract implementation and synthetic tests only. It
does not authorize real ARI, NMI, balanced accuracy, confidence intervals,
method ranking, or biological interpretation.

## Frozen sources

The label source is the historical Git-ignored
`joined_metadata_cellcounts.csv`, supplied to a later runner as an explicit
external path. Its required SHA-256 is
`e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`.
It contains 124 unique samples and 18 studies; the previously frozen five-
tissue candidate subset contains 90 samples in 15 studies. The execution runner
must abort before reading labels if this identity changes.

Clustering sources are the committed MV5-Q selected training partitions,
held-out assignments, group completions, stability summaries, queue, and
artifact manifest. All 150 analysis groups, five seeds, two algorithms, 90
sample IDs, and source hashes must validate before the label file may be read.

## Algorithms and distance roles

PAM (`pam_stability_k_v1`) is primary. Average linkage
(`hclust_average_v1`) is sensitivity only. Both reuse the PAM-selected k. No
refit, oracle k, spectral method, seed removal, or method change is permitted.

H0 and H1 remain separate confirmatory topology components. The raw H0/H1
Euclidean combination is descriptive, energy is the matched cell baseline, and
pseudobulk is context. These roles cannot be promoted from outcome results.

## Training-partition alignment endpoints

For every fold, seed, analysis group, and algorithm, calculate sample-level ARI
and NMI against tissue, study, and approach labels. These six endpoints describe
alignment of overlapping training partitions; they are not held-out
generalization estimates because LOSO training sets overlap heavily.

Write all five seed values, their mean within fold, and a leave-one-seed-out
jackknife technical SE. Report the complete fold distribution. Do not compute
p-values, treat seeds as independent observations, or pool fold-specific cluster
IDs.

ARI is the primary external-agreement metric within this secondary clustering
analysis. NMI is supportive. Study and approach agreement are technical
diagnostics, not biological success criteria.

## Held-out generalization endpoints

For tissue and approach separately, learn a label map from each frozen training
partition only: each cluster receives its plurality training label, with exact
count ties resolved by lexicographically smallest label. Apply that fixed map to
the already frozen held-out cluster assignments.

Each sample is scored only in its held-out-study fold. Average the five seed
correctness indicators within sample, then macro-average equally across tissue
for the biological endpoint and equally across approach for the technical
endpoint. The tissue and approach endpoints are reported separately and are
never combined into an overall score.

Use 2,000 deterministic percentile bootstrap replicates with study as the
resampling block. For tissue, studies are tissue-homogeneous, so resample study
blocks within tissue strata. For approach, three studies contain both approaches;
resample all 15 study blocks globally and recompute the approach-macro statistic
in each replicate. Do not split mixed studies or falsely stratify a study into
two independent blocks. Report 95% intervals and support counts. Do not compute
a p-value. If a bootstrap stratum or replicate has insufficient independent
study support, retain the point estimate and mark interval support insufficient
rather than substituting another analysis.

Held-out study prediction is prohibited because each outer-fold study label is
absent from training. ARI/NMI of pooled held-out cluster IDs is also prohibited:
canonical cluster numbers are fold-specific and not exchangeable across folds.

## Multiplicity and reporting

All clustering outcomes are secondary or sensitivity and remain outside F1–F3.
There is no clustering p-value family and no outcome-driven selection. Publish
complete endpoint tables for every representation, distance, algorithm, fold,
and seed, including failures and denominators. A favorable clustering outcome
cannot override the previously completed retrieval conclusions or promote the
descriptive composite, pseudobulk, or average-linkage sensitivity.

## Missingness and failure states

Every queued unit must finish as one of:

- `completed`;
- `missing_or_drifted_label_source`;
- `sample_axis_mismatch`;
- `incomplete_seed_set`;
- `missing_training_label_class`;
- `degenerate_label_or_partition_metric_not_identifiable`;
- `bootstrap_support_insufficient`;
- `immutable_source_hash_mismatch`;
- `software_failure`.

Failures are retained and never replaced by another fold, seed, k, algorithm,
distance, label, or favorable result.

## Prospective queue and execution gate

The immutable queue contains 2,400 units: 150 analysis groups × two algorithms
× eight endpoints. A later execution sprint must first revalidate the queue,
external label hash, MV5-Q artifact hashes, complete axes, and zero prior outcome
counters. It must write atomic status/output pairs and pass deterministic repeat,
immutable resume, public safety, and independent metric oracles.

## Stop boundary

MV5-R stops at a committed prefreeze. It does not compute or open any real
clustering outcome. Method/representation/component/fold/tissue selection,
spectral promotion, robustness execution, gene topology, fusion, new data,
optimization/Rust, package-default changes, manuscript claims, PDFs, reviewer
material, private-cache tracking, `example_run.r`, and pushing remain prohibited.
