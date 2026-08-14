# MV5-S prediction-locked clustering-outcome execution specification v1

Date frozen: 2026-08-10

Accepted contract base: `3c3271e`

Outcome state at freeze: not computed

## Purpose

MV5-S implements and executes the committed MV5-R 2,400-unit evaluation queue
without changing its scientific choices. Execution is permitted only after all
18 MV5-R source identities, the exact 90-sample label join, the 150 analysis
groups, five seeds, two algorithms, and zero-outcome counters validate.

This specification resolves implementation details that MV5-R intentionally
left below the scientific-contract level. These details are frozen and tested
before real ARI, NMI, or balanced accuracy is calculated.

## Exact metric definitions

Adjusted Rand index uses the Hubert-Arabie chance correction over the complete
contingency table and must agree with `mclust::adjustedRandIndex` to numerical
tolerance on independent fixtures.

Normalized mutual information uses
`aricode::NMI(first, second, variant = "max")` under aricode 1.0.3:

`NMI = MI(U,V) / max(H(U), H(V))`.

An independent contingency-table implementation must reproduce the library
value. A zero entropy denominator or other non-finite metric is retained as
`degenerate_label_or_partition_metric_not_identifiable`; it is not replaced.

Balanced accuracy is the unweighted mean of per-class correctness. This is
equivalent to the mean recall across true classes for the frozen single-label
prediction task.

## Training units

Each of the 1,800 training units evaluates all five frozen partitions for its
analysis group, algorithm, label axis, and metric. It writes one public row per
seed and one unit summary with the five-seed mean and leave-one-seed-out
jackknife technical SE. Fold-specific cluster IDs are never pooled. The full
15-fold distribution is retained; no fold is treated as an independent
biological replicate and no p-value is computed.

## Held-out units

Each of the 600 held-out units uses training labels only to assign a plurality
label to each frozen training cluster. Exact count ties use lexical label order.
That fixed map is applied to the already frozen held-out cluster assignment.
Changing a held-out true label must not change its prediction in a leakage
oracle.

The private unit artifact contains sample ID, study, true label, predicted
label, correctness, cluster, plurality count, and tie size for all five seeds.
Public seed and unit tables contain only identifiers, estimates, counts, and
status; they never contain true or predicted label values.

Across folds, correctness is first averaged over the five technical seeds for
each of the 90 samples. Tissue balanced accuracy is the unweighted mean of the
five tissue-specific sample means. Approach balanced accuracy is the unweighted
mean of the two approach-specific sample means. Class-level rows remain
private. The public table contains all 20 representation-distance-algorithm
contexts for both endpoints, without ordering, winner labels, or selection.

## Bootstrap details

Use 2,000 replicates, seed `20260810`, R's Mersenne-Twister generator,
`Inversion` normal generator, and `Rejection` sample generator. The same frozen
block-count matrix is applied to every method context within an endpoint.

For tissue, each study is tissue-homogeneous. Within each tissue, resample its
study IDs with replacement using the observed number of studies, concatenate
all five strata, and recompute the sample-weighted class means and tissue macro
average.

For approach, three studies contain both approaches. Resample all 15 study IDs
globally with replacement, preserve every selected study as an intact block,
and recompute the two approach means and their macro average. Never split a
mixed study. A replicate missing either approach is non-estimable. Report a
type-7 95% percentile interval only if at least 95% (1,900/2,000) replicates are
estimable; otherwise retain the point estimate and mark
`bootstrap_support_insufficient`.

## Atomic execution and resume

Run one worker. Each queue unit writes a private RDS artifact to a temporary
path, atomically renames it, then writes a status RDS containing the artifact
SHA-256, bytes, elapsed time, source-freeze identity, and terminal state.
Existing completed pairs are reused only when identity and hash match. A
missing half, mismatch, or stale identity aborts; it is not silently rebuilt.

The full first pass must remain below 300 seconds and 4 GiB per unit, two
aggregate worker-hours, and 1 GiB of public output. A second full pass must
reuse all 2,400 units with unchanged paths, hashes, sizes, and timestamps.

## Required validation

Before real execution:

- validate 18/18 frozen MV5-R source hashes and exact queue identity;
- pass independent ARI and max-NMI library comparisons;
- pass plurality lexical-tie and held-out-label leakage oracles;
- pass deterministic tissue-stratified and global study-block fixtures;
- confirm three mixed-approach studies and tissue-homogeneous studies;
- confirm public schemas contain neither true/predicted labels nor label values.

After execution, an implementation-independent validator reconstructs every
training metric, every held-out prediction and sample correctness value, both
bootstrap block matrices and intervals, every aggregate, and all unit hashes.
A clean deterministic production repeat and immutable resume are mandatory.

## Reporting and stop boundary

Report every prespecified context, failure, denominator, and interval. PAM
remains primary; average linkage remains sensitivity. H0/H1 roles and the
energy/pseudobulk roles remain as frozen. Do not calculate p-values, rank or
select methods, tune k, remove folds/seeds, promote a result, or alter the
retrieval conclusions.

MV5-S stops after validated complete reporting. Spectral promotion, robustness
expansion, gene topology, cell/gene fusion, new data, optimization/Rust,
package-default changes, manuscript claims, PDF/reviewer/private-cache
tracking, `example_run.r`, and pushing remain prohibited.
