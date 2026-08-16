# MV7-I descriptive outcome prefreeze specification v1

Date: 2026-08-16

## Purpose and stop boundary

This prefreeze fixes the complete descriptive metadata evaluation before any
metadata value is joined to an immutable MV7-I cluster assignment. It may
inspect metadata schemas, category counts, sample-set membership, label-closed
artifact hashes, selected k values, and partition axes. It may not inspect or
compute cluster/metadata association values.

Passing this prefreeze authorizes one descriptive outcome execution only. It
does not authorize method selection, confirmatory inference, causal technology
claims, external generalization, cross-view fusion, new distances, manuscript
claims, or release of private sample-level contingency artifacts.

## Immutable partition inputs

The only clustering input is the validated MV7-I selected-partition artifact
identified by its public SHA-256 manifest. It contains the six frozen
representations, five seeds, PAM and average-linkage partitions, and the
label-free PAM-selected k for each representation. No partition is refit,
filtered, or selected after metadata access. PAM remains primary and average
linkage remains an algorithm sensitivity. Representations and algorithms are
reported separately and never ranked by an outcome.

## Populations and endpoints

Two fixed populations are evaluated:

1. `full124_descriptive`: all 124 retained samples and all eight tissues. Its
   tissue, study, and canonical-approach endpoints are descriptive only.
2. `primary90_context_restriction`: the accepted 90 samples from five tissues
   represented across studies. These are the full-124 partitions restricted to
   the 90 sample IDs; they are not 90-sample refits and cannot replace the
   accepted blocked MV5-S analysis.

For the 90-sample restriction, tissue and study endpoints are scheduled.
Canonical approach is retained in the endpoint registry as
`structurally_not_estimable_single_class` because all 90 samples are scRNA-seq;
no association metric is computed for it.

## Metadata authority

Tissue and study come from the accepted MV7-D sample reconciliation.
Sequencing approach comes only from the MV7-E `canonical_approach` field.
Historical expression-heuristic approach values are prohibited. The metadata
tables and the partition artifact must each have one unique row per expected
sample axis before execution; sample-set agreement may be checked, but
metadata and cluster values are not joined during prefreeze.

## Metrics and seed summaries

For every scheduled population, label axis, representation, and algorithm,
compute the following separately for each of the five topology seeds:

- Hubert-Arabie adjusted Rand index; and
- normalized mutual information `MI(U,V) / max(H(U), H(V))`.

These reuse the independently tested MV5-S definitions. Degenerate or
nonfinite metrics remain explicit and are never replaced. For each five-seed
unit, report all five values plus mean, median, minimum, maximum, and
delete-one-seed jackknife technical standard error. Seeds are technical
realizations, not biological replicates.

There are five schedulable population/label endpoints, six representations,
two algorithms, and two metrics: 120 summary units and 600 seed-level metric
rows. No p-value, null permutation, multiplicity adjustment, winner, rank, or
favorable subset is authorized.

## Confounding and interpretation

- Tissue alignment is biological description, not evidence that topology is
  tissue-causal or free of study effects.
- Study alignment is a technical/cohort association, not proof of a batch
  mechanism.
- In the full 124, all six snRNA-seq samples are substantia nigra samples from
  study SRA850958. Canonical approach is therefore perfectly nested in that
  tissue/study combination. Its alignment is reportable only as a confounded
  proxy and cannot identify a technology effect.
- The three added tissues are each single-study. Their full-124 tissue
  alignment cannot establish cross-study generalization.
- The primary-90 restriction is contextual because its clusters were fit on
  all 124 samples. Accepted MV5-S blocked results remain the primary
  cross-study evidence.
- H1 alignment does not establish a biological cycle or mechanism. H0 and H1
  remain individually reportable; secondary composites cannot replace them.

## Outputs and privacy

Public outcome outputs may contain unit identifiers, population, label-axis
names, representation, algorithm, metric, seed, estimates, uncertainty,
counts, and status. They must report all scheduled units in fixed order and
must not contain per-sample metadata values, contingency cells, predicted
labels, favorable-result flags, rankings, or claims.

The exact sample/label join and long contingency tables remain private atomic
artifacts. Their hashes and row counts may be public. Deterministic repeat,
independent recomputation, immutable resume, and a label-leakage/public-schema
check are required before the results can inform figures or the claim map.

## Resource and failure policy

Execution uses one worker, a 900-second wall cap, a 2-GiB process-tree RSS cap,
atomic private artifacts, deterministic repeat, and immutable resume. It stops
on source-hash drift, axis mismatch, duplicate or missing metadata, unexpected
category support, partition mutation, missing seed/unit rows, nonfinite values
not explicitly marked degenerate, public label leakage, resource breach, or
any outcome-driven selection.
