# MV8-G common-475 paired reference sensitivity prefreeze v1

Date: 2026-08-16

## Purpose and boundary

MV8-G is the prospective, reference-only calculation that determines whether
removing the 25 features excluded by the HCA Cell Ranger 3.0.0 reference
materially changes the accepted 124-sample topology. It is performed before a
decision about downloading or reprocessing HCA raw reads. The accepted
500-gene artifacts are immutable comparators; the new objects use an explicit
`common475_v1` identity and can never pass an MV7-H 500-gene validator.

The calculation is a harmonization sensitivity, not a biological endpoint
analysis and not a proof of equivalence. Tissue, approach, study, class,
cluster, outcome, and earlier favorable or unfavorable results remain closed.

## Exact input and transform

Execution requires independently closed MV8-F recovery evidence: all 90 raw
sources and all 450 recovered primary SCT caches must pass exact identity and
resource validation. The retained 170 added-sample caches and all accepted
500-gene source, PH, and landscape artifacts must rehash against their public
ledgers. Any mismatch stops the stage.

For each fixed seed `20260805:20260809`, use the same 124 biological samples,
384 selected cells per sample, and exact selected-cell identities as the
accepted 500-gene calculation. Restrict the original ordered panel to the 475
published exact stable IDs, recompute reference-only gene centers and scales
across the equally represented 47,616 reference cells, and refit a
30-component cell PCA. Construct both typed views from each standardized
475-by-384 sample matrix:

- cells as points in the shared 30-component PCA space; and
- genes as points with correlation-derived dissimilarity across the same 384
  selected cells.

HCA expression does not fit, tune, select, or stop this transform.

## PH and landscape estimands

Compute complete Vietoris-Rips H0 and H1 for 1,240 typed views. Ripserr is
primary. Only a gene-view process that exceeds its accepted resource cap may
use the exact GUDHI fallback, under the already approved 12-GiB policy.

The landscape definition remains dissertation-aligned and unchanged:

- finite positive-persistence intervals only;
- essential H0 excluded;
- every consecutive active level included;
- H0 and H1 separate;
- exact or error-controlled squared-L2 integration;
- no fixed grid and no level cap; and
- streamed or chunked pair calculations.

The 20 within-475 groups each contain all 7,626 unordered sample pairs, for
152,520 distance rows. A second set of 20 matched-panel groups computes one
exact landscape distance between the 500- and 475-gene diagram for every
sample, seed, view, and dimension, for 2,480 rows. No cross-seed or
cross-component distance is permitted.

## Frozen comparisons

For each of cell-H0, cell-H1, gene-H0, and gene-H1, report all five seeds and
the across-seed distribution of:

1. Spearman and Kendall correlation across the 7,626 matched between-sample
   distances.
2. Normalized stress after the least-squares nonnegative scalar fit of the
   475-gene distances to the 500-gene distances.
3. Mean per-sample top-10 neighbor overlap, with sample ID as the deterministic
   tie breaker.
4. Each matched-panel landscape shift divided by the corresponding seed and
   component median nonzero 500-gene between-sample distance.
5. PAM partition ARI for every `k=2:10` and seed.
6. Agreement between independently selected, label-free panel-specific `k`,
   using the existing five-seed one-standard-error stability rule.
7. At the 500-gene selected `k`, PAM and average-linkage partition ARI between
   panel sizes.

The four components are primary and remain separate. H0/H1 composites, if
shown, are secondary descriptions only.

## Interpretation rule

`high_harmonization_stability` requires every component to have median
Spearman at least 0.95, median top-10 overlap at least 0.80, and median
fixed-500-selected-k PAM ARI at least 0.80. `material_panel_sensitivity`
applies when any component misses at least two of those three thresholds.
All other complete results are `mixed_harmonization_stability`.

High stability supports a defensible common-475 HCA replication and makes
exact-500 raw-read processing an optional strengthening analysis. Material
sensitivity recommends raw-read reprocessing before an external topology
claim. Mixed stability requires component-specific reporting and an owner
decision. This rule does not override the project owner's stated preference
for raw-read processing.

## Resource, repeat, and continuation gates

Source bundles use 1,800 seconds/8 GiB; cell PH 600 seconds/4 GiB; primary gene
PH 1,800 seconds/8 GiB; exact gene fallback 1,800 seconds/12 GiB; and each
landscape group 3,600 seconds/12 GiB. Heavy stages are serial, atomic,
resumable, and zero-retry. The private output aggregate is capped at 4 GiB.

Closure requires complete source and PH identity validation, H0 MST oracles,
balanced H0/H1 interval accounting, immutable resume, one complete seed source
repeat, a balanced PH repeat, maximum-depth landscape repeats, R and Persim
oracles sampled across interval depth, independent reconstruction of all
comparison summaries, and verification that public evidence contains no
expression values, cell barcodes, absolute private paths, or labels.

No FASTQ download, HCA expression calculation, raw-read reprocessing,
biological label access, or manuscript claim is automatically authorized.
