# MV8-E reference reconciliation and 500-versus-475 sensitivity prefreeze v1

Date: 2026-08-16

## Purpose and decision boundary

MV8-E resolves why the frozen 500-gene reference panel intersects the admitted
HCA count matrices at exactly 475 genes, then freezes the reference-only
500-versus-475 sensitivity required before choosing whether to reprocess HCA
raw reads. The project owner prefers raw-read reprocessing if it is justified,
but no FASTQ download or custom quantification reference is authorized in this
stage.

This is a methodological sensitivity. HCA expression, tissue, study, approach,
class, cluster outcomes, and previous favorable or unfavorable results may not
select a gene, tune a transform, alter a persistence calculation, or stop the
run. Passing this prefreeze authorizes only recovery of accepted private
normalization caches and a label-closed reference sensitivity calculation.

## Exact annotation reconciliation

The HCA feature ID/name axis must be compared directly with the official
Ensembl release 93 GRCh38 GTF after applying the exact Cell Ranger 3.0.0
documented gene-biotype allow-list. Reconciliation passes only if:

- the GTF SHA-256 is
  `810e3bb63bb24bd5a005b14397d69280dda34c41a612fb86a18f3c4836fce57d`;
- the unfiltered GTF contains 58,395 gene rows;
- the documented filter produces exactly 33,538 genes;
- filtered IDs and names equal the HCA axis in the same order;
- all 500 panel stable IDs occur in the unfiltered GTF;
- exactly 475 occur in the filtered reference; and
- all 25 excluded IDs remain valid Ensembl IDs but have biotypes outside the
  documented Cell Ranger 3.0.0 allow-list.

The conclusion is structural: an identifier or symbol crosswalk cannot recover
counts for features that were absent from the quantification reference. Symbol
substitution, zero imputation, duplicate aggregation, and replacement genes are
prohibited.

## Frozen common panel

`common475_v1` is the original ordered 500-gene panel restricted to exact
stable Ensembl IDs present in the HCA Cell Ranger 3.0.0 reference. It preserves
the original relative order. Neither HCA expression values nor biological
labels contribute to the restriction. The published panel ledger and its
newline-delimited feature-axis SHA-256 are binding for every later reference
and HCA calculation.

The accepted 500-gene result remains the main existing-data descriptive
analysis. The 475-gene result is a harmonization sensitivity and the candidate
external-validation object; it does not retroactively replace the 500-gene
analysis.

## Cache-recovery gate

The 170 accepted added-sample SCT caches and five 500-gene source bundles are
present and must rehash against their accepted ledgers. The 450 older
primary-sample SCT caches were intentionally private and are no longer present.
They may be reconstructed only from the retained individual sample Seurat
files using the already accepted MV5-D0 v2 pipeline, exact selected-cell seeds
`20260805` through `20260809`, one or two workers, and the frozen software
environment.

Every reconstructed sample/seed cache must equal its accepted public manifest
in logical cache key, selected-cell SHA-256, payload SHA-256, byte size, and
file SHA-256. Any mismatch stops MV8-E before a 475-gene center, scale, PCA,
typed view, PH diagram, or landscape distance is calculated. Environment drift
must be reported; it may not be normalized away or accepted approximately.

## Paired 500-versus-475 calculation

The sensitivity uses exactly 124 reference samples, five fixed selected-cell
seeds, 384 cells per sample, and 30 cell-PCA components. For each seed:

1. Read the exact accepted SCT matrices and verify their hashes.
2. Restrict the original ordered panel to `common475_v1`.
3. Recompute reference-only per-gene center and scale from all 124 equally
   represented samples. HCA values do not contribute.
4. Refit the 30-component cell PCA on the 47,616 standardized reference cells.
5. Construct paired cell and gene topology views from the same standardized
   475 by 384 sample matrix.
6. Compute complete Vietoris-Rips H0 and H1 for all 1,240 typed views.
7. Compute all 20 within-seed view/dimension distance groups: 7,626 unordered
   sample pairs per group and 152,520 component rows total.

The accepted 500-gene source, PH, and landscape artifacts are immutable
comparators and are not recalculated unless an identity validator fails. A
475-gene object is never mixed with a 500-gene PH or landscape object.

## Persistence-landscape contract

The dissertation-aligned definition is unchanged for both panel sizes:

- finite positive-persistence intervals only;
- essential H0 excluded;
- every consecutive active landscape level included;
- H0 and H1 calculated and reported separately;
- exact or error-controlled squared-L2 integration;
- no primary fixed grid and no universal level cap; and
- streamed or chunked pair calculation rather than dense landscape matrices.

Ripserr remains the primary PH engine. Exact GUDHI is a resource fallback only
under the accepted 12-GiB policy. The accepted exact Rust landscape kernel is
eligible for hash-verified private production; grouped exact Persim remains the
portable fallback and R remains the correctness oracle.

## Frozen sensitivity summaries

For each of cell-H0, cell-H1, gene-H0, and gene-H1, report all five seeds and a
complete across-seed summary for:

- Spearman and Kendall correlation of all 7,626 pairwise distances;
- scale-free normalized stress after the least-squares nonnegative scalar fit;
- mean per-sample top-10 neighbor overlap;
- matched-sample 500-to-475 landscape distance relative to the median nonzero
  500-gene between-sample distance;
- complete `k = 2:10` PAM partition adjusted Rand indices between panel sizes;
- agreement of the separately label-free selected `k`; and
- fixed-500-selected-`k` PAM and average-linkage partition agreement.

H0/H1 composites may be supplied as secondary descriptive context only. No
view, dimension, seed, or clustering solution may be selected using a
biological label.

For planning—not equivalence testing—the result is classified as:

- `high_harmonization_stability` when every primary component has median
  Spearman at least 0.95, median top-10 overlap at least 0.80, and median
  fixed-k PAM ARI at least 0.80;
- `material_panel_sensitivity` when any component misses at least two of those
  three thresholds; or
- `mixed_harmonization_stability` otherwise.

These thresholds determine the strength and wording of a harmonized-panel
validation claim; they do not prove biological equivalence and do not override
the owner's ability to choose exact-500 raw-read reprocessing.

## Resources, repeatability, and stop rules

Execution is serial at every PH/landscape group boundary with no automatic
retry. Cache recovery retains the accepted 1,800-second/8-GiB/40-GB limits.
Source bundles use 1,800-second/8-GiB caps; cell PH uses 600 seconds/4 GiB;
gene PH uses 1,800 seconds/12 GiB under the fallback policy; and landscape
groups use 3,600 seconds/12 GiB. Private outputs are atomic and resumable.

Required validation includes all source hashes, the complete PH/MST oracle,
balanced cell/gene and H0/H1 counts, exact landscape invariants, R/Persim
oracles sampled across interval depth, one complete seed source repeat, a
balanced PH repeat, a maximum landscape-group repeat, independent result
reconstruction, and immutable resume.

Any cache identity drift, nonfinite coordinate, invalid diagram, incomplete
level accounting, resource breach, label access, incomplete pair axis, or
repeat mismatch stops the stage. No HCA FASTQ download begins automatically.

## Raw-read decision after sensitivity

After the sensitivity is independently closed:

- high stability supports the 475-gene HCA analysis as a defensible
  harmonized-panel external replication; exact-500 raw reprocessing remains an
  optional strengthening analysis;
- material sensitivity recommends exact-500 HCA raw-read reprocessing before
  an external topology claim; and
- mixed stability requires component-specific reporting and an owner decision.

Any raw-read path must then receive a separate prospective contract freezing
the custom GTF/filter, genome build, quantifier and version, multimapping and
overlap policy, chemistry/read structure, expected 48 FASTQs (approximately
79.19 GiB), compute/storage caps, output identity, and comparability language.
It must explicitly address the remaining quantification-pipeline asymmetry with
the original 124-sample corpus.
