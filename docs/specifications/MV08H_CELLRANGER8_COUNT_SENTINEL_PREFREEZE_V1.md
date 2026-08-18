# MV8-H Cell Ranger 8.0.1 count-sentinel prefreeze v1

**Status:** prospectively frozen; execution not authorized

**Date:** 2026-08-18

**Parent decision:** D-092

## Purpose and selection

This contract fixes one complete HCA unit before any new count output exists.
The selection rule is label-closed: rank the eight complete units by decreasing
compressed FASTQ bytes and break a tie by original unit order. This chooses the
unique maximum, `HCA_BM_002`, with six FASTQs and 11,249,623,632 bytes. The rule
uses no expression, QC, donor attribute, study/tissue/approach label, or outcome.

Choosing the largest unit is conservative for the first resource and output
admission. It is not a claim that compressed bytes perfectly predict runtime;
it is a prospective, immutable burden proxy available without inspecting the
data scientifically.

## Exact command

```text
cellranger count \
  --id=mv08h_count_sentinel_hca_bm_002 \
  --transcriptome=<verified_custom_reference> \
  --fastqs=<verified_hca_bm_002_fastq_directory> \
  --sample=MantonBM2_HiSeq_9 \
  --chemistry=SC3Pv2 \
  --expect-cells=7000 \
  --include-introns=false \
  --create-bam=false \
  --nosecondary \
  --localcores=4 \
  --localmem=32 \
  --disable-ui
```

The scientific settings are unchanged from D-091/D-092. Four cores and 32 GiB
are a downward-only compute amendment from the earlier 16-core/64-GiB maximum.
The lower allocation may increase elapsed time but does not change chemistry,
reference, reads, cell expectation, exon policy, or feature-barcode counts.

## Immutable identities

- Cell Ranger: 8.0.1; runtime tree SHA-256
  `aafd39e293e0ba9d14dba3896a6aeda077304531a2702d26bda0c62c4688fdf3`.
- Launcher SHA-256:
  `4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6`.
- Custom reference: 19 regular files, zero symlinks, 20,765,871,518
  regular-file bytes, tree SHA-256
  `5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c`.
- FASTQs: the six file names, UUIDs, byte sizes, and SHA-256 values in
  `mv08h-count-sentinel-fastqs.csv`; total 11,249,623,632 bytes.

The launcher must recompute the six file hashes and full reference-tree hash
before execution. It refuses overwrite/resume and requires an absent target.

## Resource and monitoring contract

- exactly 4 local cores and 32 GiB local memory;
- one serialized unit only;
- no competing Cell Ranger, Martian, or STAR process at launch;
- 80-GiB absolute subprocess-tree RSS scientific-acceptance ceiling;
- 200-GiB run-workspace scientific-acceptance ceiling;
- at least 1 TiB free throughout;
- 96-hour elapsed scientific-acceptance ceiling; and
- 30-second subprocess-tree RSS, run-tree, free-space, and elapsed samples.

The 96-hour ceiling replaces the former 24-hour planning value solely to allow
the quartered core allocation to run longer. Monitoring is non-destructive. A
breach invalidates admission and preserves private logs and artifacts for owner
review; the launcher does not kill a process, delete files, retry, or advance.

## Completion and structural validation

A later authorized execution is admissible only after all of the following:

1. exit zero and Cell Ranger pipestance success;
2. complete resource receipt with no gate breach;
3. expected raw and filtered feature-barcode HDF5 and molecule-info outputs;
4. no BAM because `--create-bam=false`;
5. internally consistent sparse matrix dimensions and indices;
6. 33,563 unique reference genes and exact presence of the ordered exact-500
   and common-475 targets;
7. unique private barcodes, with no barcode or expression value published; and
8. a separate auditable closure before QC or another unit.

This structural validation may inspect private matrices only in a later
explicitly authorized closure. It must not evaluate QC thresholds, eligible
cells, PCA, PH, landscapes, clustering, labels, or biological outcomes.

## Firewall and stop conditions

Public evidence may contain only immutable file/runtime/reference identities,
matrix shapes and feature-axis identities, aggregate resources, and run status.
It must contain no absolute private path, expression/UMI value, cell barcode,
donor attribute, QC value or eligibility decision, study/tissue/approach label,
or biological outcome.

This sprint does not authorize or run `cellranger count`. Count execution,
matrix access, QC, remaining units, PCA, PH, persistence landscapes, clustering,
labels, outcomes, and deletion remain closed. A successful future sentinel can
open only a separate structural/QC review and remaining-unit decision.

The downstream landscape definition remains unchanged: cell and gene
observations are separate typed views; H0 and H1 remain separate; essential H0
is excluded; and every consecutive active landscape level is compared with
exact or error-controlled squared-L2 integration, no fixed grid, and no
universal level cap.
