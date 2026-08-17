# MV8-H Cell Ranger 8.0.1 `mkref` execution prefreeze v1

**Status:** prospectively frozen for one custom-reference build

**Date:** 2026-08-17

**Parent decision:** D-091

**Parent contract:** `MV08H_CELLRANGER8_RUNTIME_REFERENCE_AMENDMENT_PREFREEZE_V2.md`

## Purpose

This execution amendment applies the project owner's instruction to use fewer
resources if needed while another independent Cell Ranger workload is active.
It changes compute allocation only. The runtime, FASTA, GTF, gene panels,
scientific estimands, landscape definition, privacy firewall, and downstream
stop rules remain exactly as frozen in the parent contract.

The already verified Ensembl-93 archive was atomically decompressed before
execution. The resulting input is exactly 3,151,425,857 bytes with SHA-256
`78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c`.

## Downward-only resource amendment

The single authorized build uses four cores and 32 GiB instead of the parent
maximum of 16 cores and 64 GiB:

```text
cellranger mkref \
  --genome=GRCh38_Ensembl93_target_complete_33563 \
  --fasta=<verified_uncompressed_ensembl93_primary_assembly> \
  --genes=<verified_target_complete_33563_gtf> \
  --nthreads=4 \
  --memgb=32 \
  --localcores=4 \
  --localmem=32 \
  --disable-ui
```

The 50-GiB output-tree ceiling and 1-TiB free-space floor are unchanged. The
launcher samples the complete subprocess-tree RSS, run-tree bytes, and free
space every 30 seconds. Monitoring is non-destructive: it does not kill a
process, remove an artifact, or compete with another Cell Ranger run. Any
observed resource breach prevents scientific acceptance and requires owner
review.

## Input and execution gates

Before starting, the launcher must verify:

- Cell Ranger 8.0.1 and launcher SHA-256
  `4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6`;
- the uncompressed Ensembl-93 FASTA size and SHA-256 above;
- the 1,099,737,654-byte target-complete GTF and SHA-256
  `e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182`;
- at least 1 TiB free at launch; and
- absence of both the run root and final reference target.

The run is acceptable only after exit zero, the Cell Ranger completion state,
absence of a resource breach, and independent complete-tree, metadata, FASTA,
GTF, and exact-500 feature validation. Private logs, licensed runtime files,
expression data, barcodes, donor attributes, labels, and outcomes are not
published.

## Authorization boundary

This prefreeze authorizes exactly one `mkref` execution. It does not authorize
`cellranger count`, QC, matrix access, PCA, PH, persistence landscapes,
clustering, label access, biological-outcome analysis, or deletion. A successful
reference closure may open only a separate prospective one-unit count-sentinel
prefreeze.
