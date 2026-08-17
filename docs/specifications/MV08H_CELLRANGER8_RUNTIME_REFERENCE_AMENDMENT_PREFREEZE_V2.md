# MV8-H Cell Ranger 8.0.1 runtime/reference amendment prefreeze v2

**Status:** prospectively frozen; only custom-reference construction is opened

**Date:** 2026-08-17

**Parent decisions:** D-086 through D-090

**Public evidence:** `docs/audits/mv08h-cellranger8-runtime-reference-prefreeze-v2/`

## Purpose and amendment boundary

The project owner confirmed that Cell Ranger was already installed and
authorized this goal after review of the exact available versions. The local
environment contains Cell Ranger 7.1.0 and 8.0.1; the configured modern runtime
is 8.0.1. No exact Cell Ranger 3.0.x installation was found in the expected
locations.

This v2 prefreeze deliberately replaces only the prospective counter-runtime
clauses of the MV8-H v1 acquisition contract. It does not rewrite or invalidate
the historical v1 evidence. The following remain unchanged:

- the 48-file, eight-unit FASTQ cohort and its complete checksum closure;
- the Ensembl release-93 primary-assembly FASTA identity;
- the deterministic 33,563-gene target-complete GTF identity;
- the exact ordered 500-gene primary panel and ordered common-475 sensitivity;
- MV8-D QC, five deterministic 384-cell selections, and seeds;
- separate cell and gene topology, complete Vietoris–Rips H0/H1, and the
  accepted all-active-level landscape definition; and
- the label, outcome, privacy, resource, and no-deletion firewalls.

The admitted historical HCA matrices remain identified as Cell Ranger 3.0.0
outputs. They are comparators; this amendment does not relabel them as 8.0.1 or
claim exact runtime reproduction.

## Exact runtime identity

The selected runtime is the installed Cell Ranger 8.0.1 distribution. The
public evidence binds:

- the reported and embedded 8.0.1 version;
- a deterministic content digest over every regular file and symlink target in
  the installed distribution;
- the launcher, STAR, samtools, and Python identities; and
- the locally exposed `count` and `mkref` controls required by this contract.

The resulting complete-tree SHA-256 is
`aafd39e293e0ba9d14dba3896a6aeda077304531a2702d26bda0c62c4688fdf3`.

Only relative component paths and content identities are published. The local
installation path and licensed vendor payload are not committed.

Cell Ranger 8.0.1 supports `SC3Pv2`, but modern Cell Ranger includes intronic
reads by default. The prospective count command therefore fixes
`--chemistry=SC3Pv2` and `--include-introns=false`. It also fixes
`--create-bam=false` and `--nosecondary`; these omit unneeded storage-heavy
outputs without changing the feature-barcode count estimand used here.

## Reference and command contract

`cellranger mkref` is frozen as:

```text
cellranger mkref \
  --genome=GRCh38_Ensembl93_target_complete_33563 \
  --fasta=<verified uncompressed Ensembl93 primary assembly> \
  --genes=<verified target-complete GTF> \
  --nthreads=16 \
  --memgb=64
```

The later one-unit count sentinel is frozen as:

```text
cellranger count \
  --id=<unit_id> \
  --transcriptome=<independently verified custom reference> \
  --fastqs=<verified unit FASTQs> \
  --sample=MantonBM<unit>_HiSeq_9 \
  --chemistry=SC3Pv2 \
  --expect-cells=7000 \
  --include-introns=false \
  --create-bam=false \
  --nosecondary \
  --localcores=16 \
  --localmem=64
```

This prefreeze authorizes only the `mkref` command. It does not execute it and
does not authorize the count sentinel, the remaining seven units, QC, PCA, PH,
landscapes, clustering, label access, or biological outcomes.

## Separating panel and runtime effects

The modern raw-read analysis will keep two comparisons distinct:

1. Exact-500 versus common-475 matrices derived from the same Cell Ranger
   8.0.1 counts estimate the effect of the 25 restored genes without a runtime
   change.
2. Cell Ranger 8.0.1 common-475 versus the admitted historical Cell Ranger
   3.0.0 common-475 matrices quantify the combined historical-versus-modern
   preprocessing difference on an identical ordered feature panel.

The second comparison is a sensitivity analysis, not a claim that runtime is
the only possible cause of every observed difference. Its interpretation must
retain the frozen batch/cohort and historical-processing caveats.

## Resource, validation, and stop rules

- Reference construction remains limited to 16 cores, a 64-GiB request, a
  50-GiB output-tree cap, and a 1-TiB free-space floor.
- No reference artifact may be deleted without project-owner authorization.
- The complete reference tree, its embedded metadata, feature axis, byte size,
  and deterministic tree digest must pass a separately committed independent
  validator before any count sentinel is opened.
- The reference must contain all exact 500 ordered stable IDs and must bind the
  verified FASTA, GTF, runtime, and exact command identities.
- Any version, digest, command, resource, feature-axis, or firewall mismatch is
  a stop condition. It may not be handled by silent substitution or scientific
  retry.

## Landscape contract retained unchanged

Cell and gene observations remain separate typed views. H0 and H1 remain
separate. Essential H0 is excluded explicitly. Landscapes use every consecutive
active level with exact or error-controlled squared-L2 integration, no fixed
grid, and no universal level cap. The Cell Ranger amendment does not reopen or
modify any landscape definition.

## Authoritative documentation

- 10x Cell Ranger 8.0 release notes: <https://www.10xgenomics.com/support/software/cell-ranger/8.0/release-notes/cr-release-notes>
- 10x Cell Ranger command arguments: <https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments>
- 10x custom-reference overview: <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-3p-reference>
- Ensembl release-93 FASTA archive: <https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/>
