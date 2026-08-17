# MV8-H exact-500 HCA raw-read acquisition prefreeze v1

**Status:** prospectively frozen for acquisition; raw reprocessing remains closed

**Date:** 2026-08-17

**Parent decision:** D-086 and `docs/audits/mv08g-raw-read-decision-evidence-v1/`

**Public evidence:** `docs/audits/mv08h-raw-read-prefreeze-v1/`

## Purpose and authorization

The common-475 analysis is a valid secondary harmonization sensitivity, but it
is not an adequate sole basis for external cell-H1 topology because removing 25
targets materially changes cell-H1 distance ranks and nearest neighbors. MV8-H
therefore acquires the original reads for a prospective exact-500 external
analysis while keeping every scientific computation and label access closed.

MV8-H exact 48-file HCA FASTQ download authorized by the project owner on 2026-08-17.

This authorization covers the 85,034,239,918-byte (79.194307-GiB) HCA FASTQ
manifest and the small Ensembl release-93 reference inputs needed by the next
reference-validation gate. It does not authorize `cellranger mkref`,
`cellranger count`, QC execution, normalization, PCA, PH, landscapes,
clustering, label access, biological outcomes, or deletion of any material
artifact.

## Exact source cohort and read structure

- HCA project: `cc95ff89-2e68-4a08-a234-480eca21ce79`, “A single cell immune
  cell atlas of human hematopoietic system.”
- Catalog: HCA Azul `dcp60`.
- Biological units: exactly `HCA_BM_001` through `HCA_BM_008`; units are never
  pooled or treated as pseudoreplicates.
- Library construction: 10x 3′ v2; `paired_end=false` in official metadata.
- Each unit has exactly two lanes (`L002`, `L003`) and three read roles per lane
  (`I1`, `R1`, `R2`), for six files per unit and 48 files total.
- Every row is frozen by unit, name, HCA file UUID, opaque file version, exact
  byte count, SHA-256, CRC32C, DRS URI, Azul repository route, and mirror URI in
  `mv08h-fastq-manifest.csv`.
- The public unit mapping is independently joined to the previously admitted
  H5 files through private donor-document identities; the donor identities are
  never published.

The manifest was requested from the official HCA Azul API with the already
admitted project/organ/cell-type/one-donor filters and
`fileFormat=fastq.gz`. Temporary manifest and signed object URLs are private;
only stable public identities and routes are committed. HCA access and reuse
remain governed by its published data-release policy and CC BY 4.0 terms.

## Custom reference contract

The baseline HCA matrix uses the exact Cell Ranger 3.0.0/Ensembl-93
33,538-gene filtered axis. All 500 targets occur in the unfiltered Ensembl-93
GTF, but 25 were excluded by the historical Cell Ranger biotype allow-list.
The custom annotation is therefore a minimal union:

1. retain every gene admitted by the documented Cell Ranger 3.0.0 17-biotype
   allow-list;
2. add only the 25 exact stable IDs identified by MV8-E; and
3. retain all GTF feature rows belonging to that 33,563-gene union.

The resulting private 33,563-gene GTF is frozen at
SHA-256 `e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182`.
It contains 33,563 unique gene features and all 500 ordered target stable IDs.
This minimizes reference drift relative to the admitted HCA matrix while
restoring countability of the complete panel. It is not a 500-gene-only
reference, which could redirect reads from overlapping omitted genes.

The genome source is Ensembl release-93 GRCh38 primary assembly
`Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`, 881,214,682 bytes, frozen by
its immutable archive path and Ensembl BSD checksum record
`07119 860562`. Its local SHA-256 must be recorded and independently validated
before `mkref`; archive metadata alone cannot open that gate.

## Counter, QC, and topology contract

The selected comparator-preserving counter is the exact Cell Ranger 3.0.0
binary distribution with explicit `SC3Pv2` chemistry. It is currently absent
from WSL, so no substitution or newer runtime may be used silently. The frozen
commands and serial 16-core/64-GiB limits are recorded in
`mv08h-processing-contract.csv`. A future reference/runtime validation must
hash the executable distribution, embedded versions, command lines, FASTA,
custom GTF, and complete `mkref` output tree before authorizing one count
sentinel.

After counting, the analysis must reuse
`docs/specifications/mv08d-hca-qc-contract-v1.csv` without modification:

1. entry `min.features = 200` and `min.cells = 3`;
2. `500 <= nFeature_RNA <= 9000`;
3. mitochondrial percentage `<= 20` and ribosomal percentage `> 5`;
4. post-QC features detected in more than three cells; and
5. at least 384 eligible cells per unit.

Post-QC barcode IDs are sorted before five deterministic selections of exactly 384 cells
under seeds `20260805` through `20260809`. The exact ordered 500-gene
axis must be present before normalization or topology. Cell and gene topology
remain distinct typed views. Complete Vietoris–Rips H0 and H1 are retained.
Landscapes use every consecutive active level, exact or error-controlled
squared-L2 integration, no fixed grid, no level cap, separate H0 and H1, and
explicit exclusion of essential H0. Raw common-475 remains a secondary
harmonization sensitivity, not a replacement for exact-500 external topology.

## Acquisition and resource policy

The downloader is `scripts/fetch_mv08h_hca_fastq.py`. It must:

- verify the committed manifest SHA-256 and the authorization sentence above;
- preflight the complete remaining byte count, a 92,274,688,000-byte cache cap,
  and at least 1.5 TiB free after all remaining FASTQs;
- acquire one file at a time, beginning with one small sentinel;
- use at most five bounded network attempts and resume only when HTTP Range is
  explicitly honored;
- retain interrupted or identity-failing `.partial` files for review;
- validate exact size and SHA-256 before atomic publication;
- write an atomic private receipt without signed URLs or absolute paths; and
- treat hash, size, manifest, authorization, cap, and Range failures as
  non-scientific stop conditions—not reasons to alter data or retry silently.

Future `mkref` and `count` limits are prospective only: 50 GiB for the reference
tree, 80 GiB per-unit peak RSS, 200 GiB per-unit workspace, 24 hours per unit,
one serialized unit, and a 1-TiB free-space floor. No deletion is authorized by
this contract.

## Validation, firewall, and stop rules

Thirteen independent prefreeze checks must pass before the sentinel. The
sentinel must match exact size and SHA-256 before the remaining 47 files may
continue. Download closure requires all 48 exact files, all eight unit totals,
zero partial files, an independently reconstructed receipt, the manifest total,
and a cache-cap/free-space pass.

The complete download does not authorize processing. Processing remains
blocked until the local FASTA SHA-256, exact Cell Ranger 3.0.0 runtime, custom
reference tree, exact-500 feature oracle, and a separately committed count
sentinel contract independently pass. Labels and outcomes remain closed through
immutable external topology and prediction locks. Every future PH and landscape
artifact must reuse the accepted typed-view and landscape contracts; no result
may redefine their level count, grid, H0/H1 combination, or essential-class
policy.

## Authoritative public sources

- HCA project portal: <https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79>
- HCA Azul service/OpenAPI: <https://service.azul.data.humancellatlas.org/>
- Ensembl release-93 FASTA archive: <https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/>
- 10x custom-reference overview: <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-3p-reference>
