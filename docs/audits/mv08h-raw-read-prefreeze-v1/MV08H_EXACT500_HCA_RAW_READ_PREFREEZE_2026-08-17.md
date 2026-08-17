# MV8-H exact-500 HCA raw-read prefreeze closure

MV8-H is ready for the authorized data acquisition, but not for raw-read
reprocessing.

The official HCA Azul manifest contains exactly 48 FASTQs across the eight
previously admitted bone-marrow units. Each unit has two lanes and I1/R1/R2,
and all file UUIDs, versions, names, sizes, SHA-256 values, and unit assignments
reconcile exactly. The total is 85,034,239,918 bytes (79.194307 GiB).

The reference design is also prospectively fixed. Rather than discard the 25
targets or build an isolated 500-gene reference, it minimally extends the exact
historical 33,538-gene Cell Ranger 3.0.0/Ensembl-93 axis with the excluded 25,
yielding a 33,563-gene GTF. All 500 target IDs are present and the private GTF
hash is recorded. The primary-assembly FASTA remains blocked before `mkref`
until its local SHA-256 is acquired and frozen.

Independent validation passes 13/13 checks. The current drive preflight found
2,132,205,756,416 bytes free, sufficient for the exact 79.2-GiB acquisition
while preserving the 1.5-TiB reserve. Acquisition is serialized and resumable;
files are atomically published only after exact size and SHA-256 validation.
The prefreeze rebuild and custom GTF repeat byte-identically across 9/9
artifacts, the focused contract suite passes 49/49 expectations, and the
package-loaded complete suite passes 2,500 expectations with four established
skips and no failures or warnings.

The first execution gate is one small FASTQ sentinel. If it validates, the
same committed downloader may resume through all 48 files. Completion still
does not authorize Cell Ranger: the exact 3.0.0 runtime is not installed and a
separate reference/runtime/count-sentinel validation must precede processing.
QC, deterministic 384-cell selection, separate cell/gene topology, complete
VR H0/H1, and the all-active-level exact/error-controlled landscape contract
remain unchanged. Labels and biological outcomes remain closed.
