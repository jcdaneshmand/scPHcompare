# MV8-H FASTQ sentinel closure

The first MV8-H acquisition sentinel passed independently and opens only the
remaining 47-file download.

The sentinel is HCA file order 1,
`MantonBM1_HiSeq_9_S1_L002_I1_001.fastq.gz` for `HCA_BM_001`. The final private
cache file contains exactly 394,373,114 bytes and has SHA-256
`4464d4ea5000356ec74dd80f37dc69016ca13668bc8f33c8d2fcffdef31dd302`,
matching the committed official HCA manifest. It has gzip magic bytes, an exact
one-row receipt, no coexisting partial file, and was atomically published only
after the downloader's identity checks passed.

The independent cache validator repeated the size and SHA-256 checks without
using the downloader implementation. The remaining 84,639,866,804 manifest
bytes fit within the frozen cache cap and preserve more than the required 1.5
TiB free-space reserve.

The exact committed downloader may therefore resume through all 48 manifest
rows. Completion remains an acquisition gate only. `mkref`, `cellranger count`,
QC execution, normalization, PCA, PH, landscapes, clustering, label access,
biological outcomes, and deletion remain unauthorized.
