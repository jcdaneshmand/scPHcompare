# MV8-H complete FASTQ acquisition closure

The exact MV8-H HCA raw-read acquisition is complete and independently
validated. All 48 frozen manifest files are present in the private cache, span
the eight admitted bone-marrow units, and contain exactly 85,034,239,918 bytes
(79.194307 GiB). The acquisition receipt has exactly 48 rows, no partial files
remain, and the downloader ended with its required 48-file completion marker
and an empty error log.

The independent cache validator did not trust the downloader's receipt as its
identity proof. It reopened every final FASTQ, recomputed all 48 SHA-256
digests, checked exact byte sizes and gzip magic, and reconciled each result to
both the committed official HCA manifest and the corresponding receipt row.
All 48 files passed. The observed cache total equals the frozen total exactly,
and the storage check still preserves more than the required 1.5 TiB reserve.

This closes only the acquisition gate. It authorizes acquisition and local
identity validation of the frozen Ensembl release-93 primary-assembly FASTA.
It does not authorize `cellranger mkref`, `cellranger count`, QC execution,
normalization, PCA, persistent homology, landscapes, clustering, label access,
biological outcomes, or deletion. Exact Cell Ranger 3.0.0 runtime and custom
reference validation remain the next fail-closed gate.
