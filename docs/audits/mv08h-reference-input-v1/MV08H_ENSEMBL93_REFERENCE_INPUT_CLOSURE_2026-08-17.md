# MV8-H Ensembl-93 reference-input closure

The frozen Ensembl release-93 GRCh38 primary-assembly FASTA input has been
acquired from the official archive and independently validated. The compressed
archive contains exactly 881,214,682 bytes and has SHA-256
`2a27436d44f0d6350f86894fbe5edec56faa5467028879784508041562406aa0`.

The live official Ensembl `CHECKSUMS` file contains the exact prefrozen record
`07119 860562 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`. A local BSD
`sum` recomputation returns the same checksum and block count, and a complete
gzip read reaches end-of-stream without an integrity error. The evidence
bundle records the decompressed byte count as an additional repeatable stream
check. The downloaded archive and official checksum file remain private,
Git-ignored inputs; only their non-sensitive identities are published.

This closes only reference-input identity. Cell Ranger 3.0.0 is not installed
or validated, no software agreement has been accepted on the owner's behalf,
and `cellranger mkref`, `cellranger count`, QC, topology, landscapes,
clustering, labels, biological outcomes, and deletion remain unauthorized.
The next gate requires an owner-provided exact Cell Ranger 3.0.0 runtime and a
separate custom-reference validation contract.
