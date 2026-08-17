# MV8-H Cell Ranger 8.0.1 runtime/reference prefreeze

The project-owner-approved runtime amendment is prospectively frozen and
independently auditable. The installed Cell Ranger 8.0.1 distribution is bound
by complete deterministic tree SHA-256
`aafd39e293e0ba9d14dba3896a6aeda077304531a2702d26bda0c62c4688fdf3`
plus launcher, STAR, samtools, and Python hashes. The local CLI exposes every
required reference/count control.

The already accepted Ensembl release-93 FASTA, source GTF, independently
reproduced 33,563-gene target-complete GTF, exact-500 panel, and common-475
panel are rebound without changing their identities. All exact 500 stable IDs
remain present in the custom annotation.

The prospective count policy is explicitly `SC3Pv2`, exon-only
(`include-introns=false`), BAM-disabled, secondary-analysis-disabled, and
serialized under the existing compute and storage ceilings. Exact-500 versus
common-475 will be derived within the same 8.0.1 counts; the shared common-475
axis will separately compare modern output with the admitted historical 3.0.0
matrices.

This closure opens only `cellranger mkref`. It does not run `mkref` and does
not authorize any count, QC, PCA, PH, landscape, clustering, label, outcome,
or deletion operation. The custom reference must next pass a complete-tree,
feature-axis, embedded-version, resource, and immutable-input validation before
one complete-unit count sentinel can be prefrozen.

The corrected landscape definition remains unchanged: separate cell/gene and
H0/H1 views, every consecutive active landscape level, exact or
error-controlled squared-L2 integration, no fixed grid, no level cap, and
explicit exclusion of essential H0.
