# MV7-G typed-view and PH sentinel prefreeze v1

Date: 2026-08-15

MV7-G validates the frozen full-124 descriptive transform and both corrected
topology orientations before full PH or any landscape-distance calculation.
It uses the six MV7-D tissue/depth-extreme sample identities across all five
seeds. Their earlier tissue/depth labels are design provenance only; execution
uses only the frozen sample identities.

For each seed, one global standardization and 30-component PCA is fitted
transductively over all 124 samples, 384 deterministic cells per sample, and
the exact MV7-FP 500-gene panel. This yields 47,616 fit cells per seed. The
cell view is 384 cells in the shared 30-PC Euclidean space. The gene view is
500 panel genes with the explicit Pearson-correlation chord distance across
the same 384 cells. Missing features, nonfinite values, and effectively
constant genes stop execution; no imputation or sample exclusion is allowed.

The frozen workload contains five global-fit/view bundles and 60 corrected PH
records (six samples by five seeds by two views). Ripserr runs complete
Vietoris-Rips H0/H1 with field 2 and no filtration truncation. Every full-view
finite H0 death multiset must match an independently calculated MST edge
multiset. All H1 intervals must be finite and have positive persistence.

For the first seed, the first 32 ordered points of all 12 sentinel views are
checked with both Ripserr and `TDA::ripsDiag`/GUDHI: first ordered cells in the
cell view and first ordered genes in the explicit gene-distance view. H0 and
H1 must match within `1e-6` after removing any maximum-scale-capped essential
H0 representation. These 24 reduced checks are corroboration, not biological
results.

Execution is serial with one worker and zero retries. Each global-fit bundle
has a 1,800-second/8-GiB process-tree cap; cell PH has 600-second/4-GiB caps;
gene PH has 1,800-second/8-GiB caps. Aggregate primary plus representative
repeat work has a 14,400-second cap and 2-GiB private-state cap. One complete
seed bundle and its 12 PH records must repeat byte-for-byte. Writes are atomic;
ambiguous or stale resume state stops execution.

The measured conservative projection must place full typed-view and PH work
under 24 worker-hours and 8 GiB private storage to authorize MV7-H prefreeze.
No landscape, pairwise distance, clustering, label, or outcome job is allowed
in MV7-G.

The landscape contract remains unchanged and closed during this sprint:
finite positive-persistence intervals only, essential H0 excluded, all
consecutive active levels, H0 and H1 separate, exact or error-controlled
squared-L2 integration on dimension-specific support, no universal fixed grid
or level cap, and streamed or chunked computation.
