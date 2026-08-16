# MV7-FP full-124 global-core panel lock prefreeze v1

Date: 2026-08-15

MV7-FP applies the unchanged MV6-C availability and within-cache variance rule
to the complete label-closed 124-sample, five-seed SCT corpus. The input axis is
exactly 450 accepted primary caches plus 170 MV7-F caches, 620 records total,
with 384 deterministic cells per record.

The algorithm forms the exact feature intersection, excludes the accepted
mitochondrial/ribosomal/hemoglobin technical categories, requires finite
variance strictly above machine epsilon in every cache, ranks features by
descending variance within each cache, orders eligible canonical genes by the
median of 620 ranks with the existing deterministic tie breaks, and selects the
first 500. No alternative panel size, threshold, imputation, sample exclusion,
or result-dependent branch is allowed.

The calculation is explicitly transductive technical harmonization. It may
read only cache identities and SCT matrices; tissue, study, approach, endpoint,
clustering, and outcome labels remain closed. Public output may contain the
ordered panel and aggregate eligibility/stability summaries, never expression
values, per-cache variances/ranks, cell IDs, or private paths.

The ordered-panel SHA-256 uses a canonical row-name-free payload containing
only integer `panel_order`, character `feature_id`, and character `gene`, so
the digest is reproducible from the public ordered-panel fields.

The inventory has one 2,700-second and 8-GiB process cap. It must validate all
620 cache hashes and payload identities, run in one streaming pass, pass an
independent reconstruction, and reproduce deterministic outputs byte-for-byte.
Failure to yield 500 eligible unique genes, any resource breach, or any source
drift stops before PCA, typed views, PH, landscapes, labels, or outcomes.

Passing MV7-FP freezes one exact 124-derived panel and authorizes only a
separate MV7-G typed-view/PH sentinel prefreeze. The landscape contract remains
finite positive intervals, essential H0 excluded, all consecutive active
levels, H0/H1 separate, exact or error-controlled squared-L2 distance, no
universal grid or level cap, and streamed/chunked computation.
