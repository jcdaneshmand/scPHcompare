# ADR-015: Evaluate SCT cell retrieval only after an immutable prediction lock

## Status

Accepted and executed on 2026-08-08.

## Context

MV5-D5 produced 176,750 label-closed rankings for five frozen methods across 75
fold-seed bundles. Opening tissue labels earlier would permit outcome-informed
reranking, scaling, selection, or replacement. The five subsample seeds are
technical repetitions, while studies—not cells, pairs, or seeds—are the
independent generalization blocks.

## Decision

We require the accepted D5 ranking hash, 75 bundle identities, and 375/375
completion rows to pass before reading the frozen metadata. MV5-E then computes
only cross-study tissue MRR and fixed 1-NN balanced accuracy from the immutable
ranks. It averages seeds within biological samples, samples within tissues, and
the five tissues equally.

H0 and H1 are compared separately with the matched cell-energy baseline. The
raw H0/H1 combination remains descriptive and pseudobulk remains contextual.
Uncertainty uses paired tissue-stratified held-out-study bootstrap resampling;
the two F1 MRR comparisons use paired study sign flips and Holm correction.
Every null, negative, tissue-specific, failure, and non-estimable result is
retained.

## Evidence

The prediction lock matched commit `d2a1568`, ranking SHA-256 `ee925b…c599`, 75
bundles, and 375/375 completed method groups. The frozen metadata hash matched
the statistical-plan record and reconstructed 90 samples in 15 tissue-specific
study blocks.

All 2,250 expected query-method-seed observations were estimable. Macro MRR was
0.3874 for H0, 0.3832 for H1, 0.3878 for the descriptive raw composite, 0.3914
for matched energy, and 0.3958 for pseudobulk context. H0 minus energy was
−0.00396 (95% blocked-bootstrap CI −0.0748 to 0.0627); H1 minus energy was
−0.00819 (−0.1085 to 0.0569). Both Holm-adjusted p-values were 1. These are null
primary comparisons, not evidence that topology improves SCT tissue retrieval.

Results differed substantially by tissue: topology was strongest for PBMC,
whereas energy and pseudobulk were nearly perfect for testis; all methods were
weak for liver. This heterogeneity is retained as a diagnostic and is not used
to select a winner or redefine the endpoint.

An independent implementation reproduced every endpoint, denominator,
bootstrap matrix, interval, paired contrast, sign-flip null, p-value, and Holm
adjustment. Fifteen cleanly regenerated production artifacts were byte
identical. No upstream refit/rerank or clustering, integration, gene, or fusion
job occurred.

## Consequences

The corrected SCT cell-view primary retrieval result is trustworthy but does
not support a topology-superiority claim. The result is scientifically useful:
it narrows the claim space, demonstrates strong tissue dependence, and provides
a leakage-resistant baseline for deciding which prespecified analysis should
come next.

MV5-E does not authorize post hoc method tuning, selective tissue reporting,
integration claims, clustering, gene topology, fusion, new data, or manuscript
claim promotion. Any next stage requires a separate contract and goal.
