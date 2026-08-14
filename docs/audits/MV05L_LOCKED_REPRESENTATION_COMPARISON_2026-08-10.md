# MV5-L locked SCT-versus-integrated retrieval comparison audit

## Outcome

MV5-L is complete. Under the frozen paired difference-in-differences design,
integration did not improve either H0 or H1 topology's incremental cross-study
tissue-retrieval contrast relative to the matched within-representation energy
baseline. Both estimates were negative, both blocked intervals were broad and
included zero, and both Holm-adjusted p-values were 1.

This is a locked-result comparison, not a fully outcome-blind prospective
experiment. The marginal MV5-E and MV5-K aggregate outcomes were already known
when the MV5-L contract was written. The joint estimand, method-role map,
pairing, resampling, multiplicity family, and interpretation boundary were
committed in `b3f7e28` before any cross-representation sample-level join.

## Immutable source lock

Ten accepted source artifacts were hashed before endpoint rows were read: the
MV5-E and MV5-K query endpoints, prediction locks, independent validations,
deterministic-repeat records, and public artifact manifests. The source query
endpoint identities were:

| Evaluation | Query-endpoint SHA-256 |
|---|---|
| MV5-E SCT | `b65dfd0bae4889fd6297372e7ea180ec8f8465487f68e3c6a5bafd4b9920565f` |
| MV5-K integrated | `8fa9f70a64f35b87c5e408d32bd7fc4440d5df611c87cd201ff68968b959df2b` |

Every source hash remained unchanged after comparison. Neither accepted source
evaluation was modified.

## Compatibility and identity controls

The five frozen roles mapped exactly: H0, H1, raw H0/H1 composite, matched
energy, and shared pseudobulk. Each family contributed 450 paired
query-seed endpoints over the same 90 samples, 15 held-out studies, five
tissues, and five seeds. Fold, study, seed, sample, tissue, training denominator,
estimability, label-lock, refit, and reranking fields all agreed. There were no
available-case repairs or non-estimable pairs.

The pseudobulk branch was intentionally shared by MV5-E and MV5-K. All 450
endpoints were exactly identical for first same-tissue sample, first
same-tissue rank, reciprocal rank, nearest sample, nearest tissue, 1-NN
correctness, and nearest-distance tie state. Its direct difference is therefore
exactly zero for both endpoints, as required.

## Primary representation effects

The primary estimand is

`(integrated topology - integrated energy) - (SCT topology - SCT energy)`.

Positive was the prespecified direction. It means integration made topology's
contrast against its own matched geometric baseline more favorable; it does
not by itself establish topology superiority within either representation.

| Dimension | Macro-MRR DID | 95% paired blocked-bootstrap CI | Raw p | Holm p |
|---|---:|---:|---:|---:|
| H0 | -0.001195 | [-0.136724, 0.167949] | 0.9916 | 1.0000 |
| H1 | -0.032016 | [-0.114822, 0.052640] | 0.7125 | 1.0000 |

The H0 estimate is almost exactly null. The H1 estimate is more negative, but
its interval remains compatible with changes in either direction. Neither
result supports the prespecified favorable effect.

Supportive fixed-1-NN DIDs were also negative: H0 was -0.047465 (95% CI
[-0.228931, 0.149871]) and H1 was -0.075050 ([-0.254009, 0.030088]). These
supportive results had no p-values by contract.

## Direct representation diagnostics

Direct integrated-minus-SCT macro-MRR changes were:

| Family | Difference | 95% paired blocked-bootstrap CI |
|---|---:|---:|
| H0 | 0.015144 | [-0.108289, 0.173651] |
| H1 | -0.015678 | [-0.093975, 0.067482] |
| Raw composite | 0.014776 | [-0.107202, 0.172405] |
| Matched energy | 0.016339 | [-0.021241, 0.050036] |
| Shared pseudobulk | 0.000000 | [0.000000, 0.000000] |

These are secondary/descriptive because integration changes the underlying
coordinate geometry. They cannot replace the topology-minus-energy DID.

## Tissue heterogeneity

| Tissue | H0 DID MRR | H1 DID MRR |
|---|---:|---:|
| Bone marrow | 0.1867 | 0.0528 |
| Colon | 0.1223 | 0.2836 |
| Liver | 0.1004 | -0.1698 |
| PBMC | -0.2741 | -0.2895 |
| Testis | -0.1412 | -0.0372 |

The sign variation is retained as complete diagnostic evidence. Tissues were
not selected, reweighted, or tested individually, and this table is not a basis
for a tissue-specific claim.

## Implementation stop and accepted run

The first implementation run stopped before constructing any estimand because
a topology-versus-energy identity check compared data-frame row-name attributes
as well as the five scientific identity columns. No result from that attempt
was accepted. The check was corrected to compare the identity vectors
explicitly, focused regression coverage was added, and production was rerun in
a fresh directory. This correction changed no contract, endpoint, or formula.

## Independent validation and determinism

The independent validator did not call the production pairing, estimand,
aggregation, bootstrap, or randomization helpers. It reconstructed all 2,250
pairs, 6,300 query-estimand values, 1,260 sample summaries, 70 tissue summaries,
14 macro estimates, two 2,000-replicate paired study-block matrices, 14
intervals, one 9,999-row sign matrix, two null distributions, conservative
64-epsilon exceedance counts, raw p-values, and Holm adjustment.

All 11 validation categories passed; the largest CSV-versus-memory numerical
difference was `4.89e-15`. A clean second production run reproduced all 13
production artifacts byte for byte. The focused MV5-L regression file passed
27 expectations. The complete package test suite passed with only the two
existing CRAN-gated skips. A fresh prefiltered source tarball built
successfully, and `R CMD check --no-manual` completed with `Status: OK` on R
4.4.1 under Ubuntu 22.04.4.

## Boundary audit

| Operation | Count |
|---|---:|
| Endpoint recomputations | 0 |
| Post-label reranking operations | 0 |
| Method tuning operations | 0 |
| Method-selection operations | 0 |
| Tissue-selection operations | 0 |
| Clustering jobs | 0 |
| Gene-view jobs | 0 |
| Fusion jobs | 0 |
| New-data jobs | 0 |
| Optimization jobs | 0 |
| Accepted source evaluations modified | 0 |

The source PDFs, reviewer material, private caches, and `example_run.r` remain
untracked. Nothing was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve the locked paired DID comparison with explicit timing limitation |
| Correctness demonstrated? | Pass independently and deterministically |
| Computation feasible? | Yes |
| Biological interpretation permitted? | Confirmatory for the frozen existing-data representation DID only; no favorable H0/H1 effect supported |
| Next action | Audit which still-unanswered benchmark axis has highest value before authorizing clustering, technical-mixing, robustness, gene topology, or fusion |

## Claim boundary

MV5-L supports the narrow statement that, in this existing-data LOSO retrieval
benchmark, inductive integration did not measurably improve H0 or H1 topology's
incremental tissue-retrieval contrast over its matched energy baseline.

It does not show that integration is globally neutral or harmful, that PH is
generally ineffective, that tissue-specific signs will replicate, or that
clustering, technical mixing, gene topology, fusion, external validation, or a
manuscript claim is ready.
