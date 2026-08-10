# MV5-L locked SCT-versus-integrated retrieval comparison specification v1

## Document control

| Field | Value |
|---|---|
| Contract ID | `mv05l_locked_representation_comparison_v1` |
| Date | 2026-08-10 |
| Status | Executed and independently validated after pre-join freeze |
| Compared representations | accepted MV5-E SCT cell view and accepted MV5-K integrated cell view |
| Primary endpoint | paired change in topology-minus-energy macro MRR |
| Supportive endpoint | paired change in topology-minus-energy macro fixed-1-NN balanced accuracy |
| Inferential unit | held-out biological sample, blocked by study |
| Interpretation | locked-result comparison; not outcome-blind prospectively |

## 1. Purpose and timing disclosure

MV5-L compares two already accepted, immutable retrieval evaluations. It asks
whether inductive integration changes the incremental cross-study tissue
retrieval value of H0 or H1 topology relative to each representation's matched
same-coordinate energy baseline.

The marginal MV5-E and MV5-K aggregate results were known when this contract
was written. Therefore MV5-L is not represented as a fully outcome-blind
prospective analysis. This specification is frozen before reading, joining, or
calculating any sample-level cross-representation contrast. The primary
difference-in-differences estimands follow the pre-existing scientific question
and the already frozen H0/H1-versus-energy hierarchy; they are not selected
from joint results.

MV5-L does not recompute or change rankings, endpoints, methods, features,
scaling, coordinates, persistence diagrams, landscapes, component weights,
labels, eligible samples, tissues, seeds, or folds. It does not tune, select a
winner, cluster, construct gene topology, fuse views, introduce new data,
change package defaults, optimize code, or promote a manuscript claim.

## 2. Immutable input lock

The source revision before MV5-L implementation is
`c22d55540d6eb3c336359f61c0aaf0952202a92b`.

The following files and SHA-256 identities must pass before endpoint rows are
read:

| Input | SHA-256 |
|---|---|
| MV5-E query endpoints | `b65dfd0bae4889fd6297372e7ea180ec8f8465487f68e3c6a5bafd4b9920565f` |
| MV5-E prediction lock | `ed2ff4e8d5c52ddabb527853ba5ff60f81dc6017f79bbea068bdcd07939cd053` |
| MV5-E independent validation | `0f8173d57e1bef27311d8e655de867e2b694c6e39c05e7253a52361c3eb2561b` |
| MV5-E deterministic repeat | `8974196f2e28f712b3b726181b7ebdcc056b7ef8cf779bc4c53ac362e9e7bd07` |
| MV5-E artifact manifest | `2a2401ac30a52067530139fcfb4c6d642d681da3e4ada5a8270db13edf97329a` |
| MV5-K query endpoints | `8fa9f70a64f35b87c5e408d32bd7fc4440d5df611c87cd201ff68968b959df2b` |
| MV5-K prediction lock | `e19334a5adfc2fed09dfd0297958f8a837ddbf5ff6753d97ba8e4ea60696ae56` |
| MV5-K independent validation | `bd80254875753b4488b7488545f348ac7ccf478a1b9f2b89a4cc87aefeffbfd2` |
| MV5-K deterministic repeat | `296ed828e30436782b6053c8ca5a35b71e844d1bd6387ac9ef02c0b10d205f49` |
| MV5-K artifact manifest | `983b6be20b1a5858d14b4d08959f382fb09453b64f08c70ed3931e8d10c51a06` |

Both inputs must retain the same frozen metadata SHA-256
`e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`,
90 samples, 15 held-out studies, five tissues, five seeds, 2,250 estimable
query-method-seed endpoints, zero upstream refits, and zero post-label
reranking.

## 3. Frozen role mapping

| Family | SCT method | Integrated method | Analysis role |
|---|---|---|---|
| H0 | `cell_landscape_h0_v1` | `integrated_cell_landscape_h0_v1` | primary topology component |
| H1 | `cell_landscape_h1_v1` | `integrated_cell_landscape_h1_v1` | primary topology component |
| Raw composite | `cell_landscape_h0_h1_raw_euclidean_v1` | `integrated_cell_landscape_h0_h1_raw_euclidean_v1` | descriptive only |
| Energy | `cell_distribution_energy_shared_pca_v1` | `integrated_cell_distribution_energy_v1` | matched within-representation baseline |
| Pseudobulk | `pseudobulk_shared_panel_euclidean_v1` | `pseudobulk_training_standardized_panel_v1` | shared context and identity control |

The raw H0/H1 composite has no fitted component scale and cannot replace the
separate H0/H1 primary components. Pseudobulk is expected to be identical
because it was intentionally shared; exact equality of ranks, selected
neighbors, reciprocal rank, and 1-NN correctness is a required pipeline
identity control, not a biological result.

## 4. Pairing and eligibility

Rows pair exactly on fold ID, held-out study, seed, query sample ID, and query
tissue after applying the frozen role map. Training sample/study denominators,
endpoint status, and label/refit/rerank flags must agree. Every family must have
exactly 450 paired query-seed rows. No available-case repair is allowed.

If either representation is non-estimable, the pair receives a structured
disposition and is excluded from that estimand's denominator. The expected
accepted state is complete 2,250/2,250 pairing; any mismatch stops production.

## 5. Estimands and aggregation

For endpoint value `Y`, representation `R`, topology dimension `d`, and matched
energy `E`, define the query-level topology increment

`I(R,d) = Y(R,d) - Y(R,E)`.

The two primary representation effects are

`DID(d) = I(integrated,d) - I(SCT,d)`, for `d` in H0 and H1.

The primary endpoint uses reciprocal rank. Positive is the prespecified
direction, meaning integration increased topology's contrast against its own
energy baseline relative to SCT. A positive DID does not by itself mean that
topology beats energy in either representation.

The fixed-1-NN analogues are supportive. Direct integrated-minus-SCT contrasts
for H0, H1, raw composite, energy, and pseudobulk are secondary descriptive
diagnostics because they do not isolate topology's incremental contribution.

Every estimand uses the accepted aggregation order:

1. average five technical seeds within biological sample;
2. average biological samples within tissue; and
3. average the five tissue means equally.

Seeds never increase the independent sample or study count.

## 6. Uncertainty and multiplicity

All 14 estimand-endpoint combinations receive paired 95% intervals from 2,000
tissue-stratified held-out-study bootstrap replicates with seed `20260810`.
Within each tissue, studies are sampled with replacement at the original study
count; all samples from a drawn study move together; the identical draw is used
for every representation, method family, and estimand. Percentiles use R
`quantile` type 7.

The sole confirmatory family `F1_representation_did_mrr` contains exactly two
tests: H0 DID and H1 DID on macro MRR. Two-sided paired study-block sign flips
use 9,999 replicates and seed `20260811`; `(b + 1) / (B + 1)` is used and Holm
adjustment covers exactly these two p-values. Absolute null statistics within
`64 * .Machine$double.eps * max(1, abs(null), abs(observed))` of the observed
absolute statistic count as exceedances. Supportive 1-NN and all direct
contrasts receive intervals but no p-values.

## 7. Required outputs and validation

Production must emit the input lock, method map, compatibility/disposition
records, all paired query endpoints, sample/tissue/macro estimands, intervals,
primary contrasts, bootstrap/randomization audits, pseudobulk identity control,
and zero-operation summary.

An independent validator must reconstruct without calling the production
pairing, estimand, aggregation, bootstrap, or randomization helpers:

- every input hash and the exact five-family role map;
- all 2,250 paired rows and common identity/denominator fields;
- exact pseudobulk identity;
- all direct and H0/H1 DID values at query, sample, tissue, and macro levels;
- the exact bootstrap matrices, intervals, sign matrix, null distributions,
  conservative exceedance counts, raw p-values, and Holm adjustment; and
- zero recomputation, tuning, selection, clustering, gene, fusion, new-data,
  optimization, and SCT/MV5-K mutation counters.

A clean second production run must reproduce every production artifact byte for
byte. Focused and complete tests and a clean prefiltered source-package check
are required.

## 8. Interpretation boundary

MV5-L may state how integration changed the relative topology-minus-energy
retrieval contrast under the frozen existing-data benchmark. It may not claim
that integration is globally beneficial or harmful, that PH is generally
effective or ineffective, that a tissue-specific pattern is confirmed, or
that clustering, gene topology, fusion, external validity, or manuscript
readiness has been established.

## 9. Execution record

The contract was frozen in commit `b3f7e28` before the first endpoint join.
All 2,250 query-method-seed pairs aligned, all 450 shared-pseudobulk endpoints
were exactly identical, and both primary differences-in-differences were
negative with broad intervals spanning zero and Holm-adjusted p-values of 1.
Independent reconstruction passed all 11 categories, and all 13 production
artifacts repeated byte for byte. See
`docs/audits/MV05L_LOCKED_REPRESENTATION_COMPARISON_2026-08-10.md`.
