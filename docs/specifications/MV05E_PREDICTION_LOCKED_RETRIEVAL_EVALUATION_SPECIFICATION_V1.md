# MV5-E prediction-locked SCT cell retrieval evaluation specification v1

## Document control

| Field | Value |
|---|---|
| Contract ID | `mv05e_prediction_locked_retrieval_evaluation_v1` |
| Date | 2026-08-08 |
| Status | Executed and independently validated |
| Representation | SCT fold-fitted cell view only |
| Primary endpoint | `cross_study_tissue_mrr_v1` |
| Supportive endpoint | `cross_study_tissue_1nn_balanced_accuracy_v1` |
| Inferential unit | Held-out biological sample, blocked by study |
| Outcome-opening boundary | After accepted MV5-D5 prediction hashes |

## 1. Purpose and boundary

MV5-E evaluates the completed MV5-D5 held-out-query rankings without changing
them. It is the first stage allowed to join tissue labels to the corrected
90-sample SCT cell analysis. It answers only whether each frozen method
retrieves samples of the same tissue across studies.

MV5-E does not refit features, scaling, PCA, persistence diagrams, landscapes,
distances, or component weights. It does not rerank neighbors, select a method,
cluster samples, evaluate integration, construct gene topology, fuse views,
introduce new data, or promote a manuscript claim.

## 2. Prediction lock

Labels may open only after all of the following pass:

1. Git base commit is
   `d2a15680a1f2f7d0113b34872f6649ec56902df2`.
2. The MV5-D5 ranking gzip has SHA-256
   `ee925bb073847567dffae61e13cfa688e193267e4708ba3f9edd4a52ceb0c599`.
3. The public group index contains 75 unique fold-seed bundles.
4. The completion table contains 375 unique group-method rows, all completed.
5. The five-method and five-seed axes are exact.
6. Every input remains outcome-label closed and reports no prior biological
   outcome calculation.

The lock record is written before the metadata file is read. The ranking hash
is checked again after it is loaded.

## 3. Frozen labels and eligibility

The authoritative label source is the historical
`joined_metadata_cellcounts.csv` used to freeze MV-05, SHA-256
`e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`.
The tracked app metadata is not substituted for this source.

The source contains 124 large-cohort samples across 18 studies. Eligibility is
reconstructed by the prespecified rule of at least two studies per tissue:

| Tissue | Samples | Studies |
|---|---:|---:|
| Bone marrow | 31 | 3 |
| Colon | 13 | 2 |
| Liver | 6 | 2 |
| PBMC | 12 | 4 |
| Testis | 28 | 4 |
| **Total** | **90** | **15** |

No eligible study spans more than one tissue. Every held-out sample must have a
same-tissue sample in another training study. The structured dispositions are:

- `estimable`;
- `single_study_tissue_not_estimable`; and
- `training_tissue_absent_not_estimable`.

A non-estimable sample is retained in disposition tables and excluded from the
endpoint denominator; it is never scored as an error or silently removed.

## 4. Frozen method panel

The order is fixed and carries no performance ranking:

1. `cell_landscape_h0_v1` — confirmatory topology component;
2. `cell_landscape_h1_v1` — confirmatory topology component;
3. `cell_landscape_h0_h1_raw_euclidean_v1` — descriptive secondary only;
4. `cell_distribution_energy_shared_pca_v1` — matched cell baseline; and
5. `pseudobulk_shared_panel_euclidean_v1` — context baseline.

H0 and H1 remain separate. The raw composite is not promoted because it has no
training-fitted component scale.

## 5. Endpoint definitions

For each method, seed, and held-out query sample:

1. Use the immutable `neighbor_rank`; do not sort again after labels open.
2. Find the first training sample with the query tissue.
3. Record reciprocal rank `1 / first_same_tissue_rank`.
4. For supportive fixed 1-NN, use the tissue of the row whose immutable rank is
   exactly one and record correctness.
5. Retain the chosen sample IDs, tissues, tie flag, training denominator, fold,
   seed, and estimability status.

The primary macro MRR is computed in this frozen order:

1. average the five seeds within biological sample;
2. average samples within tissue; and
3. average the five tissue means equally.

The supportive balanced accuracy uses the identical order with the binary
rank-1 correctness value. Seeds are repeated technical realizations and never
increase the independent sample or study count. Because all stages are linear,
seed-specific macro summaries are also emitted as stability evidence.

## 6. Uncertainty and F1 comparisons

Absolute method intervals and paired topology-minus-energy contrasts use 2,000
bootstrap replicates with deterministic seed `20260808`.

For every replicate:

1. Within each of the five tissue strata, sample that tissue's held-out studies
   with replacement, drawing the original number of studies.
2. Move every sample from a selected study together; repeated study draws repeat
   the entire sample block.
3. Use the same draw for every method in the endpoint.
4. Recompute the sample-within-tissue mean and equal-tissue macro mean.
5. Use the 2.5th and 97.5th percentiles (`quantile` type 7).

The F1 primary family contains exactly two MRR contrasts:

- H0 minus matched energy; and
- H1 minus matched energy.

Two-sided paired study-block sign-flip inference uses 9,999 replicates, seed
`20260809`, and `(b + 1) / (B + 1)`. Holm adjustment is applied to the two raw
F1 p-values. A favorable topology result requires a positive contrast and an
interval against the matched energy baseline; a nominal comparison with the
context pseudobulk method cannot substitute for this test.

If fewer than four independent study blocks contribute, the contrast remains
descriptive, p-values are omitted, and status is
`insufficient_independent_study_blocks`.

## 7. Required validation

An independent script must reconstruct, without calling the production
endpoint helper:

- all label joins and held-out-study exclusions;
- all 2,250 query-method-seed endpoint rows;
- sample, tissue, seed, and method denominators;
- all macro MRR and balanced-accuracy estimates;
- the exact 2,000-replicate block-bootstrap matrices and percentile intervals;
- all paired effects;
- the exact 9,999-replicate sign matrix, null distributions, raw p-values, and
  Holm adjustment;
- unchanged ranking and metadata hashes; and
- zero refit, rerank, clustering, integration, gene-view, or fusion counters.

A clean second production run must reproduce all 15 output files byte for byte.

## 8. Required public outputs

The public audit set contains:

- prediction-lock and label-source records;
- sample eligibility and complete disposition counts, including zero counts;
- query-level endpoints;
- tissue-seed, seed-macro, sample, tissue, and method summaries;
- method intervals and paired contrasts;
- bootstrap and randomization identity audits;
- production summary;
- independent-validation and deterministic-repeat evidence; and
- a final artifact manifest.

Null, negative, heterogeneous, tied, failed, and non-estimable results are
retained. Outcome artifacts are explicitly marked as biological evaluation;
upstream prediction artifacts remain unchanged.

## 9. Acceptance gate

MV5-E passes only when:

- the prediction lock precedes label access;
- all 90 samples, 15 studies, five tissues, five seeds, and five methods
  reconcile exactly;
- every expected query endpoint is present or has a structured disposition;
- independent endpoint and inference reconstruction passes;
- the clean repeat is byte-identical;
- focused and complete package tests pass;
- source-package check is clean; and
- no prohibited downstream or upstream operation occurred.

Passing this gate establishes a trustworthy SCT retrieval result. It does not
establish topological superiority, integration benefit, clustering validity,
gene-view validity, fusion benefit, or manuscript-ready confirmation.
