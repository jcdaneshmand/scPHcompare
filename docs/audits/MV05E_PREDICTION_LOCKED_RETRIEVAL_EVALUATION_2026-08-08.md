# MV5-E prediction-locked SCT cell retrieval evaluation audit

## Outcome

MV5-E is complete. The first corrected full existing-data SCT cell-view
retrieval evaluation is prediction locked, independently reproduced, and
deterministic. Its primary result is null: neither H0 nor H1 demonstrates an
improvement over the matched cell-energy baseline in cross-study tissue MRR.

This result is retained without replacement, tuning, winner selection, or
manuscript promotion. Tissue heterogeneity is substantial and reported in full.

## Prediction and label boundary

| Item | Accepted identity/result |
|---|---|
| Base commit | `d2a15680a1f2f7d0113b34872f6649ec56902df2` |
| D5 ranking SHA-256 | `ee925bb073847567dffae61e13cfa688e193267e4708ba3f9edd4a52ceb0c599` |
| Fold-seed bundles | 75 / 75 |
| Method-group completions | 375 / 375 |
| Ranked rows | 176,750 |
| Frozen metadata SHA-256 | `e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0` |
| Eligible samples/studies/tissues | 90 / 15 / 5 |
| Seeds | `20260805:20260809` |

The script writes the prediction-lock record before opening the metadata. It
rehashes rankings after loading. The metadata source is the same historical
joined table used to freeze MV-05; it was not replaced with a later app table.

## Eligibility and completion

The frozen eligible strata are bone marrow (31 samples/3 studies), colon
(13/2), liver (6/2), PBMC (12/4), and testis (28/4). No study spans tissues.
Every held-out query has same-tissue training support from another study.

The evaluation contains 2,250 query-method-seed observations:

- 90 biological samples;
- five technical subsample seeds;
- five fixed methods; and
- zero non-estimable observations, failures, or exact nearest-distance ties.

The zero counts are present in the disposition artifacts rather than omitted.

## Complete method results

| Method | Role | Macro MRR | Macro 1-NN balanced accuracy |
|---|---|---:|---:|
| H0 landscape | Confirmatory topology | 0.387402 | 0.270632 |
| H1 landscape | Confirmatory topology | 0.383180 | 0.262028 |
| Raw H0/H1 Euclidean | Descriptive secondary | 0.387836 | 0.272061 |
| Cell energy | Matched baseline | 0.391367 | 0.273148 |
| Pseudobulk Euclidean | Context baseline | 0.395837 | 0.279423 |

These values are reported in the frozen registry order. The close numerical
values do not imply equivalence; uncertainty is evaluated on the paired
topology-versus-energy contrasts below.

## Primary F1 comparisons

Inference uses 2,000 paired tissue-stratified held-out-study bootstrap
replicates and 9,999 paired study sign flips. Both contrasts use all 90 samples,
15 independent study blocks, five tissues, and five completed technical seeds.

| Contrast | MRR difference | 95% blocked-bootstrap CI | Raw p | Holm p |
|---|---:|---:|---:|---:|
| H0 − energy | −0.003965 | [−0.074765, 0.062663] | 0.9557 | 1.0000 |
| H1 − energy | −0.008187 | [−0.108536, 0.056906] | 0.8989 | 1.0000 |

The favorable direction was prespecified as positive. Both point estimates are
slightly negative, both intervals include zero broadly, and neither adjusted
test supports improvement. The supportive 1-NN differences are also slightly
negative with broad intervals:

- H0 − energy: −0.002516, CI [−0.079630, 0.088096];
- H1 − energy: −0.011120, CI [−0.099686, 0.094613].

## Tissue heterogeneity

| Tissue | H0 MRR | H1 MRR | Energy MRR | Pseudobulk MRR |
|---|---:|---:|---:|---:|
| Bone marrow | 0.2187 | 0.1495 | 0.3050 | 0.2856 |
| Colon | 0.2764 | 0.3240 | 0.2363 | 0.2894 |
| Liver | 0.0269 | 0.0286 | 0.0235 | 0.0221 |
| PBMC | 0.8513 | 0.7719 | 0.3992 | 0.3856 |
| Testis | 0.5637 | 0.6418 | 0.9929 | 0.9964 |

This pattern shows why the equal-tissue macro endpoint and blocked uncertainty
were necessary. PBMC favors the topology views, testis strongly favors the
conventional baselines, colon is mixed, and liver is weak for every method.
These strata are diagnostics, not a basis for selecting a method or rewriting
the primary question.

## Independent validation

The independent validator does not call the production endpoint helper. It
reconstructed:

- all 176,750 label joins and study exclusions;
- all 2,250 first-same-tissue ranks and fixed rank-1 predictions;
- all 450 sample summaries and five-seed denominators;
- all tissue and macro summaries;
- the exact 2,000-replicate bootstrap matrices by hash;
- every interval and topology-minus-energy effect;
- the exact 9,999-replicate sign matrix and both null distributions by hash;
- raw and Holm-adjusted p-values; and
- unchanged source hashes and zero forbidden-operation counters.

All 11 independent validation categories passed. Maximum numerical differences
were below `1e-12`. A clean second production run regenerated all 15 files byte
for byte.

## Software gates

The focused MV5-E test file passes 23/23 expectations, including a structured
non-estimability regression. The complete source test suite passes 601/601
expectations with no failures, warnings, or skips. A
prefiltered source build and installed-package `R CMD check --no-manual` pass
with `Status: OK`, zero errors, zero warnings, and zero notes.

The first disposable check staging directory activated its copied renv
autoloader without the project library and therefore reported dependencies as
unavailable; it did not test the package. The corrected run disabled staging
autoload and used the existing locked library. No project dependency or lockfile
was installed or changed. The existing `renv` status notice is unchanged from
prior sprints and does not reflect an MV5-E dependency mutation.

## Boundary audit

| Operation | Count |
|---|---:|
| Upstream refits | 0 |
| Post-label reranking operations | 0 |
| Clustering jobs | 0 |
| Integration jobs | 0 |
| Gene-view jobs | 0 |
| Fusion jobs | 0 |
| New datasets | 0 |

The source PDFs, private reviewer material, private compute caches, and
`example_run.r` remain untracked. Nothing was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve the prediction-locked SCT retrieval evaluation |
| Correctness demonstrated? | Pass independently and deterministically |
| Computation feasible? | Yes |
| Biological interpretation permitted? | Confirmatory for the prespecified SCT retrieval endpoint; result is null for topology improvement |
| Next action | Narrow scope and separately specify the next prespecified stage; do not tune or select from MV5-E outcomes |

## Claim boundary

MV5-E supports the statement that corrected SCT H0/H1 topology was evaluated
fairly against a same-cell, same-coordinate energy baseline under study-held-out
retrieval and did not outperform it on the prespecified macro MRR endpoint.

MV5-E does not establish that persistent homology is generally ineffective,
that integration helps or harms topology, that any tissue-specific pattern will
replicate, or that clustering, gene topology, or fusion is valid. Those require
separate, prospectively bounded evidence.
