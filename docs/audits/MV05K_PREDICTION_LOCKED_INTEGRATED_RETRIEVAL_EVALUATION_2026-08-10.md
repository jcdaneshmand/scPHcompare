# MV5-K prediction-locked integrated cell retrieval evaluation audit

## Outcome

MV5-K is complete. The integrated cell retrieval evaluation was locked before
label access, independently reconstructed, and repeated byte for byte. Neither
integrated H0 nor integrated H1 demonstrates improvement over the matched
integrated-coordinate energy baseline on the prespecified cross-study tissue
macro mean reciprocal rank endpoint.

This sprint did not read or compare the accepted SCT outcome. The integrated
result is retained without replacement, tuning, winner selection, selective
tissue reporting, or manuscript claim promotion.

## Prediction and label boundary

| Item | Accepted identity/result |
|---|---|
| Base commit | `6836d1cec28bf23f2b920c3e88bffb2e12385e24` |
| MV5-J ranking SHA-256 | `4588902bce89a04cae0c7676b4f21f81e83013a29120ca2a4b39f3ffacfb677e` |
| MV5-J completion SHA-256 | `971e6902c1b6c2b8d39c74a9d07abc9c6d05c99c7b47bb8bd10e4a2209a08ed9` |
| MV5-J group-index SHA-256 | `0fd3e05a312079cebd4e49f4d0f7326cda037fd460772176385d6ea9e0e023ab` |
| MV5-J method-registry SHA-256 | `becb8451948e421916915fa231a50a83adb22fbd62cc10cadb3a6d934fd7de5d` |
| MV5-J scale-disposition SHA-256 | `9f6f93e73c072088e9917d395dd2f98005bfeb9d4ec6ec6dec33c0d50333e8bf` |
| MV5-J assembly-summary SHA-256 | `aabfdb0338826451d08efb27ebacea4eb6b56cad5fc30bd6f9159f60baeff68c` |
| Fold-seed bundles | 75 / 75 |
| Method-group completions | 375 / 375 |
| Frozen metadata SHA-256 | `e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0` |
| Eligible samples/studies/tissues | 90 / 15 / 5 |

The lock passed and was written before the metadata was read. Rankings were
rehashed after loading. The method registry and all 225 scale-disposition rows
confirmed the frozen five-method axis and zero held-out-derived scale fitting.

## Eligibility and completion

The same frozen strata as MV5-E were used: bone marrow (31 samples/3 studies),
colon (13/2), liver (6/2), PBMC (12/4), and testis (28/4). Every held-out query
had same-tissue support in another training study.

All 2,250 expected query-method-seed observations were estimable: 90 biological
samples, five technical seeds, and five fixed methods. There were zero failures,
non-estimable observations, or exact nearest-distance ties. Zero-count
dispositions remain explicit in the public tables.

## Integrated method results

| Method | Role | Macro MRR | Macro 1-NN balanced accuracy |
|---|---|---:|---:|
| Integrated H0 landscape | Confirmatory topology | 0.402546 | 0.262062 |
| Integrated H1 landscape | Confirmatory topology | 0.367502 | 0.225873 |
| Raw integrated H0/H1 Euclidean | Descriptive secondary | 0.402612 | 0.260772 |
| Integrated coordinate energy | Matched baseline | 0.407705 | 0.312043 |
| Training-standardized pseudobulk | Context baseline | 0.395837 | 0.279423 |

The raw composite remains descriptive and pseudobulk remains contextual. The
primary comparisons are the separate H0 and H1 contrasts against matched
integrated-coordinate energy.

## Primary F1 comparisons

Inference used 2,000 paired tissue-stratified held-out-study bootstrap
replicates and 9,999 paired study sign flips over all 90 samples, 15 study
blocks, five tissues, and five completed technical seeds.

| Contrast | MRR difference | 95% blocked-bootstrap CI | Raw p | Holm p |
|---|---:|---:|---:|---:|
| Integrated H0 - integrated energy | -0.005159 | [-0.140977, 0.164341] | 0.9625 | 1.0000 |
| Integrated H1 - integrated energy | -0.040203 | [-0.165125, 0.062618] | 0.6395 | 1.0000 |

Both favorable directions were prespecified as positive. Both point estimates
are negative, both intervals include zero broadly, and neither adjusted test
supports topology improvement. Supportive 1-NN differences are also negative:
H0 minus energy is -0.049981 (CI [-0.220025, 0.145265]) and H1 minus energy is
-0.086170 (CI [-0.263772, 0.024701]).

## Tissue heterogeneity

| Tissue | H0 MRR | H1 MRR | Energy MRR | Pseudobulk MRR |
|---|---:|---:|---:|---:|
| Bone marrow | 0.5473 | 0.3443 | 0.4469 | 0.2856 |
| Colon | 0.1860 | 0.3949 | 0.0235 | 0.2894 |
| Liver | 0.2957 | 0.0272 | 0.1919 | 0.0221 |
| PBMC | 0.5541 | 0.4594 | 0.3761 | 0.3856 |
| Testis | 0.4297 | 0.6117 | 1.0000 | 0.9964 |

This heterogeneity is complete diagnostic evidence, not a basis for selecting a
method, tissue, or revised endpoint.

## Numerical-boundary correction

The first independent validation rejected the initial output because two
sign-flip null statistics lay at the observed absolute-statistic boundary. CSV
roundoff changed the exceedance count by 2/10,000 even though the null hashes
were identical. No result from that pass was accepted.

The implementation and specification were corrected with one uniform,
conservative policy: values within 64 machine epsilons of the two-sided
boundary count as exceedances. Focused regression tests pass, the outputs were
regenerated from scratch, and the independent validator now exactly reproduces
the exceedance counts, p-values, and Holm adjustment. This correction is
representation- and outcome-agnostic and does not change ranks, endpoint
estimates, bootstrap matrices, or intervals.

## Independent validation and determinism

The independent validator did not call the production endpoint helper. It
reconstructed all 176,750 label joins, 2,250 endpoints, 450 sample summaries,
tissue and macro denominators, exact 2,000-replicate bootstrap matrices, all
intervals and paired effects, the exact 9,999-replicate sign matrix and null
distributions, conservative exceedance counts, raw p-values, and Holm values.

All 11 validation categories passed with maximum numerical difference below
`5e-16`. A clean second production run reproduced all 15 public artifacts byte
for byte.

The focused MV5-K regression file passed all 25 expectations. The complete
package test suite passed with only the two existing CRAN-gated skips. A fresh,
prefiltered source tarball built successfully, and `R CMD check --no-manual`
finished with `Status: OK` on R 4.4.1 under Ubuntu 22.04.4.

## Boundary audit

| Operation | Count |
|---|---:|
| Upstream refits | 0 |
| Post-label reranking operations | 0 |
| Clustering jobs | 0 |
| Integration jobs | 0 |
| Gene-topology jobs | 0 |
| Fusion jobs | 0 |
| New-data jobs | 0 |
| SCT outcome files read | 0 |

The source PDFs, reviewer material, private caches, and `example_run.r` remain
untracked. Nothing was pushed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve the prediction-locked integrated retrieval evaluation |
| Correctness demonstrated? | Pass independently and deterministically |
| Computation feasible? | Yes |
| Biological interpretation permitted? | Confirmatory for this integrated retrieval endpoint; result is null for topology improvement over matched integrated energy |
| Next action | Separately specify comparison of already locked SCT and integrated results; do not tune or select from either outcome |

## Claim boundary

MV5-K supports the statement that integrated H0/H1 topology was evaluated
fairly against a same-cell, same-coordinate integrated energy baseline under
study-held-out retrieval and did not outperform it on the prespecified macro
MRR endpoint.

MV5-K alone does not establish whether integration helps or harms relative to
SCT, that persistent homology is generally ineffective, that tissue-specific
patterns will replicate, or that clustering, gene topology, or fusion is valid.
Those questions require separate, prospectively bounded evidence.
