# MV6-H blocked fusion outcome execution audit

## Disposition

`fusion_rejected_report_cell_and_gene_views_separately`

The exact-commit, prediction-locked MV6-H blocked outcome evaluation completed
and independently validates 15/15 categories. Equal-weight cell/gene fusion did
not satisfy the prospective requirement to improve MRR over both component
views. G-MV6 therefore closes with a negative fusion disposition: retain the
cell and gene topology views separately and do not admit advanced fusion.

## Durable boundary

- prediction-lock commit:
  `42dc894a6e69a09e8c7bdd3bf0c2f46ac7b9f865`;
- prediction-manifest root:
  `c752408f064800d4539b475455e5bca7e9fcb010019e90a1fc261824c865f4fd`;
- locked rankings: 318,150 rows across 75 groups;
- outcome rows: 4,050 query-method records;
- authoritative metadata:
  `joined_metadata_cellcounts.csv`, SHA-256
  `e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`.

The runner verified and loaded every prediction, then wrote
`mv06h-label-open-receipt.csv` before hashing or parsing metadata. All 75
ranking hashes and the metadata hash matched again after evaluation.

## Primary results

The primary equal-weight fusion MRR was `0.3962926`, compared with `0.4064869`
for the cell composite and `0.3694296` for the gene composite.

| Prespecified contrast | Estimate | 95% blocked-bootstrap interval | Raw p | Holm p |
|---|---:|---:|---:|---:|
| `fusion_gene_weight_050 - cell_composite` | -0.0101943 | [-0.0648421, 0.0360127] | 0.835 | 1.000 |
| `fusion_gene_weight_050 - gene_composite` | 0.0268630 | [-0.0847167, 0.0810076] | 0.561 | 1.000 |

Both contrasts were required to be positive. The fusion-versus-cell contrast
was negative, both intervals included zero, and both Holm-adjusted p-values
were 1. The prespecified fusion-benefit condition therefore fails.

The supportive 1-NN balanced accuracies were `0.3108189` for equal-weight
fusion, `0.3037623` for the cell composite, and `0.2626728` for the gene
composite. This supportive difference cannot override the primary MRR rule.

The descriptive weight `fusion_gene_weight_025` had the largest observed macro
MRR (`0.4179082`). It was explicitly frozen as a sensitivity rather than a
selectable method and is not promoted after label opening.

## Heterogeneity and interpretation

The tissue summaries show material heterogeneity. Equal fusion exceeded the
cell composite in bone marrow, liver, PBMC, and testis but was markedly lower
in colon; the gene composite was strongest in bone marrow and remained weak in
colon. These patterns motivate confounding/robustness synthesis but do not
justify tissue-specific post hoc weights.

This result does not invalidate either topology view. It shows that the fixed
equal-weight combination does not add reliable blocked retrieval value over
the stronger cell composite in this existing-data corpus. The fixed global
gene panel remains technically transductive at the label-free panel-selection
level, so none of these results are external validation.

## Independent validation and repeat

The standalone validator did not call the outcome engine. It independently:

- rejoined 90 samples to 15 held-out studies and five tissues;
- reconstructed all 4,050 MRR/1-NN records from immutable ranks;
- averaged five seeds within each of 810 sample-method records;
- reconstructed 45 tissue-method and nine equal-tissue macro summaries;
- reproduced both 2,000-replicate bootstrap matrices by hash;
- reproduced all 18 method intervals;
- reproduced the 9,999-by-15 sign matrix and both null distributions by hash;
- reproduced both raw and Holm-adjusted p-values;
- confirmed zero refit, reranking, method-selection, advanced-fusion, or
  clustering operations.

It passed 15/15 categories. A second exact-commit outcome run reproduced all
13 production artifacts byte-for-byte.

Evidence identities:

- primary contrasts: `2d529cf83b26d8cf3411998cc2f9b6f6a3d350e7a1258572e312f3ddf0f03437`;
- decision: `37aa0d659669c6cf0306b02cc045e8f6e82498fb28a2fe7b65f3ea364fde7ab6`;
- independent validation:
  `90ffc89cd9cc438b2f435671a8cdaa896ad32c82b7aac24e1856fa7302564b7e`;
- repeat ledger:
  `fd176f259fb1fe98a23a240306d90f082ccd27032e394c69cd8531ad48d98c87`.

## G-MV6 decision and next gate

MV-06 is complete with fusion rejected. Do not run Similarity Network Fusion,
multikernel learning, learned weights, tissue-specific weights, or fusion
clustering. Advance to MV-07 robustness/confounding synthesis with the cell and
gene views reported separately, the complete negative fusion result retained,
and the descriptive weight grid presented without winner selection.
