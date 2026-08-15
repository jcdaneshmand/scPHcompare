# MV6-G complete-corpus fusion prefreeze specification v1

## Document control

| Field | Frozen value |
|---|---|
| Contract | `mv06g_complete_fusion_prefreeze_v1` |
| Date | 2026-08-14 |
| Accepted input | Complete MV6-F at `bba0b11` |
| Samples / studies / tissues | 90 / 15 / 5 |
| Folds / seeds | 15 LOSO folds / 5 fixed seeds |
| Views | `cell_topology_v1`, `gene_topology_v1` |
| Dimensions | H0 and H1 retained separately |
| Outcome state | Closed during all scale/ranking work |
| Current authorization | One maximum-group training-scale/ranking sentinel only |

## 1. Purpose and scientific limitation

MV6-G asks whether equal-weight fusion of independently valid cell and gene
topology improves cross-study tissue retrieval over both view-specific
composites. It extends the interpretable arithmetic fusion tested in MV6-A to
the complete 90-sample blocked corpus. It does not select a weight, learn a
kernel, run Similarity Network Fusion, cluster samples, or establish a claim
before the prediction and evaluation gates pass.

The fixed 500-gene panel was selected using label-free presence and variance
rules across all accepted caches. It is therefore explicitly transductive at
the technical panel-definition level, although every source transform and
component scale is fold-training-only. Any eventual biological interpretation
must retain that limitation; MV6-G is not external validation.

## 2. Immutable prediction inputs

The only prediction inputs are the 75 accepted MV6-F group directories, the
queue root `f5471633…10bb5`, implementation root `5a1258e8…8d292`, accepted
Rust library `51d3fca4…160d`, complete inventory, and 376-row immutable-resume
ledger. Each group contains 90 matched cell/gene diagram records with explicit
training/query roles and the accepted 35,350 query-to-training biological
pairs across the corpus. Tissue and approach fields are prohibited.

## 3. Required training-only component scales

Fusion cannot use raw components because cell/gene and H0/H1 scales differ.
For every fold, seed, and component, compute exact all-active-level landscape
distances for every unordered pair of training samples and fit

`s_j = median{D_j[a,b] : a < b; a,b are training samples}`.

Each scale must be finite and greater than `sqrt(.Machine$double.eps)`. Apply
`Z_j = D_j / s_j` to the already accepted query-to-training component rows.
No held-out query distance, label, clipping, centering, or outcome may affect a
scale. The additional frozen workload is 262,675 training biological pairs,
1,050,700 component rows, and 300 scales. PH is not rerun.

## 4. Frozen method panel and ranking

Within-view composites are

- `C = 0.5 Z_cell,H0 + 0.5 Z_cell,H1`; and
- `G = 0.5 Z_gene,H0 + 0.5 Z_gene,H1`.

Fusion is `F_w = (1-w)C + wG` for the complete fixed grid
`w in {0, 0.25, 0.5, 0.75, 1}`. The nine immutable ranking methods are the four
normalized components, `C`, `F_0.25`, `F_0.5`, `F_0.75`, and `G`. `F_0.5` is
the sole primary fusion; `C` and `G` are its two required comparators;
intermediate weights and individual H0/H1 components are descriptive
sensitivities and cannot replace or tune the primary.

For every fold, seed, method, and held-out query, rank training samples by
ascending distance and canonical training sample ID for exact ties. The frozen
output contains 318,150 label-closed ranking rows. All scales, rankings,
completion records, hashes, repeat evidence, and resume evidence must be
accepted before the metadata file is read.

## 5. Prediction staging and resources

Stage one computes only the previously accepted maximum group, containing
2,080 training pairs, 8,320 component rows, four scales, 1,625 query pairs, and
14,625 ranking rows. It runs with one worker, 1,800-second and 12-GiB
process-tree limits, a 5-GiB private-storage ceiling, and no retry. Full
production remains closed until exact R/Persim oracles, deterministic repeat,
immutable resume, and a measured projection below 12 worker-hours pass.

Full production, if separately admitted, remains serial and stops on the first
resource, numerical, identity, or atomic-publication failure. It cannot open
labels or launch evaluation.

## 6. Prediction lock and label firewall

The outcome stage may begin only after a committed prediction manifest binds
all source, scale, ranking, implementation, and completion hashes. The
authoritative metadata SHA-256 is
`e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`.
The outcome code must write its lock record before reading that file and must
rehash predictions afterward. Post-label scaling, reranking, weight selection,
method replacement, tissue selection, and endpoint selection are prohibited.

## 7. Frozen endpoints and inference

The primary endpoint is cross-study tissue mean reciprocal rank. Fixed 1-NN
balanced accuracy is supportive. Aggregation matches MV5-E: average seeds
within biological sample, samples within tissue, then the five tissue means
equally. The inferential unit is the held-out biological sample, blocked by
study; seeds never increase the independent count.

The F1 family contains exactly two MRR contrasts:

1. `F_0.5 - C`; and
2. `F_0.5 - G`.

Fusion benefit requires both contrasts to be positive. Uncertainty uses 2,000
tissue-stratified held-out-study bootstrap replicates with seed `20260815`.
Two-sided paired held-out-study sign flips use 9,999 replicates with seed
`20260816`; raw p-values receive Holm adjustment across the two F1 contrasts.
Fewer than four independent study blocks yields descriptive status only.
Intermediate weights and individual components receive estimates and
intervals but no confirmatory selection or substitution.

## 8. Validation and stop boundary

Prediction production requires independent reconstruction of training-pair
membership, scales, all nine formulas, ranks, ties, identities, source
immutability, deterministic repeat, and zero-rebuild resume. Outcome production
requires independent label joins, endpoints, aggregation, bootstrap matrices,
sign matrices, contrasts, multiplicity, prediction hashes, and byte repeat.

The current prefreeze authorizes only the maximum label-closed training-scale
and ranking sentinel. Full training-scale production, label opening, outcome
evaluation, clustering, advanced fusion, new data, package defaults, release,
and claims remain closed until their separate gates pass.
