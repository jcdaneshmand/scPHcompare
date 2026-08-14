# MV5-AC cosine robustness-outcome prefreeze audit

Date: 2026-08-11

Accepted base: `3fa96fa`

Status: complete; later retrieval execution eligible but not authorized here

## Outcome

MV5-AC proves that the accepted cosine-chord calculation can be compared
one-to-one with the accepted 384-cell/30-PC Euclidean retrieval benchmark. It
freezes the complete retrieval analysis before any cosine ranking or label
value access and formally excludes clustering from this artifact set.

No distance result was compared, no cosine rank was constructed, no tissue
value was read, no accepted endpoint row was opened, and no biological outcome
or method selection was calculated.

## Bound scope

| Item | Frozen value |
|---|---:|
| Source identities | 187 |
| Private cosine group manifests | 150 |
| Representations | 2 |
| Fold-seed groups | 150 |
| Held-out/training biological pairs | 70,700 |
| Pair-by-method rows | 282,800 |
| Later query-by-method rows | 3,600 |
| Later long endpoint rows | 7,200 |
| Method pairing axes | 8 |
| Fixed endpoints | 2 |
| Fixed estimands | 24 |
| Confirmatory MRR DID tests | 4 |

All eight SCT/integrated H0/H1/raw-composite/energy axes contain 35,350 rows
per representation-family and match exactly on representation, fold, seed,
query, training reference, and mapped method. Missing and excess counts are
zero. The matching SHA-256 axis identities are recorded in
`mv05ac-prediction-axis-compatibility-2026-08-11.csv`.

## Frozen scientific contrast

Cosine chord changes only point geometry. Both arms use 384 cells and the same
30 fold-compatible coordinates. The baseline uses Euclidean distance on the
coordinate rows; cosine chord first row-unit-normalizes each nonzero cell
vector and then applies Euclidean chord distance.

H0 and H1 remain separate exact all-active-level landscape components. The
unscaled raw H0/H1 composite remains descriptive. Energy is recomputed in each
arm's own geometry and is the matched non-topological baseline.

The registry contains 16 direct cosine-minus-Euclidean endpoint changes and
eight topology-increment difference-in-differences. Four H0/H1 MRR DIDs form
the only multiplicity-adjusted family. No equivalence or noninferiority claim
is authorized.

## Prediction and label firewall

Future ranks are ordered independently within representation, fold, seed,
method, and query by ascending immutable distance, then ascending canonical
training sample ID for exact ties. All 282,800 rankings must be durably locked
before tissue access. Nonfinite distances, denominator drift, or any post-label
reranking abort the execution.

The historical metadata SHA is unchanged. Only sample and study keys were read
in this sprint: 124 source rows, 18 source studies, 90 analysis samples, and 15
held-out studies match exactly. Tissue and approach values were skipped.
Future execution may read only `Tissue.x`; all other outcomes remain forbidden.

## Aggregation and inference

The fixed aggregation averages five technical seeds within sample, samples
within tissue, then the five tissues equally. Biological sample is the
inferential unit and held-out study within tissue is the resampling block.

Every estimand receives a paired 2,000-replicate tissue-stratified study-block
bootstrap interval (seed `20260814`). Only the four MRR DIDs receive two-sided
9,999-replicate paired study-block sign-flip p-values (seed `20260815`) and
Holm adjustment. Complete null, negative, heterogeneous, failed, and
non-estimable reporting is mandatory.

## Clustering disposition

MV5-AB has 70,700 held-out/training pairs and zero within-training cosine
pairs. Compatible clustering would require 262,675 within-training biological
pairs per representation, 525,350 total before component expansion, followed
by frozen training partitions and held-out assignments. Cosine clustering is
therefore not identifiable from MV5-AB; neither imputation nor reuse of an
incompatible Euclidean partition is allowed.

## Validation and reproducibility

An implementation independent of the MV5-AC helper passed 11 categories:

- all 36 tracked source files and the external label identity match;
- all 150 private group manifest hashes and identities match;
- all 282,800 private axes are unique and exactly pair the baseline;
- the 90-sample structural label join passes without label values;
- all method, endpoint, estimand, prediction-lock, decision, and clustering
  contracts have the frozen cardinalities and states; and
- public artifacts contain no cosine rankings, labels, or outcomes.

A clean second assembly reproduced all 14 generated contract ledgers byte for
byte. The independent-validation and repeat ledgers are public. Private
cosine results and the metadata source remain ignored and untracked.

## Decision and next action

Retrieval evaluation is scientifically identifiable. Clustering is not.
MV5-AC approves only a later, separately committed prediction-locked cosine
retrieval execution. The current 150-row queue remains unauthorized, and both
nested-cell configurations remain closed.

The next sprint should implement the exact frozen ranking/outcome engine, bind
all 187 source identities, commit the 282,800-row cosine prediction lock before
tissue access, then execute and independently reconstruct all 7,200 endpoint
rows and 24 estimands. It must stop if any frozen identity or pairing changes.
