# MV5-AG nested-192 robustness-outcome prefreeze specification v1

Date frozen: 2026-08-11

Accepted nested-192 calculation base: `3286e9f`

Outcome state: not computed

## Purpose

MV5-AG prospectively freezes a retrieval-only comparison between the accepted
384-cell, 30-coordinate Euclidean benchmark and the accepted nested-192-cell,
30-coordinate Euclidean calculation. The scientific question is whether
reducing deterministic cell representation depth from 384 to 192 changes
cross-study tissue retrieval or changes topology's increment relative to
energy divergence calculated from the same selected coordinates.

This is a one-factor-at-a-time cell-depth sensitivity. Coordinate count,
point metric, fold-fitted representation, frozen 384-cell realization, fold,
seed, query, reference, landscape definition, and method family remain fixed.
The 192 cells are the first 192 elements of the frozen sample/seed/cell-ID
SHA-256 order and are contained in the still-closed 256-cell set. SCT and
inductive-integrated views remain separate strata.

The contract is frozen before any nested192 ranking, tissue-label join, endpoint,
estimate, comparison, method selection, or clustering calculation. Existing
Euclidean outcomes are known historical evidence, so complete fixed-panel
reporting and immutable identities—not a claim of total blindness—protect this
analysis from favorable post-hoc selection.

## Accepted scientific objects

The accepted Euclidean and nested192 artifacts share exactly:

- 90 biological samples in 15 leave-one-study-out folds;
- five technical cell-subsample seeds (`20260805` through `20260809`);
- the same frozen 384-cell realization, with 384 baseline cells versus its
  deterministic nested 192-cell subset, and 30 fold-compatible coordinates;
- 75 SCT and 75 inductive-integrated fold-seed groups;
- 70,700 directed held-out-query/training-reference biological pairs;
- 282,800 rows across four paired method families;
- exact all-active-level persistence-landscape L2 distances for H0 and H1,
  retained separately;
- the unscaled raw `sqrt(H0^2 + H1^2)` composite as descriptive only; and
- energy divergence as the matched same-geometry non-topological baseline.

The sole changed factor is cell count:

| Configuration | Deterministic cell selection | Point distance |
|---|---|---|
| Accepted baseline | all frozen 384 cells | Euclidean |
| Nested 192 | first 192 cells in frozen SHA-256 order | Euclidean |

No coordinate row is normalized, rescaled, refit, or reselected after the
prefreeze. Nonfinite coordinates, selection-order drift, or failure of
192-within-256 inclusion causes an abort before label access.

## Exact method pairing

Each representation has four paired families: H0, H1, raw H0/H1 composite,
and energy. The generic nested192 method IDs are mapped to the accepted integrated
baseline names only for identity pairing. No value is changed by this mapping.

All eight representation/family axes must match one-to-one on representation,
fold, seed, held-out query sample, training reference sample, and mapped method.
Duplicate, missing, or excess keys stop the sprint. Pairability establishes a
valid paired design; it does not compare distances, ranks, or outcomes.

The coordinate source must also match exactly in every fold/seed group. For
SCT, the nested-192 source hash must equal the accepted MV5-D1 384-cell private
cache hash and its fold cache key must equal the MV5-D5 baseline key. For the
inductive-integrated view, it must equal the MV5-J-bound MV5-G 384-cell private
file hash. The private MV5-AF source identity must independently repeat the
same hash. Thus the 192 cloud is a deterministic subset of the exact accepted
384-cell source, rather than merely a cloud with compatible dimensions.

## Prediction lock and label firewall

A later runner must first verify every frozen source and then rank all and only
the training references independently within representation, fold, seed,
method, and query. Ordering is ascending immutable distance; an exact distance
tie is broken by ascending canonical training sample ID using radix ordering.
Any nonfinite distance aborts before a lock is published.

The complete 282,800-row ranking artifact and its manifest must be durably
committed before tissue access. No reranking, rescaling, refitting, method
removal, denominator change, or tie-policy change is allowed afterward.

The historical metadata source remains untracked and is bound by SHA-256
`e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`.
This prefreeze reads only sample and study keys and proves the 90-sample,
15-held-out-study join. It skips all label values. At execution, only
`Tissue.x` may be read after the prediction lock; `Approach.x` and every other
outcome are prohibited. Accepted Euclidean endpoint rows are hash-bound but
must not be opened until the nested192 prediction lock is durable.

## Endpoints and complete reporting

Two endpoints are inherited unchanged:

- primary: `cross_study_tissue_mrr_v1`;
- supportive: `cross_study_tissue_1nn_balanced_accuracy_v1`.

For each query, reciprocal rank uses the first same-tissue training reference;
fixed rank one yields the supportive correctness indicator. Expected later
output is 3,600 query-method rows and 7,200 query-method-endpoint rows.

Every expected row must contain an estimate or a structured disposition.
Allowed non-estimable states include a single-study tissue, an absent training
tissue, prediction-axis mismatch, immutable-source mismatch, and software
failure. Available-case repair, replacement folds/seeds, and omission of null,
negative, heterogeneous, tied, failed, or non-estimable results are prohibited.

## Frozen estimands

For endpoint value `Y`, representation `R`, and family `m`, define:

`DIRECT(R,m,Y) = Y(R,m,nested192_cells) - Y(R,m,384_cells)`

for two representations, four families, and two endpoints: 16 estimands.

For topology dimension `d` in H0/H1, define:

`DID(R,d,Y) = [Y(R,d,nested192)-Y(R,energy,nested192)] - [Y(R,d,384cell)-Y(R,energy,384cell)]`

for two representations, two topology dimensions, and two endpoints: eight
estimands. The complete registry therefore has 24 rows.

Positive direct values mean the named retrieval quantity increased under
the nested-192 setting. Positive DID values mean topology's increment relative to its
matched energy baseline increased; this alone does not establish that topology
outperforms energy in either geometry.

The only confirmatory family comprises the four H0/H1 DID estimands on MRR.
Direct MRR changes are secondary sensitivity diagnostics and all fixed-1-NN
changes are supportive. No equivalence or noninferiority margin is justified,
so nonsignificance cannot be called equivalence, invariance, or robustness.

## Aggregation and uncertainty

Aggregation is fixed in this order:

1. average five technical seeds within biological sample;
2. average biological samples within tissue;
3. average the five tissue means equally.

The biological sample is the inferential unit. Held-out study within tissue is
the resampling block; technical seeds never inflate sample size.

All 24 estimands receive paired 95% percentile intervals from 2,000
tissue-stratified held-out-study bootstrap replicates using seed `20260814`.
Within each tissue, studies are sampled with replacement at the observed count;
all samples/seeds from a selected study move together, and the identical draw
matrix is applied to both geometries and every method. R type-7 percentiles are
used.

Only the four primary MRR DID tests receive p-values. Two-sided paired
study-block sign flips use 9,999 replicates and seed `20260815`, followed by
Holm adjustment across exactly four tests. The finite-replicate correction is
`(b + 1) / (B + 1)` and tolerance-equivalent absolute null statistics count as
exceedances. Other p-values are prohibited. Fewer than four contributing study
blocks yields descriptive output with
`insufficient_independent_study_blocks` and no interval or p-value.

## Clustering identifiability

MV5-AF contains directed held-out-query/training-reference pairs only. It has
zero within-training nested192 distances, while compatible clustering requires
262,675 within-training biological pairs per representation (525,350 across
both representations before component expansion), training matrices, frozen k
selection, partitions, and held-out assignments.

Therefore nested192 retrieval sensitivity is identifiable, but nested192 clustering
is not identifiable from the accepted artifacts. Imputation, symmetrizing the
directed rows, or reusing an incompatible Euclidean partition is prohibited.
No new clustering calculation is authorized in MV5-AG.

## Atomicity, validation, repeat, resume, and resources

The evaluation queue has exactly 150 groups and remains unauthorized in this
sprint. A later execution must publish immutable prediction/result/status
pairs atomically and reuse a completed unit only when all source, algorithm,
label, and contract identities match. Partial or stale units are quarantined.

Independent validation must reconstruct source hashes, all 150 group
identities, all eight pairing axes, canonical ranks/ties, the structural label
join, all 7,200 endpoint rows, all 24 estimands, aggregation, bootstrap draws,
sign-flip nulls, exceedance counts, and Holm adjustment without production
scientific helpers. Four minimum/maximum-pair groups at seed `20260805` form
the prospective private repeat. A clean public assembly must be byte-identical;
a full resume must preserve paths, hashes, sizes, and timestamps.

Execution limits are one worker, 300 seconds and 4 GiB process-tree RSS per
group, two aggregate worker-hours, and 1 GiB public output. A breached limit
blocks acceptance and does not authorize a method change.

## Decision and boundary

MV5-AG may approve only a later, separately committed, prediction-locked
nested-192 retrieval evaluation. It does not execute or authorize rankings,
labels, endpoints, nested 256, clustering, gene topology,
cell/gene fusion, new data, Rust optimization, package defaults, manuscript
claims, pushing, confidential-material tracking, PDF tracking, or
`example_run.r`.

Any accepted hash drift, pairing failure, leakage, post-hoc contract change,
unfreezable estimator, or inability to specify atomic repeat/resume/validation
rules is a major stop condition.

The public acceptance-criteria ledger must make the eight prefreeze gates
explicit: prospective commit identity, complete source binding, exact
coordinate-source nesting, complete prediction-axis pairing, structural label
firewall, fixed analysis panel, clustering disposition, and zero-execution
boundary. All eight must pass before a later runner can be authorized.
