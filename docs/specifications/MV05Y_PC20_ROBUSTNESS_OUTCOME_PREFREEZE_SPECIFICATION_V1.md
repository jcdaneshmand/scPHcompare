# MV5-Y PC20 robustness-outcome prefreeze specification v1

Date frozen: 2026-08-11

Accepted PC20 calculation base: `f69c6e8`

Outcome state: not computed

## Purpose and scientific question

MV5-Y prospectively defines how the accepted 20-coordinate sensitivity
calculation will be compared with the accepted 30-coordinate SCT and inductive
integrated retrieval benchmarks. The narrow question is whether truncating the
same 384-cell Euclidean coordinate views from 30 to 20 coordinates changes
cross-study tissue retrieval or changes topology's increment relative to the
matched same-coordinate energy baseline.

This sprint freezes sources, pairing, endpoints, estimands, aggregation,
uncertainty, multiplicity, reporting, validation, and execution guards. It does
not sort a PC20 distance, calculate a neighbor rank, open a label against a
PC20 prediction, calculate an endpoint, compare a result, select a method, or
authorize execution.

The PC20 configuration was selected by the outcome-free MV5-T gate and executed
without labels in MV5-X. Accepted 30-coordinate marginal results existed before
this specification, so this contract is prospective for every PC20 result and
paired PC20 contrast but is not described as wholly blind to the earlier
baseline literature. Complete fixed-panel reporting and immutable source hashes
prevent the known baseline from selecting a favorable PC20 method or endpoint.

The dissertation, paper, and confidential reviewer material were consulted
privately only to recover generalized requirements: retain separate H0 and H1,
use comparable projections and matched baselines, separate biological from
technical interpretations, report sensitivity and failures completely, and
make method details reconstructable. No private document or reviewer wording is
copied into this contract or Git.

## Exact compatibility result

The accepted baseline and PC20 artifacts share exactly:

- 90 biological samples in 15 leave-one-study-out folds;
- five fixed technical cell-subsample seeds (`20260805` through `20260809`);
- 384 matched cells per sample realization;
- the held-out-query to training-reference direction;
- 35,350 biological query/reference pairs per representation;
- separate raw exact all-active-level H0 and H1 landscape distances;
- the unscaled raw H0/H1 Euclidean composite as descriptive output; and
- an energy baseline calculated in the same coordinate system as its topology.

All eight representation/family axes match exactly on fold, seed, query sample,
training sample, and mapped method role. This establishes scientific pairing;
it does not imply equal distances or equal outcomes.

The coordinate comparison is one-factor-at-a-time within representation:

| Representation | Accepted baseline | Sensitivity |
|---|---|---|
| SCT fold-fitted | 384 cells, 30 training-fit PCs, Euclidean | same realization and fit, first 20 PCs |
| Inductive integrated | 384 cells, 30 inductive integrated coordinates, Euclidean | same realization and fit, first 20 coordinates |

SCT and integrated results remain separate strata. MV5-Y does not combine them
into a single average or treat representation as a technical replicate.

## Frozen method map

Each representation contains four paired families:

1. H0 exact landscape distance — confirmatory topology component;
2. H1 exact landscape distance — confirmatory topology component;
3. raw `sqrt(H0^2 + H1^2)` — descriptive uncalibrated composite; and
4. energy divergence — matched same-cell, same-coordinate baseline.

The accepted SCT method IDs are reused directly. The generic PC20 method IDs
for the inductive integrated rows are mapped to the corresponding accepted
`integrated_*` baseline IDs before pairing. That name mapping changes no value.

Pseudobulk is not in the MV5-X PC20 calculation because PCA truncation cannot
change the frozen 500-gene pseudobulk baseline. It remains an identity/context
control in the accepted evaluations but is not invented as a duplicated PC20
row. H0 and H1 remain separate; the raw composite cannot override them.

## Endpoints, eligibility, and missingness

The endpoints are inherited unchanged from MV5-E and MV5-K:

- primary: `cross_study_tissue_mrr_v1`; and
- supportive: `cross_study_tissue_1nn_balanced_accuracy_v1`.

At execution, immutable PC20 distances will be ordered within fold, seed,
method, and query by ascending distance, with exact ties broken by canonical
training sample ID. Labels may be joined only after all source and ranking
identities are written and validated. For each query, the first same-tissue
rank yields reciprocal rank; fixed rank one yields the supportive correctness
indicator. No post-label reranking is permitted.

The accepted eligibility states are retained:

- `estimable`;
- `single_study_tissue_not_estimable`;
- `training_tissue_absent_not_estimable`;
- `baseline_or_pc20_prediction_axis_mismatch`;
- `immutable_source_hash_mismatch`; and
- `software_failure`.

Every expected query/method/endpoint row must exist with an estimate or a
structured disposition. Available-case repair, replacement folds, replacement
seeds, and favorable-case omission are prohibited.

## Estimands and directionality

For endpoint value `Y`, representation `R`, method family `m`, and coordinate
count `p`, define the aggregated method value `Y(R,m,p)`.

Sixteen direct sensitivity estimands are frozen:

`DIRECT(R,m,Y) = Y(R,m,20) - Y(R,m,30)`

for two representations, four families, and two endpoints.

Eight topology-increment difference-in-differences estimands are frozen:

`DID(R,d,Y) = [Y(R,d,20) - Y(R,energy,20)] - [Y(R,d,30) - Y(R,energy,30)]`

for two representations, `d` in H0/H1, and two endpoints.

Positive means the named quantity increased at PC20, negative means it
decreased, and zero means it was unchanged. For DID, positive means topology's
increment relative to matched energy increased at PC20; it does not by itself
mean topology exceeds energy at either coordinate count.

The four H0/H1 DID estimands on MRR are the only confirmatory sensitivity
family. Direct MRR changes are secondary diagnostics. All fixed-1-NN changes
are supportive. The complete registry has exactly 24 rows.

No scientifically justified equivalence or noninferiority margin has been
established. Therefore MV5-Y must not turn a nonsignificant difference into an
equivalence, invariance, noninferiority, or robustness-success claim. It reports
signed changes, uncertainty, and full heterogeneity; a later claim gate may
define a practical margin only before seeing the relevant results.

## Aggregation and inferential unit

Every endpoint and contrast uses the accepted order:

1. average five technical seeds within biological sample;
2. average biological samples within tissue; and
3. average the five tissue means equally.

The biological sample is the inferential unit and held-out study is the
resampling block. Seeds represent repeated cell realizations and never increase
the sample or study count. Fold-specific values are retained for heterogeneity
diagnostics; overlapping training folds are not treated as independent
biological replicates.

## Uncertainty and multiplicity

Every one of the 24 estimands receives a paired 95% percentile interval from
2,000 tissue-stratified held-out-study bootstrap replicates using seed
`20260812`. Within each tissue, studies are sampled with replacement at the
observed study count; every sample and all five seeds from a selected study
move together. The identical block-count matrix is applied to PC20, baseline,
every method, and every estimand. Percentiles use R `quantile` type 7.

The sole p-value family is `F1_pc20_topology_increment_mrr_four_tests`: SCT-H0,
SCT-H1, integrated-H0, and integrated-H1 DID on MRR. Two-sided paired
study-block sign flips use 9,999 replicates and seed `20260813`, with
`(b + 1) / (B + 1)`. Absolute null statistics within
`64 * .Machine$double.eps * max(1, abs(null), abs(observed))` of the observed
absolute statistic count as exceedances. Holm adjustment covers exactly the
four tests. No other p-values are authorized.

If fewer than four independent study blocks contribute to an estimand, its
point estimate remains descriptive, inferential values are omitted, and the
status is `insufficient_independent_study_blocks`.

## Label and accepted-outcome firewall

The external historical metadata is SHA-256 bound but remains untracked. This
prefreeze reads only sample and study keys after verifying its hash; tissue,
approach, and other label-value columns are skipped. It verifies the exact
90-sample/15-fold key join without joining a label to a PC20 distance.

At later execution, only tissue is permitted for these endpoints. Approach and
all other outcomes remain prohibited. Accepted MV5-E and MV5-K query-endpoint
files are hash-bound but not opened in this prefreeze. A later runner may read
them only after the PC20 ranking lock is durable, and must verify their hashes,
pair keys, denominators, no-refit flags, and no-rerank flags before paired
comparison.

## Queue and atomic execution contract

The queue contains 150 groups: 75 SCT and 75 inductive integrated groups. Each
group binds the PC20 private artifact manifest, the accepted baseline group
identity, fold, seed, representation, query/training counts, 4 method families,
and 2 endpoints. Expected later output is 3,600 query/method rows or 7,200
query/method/endpoint rows before aggregation. Every queue row has
`execution_authorized=FALSE` in this sprint.

A later execution must write atomic private result/status pairs. Existing
completed units are reused only when every source, algorithm, label, and
contract identity matches. Missing halves, stale identities, or hash mismatch
cause refusal rather than silent rebuild or overwrite.

## Validation, repeat, resume, and resources

Before labels open, a later execution must independently verify all source
hashes, all 150 group identities, all eight exact prediction axes, canonical
ranking/ties, and the external label identity. After execution, an independent
implementation must reconstruct all 7,200 endpoint rows, all 24 estimands,
aggregation, the paired bootstrap matrices and intervals, sign-flip nulls,
exceedance counts, and Holm adjustment without calling production helpers.

Minimum/maximum-pair groups for each representation at seed `20260805` are the
four prospective private repeats. A clean full public assembly must reproduce
all public artifacts byte for byte. A completed 150-group resume must reuse all
private result/status pairs with unchanged paths, hashes, sizes, and timestamps.

Execution limits are one worker, 300 seconds and 4 GiB process-tree RSS per
group, two aggregate worker-hours, and 1 GiB public output. A guard breach
blocks acceptance pending audit; it does not authorize a changed method.

## Why clustering is excluded

MV5-X contains only directed held-out-query to training-reference distances.
It contains zero within-training PC20 pairs and therefore cannot reconstruct
the training distance matrices, frozen k selection, training partitions, or
held-out cluster assignments required by MV5-R/S. The missing scope is 262,675
within-training biological pairs per representation (525,350 total before
method/component expansion).

Consequently PC20 retrieval sensitivity is identifiable from MV5-X, but PC20
clustering sensitivity is not. Reusing directed retrieval rows as a full
matrix, imputing missing entries, or comparing incompatible cluster artifacts
is prohibited. Clustering would require a separately justified, label-closed
PC20 matrix-and-clustering calculation sprint before another outcome prefreeze.

## Reporting and interpretation boundary

All 24 estimands must be published in registry order, including null, negative,
heterogeneous, tied, failed, and non-estimable rows. No winner-only display,
outcome-driven method promotion, seed removal, tissue removal, or selective
configuration continuation is permitted.

MV5-Y can later support a bounded statement about sensitivity of the existing
cross-study tissue-retrieval benchmark to 20 versus 30 coordinates. It cannot
establish equivalence, general PH superiority, clustering robustness,
biological mechanism, external validity, gene topology, fusion benefit, or
manuscript readiness.

## Stop boundary and decision

The compatible retrieval contract is approved for a later, separately bound
prediction-locked execution. Execution is not authorized now. PC20 clustering
from MV5-X is rejected as non-identifiable. The other three robustness
configurations, spectral methods, gene topology, cell/gene fusion, new data,
optimization/Rust, package defaults, manuscript claims, PDF/reviewer tracking,
`example_run.r`, pushing, and author-credit decisions remain outside this
sprint.
