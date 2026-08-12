# MV5-AL prediction-locked nested-256 retrieval-robustness execution specification v1

Date frozen: 2026-08-11

Contract base: `9bcf650a1f4b5910403cc32ecad41c5bccc4a9d1`

Outcome state at engine freeze: not computed

## Authorized analysis

MV5-AL implements exactly the accepted MV5-AK retrieval-only contract. It
first constructs canonical label-closed nested-256 rankings for all 150 groups,
independently reconstructs them, and commits their identities durably. Tissue
labels and accepted MV5-E/K endpoint rows may be read only after that commit.

The analysis compares the accepted 384-cell calculation with its deterministic
nested-256-cell subset while retaining all 30 coordinates and Euclidean
geometry. Its first 192 cells are the exact accepted nested-192 prefix, and all
256 are a strict subset of the same 384-cell realization. SCT and
inductive-integrated representations remain separate. H0 and H1 remain
separate exact all-active-consecutive-level landscape components; the raw
composite is descriptive; energy is the matched same-coordinate baseline.

## Prediction-lock phase

All 196 MV5-AK sources and all 150 MV5-AJ group manifests must revalidate
before ranking. Each group publishes an atomic private ranking/status pair.
Existing pairs are reusable only when queue, source-freeze, implementation,
and artifact identities match exactly. Partial or stale pairs abort.

Rankings contain all and only training references for each representation,
fold, seed, method, and held-out query. Ordering is ascending immutable
distance and then ascending canonical training sample ID using radix order.
Exact tie sizes are retained. Nonfinite or negative distances, duplicate keys,
missing references, denominator drift, or post-label ranking abort execution.

The complete label-closed lock comprises 150 groups and 282,800 ranking rows.
An independent implementation must reconstruct every distance identity, rank,
and exact tie size from MV5-AJ method rows without calling production helpers.
The lock, independent evidence, and source/implementation identities must be
committed locally before any tissue value or accepted endpoint row is opened.

## Outcome phase

Only after the prediction-lock commit may the runner revalidate all sources,
parse `Tissue.x`, and open the hash-bound accepted Euclidean query endpoints.
`Approach.x` and all other outcome columns remain unread. No refit, rescaling,
reranking, method removal, denominator change, or held-out-derived selection is
allowed after label access.

The complete axes are:

- 3,600 nested-256 query-method outcomes;
- 7,200 long query-method-endpoint rows;
- 10,800 query estimand rows;
- 2,160 sample estimands;
- 120 tissue estimands;
- 24 macro estimands and 24 intervals; and
- exactly four Holm-adjusted primary MRR DID tests.

The two endpoints, 16 direct nested-256-minus-384 cell-depth changes, and eight
topology-increment DIDs are exactly those in the MV5-AK registries. Every row
must be reported or receive its frozen structured failure disposition.

## Aggregation and inference

Five technical seeds are averaged within biological sample, samples are
averaged within tissue, and the five tissues are weighted equally. Biological
sample is the inferential unit; held-out study within tissue is the resampling
block.

All 24 estimands receive paired 95% type-7 percentile intervals from exactly
2,000 tissue-stratified held-out-study bootstrap replicates using seed
`20260814`. The same block-count matrix applies to both cell depths and every
method/estimand.

Only the four H0/H1 MRR topology-increment DIDs receive p-values. Exactly 9,999
paired study-block sign flips use seed `20260815`; the two-sided finite-sample
p-value is `(b + 1) / 10000`, numerical-boundary ties count as exceedances, and
Holm adjustment covers exactly those four tests.

No equivalence or noninferiority margin exists. Null or adjusted-nonsignificant
results must not be described as equivalence, invariance, robust success,
superiority, or a default-setting result.

## Independent validation and reproducibility

Without production scientific helpers, validation must reconstruct every
ranking, tissue endpoint, Euclidean pairing, query/sample/tissue/macro
aggregation, deterministic bootstrap and sign matrix, interval, raw p-value,
Holm value, and private hash.

The four prospective private repeats are the minimum- and maximum-pair fold for
each representation at seed `20260805`. A clean full public outcome assembly
must be byte-identical. Full prediction and outcome resumes must preserve all
private artifact/status paths, hashes, sizes, and timestamps.

Limits are one worker, 300 seconds and 4 GiB process-tree RSS per group, two
aggregate worker-hours for each phase, and 1 GiB public output. A breach blocks
acceptance; it does not authorize changing the method.

## Stop boundary

MV5-AL excludes nested-256 clustering and within-training pairs, alternate
configurations, coordinate/PH/landscape/energy recalculation, approach or other
labels, method/configuration/tissue/seed selection, equivalence claims,
gene/fusion/new-data/optimization/Rust/default/manuscript-claim work, private
source/PDF/reviewer tracking, pushing, and `example_run.r`.

Any source drift, rank mismatch, premature label access, endpoint-pairing
failure, inference mismatch, atomicity failure, non-reproducible repeat/resume,
or estimator-invalidating scientific issue is a major stop condition.
