# MV5-AK nested-256 robustness-outcome prefreeze

Date: 2026-08-11

Accepted calculation base: `5e88fc149a2a9ef88de8af89372f7ca793b0cf41`

Prospective engine: `5c882d513bc8547214c49cc7338f0366547b78be`

## Decision

MV5-AK passes and approves only a later, separately committed,
prediction-locked retrieval evaluation comparing the accepted 384-cell and
nested-256-cell Euclidean calculations. It does not execute that evaluation.
Ranks, tissue-label joins, endpoints, estimates, p-values, method selection,
clustering, alternate configurations, and new calculations remain zero.

## Exact scientific comparison

The contrast changes only deterministic cell depth: 384 cells versus the first
256 cells in the accepted sample/seed/cell-ID SHA-256 order. The first 192 of
those 256 cells are byte-identical to the accepted nested-192 selection, and
all 256 are a strict subset of the same frozen 384-cell realization. Both sides
retain the same 30 training-fit coordinates, Euclidean point metric, fold,
seed, query, reference, representation, exact all-active-consecutive-level
landscape definition, and matched energy baseline. SCT and
inductive-integrated representations remain separate strata; H0 and H1 remain
separate methods.

Every one of the 150 groups independently cross-links its nested source to the
exact accepted 384-cell coordinate source. SCT uses the MV5-D1 private-cache
hash and matching MV5-D5 fold key. Integrated groups use the MV5-J-bound MV5-G
private-file hash. The private MV5-AJ source identity repeats the same source
hash. The independent MV5-AJ geometry evidence additionally proves the exact
192 prefix and 256-of-384 subset for every group. All 150 cross-links and both
nesting identities pass.

## Pairing, prediction lock, and firewall

All eight representation/method-family axes pair one-to-one over the complete
fold, seed, held-out query, training reference, and method key. There are
282,800 baseline rows and 282,800 nested-256 rows, with zero missing, excess,
or duplicate keys. Distances were not compared and ranks were not computed.

A later runner must rank ascending immutable distance and break exact ties by
ascending canonical training sample ID using radix order. Its complete
282,800-row prediction artifact must be durably committed before tissue access.
The current prefreeze reads only `orig.ident` and `SRA` from the hash-bound
external metadata, proving 124 unique source keys, 90 analysis samples, and 15
held-out-study mappings. `Tissue.x`, `Approach.x`, accepted endpoint values,
the completed nested-192 outcome panel, and all other outcomes remain unread.

## Frozen analysis

The fixed panel inherits the accepted endpoint and reporting structure without
selection: two endpoints, 16 direct nested-256-minus-384 cell-depth changes,
and eight topology-increment difference-in-differences estimands. All 24 rows
must be reported. The four H0/H1 MRR DIDs form the sole confirmatory family and
receive Holm adjustment. Paired uncertainty uses 2,000 tissue-stratified
held-out-study bootstrap draws; the four primary tests use 9,999 paired
study-block sign flips. Technical seeds are averaged within sample and never
treated as independent observations. No equivalence, invariance,
noninferiority, superiority, default, or manuscript conclusion is authorized.

## Clustering disposition

Retrieval is identifiable, but clustering is not. MV5-AJ contains directed
held-out-query-to-training-reference rows and zero nested-256 within-training
distances. A compatible two-representation clustering comparison would require
525,350 missing within-training biological pairs before component expansion.
Imputation, directed-row symmetrization, and reuse of an incompatible 384-cell
partition are prohibited.

## Validation and reproducibility

The production prefreeze emits 15 ledgers over 196 bound sources: 45 committed
sources, 150 private group manifests, and the external metadata identity. It
covers 150 groups, eight exact method axes, 24 estimands, and eight explicit
acceptance criteria. An independent implementation passes 12 categories,
including live source and private-manifest hashes, complete axis
reconstruction, explicit 192/256/384 nesting, the structural metadata-key
firewall, prediction-lock semantics, clustering disposition, and zero-outcome
closure.

Two clean assemblies at `5c882d5` reproduce all 15 ledgers byte-for-byte. The
two independent validators also reproduce byte-for-byte, for 16/16 exact
comparisons. A later execution must additionally create atomic
prediction/result units and preserve complete path/hash/size/timestamp
snapshots through validation-only resume.

## Next boundary

The next sprint is MV5-AL: implement and commit the exact prediction runner,
construct and independently validate all nested-256 ranks, and durably commit
the complete prediction lock before opening `Tissue.x`. Outcome execution must
remain a later atomic phase of that sprint. Clustering, other configurations,
gene/fusion, new data, Rust/default changes, manuscript claims, confidential
materials, PDFs, `example_run.r`, and pushing remain outside scope.
