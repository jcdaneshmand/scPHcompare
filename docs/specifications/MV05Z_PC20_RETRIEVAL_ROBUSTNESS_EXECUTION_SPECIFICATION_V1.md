# MV5-Z prediction-locked PC20 retrieval-robustness execution specification v1

Contract base: `111ef82`

Outcome state at engine freeze: not computed

MV5-Z implements and executes exactly the accepted MV5-Y retrieval-only
contract. It first constructs canonical label-closed PC20 rankings for all 150
groups and makes their hashes and completion records durable. Tissue labels and
accepted MV5-E/K endpoint rows may be read only after that prediction lock
passes. The runner then evaluates all four fixed method families and both fixed
endpoints, pairs them to the immutable 30-coordinate baseline, and reports all
24 prespecified estimands.

Every private group publishes an atomic artifact/status pair. Existing pairs
are reusable only when the queue, implementation, source-freeze, and artifact
hash identities are exact. Partial or stale pairs abort. Canonical rankings sort
by distance and then training sample ID; exact tie sizes are retained.

The required complete axes are 150 groups, 282,800 ranking rows, 3,600
query/method outcomes, 7,200 long query endpoints, 10,800 query estimand rows,
2,160 sample estimands, 120 tissue estimands, 24 macro estimands and intervals,
and exactly four Holm-adjusted MRR DID tests. Bootstrap and randomization
identities, formulas, seeds, blocking, directionality, and missingness are
unchanged from MV5-Y.

An implementation-independent validator must reconstruct every ranking from
MV5-X distances, every endpoint from immutable ranks and tissue labels, all
baseline pairings and estimands, both deterministic random matrices, every
interval and p-value, and every artifact hash without calling MV5-Z production
scientific helpers. Four frozen group repeats, a clean full public repeat, and a
full immutable resume are mandatory.

The stop boundary excludes PC20 clustering, the other three robustness
configurations, refitting or recalculation of PCA/PH/landscapes/energy,
post-label reranking, method or subgroup selection, approach outcomes,
equivalence claims, gene/fusion/new-data/optimization/Rust/default/claim work,
private-source tracking, pushing, and `example_run.r`.

## Execution record

The engine was frozen in `3756f7e`; the independently validated prediction lock
was committed in `c16f2b2` before tissue access. All 150 groups, 7,200 endpoint
rows, 24 estimands, 24 intervals, and four Holm-adjusted primary tests completed.
Independent reconstruction passed 15 categories, all 150 private artifacts and
16 deterministic public files repeated byte-for-byte, and all 300 private
result/status files remained unchanged through a full resume. See
`docs/audits/MV05Z_PC20_RETRIEVAL_ROBUSTNESS_EXECUTION_2026-08-11.md`.
