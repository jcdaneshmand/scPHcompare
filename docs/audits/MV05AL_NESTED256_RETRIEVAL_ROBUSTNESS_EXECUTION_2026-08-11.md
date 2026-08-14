# MV5-AL prediction-locked nested-256 retrieval robustness execution

Date: 2026-08-11

Accepted prefreeze: `9bcf650`

Prospective engine: `fdd8ce1`

Durable prediction lock: `d889838`

## Result

MV5-AL completed the frozen retrieval-only nested-256 cell-depth sensitivity.
The label-closed phase produced and independently reconstructed all 282,800
ranking rows across 150 groups. The complete prediction lock was committed at
`d889838` before any tissue value or accepted baseline endpoint row was parsed.

Only then did the outcome phase read `Tissue.x`. `Approach.x` and every other
label remained unread. The execution produced 3,600 query-method outcomes,
7,200 query-endpoint rows, 10,800 query estimands, 2,160 sample estimands, 120
tissue estimands, 24 macro estimands and intervals, and exactly four primary
MRR tests. No method, representation, tissue, seed, or result was selected or
omitted.

## Primary result

The four predeclared topology-increment MRR changes for nested 256 minus the
accepted 384-cell calculation are:

| Representation | Dimension | Estimate | 95% blocked bootstrap interval | Raw p | Holm p |
|---|---:|---:|---:|---:|---:|
| Inductive integrated | H0 | -0.00951 | [-0.03399, 0.06169] | 0.7198 | 0.7758 |
| Inductive integrated | H1 | -0.04579 | [-0.06398, -0.00569] | 0.0149 | 0.0596 |
| SCT whole | H0 | 0.00808 | [-0.00340, 0.03157] | 0.3140 | 0.7758 |
| SCT whole | H1 | -0.01386 | [-0.04671, 0.00354] | 0.2586 | 0.7758 |

No primary test passes the frozen four-test Holm family. The integrated H1
estimate is directionally negative, its blocked bootstrap interval excludes
zero, and its raw randomization p-value is 0.0149, but its multiplicity-adjusted
p-value is 0.0596. It must therefore be reported as suggestive sensitivity
evidence, not a confirmatory detection. This analysis does not establish
equivalence, invariance, noninferiority, or a preferred default.

Direct MRR changes are also retained without selection: integrated
energy/H0/H1/raw are -0.00344/-0.01295/-0.04923/-0.01446, while SCT
energy/H0/H1/raw are -0.00129/+0.00680/-0.01515/+0.00981. All 24 fixed
estimands remain available in the public ledgers.

## Prediction firewall and validation

All 196 MV5-AK sources revalidated before ranking. Every rank used ascending
immutable distance and canonical training-sample ID radix tie-breaking. The
independent implementation reconstructed all 282,800 distances, ranks, tie
sizes, and 150 private artifact/status identities. Ten prediction-validation
categories pass.

Two clean prediction builds reproduced the complete ranking gzip, six stable
public ledgers, and all 150 private ranking payloads. A full resume reused all
150 groups and preserved all 300 private paths, hashes, byte sizes, and
timestamps. Every prediction and repeat row records labels closed and outcomes
absent.

The outcome validator independently reconstructed source and label boundaries,
all tissue endpoints, accepted 384-cell pairings, query/sample/tissue/macro
aggregations, bootstrap count and estimate matrices, sign matrix, four null
distributions, exceedance counts, raw p-values, and Holm adjustment. All 17
validation categories pass.

Two clean outcome builds reproduce 17 deterministic public/validation ledgers,
all 150 private outcome payloads, and the inference matrix exactly. Resource
ledgers differ only in measured runtime/RSS. A full outcome resume preserves
all 301 private paths, hashes, byte sizes, and timestamps and reuses the
inference matrix.

## Resources and recovery audit

The accepted prediction phase used 29.867 aggregate group-seconds, 1.193
seconds maximum per group, 905,879,552 bytes maximum RSS, 12,113,200 bytes of
private ranking payloads, and a 16,169,139-byte public ranking file. The
accepted outcome phase used 9.831 aggregate group-seconds, 0.681 seconds
maximum per group, 1,081,872,384 bytes maximum RSS, and 140,568 bytes of private
outcome payloads. All frozen caps pass.

After an execution interruption, the recovery audit confirmed the committed
engine HEAD, clean canonical scope, absent public outcomes, closed labels, and
two complete private prediction builds. It found two bookkeeping defects only
in an unexecuted local repeat helper: an unsupported `list.files()` argument
and path-bearing vector names in equality checks. Both were corrected before
the repeat ledger was created. No production code, prediction, outcome,
estimand, seed, or published artifact was affected.

## Boundary and next decision

MV5-AL does not authorize clustering, another configuration, gene/fusion
analysis, new data, Rust/default changes, or manuscript claims. The nested-256
robustness chain is now complete. The next sprint should be a
selection-resistant synthesis gate over all four complete configurations
(PC20, cosine, nested 192, and nested 256), with the entire fixed result panels
bound before choosing any further calculation or analysis. It must not select
the integrated-H1 near-threshold result, or any favorable representation,
dimension, endpoint, tissue, interval, or p-value.
