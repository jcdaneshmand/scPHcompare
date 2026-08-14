# MV5-AH prediction-locked nested-192 retrieval robustness execution

Date: 2026-08-11

Accepted prefreeze: `6f36741`

Prospective engine: `d9334bc`

Durable prediction lock: `1a197a8`

Immutable-resume correction: `41bc7c7`

## Result

MV5-AH completed the frozen retrieval-only cell-depth sensitivity. The
label-closed phase produced and independently reconstructed all 282,800
nested-192 ranking rows across 150 groups. That complete prediction lock was
committed at `1a197a8` before any tissue value or accepted baseline endpoint
row was parsed.

Only then did the outcome phase read `Tissue.x`. `Approach.x` and every other
label remained unread. The execution produced 3,600 query-method outcomes,
7,200 query-endpoint rows, 10,800 query estimands, 2,160 sample estimands, 120
tissue estimands, 24 macro estimands and intervals, and exactly four primary
MRR tests. No method, representation, tissue, seed, or result was selected or
omitted.

## Primary result

The four predeclared topology-increment MRR changes for nested 192 minus the
accepted 384-cell calculation are:

| Representation | Dimension | Estimate | 95% blocked bootstrap interval | Raw p | Holm p |
|---|---:|---:|---:|---:|---:|
| Inductive integrated | H0 | -0.01080 | [-0.03676, 0.05266] | 0.7511 | 1.0000 |
| Inductive integrated | H1 | -0.01115 | [-0.04821, 0.05495] | 0.7139 | 1.0000 |
| SCT whole | H0 | 0.00131 | [-0.01601, 0.02958] | 0.9225 | 1.0000 |
| SCT whole | H1 | -0.01366 | [-0.04051, 0.01803] | 0.5091 | 1.0000 |

No primary test detects a cell-depth change in topology's increment over the
matched energy baseline. This is sensitivity evidence for the current data and
benchmark, not proof of equivalence, invariance, noninferiority, or a preferred
default. Direct MRR changes are likewise small: integrated energy/H0/H1/raw are
-0.00495/-0.01576/-0.01610/-0.01316, while SCT energy/H0/H1/raw are
-0.00077/+0.00054/-0.01443/+0.00126. All 24 fixed estimands remain reported in
the public ledgers.

## Prediction firewall and validation

All 188 MV5-AG sources revalidated before ranking. Every rank used ascending
immutable distance and canonical training sample ID radix tie-breaking; the
independent implementation found zero tied rows and reconstructed all ranks,
distances, and 150 private hashes. Two clean prediction builds reproduced the
complete ranking gzip, prediction lock, six deterministic public ledgers, and
all 150 private ranking payloads. A 150-group resume preserved all 300 private
paths, hashes, sizes, and timestamps.

The outcome validator independently reconstructed all tissue endpoints,
accepted 384-cell pairings, query/sample/tissue/macro aggregations, the exact
bootstrap count and estimate matrices, sign matrix, four null distributions,
exceedance counts, raw p-values, and Holm adjustment. All 15 categories pass.

## Audited operational correction

The first diagnostic outcome run exposed a resume-only defect: an already
identical inference-matrix RDS would have been rewritten, changing its
timestamp. No estimand, ranking, random seed, formula, or result was changed.
The runner was corrected and prospectively committed at `41bc7c7`; accepted
outcomes were then regenerated in a fresh private root. Two corrected clean
runs reproduce 17 deterministic public ledgers, all 150 private outcome
payloads, and the inference matrices exactly. The resource ledger differs only
by observed timing. The corrected full resume preserves all 301 private paths,
hashes, sizes, and timestamps.

## Resources and boundary

The accepted prediction phase used 25.838 aggregate group-seconds, 1.174
seconds maximum per group, 897,257,472 bytes maximum RSS, and 12,097,200 bytes
of private ranking payloads. The accepted outcome phase used 9.631 aggregate
group-seconds, 0.629 seconds maximum per group, 999,366,656 bytes maximum RSS,
and 140,548 bytes of private outcome payloads. All frozen caps pass.

MV5-AH does not authorize clustering, nested 256, gene/fusion analysis, new
data, Rust/default changes, or manuscript claims. The next sprint must be a
selection-resistant continuation gate that binds the complete PC20, cosine,
and nested-192 result panels before deciding whether the still-closed nested
256 calculation is justified. It may not select a favorable representation,
dimension, endpoint, tissue, estimate, interval, or p-value.
