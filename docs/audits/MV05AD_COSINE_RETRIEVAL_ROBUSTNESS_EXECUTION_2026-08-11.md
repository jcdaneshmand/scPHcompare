# MV5-AD prediction-locked cosine retrieval-robustness execution audit

Date: 2026-08-11

Accepted contract: MV5-AC commit `260ac6d`

Engine freeze: `c22d667`

Prediction-lock commit: `2e2f9f3`

Status: complete

## Executive result

The fixed cosine-chord sensitivity analysis completed exactly as prespecified.
All four primary changes in topology's MRR increment relative to matched energy
were negative. SCT H0 and H1 survived the four-test Holm adjustment; the two
integrated contrasts did not.

| Representation | Component | Cosine-minus-Euclidean topology-increment MRR | 95% paired block CI | Raw p | Holm p |
|---|---|---:|---:|---:|---:|
| Inductive integrated | H0 | -0.07216 | [-0.21895, 0.00552] | 0.2734 | 0.2734 |
| Inductive integrated | H1 | -0.10106 | [-0.21066, -0.02869] | 0.0526 | 0.1052 |
| SCT whole | H0 | -0.11833 | [-0.20089, -0.06917] | 0.0070 | 0.0224 |
| SCT whole | H1 | -0.10088 | [-0.15610, -0.07495] | 0.0056 | 0.0224 |

These are geometry-sensitivity estimates. A negative DID means that replacing
Euclidean point geometry with cosine chord reduced topology's retrieval
increment relative to energy. It is not an equivalence, universal Euclidean
superiority, biological mechanism, or default-setting result.

## Prediction firewall

The prospective engine was committed at `c22d667`. It revalidated all 187
MV5-AC sources and then constructed all 282,800 label-closed cosine ranking
rows in 150 groups. Ranking was ascending immutable distance, then canonical
training sample ID in radix order. There were zero exactly tied rows in this
dataset, but the tie policy and tie sizes were still independently verified.

The independent prediction validator reconstructed all 282,800 distances,
ranks, and tie sizes without production helpers and rehashed all 150 private
artifact/status pairs. Labels and outcomes remained zero. Four prospectively
selected private ranking artifacts repeated byte-for-byte, and all 300 private
ranking unit files retained identical paths, hashes, sizes, and timestamps on
full resume.

The prediction lock and its validation were committed at `2e2f9f3`. Tissue
values and accepted Euclidean endpoint rows were not opened before this commit.

## Outcome execution

After the durable lock only, the runner revalidated all 187 sources, parsed
only `Tissue.x`, and evaluated the fixed retrieval endpoints. `Approach.x` and
all other label columns remained unread. No coordinate, PH, landscape, energy,
or ranking value was refit or recomputed.

The complete output contains:

| Object | Rows |
|---|---:|
| Cosine query-method outcomes | 3,600 |
| Long query-method endpoints | 7,200 |
| Paired cosine/Euclidean query-method rows | 3,600 |
| Query estimands | 10,800 |
| Sample estimands | 2,160 |
| Tissue estimands | 120 |
| Macro estimands and intervals | 24 each |
| Primary tests | 4 |

All expected rows were estimable. Five technical seeds were averaged within
sample, samples within tissue, and five tissues equally. The 2,000-replicate
paired tissue-stratified held-out-study bootstrap used seed `20260814`.
Exactly four MRR H0/H1 DIDs used 9,999 paired study-block sign flips with seed
`20260815` and Holm adjustment.

## Direct geometry changes

The primary DID pattern is explained by the direct MRR changes. Under cosine
chord, integrated energy changed by +0.00124 while integrated H0/H1 changed by
-0.07092/-0.09982. SCT energy changed by +0.01208 while SCT H0/H1 changed by
-0.10625/-0.08880. Thus the negative DIDs are not created by a shared decline
of both topology and its matched energy baseline.

This supports a bounded conclusion: radial magnitude in the accepted cell
coordinate systems carries retrieval information used by the topological
distances, particularly in SCT. The result does not isolate which biological
or technical feature supplies that magnitude information and does not license
an outcome-driven method replacement.

## Independent validation and reproducibility

The helper-independent outcome validator passed 15 categories and reconstructed:

- every canonical ranking and tissue-only label join;
- all 3,600 cosine endpoint objects and 3,600 accepted-baseline pairings;
- all query, sample, tissue, and macro estimands;
- the exact 2,000 by 15 bootstrap count matrix and 2,000 by 24 replicate matrix;
- all 24 type-7 intervals;
- the exact 9,999 by 15 sign matrix, four null distributions, exceedance counts,
  raw p-values, and Holm values; and
- all 150 private outcome artifact/status identities.

Four prospective private outcome artifacts, 16 deterministic runner outputs,
the private inference matrices, and the independent validation ledger repeated
byte-for-byte (22 recorded comparisons). All 300 private outcome unit files
retained identical paths, hashes, sizes, and timestamps through full resume.

Prediction ranking used 28.589 aggregate group-seconds, at most 1.257 seconds
per group and 910,077,952 bytes RSS. Outcome evaluation used 10.209 aggregate
group-seconds, at most 0.688 seconds per group and 1,010,532,352 bytes RSS. All
caps passed.

## Scientific boundary and next action

The accepted evidence argues against cosine chord as a replacement geometry
for the current topological retrieval benchmark, especially for SCT, while
remaining a valuable negative sensitivity result. It does not prove that
Euclidean geometry is globally optimal, that PH is superior to matched
baselines, or that the result transfers to clustering, gene topology, fusion,
or new data.

Cosine clustering remains non-identifiable because MV5-AB contains no
within-training cosine distances. Both nested-cell configurations remain
unexecuted. The next sprint must be a selection-resistant continuation gate
that binds the complete cosine result and canonical MV5-V order before deciding
whether nested 192 cells is justified. No subgroup, tissue, representation,
component, interval, or p-value may select that continuation.
