# MV5-Z prediction-locked PC20 retrieval-robustness execution audit

Date: 2026-08-11

Engine commit: `3756f7e`

Committed prediction lock: `c16f2b2`

## Outcome

MV5-Z completed the exact retrieval-only MV5-Y contract. The prediction stage
first produced and independently reconstructed 282,800 canonical PC20 ranking
rows across all 150 groups. The ranking SHA-256 is
`b73de06567e41980ae0320570cea90fcc37fc425c4d2fd2fad7d9ce95af96fc4`.
All 178 frozen sources passed immediately before the lock and again immediately
before label access. The lock was committed before tissue values or accepted
baseline endpoint rows were opened.

Only sample ID, study, and tissue were parsed. Approach and every other outcome
column were skipped. Labels were never used for PCA, PH, landscapes, energy,
distance ordering, tie handling, method selection, or any other fit.

## Complete execution axes

- 150/150 PC20 outcome groups completed;
- 3,600 query/method outcomes;
- 7,200 long query/method/endpoint rows;
- 3,600 exact baseline/PC20 query-method pairs;
- 10,800 query-estimand rows;
- 2,160 sample estimands, each averaging exactly five technical seeds;
- 120 tissue estimands;
- 24 macro estimands and paired 95% intervals; and
- exactly four MRR topology-increment DID sign-flip tests with Holm adjustment.

No missing, tied-nearest, failed, or non-estimable query row was silently
removed. All fixed registry rows, including null and negative changes, are
published.

## Primary PC20 sensitivity family

The signed estimand is the PC20-minus-PC30 change in topology's increment over
its matched energy baseline. Positive means the topology increment increased at
PC20; negative means it decreased.

| Representation | Component | Estimate | Paired 95% interval | Raw p | Holm p |
|---|---|---:|---:|---:|---:|
| Inductive integrated | H0 | -0.02036 | [-0.08386, 0.02537] | 0.5893 | 1.000 |
| Inductive integrated | H1 | -0.06444 | [-0.08203, -0.03374] | 0.0880 | 0.352 |
| SCT | H0 | -0.02567 | [-0.05116, -0.00042] | 0.1046 | 0.352 |
| SCT | H1 | 0.00537 | [-0.03537, 0.02673] | 0.8102 | 1.000 |

The blocked bootstrap and sign-flip procedures answer related but nonidentical
questions and can differ in small, heterogeneous blocked designs. Two paired
bootstrap intervals lie below zero, but none of the four prespecified
Holm-adjusted sign-flip tests is significant. The correct bounded conclusion is
that PC20 sensitivity is heterogeneous and includes evidence of degradation for
some topology/representation combinations; the analysis does not support a
claim that PC20 is equivalent, uniformly robust, superior, or a better default.

For context, the direct integrated H1 MRR change is -0.05917 with paired 95%
interval [-0.07738, -0.03527]. Direct changes are secondary and receive no
p-value under the frozen multiplicity contract.

## Independent validation

A standalone validator that does not call MV5-Z production scientific helpers
passed 15 categories. It independently reconstructed:

- all 282,800 canonical ranking rows and exact tie sizes;
- the 90-sample, 15-study, five-tissue join;
- all 3,600 query/method outcomes;
- all accepted baseline pairings and 10,800 query estimands;
- seed, sample, tissue, and macro aggregation;
- the exact 2,000 x 15 bootstrap count matrix and 2,000 x 24 estimate matrix;
- all 24 type-7 percentile intervals;
- the exact 9,999 x 15 sign matrix and four null distributions;
- all raw and Holm-adjusted p-values; and
- all 150 private outcome artifact hashes and status identities.

## Determinism, resume, and resources

All 150 private outcome artifacts are byte-identical in a clean fresh-root
repeat. All 16 deterministic public contract/result files are byte-identical;
the separately measured resource file is intentionally excluded from the clean
repeat comparison. The four prospectively named minimum/maximum-pair group
repeats pass explicitly.

A full same-root resume reused all 150 units. All 300 private result/status paths,
hashes, sizes, and timestamps were unchanged, and all 17 runner-produced public
files were byte-identical.

First-pass outcome-group work used 9.967 aggregate seconds, 0.652 seconds
maximum per group, and 1,012,342,784 bytes maximum process RSS. All frozen caps
passed.

## Boundary and next action

MV5-Z does not calculate PC20 clustering: MV5-X still lacks the required
within-training distance matrices. It does not authorize the other three
robustness configurations, select a method, set an equivalence margin, or alter
package defaults. Gene topology, cell/gene fusion, new data, optimization/Rust,
and manuscript claims remain closed.

The next reasonable sprint is an outcome-informed robustness continuation gate.
It should reconcile the heterogeneous PC20 result with the selection-resistant
MV5-T ordering and decide prospectively whether executing the next already
frozen one-factor-at-a-time configuration is scientifically useful and
computationally justified. It must not choose the next axis or method because
of a favorable PC20 subgroup.
