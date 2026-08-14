# MV5-S prediction-locked clustering-outcome execution

Date: 2026-08-10

Accepted contract base: `3c3271e`

Execution-source freeze: `4f7f73c`

Outcome status: complete, independently validated, secondary

## 1. Result

MV5-S executed the complete committed MV5-R contract without refitting,
reselecting `k`, changing algorithms, removing folds or seeds, calculating a
p-value, or selecting a method from outcomes. All 2,400 queued units completed:
1,800 training-alignment units and 600 held-out prediction units.

The public evidence contains 9,000 training seed metrics, 1,800 training unit
summaries, 120 descriptive training contexts, 3,000 held-out seed summaries,
600 held-out unit summaries, and 40 prespecified held-out macro contexts. The
private root retains 18,000 sample-seed held-out predictions and 3,600
sample-context summaries with true/predicted labels. Label values never enter a
public artifact.

## 2. Prediction lock and implementation

The prospective engine was committed as `c3f8da0` before any real outcome was
computed. A separate execution freeze at `4f7f73c` binds ten engine/source files
byte-for-byte. The pre-outcome focused gate passed 28/28 checks and the full
repository suite passed 831 checks with two expected CRAN-only skips.

PAM with the already selected stability `k` remains primary; average linkage at
the same `k` remains sensitivity. Training ARI and NMI use the exact frozen
partitions. NMI is max-entropy-normalized under aricode 1.0.3. Held-out labels
are predicted only from training-cluster plurality, with lexical ties.

Tissue uncertainty uses 2,000 tissue-stratified study-block replicates.
Approach uncertainty uses 2,000 global study-block replicates because three
studies contain both approaches; mixed studies remain intact. All 2,000 tissue
replicates and 1,982/2,000 approach replicates are estimable in every context,
above the frozen 95% threshold.

## 3. Held-out tissue balanced accuracy

The table is in fixed representation/distance order, not performance order.
Each cell is estimate `[95% study-block percentile interval]`.

| Representation | Distance | Average linkage sensitivity | PAM primary |
|---|---|---:|---:|
| Integrated | Energy | 0.0910 [0.0400, 0.3107] | 0.2578 [0.2247, 0.4000] |
| Integrated | raw H0/H1 | 0.0191 [0.0000, 0.2062] | 0.0387 [0.0000, 0.2000] |
| Integrated | H0 | 0.0191 [0.0000, 0.2062] | 0.0387 [0.0000, 0.2000] |
| Integrated | H1 | 0.0555 [0.0047, 0.2407] | 0.0482 [0.0195, 0.1533] |
| Integrated | Pseudobulk | 0.0236 [0.0013, 0.2092] | 0.2263 [0.1449, 0.4872] |
| SCT | Energy | 0.0358 [0.0077, 0.2062] | 0.1631 [0.0813, 0.3800] |
| SCT | raw H0/H1 | 0.0814 [0.0556, 0.2646] | 0.0106 [0.0000, 0.1225] |
| SCT | H0 | 0.0783 [0.0556, 0.2615] | 0.0106 [0.0000, 0.1225] |
| SCT | H1 | 0.0103 [0.0000, 0.0800] | 0.0182 [0.0048, 0.0893] |
| SCT | Pseudobulk | 0.0236 [0.0013, 0.2092] | 0.2263 [0.1449, 0.4872] |

With five equally weighted tissue classes, 0.2 is the balanced-accuracy value
of a class-independent chance reference. MV5-S did not perform a test against
that value. Estimates span 0.0103 to 0.2578 and intervals are generally broad.
This is evidence of limited held-out tissue generalization under the frozen
cluster-to-label rule, not evidence for selecting or rejecting an individual
representation or distance.

The integrated H0 and raw H0/H1 rows are numerically identical in both
algorithms; SCT PAM is also identical for those two distances. This is
consistent with the already documented H0 dominance of the unweighted raw
combination, but MV5-S does not promote that descriptive recurrence into a new
claim.

## 4. Held-out approach balanced accuracy

Approach balanced accuracy is 0.4925 to 0.5000 across all 20 contexts. All PAM
contexts are 0.4950 or 0.5000; average-linkage contexts are 0.4925 to 0.5000.
Intervals stay close to 0.5. With two equally weighted approach classes, 0.5 is
the class-independent chance reference. No hypothesis test was performed. The
prespecified technical endpoint therefore shows essentially chance-level
held-out approach discrimination under the frozen mapping.

## 5. Descriptive training alignment

Training folds overlap heavily, so these are descriptive distributions rather
than independent generalization evidence. Across the 20 fixed contexts per
metric:

| Label | Metric | Minimum fold mean | Median fold mean | Maximum fold mean |
|---|---|---:|---:|---:|
| Tissue | ARI | -0.0066 | 0.2398 | 0.4796 |
| Tissue | NMI-max | 0.0295 | 0.2697 | 0.5395 |
| Study | ARI | -0.0176 | 0.1907 | 0.3815 |
| Study | NMI-max | 0.0240 | 0.2338 | 0.4676 |
| Approach | ARI | -0.0339 | 0.0988 | 0.1976 |
| Approach | NMI-max | 0.0052 | 0.0539 | 0.1078 |

The gap between some moderate training alignment values and weak held-out
tissue performance illustrates why the LOSO prediction endpoint was necessary.
It is not a post-hoc basis for method selection.

## 6. Validation, repeat, and resume

The separately committed validator passed 12 categories:

- 2,400/2,400 unit identities and artifact hashes;
- all 9,000 ARI/NMI values reconstructed with mclust/aricode;
- all 18,000 held-out predictions reconstructed from training labels only;
- all seed, unit, sample, context, macro, and interval aggregations;
- both 2,000-replicate block-count matrices and every replicate-estimate hash;
- all public schemas and label-value scans.

A clean production repeat reproduced 2,400/2,400 private outcome artifacts,
8/8 outcome-bearing public tables, and the private sample summary byte-for-
byte. The production resume reused 2,400/2,400 units with zero execution; all
4,800 artifact/status paths, hashes, sizes, and timestamps remained unchanged.

Final staged-evidence verification parses 5/5 MV5-S scripts and passes 831
repository checks with zero failures/errors/warnings and two expected CRAN-only
skips. The isolated standard package check is `Status: OK`. Its candidate
source tarball SHA-256 is
`6e7f85036dff85ac8ee13c1c34cd121be6b4923ed96b8ad89b629c290c2cf667`;
the check-log SHA-256 is
`3087e3407ce0a5cf4f676940389270d283302e281f8e5e6071d403a331ca1bd9`.
Those hashes precede only this verification paragraph and the resulting
manifest refresh; no R code, result, or test changed afterward.

## 7. Resources and public safety

The first pass used 51.888 aggregate unit-seconds, at most 0.570 seconds per
unit, and at most 740,687,872 bytes process RSS. Private unit outcomes occupy
2,211,948 bytes. Runner-measured public outputs occupy 6,532,590 bytes. These
are below every frozen cap. The complete runner wall time was 325 seconds,
including package loading, 4,800 atomic writes, aggregation, bootstrap, and
public output generation.

The external label file remains untracked. Public tables contain no
`true_label` or `predicted_label` column and none of the five tissue or two
approach label values. `example_run.r`, PDFs, reviewer material, and private
production/repeat roots remain untracked.

## 8. Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | approve |
| Correctness demonstrated? | pass |
| Computation feasible? | yes |
| Biological interpretation permitted? | secondary, complete-reporting only |
| Next action | freeze an outcome-aware but selection-resistant robustness/gap gate before any robustness execution |

MV5-S completes the originally selected clustering benchmark action. The next
sprint should not tune methods around these values. It should audit the already
planned robustness families, identify which perturbations remain scientifically
identifiable and affordable, and freeze any execution without choosing a
method/representation from MV5-S outcomes. Spectral promotion, gene topology,
cell/gene fusion, new data, optimization/Rust, package-default changes, and
manuscript claims remain closed.
