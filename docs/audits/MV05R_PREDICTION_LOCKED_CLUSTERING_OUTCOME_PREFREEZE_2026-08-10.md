# MV5-R prediction-locked clustering-outcome prefreeze

Date: 2026-08-10

Accepted base: `f16321c` (`Complete MV5-Q clustering artifact production`)

Execution state: prospective contract complete; real outcomes not computed

## 1. Result

MV5-R freezes a label-safe evaluation contract for the immutable MV5-Q
clustering artifacts. The contract binds 18 sources, including the external
metadata file by basename and SHA-256 without copying it into the repository;
defines eight endpoints for two frozen clustering algorithms; and creates an
immutable 2,400-unit execution queue spanning all 150 analysis groups.

This sprint inspected metadata structure and verified the exact candidate
sample join, but it did not calculate, expose, rank, or interpret any ARI, NMI,
balanced-accuracy, biological, or technical outcome.

## 2. Frozen external label source

| Property | Frozen value |
|---|---|
| Basename | `joined_metadata_cellcounts.csv` |
| SHA-256 | `e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0` |
| Source rows / studies | 124 / 18 |
| MV5 candidate samples / studies | 90 / 15 |
| Tissue classes | 5 |
| Approach classes | 2 |
| Studies containing both approaches | 3 |
| Exact MV5-Q sample-axis match | yes |
| Copied or tracked | no |

The design audit records only counts, identities, hashes, and structural
properties needed to prove evaluability. Label values are not copied into any
public MV5-R artifact.

## 3. Why training and held-out endpoints differ

LOSO training partitions overlap across folds. Training ARI/NMI therefore
describe foldwise alignment and are not independent replications. They are
reported for tissue, study, and approach without pooling fold-specific cluster
IDs or making an independence claim.

Held-out cluster IDs have meaning only through their frozen training cluster.
For each seed, algorithm, and analysis group, MV5-R will learn a cluster-to-
label map from training samples only. The predicted label is the plurality
training label in the assigned cluster; exact count ties use lexicographic
label order. The held-out label never participates in the map. Tissue and
approach are evaluable this way. Held-out study prediction is prohibited
because the LOSO study is absent from training and cannot be learned as a
cluster label.

## 4. Endpoint registry

PAM with the already selected stability-based `k` is primary. Average linkage
at that same frozen `k` is a sensitivity analysis. Refit, oracle `k`, and
outcome-driven tuning are prohibited.

| Scope | Label | Metrics | Role |
|---|---|---|---|
| Overlapping training partition | tissue | ARI, NMI | biological alignment |
| Overlapping training partition | study | ARI, NMI | technical alignment |
| Overlapping training partition | approach | ARI, NMI | technical alignment |
| Held-out prediction from training plurality | tissue | balanced accuracy | biological generalization |
| Held-out prediction from training plurality | approach | balanced accuracy | technical generalization |

Training endpoints retain seed rows, a seed mean within each fold, and the
complete fold distribution. They receive descriptive seed-jackknife
uncertainty only and no p-values. Held-out correctness is first averaged over
the five technical seeds within each sample, then macro-averaged equally over
tissue or approach. Its uncertainty is a 2,000-replicate study-block bootstrap.
Tissue-homogeneous studies are resampled within tissue strata. Because three
studies contain both approaches, the approach analysis resamples all study
blocks globally and recomputes the approach macro-average; mixed studies are
never split or treated as two independent blocks. Non-estimable support is
reported rather than silently replaced.

All endpoints are secondary or supportive and belong to a complete-reporting
family with no p-values and no multiplicity-driven selection. The evaluation
cannot promote a representation, distance, algorithm, component, fold, tissue,
or approach into a winner.

## 5. Immutable execution queue

The queue is the Cartesian product of:

- 150 frozen MV5-Q analysis groups;
- 2 algorithms (`pam_stability_k_v1`, `hclust_average_v1`); and
- 8 endpoints.

This produces exactly 2,400 unique evaluation units: 1,200 per algorithm and
300 per endpoint. Each unit identity includes the analysis group, algorithm,
endpoint, and complete source-freeze identity. Every queue row states
`outcomes_computed=FALSE`, `evaluation_executed=FALSE`, and
`method_selection_executed=FALSE`.

## 6. Validation and abort contract

The prospective validation plan requires external-label hash/schema checks,
all MV5-Q source hashes, complete group/seed/algorithm axes, exact 90-sample
join, independent synthetic ARI/NMI oracles, plurality-tie and training-only
mapping oracles, deterministic study-block bootstrap, immutable resume, and
public-safety checks.

Ten abort rules stop execution for source drift, incomplete axes, join or
missingness errors, upstream refit/reselection, cross-fold cluster pooling,
held-out-label leakage, seed pseudoreplication, stale partial artifacts,
failed independent/repeat/resume validation, or public leakage. The frozen
envelope is one worker, 300 seconds and 4 GiB per unit, two aggregate worker-
hours, 2,000 bootstrap replicates at seed `20260810`, and 1 GiB of public
outputs.

## 7. Evidence produced

- `mv05r-source-freeze-2026-08-10.csv`: 18 immutable source identities.
- `mv05r-external-label-design-audit-2026-08-10.csv`: structural join audit.
- `mv05r-algorithm-registry-2026-08-10.csv`: two fixed algorithm roles.
- `mv05r-endpoint-registry-2026-08-10.csv`: eight frozen endpoints.
- `mv05r-evaluation-queue-2026-08-10.csv`: 2,400 prospective units.
- `mv05r-validation-plan-2026-08-10.csv`: 12 required validations.
- `mv05r-abort-rules-2026-08-10.csv`: 10 hard-stop conditions.
- `mv05r-resource-envelope-2026-08-10.csv`: execution and storage caps.

The helper module enforces a pre-outcome firewall against outcome-result
columns or any executed/computed flag. Final acceptance evidence is:

- focused MV5-R tests: 18 passed, zero failed/errors/warnings/skips;
- full repository suite: 803 passed, zero failed/errors/warnings, two expected
  CRAN-only skips;
- MV5-R module and builder parse: 2/2;
- frozen sources: 18/18 present with exact SHA-256 identities;
- clean regeneration: 8/8 public CSV artifacts byte-identical;
- public label-value scan: none of the five tissue values present;
- exact staged-tree source tarball SHA-256:
  `cebf70075e0fe654a78ec3654565f39a218098a72b04baf455a669311904214c`;
- standard `R CMD check --no-manual`: `Status: OK`, log SHA-256
  `1c7485cdcc4f849a914b64fd85ebf5d82b22d61ddb5e4047fc9fd764d8c0f852`;
- stricter `--as-cran` check: code, install, examples, and tests pass; final
  status is one warning and two notes from the pre-existing non-mainstream
  `SeuratDisk` dependency, 41 imports, and unverifiable system time. Its log
  SHA-256 is
  `d2106865c5ad377593a71ac0b28aea35080ac6e723a40c2a18cd056d0fb2d91d`.

## 8. Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | approve |
| Correctness demonstrated? | pass for prospective construction; execution validation remains future work |
| Computation feasible? | yes, within frozen envelope |
| Biological interpretation permitted? | prohibited in MV5-R |
| Next action | commit this prefreeze, then separately execute MV5-S without changing it |

MV5-R stops before real evaluation. The next sprint may execute only the
frozen 2,400-unit contract after revalidating every source identity. Robustness
expansion, spectral promotion, gene topology, cell/gene fusion, new data,
optimization/Rust, package-default changes, and manuscript claims remain
closed.
