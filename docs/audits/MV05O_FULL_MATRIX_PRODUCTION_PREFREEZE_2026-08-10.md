# MV5-O complete-training-matrix production prefreeze audit

## Decision summary

| Question | Result |
|---|---|
| Exact MV5-N scope preserved? | **Yes**: 150 groups and accepted group/chunk identity roots |
| Landscape definition preserved? | **Yes**: separate exact all-active-level H0/H1, no cap/grid |
| Executable implementation frozen? | **Yes**: stager, group-cached landscape runner, baseline runner |
| Queue complete? | **Yes**: 4,565 units / 1,838,725 values |
| Resource envelope | **Pass**: 16.117 h and 1.278 GB versus 21.6 h / 10 GiB caps |
| Independent validation | **Pass**: 19/19 prefreeze checks |
| Deterministic repeat | **Pass**: 7/7 public artifacts byte-identical |
| Real runner fixture | **Pass**: 10/10 plus immutable 2-artifact resume |
| Canonical verification | **Pass**: 17/17 focused; 776/776 full; clean package `Status: OK` |
| Production or clustering executed? | **No**: both counters remain zero |
| Labels/outcomes opened? | **No** |

## 1. What was audited

The audit reconciled the MV5-N exact pair generator and group/chunk roots, the
bounded Persim critical-pair engine, matched energy/pseudobulk formulas, and
the immutable queue/status patterns used by MV5-D3/D4/D5 and MV5-G/H/I/J.

The bounded admission runner could not be promoted unchanged: its fixed
12×32-row boundary is correct for profiling but not full production, and
process-per-chunk execution would lose the caching assumed by the MV5-N
projection. MV5-O therefore adds a production group stager, a group-cached
landscape runner with chunk-level atomic outputs, and a generalized baseline
runner. Production has not been launched.

## 2. Frozen identity and queue results

The 18-artifact source root is
`541e7d3aa8acce5d512bbb4819c034735eef47387e91a63abccfa259f53d6de1`.
Four implementation hashes are recorded in the public source-freeze and
summary artifacts.

The exact queues contain:

- 150 fold-seed-representation groups;
- 4,340 landscape chunks and 1,050,700 H0/H1 rows;
- 150 energy groups and 525,350 pairs;
- 75 shared pseudobulk groups and 262,675 pairs; and
- 4,565 resumable units and 1,838,725 total values.

Every group/chunk ID is unique. All sources and current file hashes validate.
Non-materialized exact pair rows must be regenerated from the accepted
generator and reproduce their public roots before computation.

## 3. Resource and storage correction

The accepted worker projection is unchanged at 16.117047 hours including its
10% reserve. Execution permits at most two group workers, 900 seconds per
group/unit, 4 GiB per process tree, 21.6 aggregate worker-hours, and 10 GiB
additional private storage.

MV5-O explicitly adds a conservative baseline-storage estimate that MV5-N's
landscape output projection did not need to express:

| Component | Bytes |
|---|---:|
| Landscapes | 618,845,884 |
| Energy + shared pseudobulk at 512 bytes/value | 403,468,800 |
| 25% status/validation/storage reserve | 255,578,671 |
| **Total** | **1,277,893,355** |

## 4. Validation evidence

The independent validator reconstructed 19 categories: source existence and
hashes, one source root, exact group/unit/value counts, caps, unique IDs,
validation/abort registries, label closure, zero outcomes, zero clustering,
zero production, and summary arithmetic. All 19 pass.

Two clean prefreeze builds produced seven byte-identical artifacts. The frozen
validation plan names 12 exact R oracle requests, full output repeats for the
maximum group in each representation (33 units each), and immutable resume
over all 4,565 units.

The production landscape runner was also exercised on 32 corrected real
diagram pairs:

- 32/32 distances match the accepted MV5-N engine within `1e-12`;
- exact/all-level/no-cap flags all pass;
- one independent R exact oracle passes;
- clean output is byte-identical; and
- an in-place rerun reused the chunk while output/status hashes, sizes, and
  timestamps remained unchanged.

The production baseline runner retains the already validated MV5-N energy
V-statistic and pseudobulk Euclidean helpers, adds group/source/hash/status
guards, and parses without error. No real full-group baseline was executed.

## 5. Abort and recovery

Ten non-retrying abort rules cover every source/identity/label/resource/
partial-artifact/correctness boundary. Completed immutable artifacts are
preserved; partial state is quarantined. A correctness or cap failure revokes
authorization and requires a corrective sprint. It cannot silently narrow H1,
levels, samples, folds, seeds, or validation.

## 6. Decision and next action

MV5-O prefreeze passes technically. Focused MV5-O tests passed 17/17; the full
source suite passed 776/776 with zero failures, warnings, or skips; ten R
sources and the Python runner parsed; and a clean Git-archive source package
passed `R CMD check --no-manual` with `Status: OK`. The source tarball SHA-256
is `0235bf2987b24037d11c7c16c316d4329967471bf17bdb0efa9c0f96c1c19026`;
the check-log SHA-256 is
`c26efba96e247138357a8666e69c18da6a3b2dd94e46923a8b5f2022ead978ca`.
The local prefreeze commit may authorize the next goal to execute only the
frozen label-closed distance queues.

The next goal must stop again before fitting production clusters or opening
labels. Later descriptive training-partition alignment must remain separate
from inductive held-out generalization; neither may retune the pipeline.
