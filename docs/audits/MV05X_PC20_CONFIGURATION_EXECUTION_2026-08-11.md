# MV5-X PC20 configuration execution

Date: 2026-08-11

Status: complete; calculation accepted, outcome evaluation still unauthorized

Accepted predecessor: MV5-W evidence commit `eb8540e`

Execution engine: `1b8e25782706ae82161db91d7914fd01b8f0d39c`

Binding commit: `e2bf8356e71904b4d80496ca06eba4246e216498`

Execution HEAD: `bbdeac17370b20c9a479f18a0cc7e3171e2d186b`

## Decision

The complete `cells384_pc20_euclidean_v1` calculation is accepted as an
auditable, label-closed robustness input. Exactly 150 prospectively authorized
groups completed: 75 `inductive_integrated` and 75 `sct_whole` groups over 15
folds and five seeds. All 450 rows belonging to the other three MV5-V
configurations remained unauthorized and unexecuted.

This decision accepts calculation integrity and configuration-scale
feasibility only. It does not authorize retrieval, clustering, robustness, or
biological outcome evaluation; it does not rank or compare methods; and it
does not authorize another robustness configuration.

## Completed calculation scope

| Quantity | Validated total |
|---|---:|
| Groups | 150 |
| Private files | 1,650 |
| Manifest-declared artifacts | 1,350 |
| Sample views with PH/MST validation | 13,500 |
| Landscape summaries | 27,000 |
| Heldout-to-training biological pairs | 70,700 |
| Exact H0/H1 landscape rows | 141,400 |
| Energy-distance rows | 70,700 |
| Four-method rows | 282,800 |

Every landscape row retains the revised dissertation-aligned contract: H0 and
H1 are separate, all consecutive active levels are included, no fixed level
cap is applied, and integration uses exact linear critical-pair segments.
Landscape requests remain in deterministic subchunks of at most 250 rows.

## Execution evidence

- All 150 groups completed on the first pass; there were no failed, skipped,
  timeout, RSS, storage, or manual-retry dispositions.
- Total measured group time was 11,163.624 seconds (3.101 worker-hours).
- Maximum group time was 185.348 seconds versus the 600-second cap.
- Maximum process-tree RSS was 638,365,696 bytes versus the 4-GiB cap.
- Final private storage was 2,487,457,825 bytes versus the 16-GiB cap.
- Exactly one worker was used throughout.
- The result payloads remain under ignored `tmp/mv05x/` storage; only
  label-closed provenance and aggregate evidence are published.

## Independent validation

The independent validator passed 15 categories. It rehashed every declared
artifact and independently checked:

- the 150/450 authorization boundary and bound execution identity;
- all private source hashes, group statuses, manifests, and file sizes;
- 13,500 view records, including point/coordinate counts, H0 interval
  cardinality, essential-H0 policy, finite-interval reconciliation, and every
  stored H0 MST oracle result;
- complete heldout-to-training pair coverage for 70,700 biological pairs;
- exact/all-active H0/H1 flags, squared-distance consistency, subchunk limits,
  and all landscape/energy cardinalities;
- reconstruction of all 282,800 four-method distances from their H0, H1, and
  energy inputs within a strict numerical tolerance;
- aggregate resource caps and label/outcome closure.

The validation step did not compute ranks, winners, retrieval correctness,
cluster agreement, p-values, or biological summaries.

## Determinism and resume

Two groups were marked prospectively for clean repeat: one per representation.
All eight deterministic artifacts per group reproduced byte-for-byte, for
16/16 exact artifact matches.

The full completed tree was hashed before and after a 150-group
`validate_resume` run. Every group returned `reused_validated`, and all 1,650
relative paths, byte counts, and SHA-256 values were unchanged exactly.

The original single-smoke snapshot helper correctly rejected the 150-group
tree. A dedicated MV5-X helper now enforces exactly 150 published directories
and 1,650 files; it is validation tooling and is not part of the frozen
calculation implementation hash.

## Verification

- Focused MV5-X contract tests: 12 expectations passed.
- Complete repository suite: passed with the two established CRAN skips.
- Independent production validator: passed for 150 groups and 70,700 pairs,
  with labels opened = 0 and outcomes computed = 0.

## Next action

The next sprint may only prefreeze the PC20 robustness-outcome evaluation. It
must define estimands, aggregation across seeds/folds, method comparisons,
uncertainty, missingness, multiplicity/reporting, label access, and independent
validation before opening labels or computing an outcome. The other three
configurations, gene topology, fusion, new data, optimization/Rust, package
defaults, and manuscript claims remain outside that authorization.
