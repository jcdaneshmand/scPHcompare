# MV5-D4 exact cell landscape-distance audit

| Field | Result |
|---|---|
| Date | 2026-08-07 |
| Input | 6,750 validated MV5-D3 H0/H1 records |
| Scope | 35,350 query-to-training biological pairs |
| Dimension rows | 70,700: equal H0 and H1 |
| Chunks / groups | 360 / 75 |
| Outcome-label state | Closed |
| Result | MV5-D4 complete; stop before endpoints and outcomes |

## Outcome

The complete corrected SCT cell landscape-distance layer now exists under the
owner-approved dissertation-aligned definition. Every eligible held-out query
was compared only with training references from its fold and seed. All outputs
were created and independently validated without tissue, approach, or outcome
labels.

## Scope and interval boundary

| Measurement | Result |
|---|---:|
| Biological query-to-training pairs | 35,350 |
| Separate H0 rows | 35,350 |
| Separate H1 rows | 35,350 |
| Finite positive-persistence intervals staged | 4,682,085 |
| Finite H0 intervals | 2,585,250 |
| Finite H1 intervals | 2,096,835 |
| Essential H0 intervals excluded | 6,750 |

The full request manifest is retained as a 13 MB gzip rather than an 88.7 MB
plain CSV. Its decompressed SHA-256 is
`6b4bb2967669db6e510db7e759157d4c1d05cb3e9d2f4a00030b462678d897b0`,
which is the identity bound into every accepted request and chunk.

## Admission and exact oracle

The first admission group contained 850 rows in four chunks; the maximum-scope
group contained 3,250 rows in 14 chunks. Both passed complete structural and
resource validation. One full H0 and H1 pair from each group was recomputed by
the independent R exact breakpoint-stream oracle.

| Dimension checks | Maximum absolute difference | Result |
|---|---:|---|
| 2 H0 | `1.4210854715202e-14` | Pass |
| 2 H1 | `4.88498130835069e-15` | Pass |

No reduced diagram, fixed grid, or capped landscape depth entered these checks.

## Accepted production evidence

| Measurement | Result |
|---|---:|
| Completed / failed groups | 75 / 0 |
| Completed / failed chunks | 360 / 0 |
| Completed / failed dimension rows | 70,700 / 0 |
| Group worker time | 4,194.152 s (1.165 h) |
| Exact pair-operation time | 2,613.887 s |
| Median / maximum group elapsed | 49.958 / 156.214 s |
| Peak process-tree RSS | 366,899,200 B (349.9 MiB) |
| Private output/status storage | 65,057,103 B (62.0 MiB) |
| Accepted implementation SHA-256 | `1c40f910ece22cff9a5243b7ddfb770739e4f6ec1a3711365ea3cde506973a3e` |

Every result is finite, non-negative, exact, all-active-level, and matched to
its source record and diagram hashes. Every stored output hash was recomputed.
H0 and H1 were paired into 35,350 public compressed component rows with their
secondary combined distance and H1 squared-distance fraction.

## Determinism and correction audit

The first complete pass embedded `pair_seconds` in scientific CSVs. That made
file bytes nondeterministic even though the distances were correct. It was
rejected, preserved privately, and replaced by v3, which stores timing only in
status/resource records.

| Check | Result |
|---|---:|
| Matched superseded/accepted request IDs | 70,700 / 70,700 |
| Scientific fields identical across versions | Yes |
| Accepted scientific files contain timing | No |
| Fresh complete-group repeat rows | 850 / 850 |
| Byte-identical repeated scientific files | 4 / 4 |
| Immutable admission resume files | 8 / 8 |

The monitor's first attempt also failed safely before producing distances
because path normalization resolved the virtual-environment Python symlink to
system Python, where Persim was absent. The two zero-row failures are retained
as public evidence. The corrected launcher preserves the venv executable path;
the subsequent v2 and accepted v3 queues completed fully.

## Measured primary resource decision

| Component | Worker-hours |
|---|---:|
| SCT normalization | 2.562 |
| Training-only cell coordinates | 2.376 |
| Complete cell PH | 1.047 |
| Complete exact cell landscape distances | 1.165 |
| Total measured cell-primary precomputation | 7.150 |
| Planning cap with reserve | 21.600 |
| Remaining margin | 14.450 |

The former 3.572-hour landscape projection has been replaced, not added. All
primary cell-topology precomputation terms are now measured.

## Privacy, checks, and scope

- Private interval bundles, chunks, statuses, repeats, superseded results, and
  logs remain under ignored `tmp/mv05d4` paths.
- PDFs, reviewer correspondence, and `example_run.r` remain untracked.
- No dependency or lockfile was changed and nothing was pushed.
- Clustering, integration, gene-view, and biological-outcome counters are zero.
- Outcome labels remained closed.
- The complete current-source test suite passed with only its two intentional
  CRAN skips.
- The isolated source-package check under R 4.4.1 reported `Status: OK` with
  zero errors, warnings, or notes. A first wrapper timed out after 20 minutes;
  a clean rerun from the same built tarball completed in 302.1 seconds.

## Decision table

| Question | Disposition |
|---|---|
| Dissertation-aligned landscape definition implemented? | Yes |
| H0/H1 components retained separately? | Yes; combination secondary only |
| Complete eligible pair scope? | Yes: 35,350 pairs / 70,700 rows |
| Correctness and determinism demonstrated? | Yes: R oracles, all-record hashes, exact repeat, resume |
| Resource gate passed? | Yes: 7.150 h measured total, 14.450 h margin |
| Biological interpretation permitted? | Prohibited |
| Next action | Label-closed assembly of training-fitted topological and matched-baseline retrieval bundles |
