# MV5-X PC20 configuration execution prefreeze

Date: 2026-08-10

Status: execution contract frozen; production not yet started

Accepted predecessor: MV5-W evidence commit `eb8540e`

Execution-engine commit: `1b8e25782706ae82161db91d7914fd01b8f0d39c`

## Decision

MV5-X authorizes exactly one bounded, label-closed robustness configuration:
`cells384_pc20_euclidean_v1`. It contains 150 groups: 75
`inductive_integrated` and 75 `sct_whole` groups spanning 15 folds and five
seeds. The other 450 MV5-V groups remain unauthorized.

This authorization permits calculation and validation of private distance
artifacts only. It does not permit opening labels or computing, inspecting,
ranking, comparing, or reporting retrieval, clustering, robustness, or
biological outcomes.

## Prospectively bound identities

| Identity | Frozen value |
|---|---|
| Prospective engine HEAD | `1b8e25782706ae82161db91d7914fd01b8f0d39c` |
| Implementation SHA-256 | `09f24619f513512537a4b7aac4f8c91565bef176822e84fd2da493dec6704b98` |
| Exact-landscape script SHA-256 | `1563b5576cff06e1beab2006917ffc4ab9120188da01c2df8105d6c1828c8f60` |
| Queue SHA-256 | `57df129cf519b60c78d3bb1a602ab3e97670f5ce8c8447e6f44f6da77a389c56` |
| Scope SHA-256 | `15fe24d6f096bc97f683b7a849246ccbd362c55bc249f77de6e992493fbd1663` |
| Source-freeze SHA-256 | `4f3322c0c8f583d6df5e164bfd2648ae21d6ba55eb6e70a8816b84cc907064a8` |
| Coordinate sources | 150 private, content-addressed records |
| Workers | exactly 1 |

The source-freeze ledger also binds the Python executable and its reported
Python, Persim, NumPy, and SciPy versions. Private locators are represented by
opaque provenance locators; private coordinate payloads are not tracked.
The monitor additionally requires the exact ledger-bearing execution HEAD at
launch and records that HEAD in every completed group status; this later commit
does not alter any file in the implementation hash.

## Frozen calculation scope

| Quantity | Authorized total |
|---|---:|
| Groups | 150 |
| Views | 13,500 |
| Heldout-to-training biological pairs | 70,700 |
| Exact H0/H1 landscape rows | 141,400 |
| Landscape subchunks (at most 250 rows each) | 720 |
| Energy-distance rows | 70,700 |
| Four-method assembly rows | 282,800 |

Each group retains the accepted MV5-W semantics: 90 complete sample views,
H0 and H1 computed separately, exact all-active-level landscape integration,
no fixed landscape-level cap, deterministic subchunks of at most 250 requests,
matched energy distance, and four label-closed method rows per biological
pair. No dense landscape matrix is authorized.

## Resource and publication controls

- Per group: 600 seconds and 4 GiB process-tree RSS.
- Whole configuration: 30 elapsed worker-hours and 16 GiB private storage.
- Atomic publication: a group becomes reusable only after its status,
  manifest, hashes, and nine declared artifacts validate.
- Checkpointing: the one-worker monitor publishes an atomic resource checkpoint
  after every group and accepts only a valid completed prefix on restart.
- Fail closed on stale Git/runtime/source/queue identity, authorization leakage,
  partial or hash-invalid output, cap breach, or any label/outcome exposure.

## Prefreeze validation

- The focused MV5-X contract suite passed 12 expectations.
- The complete repository suite passed after the engine change.
- Independent ledger inspection found 600 queue rows, exactly 150 authorized
  PC20 rows, 75 rows per representation, 450 unauthorized rows, 150 private
  coordinate identities, and zero label/outcome authorization.
- Production results, resource rows, repeats, resumes, and outcomes are all zero
  at this prefreeze boundary.

## Next action

Commit these immutable binding ledgers, then run the bound queue with the
dedicated one-worker MV5-X monitor. Validation and repeat/resume checks remain
label closed. No other configuration may begin during this sprint.
