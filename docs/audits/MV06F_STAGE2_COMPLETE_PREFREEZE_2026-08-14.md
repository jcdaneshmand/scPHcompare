# MV6-F complete production and immutable-resume prefreeze

| Field | Frozen value |
|---|---|
| Date | 2026-08-14 |
| Serial policy commit | `d280453` |
| Queue / implementation / Rust | `f5471633…10bb5` / `5a1258e8…8d292` / `51d3fca4…160d` |
| Complete groups | 75/75 |
| Biological pairs / component rows | 35,350 / 141,400 |
| Component balance | 35,350 each for cell-H0, cell-H1, gene-H0, gene-H1 |
| Actual stage-two work | 21,538.531 s (5.983 h) |
| Peak process-tree RSS | 9,575,215,104 B |
| Groups above former 8-GiB cap | 11 |
| Retained private state | 624,237,551 B |
| Canonical monitor | 74/74 validated reuse; zero failures |
| Independent validation | 6/6 categories pass |
| Resume-checker SHA-256 | `eb4beb16095707db107594bf9dc48b9a87ceaf2da1a832924779306642581aa5` |
| Scientific firewall | Labels closed; fusion, clustering, outcomes, and claims remain zero/closed |

## Production result

The separately committed serial 12-GiB driver completed all 74 authorized
stage-two rows. It reused the admitted exact-group diagnostic and computed the
remaining 73 groups without retry. Every new group published a complete atomic
directory and a successful per-group resource metric. No group exceeded 12
GiB or 1,800 seconds; 11 valid groups did exceed the superseded 8-GiB ceiling,
confirming that the resource amendment prevented false failures without
changing H1 or the scientific workload.

The original frozen monitor then invoked every completed group through its
resume path and produced 74/74 `reused_validated` canonical rows. Its metrics
file SHA-256 is
`b45ad98cf413b694721729b32822156a8a8f67593cf19228e6b0ff5ba56b4101`.
The independent validator separately checked all 75 group directories,
reconstructed every pair and component total, required globally unique pair
identities and balanced cell/gene H0/H1 rows, and passed 6/6 categories. The
complete inventory SHA-256 is
`43c4af75cf00ecf8b38309b27415bace9f762fbbc8573ca1961a8a9302488b41`.

## Prefrozen immutable-resume gate

Before any fusion prefreeze, the committed checker must snapshot SHA-256,
byte size, and modification time for the five scientific artifacts in each of
75 group directories (375 files) plus the canonical resource metrics. It then
reruns the original frozen monitor against the complete corpus and requires all
three identities to remain unchanged for every file. The run must contain no
partial directory and must exit zero. Evidence is written only after the
before/after comparison.

This checker is a validation utility outside the already fixed scientific
implementation root. It cannot change the queue, diagrams, H0/H1 landscape
definition, Rust binary, labels, or downstream scope. Fusion, clustering,
outcomes, defaults, release, and claims remain closed until the gate passes.

## Resume result

The committed checker reran the original monitor successfully and compared 376
rows: 375 group artifacts plus the canonical resource metrics. SHA-256, byte
size, and modification time were unchanged for 376/376 rows. No group rebuilt,
no partial directory appeared, and labels/downstream counters remained closed
and zero. The resume evidence SHA-256 is
`2742cca6c6abdb0ef41dcea255824d807e2eb630b59344b469ea4eef91d701bb`.

MV6-F is therefore complete. Only a separately committed, label-closed MV6-G
fusion prefreeze may proceed next; outcome evaluation remains closed.
