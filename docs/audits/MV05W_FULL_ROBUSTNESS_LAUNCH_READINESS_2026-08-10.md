# MV5-W full-robustness launch readiness

Date: 2026-08-10

Status: complete

Accepted MV5-V base: `64b6de0`

Prospective engine commit: `383dfd8`

Bound smoke queue commit: `5594a22`

Engine SHA-256: `a661e4470c4d4b56d8bb6862dbb373822f7b2738fb6416a01b2b2f5e73f1e8c5`

Labels opened: no

Outcomes computed: no

## Decision

MV5-W passes launch readiness for exactly the first MV5-V configuration:
384 cells, the first 20 accepted PCs, and Euclidean geometry. A later sprint
may prospectively bind and execute its complete 150 groups (75 per
representation). The other three configurations, labels, retrieval evaluation,
robustness comparisons, rankings, and claims remain unauthorized.

## Engine implemented

The dedicated runner reuses the accepted coordinate transforms, `ripserr` H0/H1
PH, independent H0 MST oracle, Persim exact critical-pair landscapes, and
matched energy implementation. It adds:

- complete accepted heldout-query-to-training pair construction;
- deterministic H0/H1 subchunks capped at 250 rows;
- exact all-active/no-cap landscape validation;
- four label-closed method rows per biological pair (H0, H1, descriptive raw
  H0/H1 Euclidean combination, and energy);
- same-filesystem atomic group publication;
- hash-valid reuse and hard failure for stale/partial publication;
- a one-worker monitor with 600-second/4-GiB group guards.

## Prospectively selected real smoke

Only queue row 1 was authorized before execution:

- fold `large_loso_v1:SRA550660`;
- seed `20260805`;
- representation `inductive_integrated`;
- configuration `cells384_pc20_euclidean_v1`;
- 90 sample views;
- 425 directed biological pairs;
- 850 H0/H1 exact-landscape rows in four subchunks;
- 425 energy rows;
- 1,700 assembled method rows.

It completed in 70.085 seconds with 541,970,432 bytes peak process-tree RSS
and 15,559,362 bytes of private output. All observed values are far below the
600-second, 4-GiB RSS, and 4-GiB smoke-storage caps. Labels and outcomes are
zero.

## Independent validation

All 13 categories pass: bound identity, exact group/view/pair axes,
subchunk coverage, all-view H0 MST, analytic exact landscape, independent
energy, four-method reconstruction, artifact hashes, clean repeat, immutable
resume, resources, and label safety.

The clean repeat matches all eight deterministic scientific artifacts
byte-for-byte. A validation-only resume reports `reused_validated`; all 11
private paths, hashes, byte sizes, and modification timestamps are unchanged.

The validator initially used bitwise `identical()` on CSV-round-tripped method
distances. Column types/order and values agreed under numerical comparison, so
the validator was corrected to require exact identity for nonnumeric fields and
absolute distance error at most `1e-12`. Production/repeat artifacts were not
changed. A separate count correction records 150 groups for a full
two-representation configuration rather than 75 per representation.

## Public evidence

- `mv05w-source-freeze-2026-08-10.csv`
- `mv05w-smoke-queue-2026-08-10.csv`
- `mv05w-prefreeze-summary-2026-08-10.csv`
- `mv05w-smoke-resources-2026-08-10.csv`
- `mv05w-independent-validation-2026-08-10.csv`
- `mv05w-deterministic-repeat-2026-08-10.csv`
- `mv05w-resume-validation-2026-08-10.csv`
- `mv05w-launch-decision-2026-08-10.csv`

All private smoke/repeat/resume outputs, coordinate caches, labels, PDFs,
reviewer material, and `example_run.r` remain untracked.

## Next action

Prospectively bind MV5-X to the 150 PC20 groups only, run one worker, preserve
configuration-level stop decisions, and independently validate all groups,
complete pair/subchunk axes, resources, selected repeats, and immutable resume
before any labels or robustness outcomes are opened. Do not combine this with
the other three configurations.
