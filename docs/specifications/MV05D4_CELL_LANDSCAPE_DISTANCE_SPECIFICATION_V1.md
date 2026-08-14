# MV5-D4 cell landscape-distance specification v1

| Field | Frozen value |
|---|---|
| Date | 2026-08-07 |
| Stage | MV5-D4 |
| Input | 6,750 independently validated MV5-D3 cell H0/H1 records |
| Biological pair scope | 35,350 held-out-query to training-reference pairs |
| Dimension rows | 70,700: 35,350 H0 and 35,350 H1 |
| Chunks | 360, at most 250 dimension rows each |
| Execution unit | One fold-seed group, containing dimension-specific chunks |
| Maximum heavy workers | 2 |
| Outcome-label state | Closed |
| Required stop | Before clustering, outcomes, integration, gene views, or fusion |

## Purpose

MV5-D4 converts the complete corrected cell persistence-diagram cache into the
exact sample-level topological distances required by the primary LOSO
retrieval design. It implements the project-owner-approved dissertation-aligned
landscape definition without opening labels or constructing a full fold matrix.

## Pair scope

Within every fold and seed, each held-out query sample is paired with every
training-reference sample and with no other sample. The full manifest binds
each H0 or H1 request to both MV5-D3 record keys, diagram SHA-256 values,
serialized-result SHA-256 values, fold, seed, roles, and the landscape
definition. It contains no tissue, approach, or outcome field.

This scope supports the frozen primary sample-retrieval endpoint. It does not
support full-matrix clustering or within-held-out-study pair contrasts; those
tasks remain deferred rather than silently approximated.

## Landscape definition

Every dimension row uses:

- separate H0 or H1 construction and reporting;
- only finite positive-persistence intervals;
- explicit exclusion of the one essential H0 interval per diagram;
- every consecutive active landscape level, with no universal level cap;
- Persim 0.3.8 exact critical-pair construction;
- project-controlled exact integration of every linear difference segment;
- no fixed uniform grid and an exact zero numerical-error estimate; and
- the R exact breakpoint-stream implementation as the independent oracle.

For each biological pair, the public component artifact retains H0 distance,
H1 distance, secondary `sqrt(H0^2 + H1^2)`, and the H1 fraction of squared
combined distance. Neither the combination nor any later normalization may
replace the primary component reporting.

## Immutable identity and staging

The full request manifest contains 70,700 unique identities and is retained in
gzip form for repository hygiene. Decompression reproduces the exact CSV
SHA-256 bound into every private request and result:
`6b4bb2967669db6e510db7e759157d4c1d05cb3e9d2f4a00030b462678d897b0`.

Each of 75 private group input bundles contains:

- its public request subset;
- finite positive-persistence intervals extracted from 90 validated MV5-D3
  records;
- source-result-set, request, interval, and manifest hashes; and
- a checkpoint proving that exactly 90 essential H0 intervals were excluded.

Existing partial, stale, mismatched, or implementation-incompatible inputs,
outputs, or status records may not be overwritten or mixed.

## Execution and guards

| Guard | Frozen value |
|---|---:|
| Requests per chunk | 250 maximum |
| Per-chunk elapsed cap | 900 seconds |
| Per-group monitored elapsed cap | 900 seconds |
| Per-group process-tree RSS cap | 4 GiB |
| Whole-stage worker-time cap | 14,400 seconds |
| Concurrent heavy groups | 2 maximum |

The monitor validates every chunk output hash, request count, implementation
identity, finite distance, exact/all-level flags, and downstream-zero counters
before recording a completed group. On a guard event it stops launching new
groups and allows active WSL processes to finish naturally.

## Correctness and determinism gates

Completion requires:

1. two admission groups spanning a representative and maximum pair count;
2. complete H0 and H1 comparisons from both admission groups against the
   independent exact R oracle;
3. 70,700/70,700 accepted rows and 360/360 chunks independently revalidated;
4. exact H0/H1 pairing into 35,350 secondary component records;
5. one fresh complete group repeat whose four scientific CSVs are byte-identical;
6. immutable resume with no output/status changes;
7. timing stored only in separate resource/status artifacts, never in immutable
   scientific distance CSVs; and
8. zero clustering, integration, gene-view, and outcome counters.

## Correction rule

An initial complete implementation embedded per-pair timing in scientific CSVs.
Its 70,700 scientific values were correct but its files could not be byte-
deterministic. Those outputs are superseded. The accepted v3 implementation
moves timing into status records; all scientific fields are identical across
all 70,700 matched request IDs.

## Completion gate

MV5-D4 completion authorizes a separate label-closed bundle-assembly sprint for
training-fitted component scaling and matched baselines. It does not authorize
opening outcome labels, clustering, integration, gene topology, fusion, or
biological claims.
