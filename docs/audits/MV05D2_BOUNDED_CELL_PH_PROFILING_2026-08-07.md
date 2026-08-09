# MV5-D2 bounded cell-PH profiling audit

| Field | Result |
|---|---|
| Date | 2026-08-07 |
| Scope | 30 representative SCT cell views |
| Coverage | 15 folds, 5 seeds, 15 held-out, 15 training |
| Mapping strata | 9 mapped held-out, 6 unmapped held-out |
| Outcome-label state | Closed |
| Result | MV5-D2 complete; stop before full PH and landscapes |

## Outcome

MV5-D2 closed the remaining cell-PH feasibility gap without crossing into
production or biological analysis. All 30 selected 384-cell by 30-PC views
produced complete typed H0/H1 diagrams. The pilot launched no landscape,
distance, clustering, integration, gene-view, or outcome job.

## Representative selection

The public manifest was built from the identities of the 75 independently
validated MV5-D1 caches. A capacity-constrained deterministic assignment places
exactly three held-out jobs in each seed while requiring a mapped view from
every fold where training-schema mapping occurs. A separately balanced
training-control schedule contributes three more jobs per seed.

| Stratum | Jobs |
|---|---:|
| Held-out, training-schema mapped | 9 |
| Held-out, no missing training features | 6 |
| Training control | 15 |

No tissue, assay-approach, or outcome field entered the selector.

## Diagram profile and correctness

Every diagram contains exactly 384 H0 intervals: 383 finite component merges
and one explicit essential interval with infinite death. The pilot retained
10,322 finite H1 intervals, ranging from 171 to 553 per view.

| Validation | Result |
|---|---:|
| Full-view H0 versus Euclidean MST | 30/30 pass |
| Maximum recorded H0 absolute error | 0 |
| Independent exact repeats | 5/5 pass |
| Byte-identical repeat files | 5/5 |
| Reduced Ripserr/GUDHI comparisons | 10/10 pass |
| Maximum recorded cross-engine absolute error | 0 |

GUDHI reports the essential H0 class at the supplied maximum filtration scale,
whereas Ripserr omits it. The cross-engine validator explicitly removed that
one capped GUDHI record before comparing the 31 finite H0 merges for each
32-cell subset. This normalization was recorded rather than hidden.

## Resources

| Measurement | Result |
|---|---:|
| Completed / failed jobs | 30 / 0 |
| Worker-seconds | 110.554 |
| Median elapsed per job | 3.601 s |
| Maximum elapsed per job | 5.433 s |
| Maximum process-tree RSS | 239,341,568 B (228.3 MiB) |
| Private pilot diagram storage | 907,330 B |
| Per-job time / RSS guards | 600 s / 4 GiB |
| PH implementation SHA-256 | `711d3c0e896bb09ddba970c08296e6d7901a05d730405dac9d3f7442814d5dd3` |

Measured time includes process startup, source-cache loading and validation,
typed PH, the independent MST oracle, and atomic serialization. It is therefore
a conservative per-view basis for the current one-view-per-process production
design.

## Full cell-PH projection

| Assumption | PH worker-hours | PH storage |
|---|---:|---:|
| Observed median | 6.752 | 205,196,625 B |
| Observed P90 | 7.135 | 230,971,500 B |
| Observed maximum | 10.187 | 253,509,750 B |

After adding measured MV5-D0 normalization, measured MV5-D1 coordinates, and
the existing exact-landscape projection, total SCT cell-primary estimates are
15.262, 15.645, and 18.697 worker-hours. Even the observed-maximum scenario is
2.903 hours below the 21.6-hour planning cap. These are feasibility projections,
not completed full-run measurements.

## Landscape definition remains unchanged

This sprint stops at diagrams. For the later landscape stage, H0 and H1 remain
separate; the essential infinite H0 interval is excluded; all finite positive-
persistence intervals contribute through every active consecutive landscape
level; and L2 integration is exact or error-controlled on dimension-specific
support. No universal level cap and no fixed uniform grid were introduced.

## Privacy and scope

- Private diagrams and execution logs are under ignored `tmp/mv05d2` paths.
- PDFs, reviewer correspondence, and `example_run.r` remain untracked.
- No dependency or lockfile was changed.
- Nothing was pushed.
- Full production PH jobs launched: 0.
- Biological outcomes computed: 0.
- The full source-loaded suite passed 550/550 expectations.
- The isolated staged source-package check reported `Status: OK` under R 4.4.1.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve for SCT cell H0/H1 diagrams |
| Correctness demonstrated? | Pass: MST, exact repeat, and cross-engine checks |
| Computation feasible? | Yes under all three bounded projection scenarios |
| Biological interpretation permitted? | Prohibited |
| Next action | Separately authorize immutable full 6,750-view cell-PH caching; stop again before landscapes |
