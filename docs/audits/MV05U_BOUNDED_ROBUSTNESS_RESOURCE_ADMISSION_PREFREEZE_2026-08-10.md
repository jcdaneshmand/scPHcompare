# MV5-U bounded robustness resource-admission prefreeze

Date: 2026-08-10

Accepted MV5-T base: `47235d8`

Prospective engine commit: `a94e9d3`

Implementation SHA-256: `8b5110606f74308da5f78b57735fce1dd1751e5869c05a1d919308ef46c76366`

Admission units executed: zero

Outcomes computed: zero

## Decision

The MV5-U execution engine and its exact 24-unit queue are ready for bounded,
label-closed resource admission. No full robustness grid is authorized.

## Frozen scope

The queue is exactly three minimum/median/maximum folds by two representations
by four one-factor-at-a-time coordinate configurations at seed `20260805`.
Every unit must construct all 90 views and full H0/H1 diagrams, validate the H0
MST death multiset, stage exact all-active landscapes, and calculate exact H0
and H1 landscape distances plus matched energy on a frozen 32-pair coverage
set. No label, retrieval score, cluster, p-value, rank, or outcome is available
to the engine.

## Source and runtime binding

The source freeze contains 170 identities: 19 committed implementation,
specification, test, and accepted-gate artifacts; all 150 accepted private SCT
and integrated coordinate hashes; and one private Python executable hash. No
absolute private path or private payload is tracked.

The isolated ignored runtime is bound by executable SHA-256
`afee1865d29d02077b338580025d90aaeb3ed48189abe8c2456879baf659d69e`
and reports CPython 3.10.20, Persim 0.3.8, NumPy 2.2.6, and SciPy 1.15.3. Its
analytic single-square H1 landscape oracle passes.

## Validation before execution

- focused MV5-U contract tests: 23/23 passed;
- full repository tests: 916/916 passed with no failures, warnings, or skips;
- all seven R source/test files parse;
- Python source parses;
- analytic exact-landscape oracle passes;
- queue has 24 unique unopened units;
- all queue rows retain `labels_opened=FALSE`, `outcomes_computed=FALSE`, and
  `admission_executed=FALSE`;
- five public prefreeze CSVs contain no tracked label values.

## Resource and stop boundary

One worker is mandatory. Each unit is capped at 600 seconds and 4 GiB
process-tree RSS; total worker time is capped at two hours and new private
storage at 2 GiB. Atomic unit directories are immutable. Any source, semantic,
oracle, resource, partial-artifact, repeat, resume, or label-safety failure
aborts the queue while preserving completed units.

Even a passing admission records `full_robustness_authorized=FALSE`. The only
possible positive continuation is a later prospective streamed full-robustness
gate based on measured resources.
