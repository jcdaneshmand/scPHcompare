# Provenance instrumentation sprint

Date: 2026-08-04

Scope: sample flow, persistent-homology attempt provenance, reconciliation, and small-input boundary handling

Scientific behavior: existing QC thresholds, `MIN_CELLS = 250`, PH dimensions, PH threshold defaults, normalization, and integration settings are unchanged

## Outcome

The pipeline now creates an auditable link from every loaded sample to its post-QC disposition and from every PH-eligible sample to its individual PH execution attempt. Numeric batch IDs remain for compatibility, but they are no longer the only identity retained for PH work.

## Implementation

### Sample flow

`R/provenance_utils.R` defines dependency-light constructors, updates, CSV writing, set reconciliation, retry summarization, and nearest-neighbor boundary helpers.

`process_datasets_PH()` now:

1. initializes one sample-flow row per loaded input matrix;
2. records loaded dimensions and a source-reported count when a recognized metadata field exists;
3. records pre-QC and post-QC cell counts;
4. assigns `eligible_for_ph` or `excluded_before_ph` using the unchanged `post_qc_cells >= MIN_CELLS` rule;
5. writes `provenance/sample_flow_<dataset_tag>.csv` before PH;
6. asserts that input IDs partition exactly into excluded and eligible IDs; and
7. returns the flow, paths, eligible IDs, and excluded IDs in the additive `ph_results$provenance` element.

The former `discard()` predicate was replaced by an equivalent logical keep mask so that the excluded and retained identities are available before subsetting. Sample order and the threshold comparison are preserved.

### PH attempts

The PH dispatch chain now carries cohort, representation, and stable sample ID through:

- `compute_ph_batch()`;
- `process_expression_list_with_monitoring()`;
- `process_batch_of_datasets()`; and
- `process_and_monitor()`.

Each worker returns a one-row attempt record to the parent. The parent appends a complete batch serially, avoiding concurrent writes from parallel workers. The attempt CSV contains cohort, representation, sample ID, legacy job ID, attempt number, threshold, exact input dimensions, timestamps, elapsed time, subprocess exit status, timeout state, PD-file state, controlled status, and bounded diagnostics.

Controlled statuses are `completed`, `ph_timeout`, `ph_child_error`, `ph_empty_output`, and `parent_error`. Retry history is retained rather than overwritten. A helper selects final status per cohort/representation/sample without discarding the earlier attempts.

The PH child now exits nonzero after a caught `ripserr` error. Parent code reads child stdout/stderr after termination and includes bounded diagnostics in failed attempt provenance. A stale PD path is removed before launch so an old artifact cannot make a new failed attempt appear successful.

### Reconciliation

After each representation finishes, result names are normalized to the expression-list sample IDs. Eligible samples are partitioned into completed and failed IDs and checked exactly. `strict_reconciliation = TRUE` is the default and stops the run if any eligible sample lacks a final nonempty PD. Diagnostic runs may set it to `FALSE`; failures remain explicit in the attempt log and returned batch metadata.

### Small-input boundary

The historical timeout estimator requested `k_nn = 100` unconditionally. It now uses `min(100, n_observations - 1)` and raises a clear error for fewer than two observations. This prevents the estimator itself from requesting an impossible neighbor count while preserving `k = 100` for inputs with at least 101 observations.

### Real subprocess fixture

`tests/testthat/test-ph-subprocess-provenance.R` exercises the actual production subprocess path with no dissertation data:

- a four-point numeric matrix completes through `ripserr::vietoris_rips()`, returns a nonempty PD, exits zero, and records `completed`;
- a same-sized deliberately nonnumeric matrix reaches the child, fails inside the PH call, exits nonzero, writes no PD, and records `ph_child_error` with diagnostic text;
- the two returned attempts are written to and read back from a real attempt CSV;
- a two-row sample-flow CSV is written and reconciled as one completed plus one failed PH-eligible sample; and
- all fixture artifacts live under a test-only temporary directory and are removed after the test.

To keep this test practical, `process_and_monitor()` now accepts `poll_interval`; its production default remains 60 seconds, while the fixture uses 0.05 seconds. The fixture also exposed and removed a legacy first-write warning in `update_progress_log()` without changing its CSV schema.

## Files changed

- `R/provenance_utils.R` — new internal provenance/reconciliation helpers
- `R/PH_Calculation.R` — sample-flow capture and additive return metadata
- `R/ph_utils.R` — PH completion reconciliation and provenance routing
- `R/PH_Functions.R` — stable sample identity, structured attempts, subprocess diagnostics, bounded `k`
- `tests/testthat/test-provenance.R` — focused deterministic tests
- `tests/testthat/test-ph-subprocess-provenance.R` — real success/failure PH subprocess and CSV fixture
- `man/process_datasets_PH.Rd` and `man/process_and_monitor.Rd` — updated public documentation

## Verification

Performed in the elevated Ubuntu WSL environment using R 4.4.1.

- All modified R and test files parsed successfully.
- Focused provenance suite: 24 assertions passed.
- Real PH subprocess/CSV fixture: 19 assertions, including known H0/H1 feature counts.
- Package loaded successfully through `devtools::test()`.
- Installed `processx` exposes the `read_all_output_lines()` and `read_all_error_lines()` APIs used by the implementation.
- `git diff --check` reported no whitespace errors before final documentation updates and should be rerun at handoff.

The five post-processing and unified-pipeline failures discovered by the complete suite were repaired in a bounded follow-up; see [`TEST_SUITE_REPAIR_2026-08-04.md`](TEST_SUITE_REPAIR_2026-08-04.md).

## Remaining validation

No dissertation-scale dataset was rerun. The small success/failure subprocess contract is now exercised directly; large-cohort reproduction remains a separate, resource-planned sprint.
