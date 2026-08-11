# MV5-U bounded robustness resource admission

Date: 2026-08-10

Status: complete

Accepted MV5-T base: `47235d8`

Prospective engine commit: `a94e9d3`

Bound execution queue commit: `6102b27`

Validator-only correction commit: `8bc2718`

Production implementation SHA-256: `8b5110606f74308da5f78b57735fce1dd1751e5869c05a1d919308ef46c76366`

Labels opened: no

Outcomes computed: no

Full robustness authorized: no

## Decision

MV5-U passes its bounded feasibility, semantics, determinism, and immutable-
resume gates. The corrected dissertation-aligned landscape calculation is
feasible for the frozen admission strata with substantial margins under every
resource cap.

This is not a robustness result. MV5-U calculated no retrieval endpoint,
clustering result, label association, p-value, rank, winner, or biological
outcome. The passing admission therefore does not authorize an immediate full
grid. Its only positive continuation is a separate prospective streamed full-
robustness execution gate.

## Frozen scope actually executed

The execution queue contains exactly 24 units:

- three fixed LOSO folds representing 65, 86, and 89 training samples;
- two representations, `sct_whole` and `inductive_integrated`;
- four one-factor-at-a-time configurations at seed `20260805`:
  384 cells/20 PCs/Euclidean, 384 cells/30 PCs/cosine chord,
  192 cells/30 PCs/Euclidean, and 256 cells/30 PCs/Euclidean;
- 90 sample views per unit;
- H0 and H1 persistent homology for every view;
- 16 frozen training-training and 16 heldout-training pair requests per unit;
- exact landscape and matched energy work only for those admission pairs.

The 192- and 256-cell configurations use the frozen deterministic nested
point order. The 384-cell PC20 and cosine configurations preserve the accepted
row identity and order. Cosine geometry is Euclidean chord distance after
row-unit normalization.

## Landscape definition enforced

MV5-U uses the revised dissertation-aligned definition:

1. H0 and H1 are calculated and compared separately.
2. The essential H0 interval is excluded.
3. Every active consecutive landscape level is included; there is no universal
   level cap.
4. There is no fixed uniform integration grid.
5. L2 distances are integrated exactly over critical pairs of the piecewise-
   linear landscapes.
6. Computation is grouped and streamed; enormous dense landscape matrices are
   not retained.

The public artifacts record 4,320 landscape summaries and 1,536 landscape-
distance rows. Every distance is flagged exact and all-active, every level-cap
flag is false, and the largest observed discrepancy between `distance^2` and
the recorded squared distance is below `4e-12` against a `1e-10` tolerance.

## Execution evidence

| Measure | Observed | Contract limit | Result |
|---|---:|---:|---|
| Admission units | 24/24 | exactly 24 | pass |
| Sample views | 2,160 | exactly 2,160 | pass |
| Finite H0/H1 intervals | 1,069,189 | valid finite intervals | pass |
| Landscape summaries | 4,320 | exactly 4,320 | pass |
| Landscape pair rows | 1,536 | exactly 1,536 | pass |
| Energy pair rows | 768 | exactly 768 | pass |
| Maximum unit time | 55.508 s | 600 s | pass |
| Total measured unit time | 895.449 s | 7,200 s | pass |
| Peak process-tree RSS | 622,227,456 B | 4,294,967,296 B | pass |
| New private production bytes | 288,635,915 B | 2,147,483,648 B | pass |
| Workers | 1 | 1 | pass |
| Labels/outcomes | 0/0 | 0/0 | pass |

The ignored one-unit smoke run completed the complete 90-view path before the
production run. The public resource decision records
`admission_complete=TRUE` and `full_robustness_authorized=FALSE`.

## Independent numerical validation

All 15 independent validation categories pass. The validator:

- rehashes all 150 private coordinate sources and the eight frozen production
  implementations;
- independently reconstructs all 2,160 transformed views;
- checks nested 192-in-256 identity, first-20-coordinate identity, cosine unit
  norms, configuration isolation, and all sample axes;
- compares every view's finite H0 death times with an independently calculated
  Euclidean minimum-spanning-tree oracle;
- validates H1 on an analytic unit square with the expected interval
  `[1, sqrt(2)]`;
- runs the analytic exact-landscape oracle;
- independently recomputes sampled energy distances;
- checks artifact hashes/cardinalities, resource caps, and public label safety.

All 24 unit-level MST and sampled-energy flags pass.

## Determinism and resume

A clean 24-unit repeat was built in a separate private root. All seven frozen
deterministic scientific artifacts in every unit match byte-for-byte: 168/168
artifact comparisons pass.

The production root was separately snapshotted before and after a validation-
only resume. All 240 private paths, hashes, byte sizes, and modification times
are identical. The resume reports 24/24 `reused_validated` units and zero
rebuilds.

## Audited validator correction

The first independent-validation attempt stopped at the landscape semantics
gate. Diagnosis showed that Python's CSV writer emitted `True`/`False`, which
R retained as character values; the validator incorrectly compared those
strings directly with logical `TRUE`/`FALSE`. Unit cardinalities, flag content,
and landscape numerical residuals were already correct.

Commit `8bc2718` adds a strict case-insensitive boolean parser that accepts only
`true` or `false`, plus tests for logical, title-case, upper-case, invalid, and
missing values. The eight-file frozen production implementation digest remains
exactly `8b511060...a23842`. Neither production nor repeat artifacts were
modified. The corrected committed validator then passed all independent gates.

## Artifact index

- `mv05u-admission-resources-2026-08-10.csv`
- `mv05u-resource-summary-2026-08-10.csv`
- `mv05u-resource-decision-2026-08-10.csv`
- `mv05u-unit-completion-2026-08-10.csv`
- `mv05u-independent-validation-2026-08-10.csv`
- `mv05u-deterministic-repeat-2026-08-10.csv`
- `mv05u-resume-validation-2026-08-10.csv`
- `mv05u-validator-correction-2026-08-10.csv`

Private coordinate caches, production/repeat units, logs, snapshots, Python
runtime, labels, PDFs, reviewer material, and `example_run.r` remain untracked.

## Final verification

The focused corrected MV5-U test file passes 26 expectations. The complete
repository suite passes with only the two established CRAN skips. `R CMD build`
and `R CMD check --no-manual` on an export of the exact staged Git index finish
with `Status: OK`.

## Next action

Prospectively freeze MV5-V as a streamed full-robustness execution gate. It
must derive its exact group/pair/chunk scope and resource reserve from MV5-U,
bind all identities before execution, define resume and partial-output rules,
and keep labels and outcomes closed. MV5-V may authorize or reject a later
full robustness calculation; it may not execute that calculation in the same
prefreeze sprint.

Spectral promotion, gene topology, cell/gene fusion, new data, optimization or
Rust work, package-default changes, and manuscript claims remain outside the
next sprint.
