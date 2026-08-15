# MV6-G first-group failure and logging rebind

## Preserved failure

The committed serial completion stopped on execution order 2, before any later
group launched. The child exited 1 after 211.699 seconds at 169,197,568 bytes
peak process-tree RSS. It breached no frozen resource cap, published no group
directory, produced no scientific CSV, kept outcomes closed, and left only an
empty atomic partial directory. The failure metric is preserved in
`mv06g-completion-attempt1-evidence/mv06g-first-group-failure.csv`; the private
partial directory remains recoverable in the ignored attempt-one quarantine.

Read-only inspection of all 180 accepted source PH records rules out an empty
homology component: every sample has 383 finite cell-H0 intervals, 499 finite
gene-H0 intervals, 165–525 finite cell-H1 intervals, and 114–3,215 finite
gene-H1 intervals. The H0/H1 landscape contract is unchanged.

## Monitor defect and correction

The scientific child wrote to process pipes that the monitor never persisted.
Consequently the nonzero exit was fail-closed but not diagnosable. The
prospective correction redirects each child's stdout and stderr to unique
private files, binds both hashes into the resource metric, and includes only
the last stderr lines in a failed metric. Existing partials, metrics, or log
paths still block rerun until explicitly quarantined. No scientific-runner
source changes.

The corrected execution root must be rebuilt and pass the original prefreeze,
repeat, focused, and full-suite gates before exactly one unchanged first-group
diagnostic rerun. Remaining groups stay closed until the cause is captured and
resolved or the unchanged group completes and independently validates.
