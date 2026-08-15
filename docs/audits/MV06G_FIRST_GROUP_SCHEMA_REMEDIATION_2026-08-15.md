# MV6-G first-group schema remediation

## Captured cause

The logging-corrected isolated rerun reproduced the first-group stop in
210.019 seconds at 167,661,568 bytes peak process-tree RSS and captured the
child error. The ranking builder received a missing expected-pair value because
the completion queue exposed `query_biological_pairs`, while the already
admitted scientific runner reads the parent MV6-F name `biological_pairs`.
The resulting guard expression evaluated `NA` and failed closed.

This is an execution-queue schema mismatch. It is not a failure of the
persistence diagrams, all-active-level landscapes, H0/H1 handling, Rust
kernel, component scaling formula, or ranking formula. No scientific group
artifact was published, and no later group, label, outcome, or fusion endpoint
ran. The public failure metric and private partial/logs are preserved in their
attempt-two evidence and ignored quarantine locations.

## Prospective correction

The completion prefreeze now adds `biological_pairs` as an exact integer alias
of `query_biological_pairs` before writing the 74-row execution queue. The
policy validator, independent validator, and live monitor all require byte-for-
byte integer equality of the two columns. A focused regression rejects a
missing or conflicting alias.

The scientific runner and scientific implementation root remain unchanged.
The execution root, queue hash, policy hash, and repeat ledger must be rebuilt
and pass all original gates before one clean first-group rerun. Later groups
remain closed until that corrected group completes and independently validates.
