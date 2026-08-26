# MV17-B null-qualification closure v1 rejection

The prospectively frozen primary and independent repeat each completed all 36
synthetic jobs. Their metric files are byte-identical at SHA-256
`de176f069bb7efd1a8ccaa3e10ecef4d5bfabcfd7e419197b6d2e49edbe3287f`,
and their status files are byte-identical at SHA-256
`bdcbb7e2a82798f69e43bb0d12ce7bfd2da401321993261927d82eeefab6aa26`.

Closure v1 passed its seven calculations but is rejected before commit because
its three public files neither bind those four production artifacts nor publish
the complete synthetic metric table. This is an audit-provenance deficiency,
not a scientific or execution failure. Production must not be rerun.

Recovery is closure-only: bind all primary/repeat bytes and hashes, enforce the
complete 36-row family/fixture/seed grid, publish complete metrics and four
family summaries, retain all frozen thresholds, and require a new exact-head
recovery prefreeze before producing closure v2.
