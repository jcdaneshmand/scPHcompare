# MV7-H landscape-stress validator portability rebind

Date: 2026-08-16

Status: accepted; clean independent-validation rerun authorized

## Second failed-closed validation attempt

The reference-engine-v3 validator completed its three numerical R oracles and
the analytical Rust fixture, then failed after 472.6 seconds at the immutable
resume check. The resumed group runner printed its usage message because the
validator's `system2()` call did not shell-quote the Rust-library path
`/mnt/e/Repositories/Jonah/PH Pipeline Repo/...`. The embedded spaces split one
logical argument into multiple command-line arguments.

The validator consequently reported only `immutable_resume` as failed. It
published no public validation directory or authorization decision. Its
private v2 oracle interval file is retained as failed-attempt evidence. The
7,626-pair primary and repeat stress artifacts were not modified, and the
remaining 19 landscape groups plus all downstream analysis remained closed.

## Correction and verification

Commit `0539857` constructs the seven runner arguments separately and applies
`shQuote()` to every argument passed through `system2()`. This is a validator
portability correction only; it does not change the production runner, Rust
library, diagrams, landscapes, pair axis, tolerances, resource gates, or
scientific estimand.

A static MV7-H contract test requires the quoted-argument path. A direct
diagnostic then invoked the corrected resume command with the actual Rust path
containing spaces. The runner reused the completed stress group successfully,
and all three production artifacts retained identical SHA-256 hashes and
modification times. The focused MV7-H suite passes 65 expectations.

A final independent-validation attempt is authorized only in a new private v3
oracle root and the still-unused public validation root. The remaining 19
landscape groups may proceed only after the validator records an explicit 8/8
pass and serial authorization.
