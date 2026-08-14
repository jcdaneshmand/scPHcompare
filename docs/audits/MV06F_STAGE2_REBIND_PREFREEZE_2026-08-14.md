# MV6-F stage-two remediated-root rebind prefreeze

| Field | Result |
|---|---|
| Date | 2026-08-14 |
| Parent stage-one commit | `93f4bef` |
| Frozen queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Parent implementation root | `599074b3cd078cf27eb4a85148eb1df2ce3f84a5bdfd3160617b80a78f78c05e` |
| Remediated implementation root | `5a1258e87eb30c367648daa4bb02d1aec1f4b40dd3799a91fbb7067e9558d292` |
| Rust SHA-256 | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Independent prefreeze validation | 6/6 categories pass |
| Focused expectations | 23/23 pass across MV6-F production and stage-two tests |
| Complete package-aware suite | Pass in 161.7 s; two established optional skips |
| Decision | `prefreeze_pass_reexecute_stage1_only` |

## Why a rebind is required

The post-merge P1 remediation changed `R/dual_view_topology.R` and
`R/mv05_resource_safe_execution.R`, both members of MV6-F's accepted
implementation inventory. The frozen runner therefore rejects the old source
root before calculation, as intended. The PCA correction also adds the fitted
standardization identity to provenance and cache identity; numerical output is
expected to remain unchanged, but artifact identity cannot be assumed.

Stage two must not mix groups produced under different implementation roots.
This sprint therefore binds all 75 eventual groups to one remediated root and
requires the maximum stage-one group to be recalculated in a new private root
before any of the remaining 74 groups run.

## Frozen orchestration

The original 75-row queue, fold/seed ordering, sample roles, fixed 500-gene
panel, and Rust binary are unchanged. The new 23-file source inventory binds
the corrected scientific dependencies, atomic group runner, stage-one monitor,
R and Persim oracle runners, immutable-resume checker, stage-two monitor,
complete-production validator, and focused tests.

Stage-two execution is serial (`maximum_workers = 1`) because concurrent
production has not been empirically admitted. Every child process is monitored
live with these unchanged conservative caps:

- 1,800 seconds and 8 GiB process-tree RSS per group;
- 12 GiB aggregate live RSS;
- 14.4 cumulative production worker-hours; and
- 10 GiB for the complete private production root.

The monitor checkpoints each completed group, validates all artifact hashes
before further admission, resumes validated groups without rewriting, refuses
partial/unexpected directories, and preserves any failed row without automatic
retry. It contains no label, fusion, clustering, or outcome path.

## Required rebind gate

The new root authorizes only one remediated maximum-group primary run and clean
repeat. Stage two remains false in the prefreeze contract until all of these
conditions pass:

1. old and remediated group directories validate under their own roots;
2. all 180 diagram matrices and scientific manifest fields are equivalent;
3. all 6,500 exact cell/gene H0/H1 landscape rows are equivalent;
4. three new-root scientific artifacts repeat byte-identically;
5. 12 balanced R landscape oracles pass;
6. 12 balanced grouped-Persim landscape oracles pass; and
7. all five group files survive a zero-rebuild resume unchanged.

Only the resulting admission CSV can open the 74-row stage-two monitor.

## Verification and command-route note

The builder reproduced the source inventory, contract, and resource plan with
identical hashes and bytes. The independent validator passed 6/6 categories,
and the focused tests passed 4/4 new expectations plus 19/19 existing MV6-F
expectations.

The first remediated root (`816186b9…a6584`) completed a primary, repeat, R
oracle, and Persim oracle, but its resume checker exposed an unquoted WSL path
containing spaces before the child runner could start. No scientific artifact
was rewritten. Because the checker belongs to the frozen inventory, that root
was invalidated rather than patched in place; its outputs were moved intact to
an ignored quarantine directory. The corrected checker is included in the
replacement root above, which must repeat the complete stage-one gate.

A bare `testthat::test_dir()` diagnostic was also attempted and failed with
missing-function errors because this legacy suite requires the package
namespace to be loaded. The supported package-aware `devtools::test()` route
then passed completely in 161.7 seconds with only the two established optional
Rust/public-audit skips. The unsupported diagnostic did not execute production
or reveal an MV6-F regression.

## Boundary

This prefreeze does not yet authorize the remaining 74 groups. It does not
authorize fusion, clustering, outcomes, public Rust adoption, binary release,
manuscript claims, new data, source PDFs, confidential reviewer material, or
`example_run.r`.

Public evidence is stored in
`docs/audits/mv06f-stage2-prefreeze-evidence/`. Caches, binaries, logs, diagrams,
and private interval staging remain excluded from Git.
