# MV6-G stale rank-block remediation

## Preserved stop

The serial driver completed nine schema-corrected groups, then stopped at
execution order 11. The failed group ran 311.165 seconds at 184,246,272 bytes
peak process-tree RSS and published no scientific directory. No later group,
label, endpoint, or outcome job launched. The nine complete directories, one
partial, ten resource metrics, and private child logs are preserved together
in the stale-index attempt quarantine/evidence locations.

## Root cause

The ranking builder created query-by-method block row indices before sorting
its result into canonical query/method/rank order. Its final invariant then
applied those stale pre-sort indices to the sorted table. For the stopped group
there are 81 training samples and nine methods: all 81 blocks pass before the
sort, while all 81 falsely fail when the interleaved pre-sort indices are used
afterward. The 6,561 ranking values and their ranks are not implicated.

The correction recomputes block indices after canonical sorting before the
final permutation check. Regression fixtures cover 72 and 81 training samples,
the multiples of nine in the supported production range. Landscape distances,
training-only scales, nine formulas, canonical order, and the tie rule are
unchanged.

Because this corrects a source bound into the scientific implementation root,
the existing nine successful directories are preserved but not reused.
Production must restart cleanly only after the corrected general runner
reproduces the accepted maximum-group sentinel byte-for-byte under the live
caps and all prefreeze/test gates pass.
