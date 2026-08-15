# MV6-G complete production validation audit

## Disposition

`complete_validation_pass_resume_pending`

The corrected-root, label-closed MV6-G prediction corpus contains the accepted
maximum-group sentinel plus all 74 clean serial-completion groups. Independent
complete-corpus validation passes 7/7 categories. Labels, outcomes, clustering,
advanced fusion, and claims remain closed until the separately committed
445-file immutable-resume gate passes.

## Frozen identities

- scientific implementation root: `8b0a1e42d9e46234edb847bcab82a2e770a47f3934d21f2a81c352c74c9cec0c`
- serial execution root: `deb03fbc3e3bcf49d7c6fcbd356256c94aae0e3dae7672ec4150b3702e36f745`
- Rust library: `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`
- complete validation queue: `77725309f87e90ce40b545030c18cd0d20ab794b6f938e0c90b597932cdbd314`

## Observed production

- 74/74 clean-restart groups completed with no scientific or resource failure;
- 20,965.105 charged worker-seconds (5.823640 hours), below the frozen
  43,200-second ceiling;
- maximum group elapsed time 327.290 seconds, below 1,800 seconds;
- peak process-tree RSS 186,503,168 bytes, below 12 GiB;
- final private state 531,251,838 bytes, below 5 GiB;
- labels remained closed and fusion-evaluation/outcome-job counters remained
  zero throughout.

The original 4.747-hour sentinel projection underestimated observed execution
by about 22.7%; the prospectively binding 12-worker-hour gate nevertheless had
substantial headroom. This projection miss must be retained in future resource
planning.

## Complete scientific validation

The independent validator accepted 75/75 groups and passed 7/7 categories:

- unique complete group inventory;
- 262,675 training biological pairs and 1,050,700 component rows;
- 35,350 query biological pairs and 318,150 ranking rows;
- exactly four training-scaled components per group;
- 74 exact-order serial resource metrics;
- closed labels and zero biological/fusion outcomes;
- no partial group state.

Evidence identities:

- canonical resource metrics: `49e2d05f157c2a0d9fdbaf3247a6aeef35f2100331acd9df9b261a3fb28beddf`
- complete group inventory: `9ab127ea07e5a83c70898b17e247948079bc6c8ba81ff0711d674dbb3060ea36`
- validation checks: `29b97f348067d9ab468339d28b84c82c3a1c684c9bb921bdd4abc0f4e79e944b`

## Bounded execution remediations

After all 74 groups completed, the original driver could not rename its final
candidate from WSL `/tmp` to the Windows-mounted worktree (`Invalid cross-device
link`). No scientific group failed or changed. The destination-atomic recovery
utility `a98ab7a2...307a2a`, outside the nine-file accepted execution root,
revalidated all 74 exact-order metrics, wrote its candidate beside the target,
and published the byte-deterministic canonical metric. Its focused regression
passed.

The first complete-validator invocation used the older stage-one status schema;
the second used the accepted corrected-root `mv06g-production-rebind-v2`
sentinel. The full workload also lacked the exact `biological_pairs` runner
alias already enforced in the 74-row production queue. Utility
`039ea626...11bec` materialized the 75-row queue with
`biological_pairs == query_biological_pairs`. These were provenance/input-schema
corrections only. A two-minute shell timeout terminated one read-only validator
attempt before evidence publication; the unchanged ten-minute rerun completed
in 126.7 seconds and passed 7/7.

## Next gate

Run the prefrozen complete-resume checker and require SHA-256, byte size, and
mtime preservation for all 445 scientific/public files before prediction-lock
or outcome-prefreeze work begins.
