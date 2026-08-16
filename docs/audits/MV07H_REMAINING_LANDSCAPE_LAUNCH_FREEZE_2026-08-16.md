# MV7-H remaining-landscape launch freeze

Date: 2026-08-16

Status: accepted; serial production authorized

## Frozen workload

The accepted stress decision authorizes exactly queue rows 2 through 20 from
`docs/audits/mv07h-fallback-prefreeze-evidence/mv07h-landscape-queue.csv`.
These 19 groups complete the five-seed by cell/gene by H0/H1 Cartesian product
after the already completed seed-20260807 gene-H1 stress group. Each group has
124 samples, all 7,626 unordered pairs, one worker, and zero retries. The full
corpus remains 152,520 view- and dimension-specific component rows.

## Execution monitor

The serial monitor is commit
`797075f934a8c901f9b485737c064c5b037c0034`; its file SHA-256 is
`14855ae77d95ab98974d4613710d02f2b9fe316308e13d45b5a6d893a3170c18`.
It verifies the frozen queue, contract, Rust library, stress authorization,
stress resource record, and completed stress directory before launching any
stage-two group.

For every group it:

- invokes the accepted atomic group runner through `processx`;
- enforces the frozen 3,600-second and 12-GiB process-tree limits;
- enforces the original 23,232.367-second total landscape budget, including
  the accepted primary stress run;
- validates all 7,626 rows, identities, hashes, exact/all-level flags, and
  H0/H1 separation before checkpointing;
- writes an atomic private resource checkpoint after each group;
- reuses only fully validated groups after interruption; and
- preserves any recorded failure and refuses an automatic retry.

The focused MV7-H suite passes 78 expectations, and the monitor parses under
the production WSL R runtime. The launch remains label-closed: dimension
combination, clustering, outcomes, and result-dependent claims are not
authorized. Private execution uses the existing stress root
`tmp/mv07h-landscape-stress-v2`, with checkpoint
`tmp/mv07h-landscape-stress-v2/stage2-resource.csv`. Public completion evidence
will be written only by a later independent validator.
