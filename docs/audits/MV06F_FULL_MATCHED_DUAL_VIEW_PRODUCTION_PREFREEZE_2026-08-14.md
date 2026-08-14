# MV6-F full matched dual-view production prefreeze

| Field | Result |
|---|---|
| Date | 2026-08-14 |
| Prospective specification commit | `3c396e3` |
| Queue groups | 75 |
| Queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Bound implementation root | `599074b3cd078cf27eb4a85148eb1df2ce3f84a5bdfd3160617b80a78f78c05e` |
| Rust library SHA-256 | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Independent validation | 10/10 categories pass |
| Focused tests | 19/19 expectations pass |
| Complete source-loaded suite | Pass; 2 established optional skips |
| Public deterministic repeat | 6/6 files byte-identical |
| Production groups executed | 0/75 |
| Labels/outcomes | Closed / zero |
| Decision | `prefreeze_pass_stage1_only` |

## Outcome

The complete matched cell/gene production scope is now explicit, reproducible,
and bound to code before calculation. The label-closed 90-sample manifest,
75-row LOSO plan, 450 accepted sample-seed caches, fixed 500-gene global panel,
and accepted Rust library all match their frozen hashes. An independent
validator reconstructs every group and workload total from those inputs.

No production group has run. The prefreeze authorizes only the separately
monitored stage-1 maximum group, followed by its clean repeat, R and Persim
oracles, and immutable-resume test. The other 74 groups remain closed until
stage 1 passes every prospective continuation gate.

## Complete workload

| Work family | Count |
|---|---:|
| Fold-seed groups | 75 |
| Typed PH jobs | 13,500: 6,750 cell and 6,750 gene |
| H0/H1 diagram-dimension records | 27,000 |
| Directed held-out-query to training-reference pairs | 35,350 |
| Separate cell-H0/cell-H1/gene-H0/gene-H1 distance rows | 141,400 |

Every group contains all 90 samples in both observation orientations. For a
group with `q` held-out and `t` training samples, it contains 180 PH jobs, 360
dimension records, `q * t` directed biological pairs, and `4 * q * t`
component rows. The independently reconstructed totals match the prospective
MV6-C workload exactly.

The label-free fail-fast schedule selects
`large_loso_v1:SRA779509` at seed `20260807` for stage 1. It has 25 held-out,
65 training, 1,625 directed biological pairs, and 6,500 component rows. It was
selected by maximum `q * t`, fold-ID tie break, and the frozen middle-seed rule;
no tissue, assay, endpoint, outcome, PH value, landscape distance, or fusion
result entered selection.

## Bound implementation

The public source inventory binds 14 accepted dependencies, scientific
contracts, reference runners, the new MV6-F helpers, the atomic group runner,
prefreeze builder, independent validator, and focused test by path, byte size,
and SHA-256. Their ordered hash map produces implementation root
`599074b3cd078cf27eb4a85148eb1df2ce3f84a5bdfd3160617b80a78f78c05e`.

The group runner performs these operations within one immutable group:

1. rehash all inputs, implementation sources, and the Rust library;
2. validate all 90 cache hashes, cache keys, and selected-cell identities;
3. fit scaling and 30-PC coordinates from training samples only;
4. construct matched cell and gene views for each sample;
5. run complete typed H0/H1 PH and retain MST audit evidence;
6. discard each view and explicit gene-distance object after PH;
7. compute all four exact all-active landscape components through Rust; and
8. validate and atomically rename one complete group directory.

A completed directory contains diagrams, a diagram manifest, component
distances, metrics, and a status binding every artifact hash and byte size.
Resume validates all identities and artifacts without rewriting. A stale,
one-sided, or partial directory fails closed and must be quarantined; there is
no automatic retry or silent fallback.

## Landscape and view definitions preserved

The new runner does not redefine the scientific objects. Cell topology still
uses 384 cells in shared training-fitted 30-PC Euclidean coordinates. Gene
topology uses the same 500 genes and 384 cells under Pearson-correlation chord
distance. They remain distinct observation spaces.

H0 and H1 remain separate. Only finite positive-persistence intervals enter
landscapes; every active consecutive level is retained; absent depth is
zero-padded; exact squared L2 is calculated without a universal level cap or
uniform grid. The Rust library accelerates this same interval-level quantity.
Grouped Persim remains the portable exact fallback and R remains the canonical
oracle, but neither is silently substituted during a Rust production group.

## Streaming, resources, and stop rules

The prior nonstreamed 8.739-GB projection was dominated by 6.961 GB of gene
views and 0.878 GB of cell views. MV6-F prohibits serializing those views.
Only compact fold identity/model evidence, typed PH records, component rows,
metrics, and status are durable. A corrected streaming projection and at least
25% validation/quarantine reserve are required before stage 1.

The frozen guards are 1,800 seconds and 8 GiB per group, 1,200 seconds and
6 GiB for stage-1 continuation, 12 GiB aggregate concurrent RSS, 14.4 complete
production worker-hours, and 10 GiB private storage. The maximum accepted
source+PH+Rust projection is 9.349 worker-hours before validation reserve.
No cap may be satisfied by dropping samples, groups, seeds, H1, active levels,
components, or validation.

Ten public abort rules cover input/implementation/Rust drift, label leakage,
axis drift, PH/MST failure, kernel/numerical failure, partial state,
oracle/cross-engine failure, repeat/resume failure, and resource/storage breach.
Every rule stops new launches and preserves already validated groups.

## Validation and deterministic evidence

The independent prefreeze validator passes all ten categories:

- frozen input and Rust hashes;
- input schemas/cardinalities;
- independent 75-group reconstruction;
- exact complete-workload totals;
- label-free stage-1 selection;
- all 14 source hashes;
- seven positive resource guards;
- ten fail-closed abort rules;
- queue/implementation/Rust contract roots with production still zero; and
- public label-firewall enforcement.

The focused R tests pass 19/19 expectations for exact workload reconstruction,
stage-1 selection, queue-root stability, axis and label failure, and directional
view/dimension-specific pair identities. They also validate a complete synthetic
180-PH-record group and prove that companion-artifact hash drift fails closed.
The builder was executed twice after
the group runner entered the source root; all six public artifacts retained
identical SHA-256 hashes and byte sizes.

The complete source-loaded suite also passes. Its two skips are the established
optional Rust-library-present fixture and public-audit-in-build exclusion, not
MV6-F failures.

The accepted original R project reports that its local `renv` state is
out-of-sync when these checks start, as in earlier MV6 runs. The bound scripts
nevertheless execute under the already accepted dependency environment, and
the focused tests/validators pass. No dependency was installed or changed.

## Gate disposition

MV6-F prefreeze completes as `prefreeze_pass_stage1_only`. The next sprint may
implement the external process-tree monitor, verify a corrected streaming
storage estimate, and launch only the frozen 6,500-row stage-1 group under one
worker. Stage 2 requires the complete stage-1 R/Persim oracle, clean-repeat,
resource, privacy, and zero-rebuild-resume dossier.

This prefreeze does not authorize the remaining 74 groups, fusion fitting,
blocked outcomes, clustering, G-MV6, new data, advanced network fusion, public
Rust defaults, binary distribution, manuscript claims, release actions,
source PDFs, confidential reviewer material, `example_run.r`, or a push.
