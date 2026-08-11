# MV5-P label-closed complete-training-distance production

## Decision summary

| Question | Result |
|---|---|
| Frozen MV5-O scope preserved? | **Yes**: 18 sources, 150 groups, and 4,565 unit identities |
| Landscape definition preserved? | **Yes**: separate exact all-active-level H0/H1, no cap or fixed grid |
| Complete distance production? | **Yes**: 4,340 landscape, 150 energy, and 75 shared-pseudobulk units |
| Matrix validation? | **Pass**: 525/525 complete, finite/nonnegative, symmetric, zero-diagonal components |
| Independent correctness? | **Pass**: 12/12 exact R oracles; maximum absolute difference `2.1316282072803e-14` |
| Deterministic repeat? | **Pass**: 66/66 frozen maximum-group outputs; supplemental SCT pseudobulk also passes |
| Immutable resume? | **Pass**: 4,565/4,565 units unchanged; zero rebuilt |
| Resource envelope? | **Pass**: `12.044379` worker-hours, `492163072` B maximum process-tree RSS, `4570070656` B private root |
| Labels, outcomes, or clustering opened? | **No** |
| Canonical verification? | **Pass**: 23/23 focused, 816/816 full, package `Status: OK` |

## 1. Scope and frozen identities

MV5-P executed only the complete-training-distance queues frozen by MV5-O at
base commit `46765a7`. An independent current-state audit rehashed all 18
source files, seven execution/contract implementations, 150 group identities,
4,340 landscape chunks, 225 baseline units, the 15-row validation plan, and
all ten non-retrying abort rules. The queues account for exactly 4,565 unique
units and 1,838,725 values.

No accepted MV5-D through MV5-O file was modified. MV5-P added bounded
orchestration and validation only. The dissertation-aligned landscape target
remained separate H0/H1 L2 distance over all active consecutive landscape
levels, evaluated exactly by critical-pair integration, with no universal
level cap and no fixed uniform grid.

## 2. Execution and recovery

The ignored production root used deterministic execution order and no more
than two group workers. Every group/unit retained a 900-second bound, every
process tree a 4-GiB RSS bound, the stage a 21.6-worker-hour bound, and private
storage a 10-GiB bound. Outputs and statuses were written atomically and could
be reused only after identity/hash validation.

The first monitor stopped after two completed groups because an orchestration
metrics replacement guard interpreted R's successful `unlink()` return value
incorrectly. No scientific output failed, no partial artifact existed, and no
completed output was changed. A direct same-filesystem atomic rename replaced
the guard. A clean two-group resource replay reproduced every training input,
landscape output, and energy output byte-for-byte; conservative replay timing
and RSS evidence were used for those two groups. Deterministic execution then
resumed at order 3 and completed all 150 groups.

## 3. Complete output accounting

The accepted production contains:

- 4,340 landscape chunks and 1,050,700 exact H0/H1 distances;
- 150 representation-specific energy units and 525,350 distances;
- 75 SCT-hosted shared pseudobulk units and 262,675 distances;
- 4,565 immutable units and 1,838,725 total values; and
- 525 complete training-distance components: 300 landscape H0/H1, 150
  energy, and 75 shared-pseudobulk matrices.

Every output/status hash and size matches its manifest. Every distance is
finite and nonnegative. Every landscape output declares exact integration,
all active levels, no level cap, and zero estimated absolute error.

## 4. Matrix and independent correctness evidence

Each of the 525 sample-distance components contains the exact expected unique
unordered training pairs, expands to a symmetric matrix, has an all-zero
diagonal, and contains only finite nonnegative entries. H0 and H1 remain
separate components; no outcome-based combination or method selection was
performed.

The 12 exact R oracle requests frozen before production span minimum,
representative, and maximum groups in both SCT and inductive-integrated
representations and both H0 and H1. All 12 agree within the frozen `1e-10`
tolerance; the maximum absolute difference is `2.1316282072803e-14`.

## 5. Repeat and immutable-resume evidence

Clean repeats of the maximum group in each representation reproduced all 32
landscape chunk outputs and the representation-specific energy output
byte-for-byte: 33/33 per group and 66/66 total. The SCT runner necessarily also
materialized its shared-pseudobulk context; that output was excluded from the
frozen 33-unit count as prespecified and was separately reported as an
additional passing byte repeat.

A before/after snapshot bound output and status hashes, sizes, and timestamps
for all 4,565 units around a completed monitor rerun. All 4,565 remained
unchanged and zero units were rebuilt.

## 6. Resource and storage reconciliation

Observed execution used `12.044379` worker-hours, with maximum measured
process-tree RSS `492163072` bytes and final private production-root
size `4570070656` bytes. All are within the 21.6-hour, 4-GiB, and 10-GiB
caps, and no group failed.

MV5-O's 1,277,893,355-byte private-storage estimate was output-focused: it
included landscape outputs, baseline bytes, and reserve but not complete
group-local repeated interval staging. The live prefix audit disclosed this
before cap pressure. Scaling all frozen group row burdens by the highest
observed bytes-per-row rate and then adding 25% reserve forecast
`6011076897` bytes, `56.0`% of the
10-GiB cap. Final observed storage remained cap-passing.
This projection miss is recorded explicitly rather than silently treating the
prefreeze estimate as accurate.

## 7. Label firewall and scientific boundary

All group, unit, validation, and public-evidence rows remain
`outcome_label_state = closed`; biological outcomes are false and clustering
job counts are zero. MV5-P did not fit PAM or average-linkage clusters, select
`k`, assign held-out samples, open biological or technical labels, calculate
ARI/NMI or other alignment/generalization outcomes, select methods or tissues,
or execute robustness, gene topology, fusion, new data, optimization, or Rust.

The result is therefore technical completion of label-closed distance inputs,
not a biological conclusion and not evidence that one representation or
method performs better.

## 8. Canonical verification and repository boundary

Focused validation passed 23/23 expectations and the complete source suite
passed 816/816 with zero failures, warnings, or skips.
A clean Git-archive package passed `R CMD check --no-manual` with
`Status: OK`; final tarball SHA-256 is
`d5627e42ff43e3a8e316311a02d5c919a4497cfb0d3704c97ad9ea0c85c54673`
and check-log SHA-256 is
`9e6ede9e4ac3fc3975ab74d280f327628363d02d7d22314b5b36ec2dc2695ab5`.

The dissertation and preprint PDFs remain ignored, confidential reviewer
material and production caches remain private, and `example_run.r` remains the
sole preserved untracked canonical file. No branch was pushed.

## 9. Decision and next action

MV5-P passes its technical completion gate. The complete label-closed training
distance inputs are now eligible for a separately scoped clustering-execution
sprint under the already frozen MV5-N contract. That next sprint must still
stop before labels and outcomes: it may fit training-only PAM, perform
label-free five-seed one-SE `k` selection, retain average linkage solely as the
frozen sensitivity, canonicalize cluster identities, and assign held-out
samples to training medoids without opening metadata outcomes.

Biological or technical outcome evaluation requires another prospective gate
after the cluster artifacts themselves pass identity, determinism, resume,
and resource validation.
