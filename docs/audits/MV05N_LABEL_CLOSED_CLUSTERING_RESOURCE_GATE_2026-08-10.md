# MV5-N label-closed sample-clustering and complete-matrix resource gate

## Decision summary

| Question | Result |
|---|---|
| Scientific contract coherent? | **Approve** training-only PAM and average-linkage sensitivity |
| Landscape definition preserved? | **Yes**: separate H0/H1, all consecutive active levels, exact critical-pair integration, no cap/grid |
| Exact missing scope reproduced? | **Yes**: 262,675 pairs / 525,350 H0-H1 rows per representation |
| Bounded real correctness | **Pass**: 384/384 rows, 12/12 R exact oracles, 12/12 byte repeats |
| Matched baselines | **Pass**: 384/384 rows; 96/96 cross-representation pseudobulk identities; byte repeat exact |
| Combined resource feasibility | **Yes**: 16.117 conservative worker-hours including 10% reserve, below 21.6-hour cap |
| Biological interpretation | **Prohibited**: labels remained closed; no clustering outcome was run |
| Full production authorized here? | **No**: requires a separate prospective execution authorization |

## 1. What changed

MV5-N converts the original sample-level clustering idea into a strictly
inductive design. Every partition is trained only on the samples outside one
held-out study. Candidate `k=2:min(10,n-1)` is selected from five-seed
partition stability with the previously frozen one-SE rule. Held-out samples
are assigned to frozen training medoids, never included in fitting.

PAM is primary. Average linkage is the sole eligible sensitivity and receives
the PAM-selected `k`; its held-out rule is minimum mean distance to a frozen
training cluster. Both methods use canonical member-based cluster labels and
fully specified tie-breaking. Spectral clustering remains ineligible pending a
separate affinity/out-of-sample implementation gate. Ward and direct k-means
on distance matrices remain excluded.

No tissue, approach, class, label, or outcome column is accepted by the
pair-generation, clustering, selection, or held-out-assignment helpers.

## 2. Landscape contract

The implementation preserves the revised dissertation-aligned definition:

- corrected cells-as-observations diagrams from MV5-D3 and MV5-H;
- finite positive-persistence intervals;
- essential H0 excluded before landscape construction;
- H0 and H1 separate;
- every consecutive active landscape level;
- exact piecewise-linear integration over critical-pair segments; and
- no universal level cap or uniform grid.

The raw `sqrt(H0^2 + H1^2)` matrix remains descriptive only. H0 and H1 are not
collapsed for primary evaluation.

## 3. Exact pair inventory

All pair identities were instantiated, validated, and digested group by group,
then non-admitted wide rows were discarded. This avoids committing or
maximum-compressing approximately one million generated request rows while
retaining exact reproducibility.

| Representation | Groups | Training pairs | H0 rows | H1 rows | Chunks |
|---|---:|---:|---:|---:|---:|
| SCT whole | 75 | 262,675 | 262,675 | 262,675 | 2,170 |
| Inductive integrated | 75 | 262,675 | 262,675 | 262,675 | 2,170 |

The public global, group, and chunk inventories bind source PH manifests,
view-resource manifests, canonical sample pairs, record keys, diagram/file
hashes, dimensions, and request-set hashes. A later authorized queue can stream
the exact identities from the same generator.

Identity generation on the Windows fallback runtime took about 325 seconds at
roughly 170 MB working set. An initial row-wise/wide maximum-compression design
was terminated after profiling showed pathological I/O overhead. The accepted
generator vectorizes pairs and publishes compact digests; this engineering
change does not alter any request identity or scientific calculation.

## 4. Bounded real admission

The profile selection used no outcomes:

| Profile | Fold | Training samples | Seed |
|---|---|---:|---:|
| Minimum | `large_loso_v1:SRA779509` | 65 | 20260805 |
| Representative (median size) | `large_loso_v1:SRA703206` | 86 | 20260805 |
| Maximum | `large_loso_v1:SRA713577` | 89 | 20260805 |

For each profile, representation, and H0/H1 dimension, 32 equally spaced
canonical request ranks were admitted: 12 groups and 384 rows total. Staging
validated every accepted PH-result file and diagram hash and extracted 197,254
finite interval rows.

The accepted `persim` 0.3.8 exact critical-pair source was reused from the
existing local project environment because the Ubuntu distribution became
temporarily unavailable during the sprint. The pure-Python exact landscape
modules ran under the bundled Windows Python; no third-party source was added
to Git. This fallback is scientifically cross-checked rather than assumed:

- 384/384 rows completed;
- H0/H1 each contributed 192 rows;
- all outputs report exact/all-active-level/no-cap status;
- maximum group elapsed time was 10.270 seconds;
- summed pair-operation time was 11.386 seconds;
- observed peak process RSS was 556,314,624 bytes (530.5 MiB);
- 12/12 independently implemented R exact breakpoint-stream oracles passed;
- maximum Python-versus-R absolute difference was `2.8422e-14`;
- 12/12 clean-repeat output files were byte-identical; and
- a resume rerun reused 12/12 groups while all 24 output/status hashes, sizes,
  and timestamps remained unchanged.

## 5. Matched baseline admission

The same 96 unique admitted training pairs were evaluated in each
representation for cell-distribution energy and training-standardized
pseudobulk, producing 384 baseline rows. All distances were finite and
nonnegative. The pair artifact repeated byte for byte. All 96 pseudobulk pairs
were exactly identical between the SCT and integrated strata, as required for
the shared context control.

Synthetic tests independently preserve the accepted empirical V-statistic
energy formula and pseudobulk Euclidean formula. No outcome label was used.

## 6. Conservative full-resource projection

For landscapes, each representation/dimension uses the maximum observed
operation-seconds-per-row across minimum, representative, and maximum profiles,
plus the maximum observed fixed group overhead times all 75 groups. Baselines
use the maximum profile seconds per pair. Pseudobulk is computed once and
reused because its values are representation-identical. A further 10% reserve
covers validation, I/O, and slow-tail variation.

| Component | Projected worker-hours |
|---|---:|
| SCT H0/H1 landscapes | 7.772 |
| Integrated H0/H1 landscapes | 2.890 |
| SCT energy | 1.961 |
| Integrated energy | 2.007 |
| Shared pseudobulk | 0.023 |
| 10% reserve | 1.465 |
| **Total** | **16.117** |

The total is below the 21.6-worker-hour planning cap. Projected landscape
output is about 619 MB; all admission groups passed 900-second and 4-GiB
guards. These measurements support feasibility, but MV5-N intentionally does
not authorize full production.

## 7. Code and fixture validation

The focused suite contains 29 passing assertions covering:

- canonical unordered pair/H0/H1 identity and chunking;
- representation-specific identity roots;
- prohibited outcome columns;
- deterministic minimum/median/maximum profile selection;
- canonical cluster labels;
- the exact five-seed candidate grid and one-SE `k` selection;
- PAM and average-linkage partitions;
- PAM medoid and average-cluster held-out tie rules;
- energy and pseudobulk baseline formulas; and
- incomplete matrices, seeds, and held-out distance rejection.

The files were copied into the canonical repository and the focused suite
passed there with 29/29 assertions. All eight MV5-N R implementation scripts
parse, and the Python worker passes syntax-tree parsing.

An initial sandboxed `wsl.exe -l -v` call incorrectly reported no registered
distributions. The required elevated read-only check showed the established
Ubuntu and Ubuntu-18.04 distributions; Ubuntu remained the default. In the
established Ubuntu R 4.4.1 environment, the complete source suite passed
759/759 assertions with zero failures, warnings, or skips.

The clean source-package check used `git archive` at full revision
`d62988860e85f6f244015413d19aa239b7d0a7da`, excluding every ignored and
untracked artifact by construction. `R CMD build --no-resave-data` produced a
237,080-byte source tarball with SHA-256
`b5ee48c0b34d5af8c6d75944cc5d815ead735494e5b260db97cef140965dc111`.
`R CMD check --no-manual` on Ubuntu 22.04.4 LTS/R 4.4.1 completed with
`Status: OK`; its check-log SHA-256 is
`97af1aad9a7d98d0c78f81dfc9c9c77b8d286ad1339397224fef31f23e9117d6`.

The first live-tree build attempt was stopped after its pre-exclusion copy
reached 28 GB because ignored scientific artifacts were being copied into the
temporary build tree. That exact temporary directory was removed, and no
result from the aborted attempt is treated as evidence. The accepted Git-only
check is both safer and a more direct test of the publishable repository.

The public artifact manifest binds 37 decision, implementation, test,
evidence, and upstream files by path, byte size, and SHA-256. A clean rebuild
of that manifest was byte-identical.

## 8. Decision and next action

The MV5-N scientific, computational, focused-test, complete-suite, and clean
source-package gates pass. Full matrix production remains closed until a
separate MV5-O execution specification binds implementation/source hashes,
queues, caps, validation sampling, repeat, resume, and stop rules before
launching any of the 1,838,725 required distance values.

MV5-O, if authorized, may produce only label-closed training distance matrices
and immutable cluster predictions. Biological/technical label opening and
cluster ARI/NMI must remain a later, separately prediction-locked sprint.
That later sprint must report descriptive alignment of frozen training
partitions separately from inductive generalization of immutable held-out
assignments; the former cannot stand in for the latter or retune the pipeline.
