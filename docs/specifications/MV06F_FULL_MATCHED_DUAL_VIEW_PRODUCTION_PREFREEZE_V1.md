# MV6-F full matched dual-view production prefreeze v1

## Document control

| Field | Frozen value |
|---|---|
| Status | Prospective execution contract; production not launched |
| Date | 2026-08-14 |
| Base revision | `78549adfce444d7bb664ba61e219d9d7e77dfb3b` |
| Parent decision | MV6-E `admit_both_with_rust_preferred_private` |
| Representation | Fixed global-core SCT cell and gene views |
| Panel | 500 genes; scientific SHA-256 `7be22cdb9056fed427c78d58be2b19e258c7c6807e6f7ac3900bd1bfa1d19eb8` |
| Production engine | Accepted Rust exact kernel, explicit private WSL only |
| Portable fallback | Grouped Persim exact critical-pair engine |
| Labels/outcomes | Closed / zero |
| Required stop | Before fusion, clustering, label opening, or biological evaluation |

## 1. Purpose and boundary

MV6-F converts the passing MV6-E acceleration gate into one immutable,
label-closed production plan for the complete matched cell/gene corpus. It
freezes all 75 fold-seed groups, training-only transforms, typed PH jobs,
landscape component rows, engine identities, streaming behavior, caps,
validation samples, clean repeat, resume, abort rules, and public/private
boundaries before the first production group is calculated.

This prefreeze does not execute a production group. It does not normalize or
combine component matrices, fit fusion weights, cluster samples, access a
biological or technical outcome, select a method, open G-MV6, or change any
package default.

## 2. Frozen upstream inputs

| Input | SHA-256 |
|---|---|
| MV5-D0 v2 cache-resource ledger | `73f757a91c202a8e38dfa746fc8816f8e272caf912c33b4b959ef55102c68308` |
| Label-closed 90-sample candidate manifest | `842c047ba821f8eca317da52504910733509fb4fddd11d6f54f7e79d9f29d0b7` |
| 75-row LOSO fold-seed plan | `50379f98cd4927c5c8cb19dbd9ca8ecc7b7b3a9af2e04eb9c8358ecb0b722c6d` |
| Public MV6-C panel CSV | `b3a5aff1a0bc01e871751fb9db0b3babfaf18835e68c5699346d8476d903d0ab` |
| MV6-D decision | `b534e0d6a582fffd7a5d2985a593c23939675d4c9179d33de8787d3115178089` |
| MV6-D worker projection | `64393e454a28a71f71f78fc52c335d5cd244cec7b7f74d7bda02a8144eeddead` |
| MV6-D storage projection | `e4bb8812a9a158c50cd53c14f26fb753a9326e14796c1fed68e2e9595b74289c` |
| MV6-E decision | `5a420f6d8bb9d441b9c6d5d87768d47fa2dc73e62f3708e4b23e94b9a16166a0` |
| MV6-E projection | `6425bac7e62dd575d5e32e8f963673b8aef3a0f169c0f0e6dea101bf5a91b446` |

The first three inputs may be read from the accepted private evidence store,
but their hashes and schemas must match the values above. The cache ledger must
resolve all 450 sample-seed cache files, and every cache must match its recorded
file hash, normalization key, selected-cell hash, and `built_atomic`
disposition. Any missing or mismatched input aborts before source construction.

The implementation freeze must bind, at minimum, these accepted scientific
sources plus every new MV6-F runner/validator by path, byte size, and SHA-256:

| Accepted source | SHA-256 |
|---|---|
| `R/mv06d_matched_profile.R` | `a8ef47fdc4d32386f79fd2411369d23b159f7c75782f00d869099c329def7f84` |
| `R/dual_view_topology.R` | `569ec50fad6d9edb2bdaa2be023a3178a86cc5098c19b4e821a0902b62060bce` |
| `scripts/run_mv06d_source_entry.R` | `e7f388efc35b47a845ff312a655176cf3d0ec127bc0e7e5b6e22039d0545ee64` |
| `scripts/run_mv06d_ph_entry.R` | `83f2c2060f82cc637add01b10b4102995939ec5c01a27c40c8c5822129fc9017` |
| `R/landscape_rust_prototype.R` | `4771ec396759dacaf36e860819b49978b37d8cc8d4224e54cc72b0077f8a5eb1` |
| `scripts/run_mv06e_rust_candidate.R` | `3a7e5331bf4063a1fed342a8c3dfeacc82783bfd2f93d22361598e2c32bd2ca6` |

## 3. Exact production scope and queue

The group axis is the Cartesian product of 15 LOSO held-out studies and the
five seeds `20260805` through `20260809`, ordered by fold ID and integer seed
using radix order. The public queue contains only technical identifiers and
counts; it contains no expression, cell identity, tissue, approach, endpoint,
outcome, or derived biological value.

| Work family | Frozen total |
|---|---:|
| Fold-seed groups | 75 |
| Samples/views per group | 90 cell + 90 gene |
| Typed PH jobs | 13,500: 6,750 cell + 6,750 gene |
| H0/H1 diagram-dimension records | 27,000 |
| Directed held-out-query to training-reference biological pairs | 35,350 |
| Separate cell-H0/cell-H1/gene-H0/gene-H1 rows | 141,400 |

For a group with `q` held-out samples and `t` training samples, `q + t = 90`,
the group must contain `180` typed PH jobs, `360` diagram-dimension records,
`q * t` directed biological pairs, and `4 * q * t` landscape component rows.
Held-out-to-held-out, training-to-training, query-to-query, reverse duplicate,
or diagonal biological requests are prohibited.

The implementation sprint must materialize and independently validate the
complete 75-row queue and a deterministic root over its scientific fields
before production authorization. Queue drift requires a new prefreeze.

## 4. Matched source and PH contract

The MV6-C panel is a fixed, explicitly transductive technical harmonization.
Conditional on that panel, each fold-seed group remains training-only:

1. validate and load the 90 accepted SCT caches for the seed;
2. extract the ordered 500-gene panel with no zero filling, substitution,
   panel shrinking, or sample exclusion;
3. fit gene centers and positive scales from pooled training cells only;
4. apply that same affine transform to training and held-out sources;
5. fit deterministic 30-component PCA on training sources only; and
6. construct matched 384-cell and 500-gene typed views for each sample.

The cell view is 384 cells in shared 30-PC Euclidean coordinates. The gene view
is the same 500 genes across the same 384 cells under Pearson-correlation chord
distance. The two observation spaces remain distinct and are never pooled into
one persistence diagram.

Every view calls the accepted typed PH route with complete Vietoris--Rips,
`max_dim = 1`, `threshold = -1`, and coefficient field 2. A view with `n`
observations must have `n - 1` finite positive H0 intervals plus exactly one
essential H0 interval. Every H1 interval must be finite and have positive
persistence. H0/H1 records retain typed source/view/cache identities and the
finite-H0 MST evidence required by the accepted MV6-D contract.

## 5. Landscape definition is immutable

The primary distance remains the dissertation-aligned definition:

- H0 and H1 are computed and stored separately;
- only finite positive-persistence intervals enter landscapes;
- every active consecutive landscape level is retained;
- missing landscape depth is zero-padded between diagrams;
- squared L2 is integrated exactly or under a recorded error certificate;
- no universal level cap or fixed uniform grid is used; and
- any later H0/H1 or cell/gene combination is derived and decomposable.

The production Rust kernel computes this exact quantity from finite intervals.
Before each group, its library must match SHA-256
`51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d`.
There is no silent engine fallback. A missing/mismatched library or nonzero
kernel status stops new launches. Switching to grouped Persim is allowed only
through a separately recorded recovery decision; grouped Persim remains the
portable scientific fallback and the R implementation remains the canonical
oracle.

## 6. Streaming and private storage

Within each group, the training transform may be held in memory, but each
sample's cell/gene view and explicit gene distance object must be discarded as
soon as its typed PH record and required audit evidence are complete. Production
must not serialize all 6,750 gene views or retain a cross-group dense landscape
matrix.

The durable private group payload contains only:

- the fold transform identity and compact model/audit fields required for
  independent reconstruction;
- 180 typed PH records (360 H0/H1 components);
- the group's canonical pair manifest or its complete reproducible identity;
- `4 * q * t` landscape component rows; and
- resource/status/hash evidence.

The MV6-D nonstreamed projection was 8,738,555,180 bytes, of which
6,961,389,750 bytes were serialized gene views and 877,925,250 bytes were cell
views. Those view payloads are prohibited here. Diagram, landscape, fold, status,
repeat, and validation artifacts must remain within a 10-GiB private-root cap.
The implementation preflight must produce a corrected streaming projection
before launch and reserve at least 25% for statuses, repeats, oracles, and
quarantined partial work.

## 7. Atomic group artifacts and immutable resume

Each group is one resumable transaction. It writes a same-parent temporary
directory containing diagrams, component distances, metrics, status, and
artifact hashes. After all schemas, counts, hashes, and label-firewall checks
pass, the temporary directory is atomically renamed to its deterministic final
group path. A status may not claim completion before every payload hash and byte
size is known.

Resume accepts a group only when its queue identity/root, source/input hashes,
implementation root, Rust hash, group output hashes/sizes/row counts, and
completion status all validate. A one-sided file, unexpected temporary
directory, stale hash, wrong row count, or mismatched implementation aborts;
production never overwrites or silently repairs it. Recovery quarantines the
partial group and requires an audited cause before retry.

Timing, RSS, and launch timestamps remain sidecars and are excluded from
scientific artifact identity. Scientific CSV/RDS serialization must be stable
under the clean-repeat environment.

## 8. Staged execution and concurrency

All BLAS/OpenMP/Rayon thread counts are one. Stage 1 runs with one worker. Stage
2 may use at most two group workers, provided aggregate process-tree RSS stays
below 12 GiB and each group remains independently atomic.

The stage-1 group is selected without labels: maximum `q * t`, then the
lexicographically first fold ID, then seed `20260807` when present (otherwise
the smallest seed). It executes once, validates, then executes into a clean
repeat root. Stage 2 is prohibited until stage 1 passes every correctness,
oracle, repeat, resource, privacy, and resume gate.

Production order for the remaining groups is descending `q * t`, then fold ID,
then seed. This is a technical fail-fast schedule, not a scientific priority.

## 9. Resource guards

| Guard | Frozen value |
|---|---:|
| One fold-seed group elapsed | 1,800 s |
| One group process-tree RSS | 8 GiB |
| Stage-1 continuation threshold | 1,200 s and 6 GiB |
| Aggregate concurrent process-tree RSS | 12 GiB |
| Complete production worker-time cap | 14.4 h |
| Private-root storage cap | 10 GiB |

MV6-D plus MV6-E maximum projections reconstruct 2.615 source hours, 1.984
cell-PH hours, 3.877 gene-PH hours, and 0.873 Rust-landscape hours: 9.349
worker-hours before validation reserve. The 14.4-hour production cap leaves a
54% margin without returning to the obsolete one-R-process-per-pair estimate.

Observed plus projected remaining work and storage are reconciled after every
completed group. A new group may launch only when both remain inside their
caps. No cap may be met by dropping samples, seeds, groups, H1, active levels,
component rows, validations, or provenance.

## 10. Correctness, oracle, repeat, and resume gates

Every group receives complete validation of input/source identities, typed-view
and PH schemas, one-essential-H0 accounting, positive finite intervals,
pair-axis membership, four-component completeness, Rust status/version/hash,
finite nonnegative distances, all-active/no-cap flags, exact group counts,
artifact hashes, label closure, resource metrics, and zero downstream counters.

Before stage 2, the maximum group must pass:

- independent regeneration of its queue and all identities;
- fresh Prim-MST checks for at least one held-out and one training sample in
  each view;
- 12 canonical R oracles: minimum, median, and maximum combined interval depth
  within each of cell-H0, cell-H1, gene-H0, and gene-H1;
- the same 12 rows through grouped Persim;
- a complete clean repeat with byte-identical scientific outputs; and
- a snapshot/resume/snapshot check with zero rebuilt or rewritten artifacts.

After all 75 groups, independent validation additionally requires:

- exactly 13,500 PH jobs, 27,000 diagram components, 35,350 biological pairs,
  and 141,400 landscape rows;
- 15 fresh fold-level source/metric audits at seed `20260807`, with one
  deterministic held-out and training sample per fold and both views checked;
- 20 R oracles stratified by view, dimension, group workload, and interval
  depth, selected by a frozen technical algorithm rather than distance value;
- at least 120 grouped-Persim cross-engine rows selected by a frozen technical
  algorithm across every fold, seed, view, and dimension;
- complete resource/storage reconciliation;
- byte-identical stage-1 repeat; and
- snapshot/resume/snapshot validation for all 75 completed group directories
  with unchanged hashes, byte sizes, and modification times.

Oracle selection may use only fold/seed/view/dimension, pair identity, group
size, and interval counts. It may not use tissue, approach, endpoint, outcome,
distance magnitude, or whether a result appears favorable.

## 11. Label firewall and public evidence

Production code and public/private result schemas accept no tissue, approach,
endpoint, outcome, class, biological label, technical label, cluster, ARI, NMI,
or fusion-selection field. Study ID is permitted only as an LOSO split
identifier. All outputs retain `outcome_label_state = closed`,
`biological_outcomes_computed = FALSE`, and zero fusion/clustering/outcome
counters.

Expression matrices, selected-cell IDs, typed views, distance matrices,
persistence diagrams/intervals, pair-level scientific distances, and complete
sample-level manifests remain under ignored private roots. Public evidence may
contain technical queue rows, aggregate interval/resource summaries, hashes,
counts, validation outcomes, projections, and decision records. It may not
contain source PDFs, reviewer text, binaries, private RDS files, or
`example_run.r`.

## 12. Abort, completion, and next authorization

Any input/implementation/Rust drift; label leakage; group-axis error; PH/MST,
oracle, repeat, or resume failure; nonfinite result; nonzero kernel status;
partial/stale artifact; or time, RSS, aggregate-work, or storage-cap breach
stops new launches and preserves already validated groups. Caps and scientific
rules cannot be relaxed after seeing production results within MV6-F.

MV6-F production completes only if all 75 groups and every final validation
gate pass. Completion authorizes only a new prediction-locked blocked-fusion
prefreeze using the immutable four component matrices. It does not itself
authorize fusion fitting, outcome access, clustering, G-MV6, new data, advanced
network fusion, public Rust adoption, binary distribution, package-default
changes, manuscript claims, releases, or pushes.
