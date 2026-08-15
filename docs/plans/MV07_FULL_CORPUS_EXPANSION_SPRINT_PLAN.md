# MV7 full-corpus expansion sprint plan

Date: 2026-08-15
Status: active; MV7-D through MV7-F complete, exact 124-derived panel lock next

## Goal

Add the 34 retained single-study-tissue samples to a corrected descriptive
cell/gene topology analysis without changing the accepted 90-sample blocked
primary result, weakening the landscape definition, or turning descriptive
coverage into unsupported cross-study generalization.

## Non-negotiable boundaries

- The accepted 90-sample LOSO benchmark remains immutable and primary.
- The 124-sample analysis is a different, explicitly transductive descriptive
  estimand; it is not a larger replication of the cross-study endpoint.
- Cell and gene views remain separate. Equal-weight fusion remains a reported
  negative result and is not reopened.
- H0 and H1 remain separate and use all-active-level exact or error-controlled
  landscapes with no universal grid or level cap.
- The 34 samples may enter only after a label-closed contract is durable.
- The three below-250 samples and the separate 25-sample validation cohort do
  not enter the 124-sample analysis.
- No new dataset is needed for this sequence.

## MV7-D — corpus reconciliation and source/SCT feasibility (complete)

Deliverables:

- 127-row sample reconciliation and estimand flags;
- eight-tissue/study summary and source/artifact coverage;
- explicit three-sample pre-PH exclusion record;
- explicit 34-sample single-study descriptive stratum;
- six prospectively selected depth-extreme source/SCT sentinels;
- 6/6 feasibility and 8/8 independent validation;
- approach-label conflict register; and
- unchanged landscape-contract ledger.

Gate result: upstream extension is feasible; expanded caches, PH, landscapes,
and outcomes remain closed.

## MV7-E — metadata provenance and descriptive estimand prefreeze (complete)

Purpose: make the 124-sample cell/gene objects mathematically comparable before
building them.

Tasks:

1. Resolve the 16 approach-label disagreements using accession/run metadata and
   original study methods. Preserve unresolved rows explicitly; do not infer
   technology from directory names or expression patterns.
2. Freeze five existing seeds (`20260805`–`20260809`) and 384 cells per sample.
3. Freeze a global descriptive cell transform fitted on all 124 samples. Mark
   it transductive and prohibit comparison with the training-only LOSO transform
   as if they were the same estimand.
4. Test the accepted 90-derived 500-gene panel for complete availability in all
   34 added samples without reading topology or labels. If complete, retain it
   for direct representation comparability. If incomplete, use the same
   deterministic global-core algorithm over all 124 and record panel overlap;
   the fallback must be selected by availability only.
5. Freeze `cell_topology_v1` as 384 cells in 30 global descriptive PCs with
   Euclidean geometry and `gene_topology_v1` as the same 500 genes with
   correlation-chord geometry across the same cells.
6. Freeze complete Vietoris–Rips H0/H1, pair identities, expected cardinalities,
   resource caps, atomic/resume rules, and a label firewall.

Acceptance: exact sample/seed/panel/fit-scope manifests; no PH; independent
validation; owner input only if the panel-availability fallback or provenance
evidence presents a genuinely different scientific choice.

Gate result: all 16 approach conflicts resolve to the public/`Approach.y`
value. The primary 90 are all scRNA-seq, so the previous mixed-approach
diagnostic is superseded as not estimable. Thirty-three added samples contain
the full accepted panel; `SRA701877_SRS3279688` lacks one feature (`KLF2`), so
the prospectively specified 124-sample global-core fallback is mandatory.
MV7-F is authorized only to create the 34 raw and 170 SCT caches. The exact
124-derived panel must be frozen in a post-MV7-F evidence commit before MV7-G.

## MV7-F — omitted-34 upstream production (complete)

Purpose: extend only the already validated raw/SCT layer.

Frozen workload:

- 34 atomic raw shards;
- 170 SCT caches (34 samples × five seeds);
- 170 deterministic 384-cell selection identities; and
- no fitted coordinates, PH, landscapes, labels, or outcomes.

Planning baseline from MV7-D: approximately 1.158 worker-hours and 1.402 GB
under simple linear scaling. Apply a prospective safety factor and reproject
from measured production before MV7-G.

Acceptance: 34/34 and 170/170 independent identity/content validation, one
byte-exact representative repeat, zero partial state, immutable resume, and an
updated measured resource projection.

Gate result: all 34 raw shards, 170 SCT caches, and 170 deterministic
selections pass. The accepted run used 6,952.050 worker-seconds, peaked at
1,697,538,048 bytes RSS, and stored 1,618,456,042 private bytes. Independent
validation passed 9/9 categories; the clean representative raw and SCT caches
were byte-identical; 204/204 cache hashes, sizes, and mtimes passed the amended
immutable gate. Two strict-comparison false negatives were corrected and
retained in separate amendment evidence with zero cache mismatch or mutation.
Authorize only the exact 124-derived global-core panel lock before MV7-G.

## MV7-G — typed-view and PH sentinel gate

Purpose: validate the newly frozen global descriptive transform and both
topology orientations before full calculation.

Use the six MV7-D depth-extreme samples across all five seeds. Build typed cell
and gene views, run complete VR H0/H1, verify H0 deaths against MSTs, cross-check
a bounded H0/H1 subset with an independent engine, repeat representative
artifacts byte-identically, and profile time/RSS/storage. Do not compute pairwise
landscape outcomes yet.

Acceptance: 60 sample-seed-view diagram records pass identity, geometry,
orientation, MST, cross-engine, repeat, and resource gates. The post-sentinel
projection must authorize or reject full PH explicitly.

## MV7-H — full 124-sample dual-view topology and landscapes

Purpose: construct the complete label-closed descriptive geometry.

Expected minimum scope under five seeds:

- 620 sample-seed states per view and 1,240 cell/gene PH diagrams total;
- 7,626 unordered sample pairs per seed;
- 38,130 biological pairs per view across five seeds;
- 76,260 view-specific pairs; and
- 152,520 H0/H1 component-distance rows.

Run PH and landscape production as separate resumable stages. Use the accepted
grouped/Rust landscape kernel only with its hash/runtime gates; keep R as the
canonical oracle. Stream pair outputs and never save dense landscape-function
matrices. Require independent counts, numerical oracles, deterministic repeats,
and immutable resume before labels are opened.

## MV7-I — descriptive topology and clustering evaluation

Purpose: use all eight tissues without overstating what the three single-study
tissues can establish.

Before label access, freeze:

- view-specific H0 and H1 matrices plus the already defined unweighted
  H0/H1 composite as descriptive only;
- seed aggregation and uncertainty summaries;
- PAM as the primary label-free clustering procedure and average linkage as
  the existing secondary sensitivity;
- a prespecified cluster-count stability rule independent of tissue labels;
- tissue/study/approach descriptive joins and the approach-provenance policy;
- complete reporting of every tissue and seed; and
- explicit prohibition of causal technology, external-generalization, or
  single-study cross-study claims.

The analysis may describe whether the added tissues occupy distinct topology
regions, how cell and gene views differ, and how cluster stability changes when
the 34 samples are included. It must not select a favorable view, tissue,
dimension, seed, or clustering method.

## MV7-J — manuscript claim map, literature update, and figures

Only after MV7-I should the result-dependent claim map become final. Link every
proposed sentence, table, and figure to its estimand, sample population, code,
artifact, uncertainty, limitation, and literature support. Regenerate the
saved persistence-diagram/landscape explanatory figure from corrected artifacts
and clearly distinguish primary 90-sample blocked findings from 124-sample
descriptive findings.

Run a current literature/reference audit in this stage and decide whether an
external dataset is necessary for a named generalization claim. New data remain
deferred unless that explicit evidence-gap gate passes.

## Optional MV7-K — below-250 threshold sensitivity

This is not part of the main sequence. If manuscript review shows that the
three exclusions materially affect a stated claim, prefreeze a common reduced
depth no greater than the smallest retained count, multiple deterministic
seeds, and a substantia-nigra-only stability/failure comparison. Treat it as a
threshold and computational-stability sensitivity, never as an addition to the
384-cell primary or descriptive cohort.

## Audit checkpoints

Each sprint requires a durable prospective specification, source hashes,
write-once public evidence, private atomic artifacts, an independent validator,
representative repeat, resume verification, measured resource reconciliation,
and an explicit authorization for the next stage. Null and failed results stay
in the record. Dr. Rouchka and Dr. Mistry remain in the credit plan; authorship
order and CRediT roles remain an author-team decision before manuscript/release.
