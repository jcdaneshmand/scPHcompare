# MV5-C one-tissue existing-data feasibility audit

| Field | Value |
|---|---|
| Date | 2026-08-06 |
| Execution contract | `mv05c_existing_data_job_v1` |
| Landscape contract | `mv05c_full_l2_exact_all_active_levels_v1` |
| Pilot | Prospectively selected six-sample, two-study existing-data tissue |
| Seeds | `20260805`–`20260809` |
| Outcome-label state | Closed |
| Biological/technical endpoints | Not computed |
| Gate | MV5-C feasibility complete; G-MV5 remains open |

## Outcome

MV5-C completed the first real-data, all-seed execution of the corrected cell
topology and the feasible portion of gene topology. The prospective selection
rule—at least two studies, every sample at least 384 cells, minimum eligible
sample count, lexical tie-break—selected liver before method execution. The
closed manifests contain sample IDs, study blocks, cell selections, folds,
seeds, methods, and fit scopes but no tissue or approach column.

All ten canonical fold/seed jobs exited successfully under the 30-minute/8-GiB
cap. The jobs used 500 training-ranked genes, 384 matched cells per sample, 30
training-fitted PCs, complete Vietoris–Rips H0/H1, and no transductive
integration fallback. The installed Seurat native SCT fallback was retained;
`glmGamPoi` was not installed and the dependency lock was not changed.

## View dispositions before labels

| Representation/view | Completed jobs | Structured unavailable jobs | Technical disposition |
|---|---:|---:|---|
| SCT cell topology | 10 | 0 | Feasible |
| SCT gene topology | 5 | 5 | Conditional on fold training diversity |
| Inductive integrated cell topology | 10 | 0 | Feasible |
| Inductive integrated gene topology | 0 | 10 | No held-out corrected gene space |

The five SCT-gene failures all occurred in the fold trained on one sample and
querying five samples. A training-only 500-gene panel contained genes that were
constant in at least one held-out sample. The correction separates cell
projection eligibility from gene-topology eligibility: constant query genes do
not invalidate training-fitted PCA projection, but they do invalidate Pearson
gene geometry. The implementation does not inspect held-out variance to alter
the training panel.

## Completed immutable artifacts

- 10 canonical job bundles and 4 retained debugging attempts, all hashed;
- 150 eligible persistence diagrams in 25 six-sample strata;
- 750 finite exact landscape distances: H0/H1 separate, every active level,
  zero numerical error estimate;
- 525 matched-baseline pair rows;
- 85 complete six-sample distance matrices;
- 30 completed held-out mappings with 39–132 anchors;
- 425 label-free held-out neighbor-ranking rows;
- 17 complete five-seed PAM stability groups with `k=2:5`, one-SE selection,
  and average-linkage sensitivity partitions; and
- no tissue/approach endpoint join, method winner, fusion, or outcome score.

The one-SE procedure selected `k=2` for 14 groups, `k=3` for two, and `k=4` for
one. These counts describe label-free partition stability only and have no
biological interpretation.

## Resources and MV5-D gate

The ten jobs used 4,760.71 worker-seconds (1.322 hours), with median 484.685
seconds, maximum 556.26 seconds, and peak RSS 2,950,356,992 bytes. Exact
landscapes completed in 136.41 internal seconds (156.37 measured wall seconds)
at 354,390,016 bytes peak RSS.

A deliberately simple lower-fidelity scaling projection for the existing
90-sample/15-study candidate set estimates:

- 75 fold/seed jobs;
- about 7,270 seconds per job under sample-linear scaling;
- about 151.46 worker-hours for representation/PH jobs;
- 1,802,250 landscape pair rows under full within-fold matrices; and
- about 91.05 additional pair-linear landscape hours.

The projection fails both current resource gates. It is not a promise of exact
runtime; it is sufficient evidence that naïve MV5-D execution is prohibited.
Normalization/model reuse, chunking, streaming, and a scientifically justified
query-to-training pair scope must be specified and re-profiled first.

## Verification and boundaries

- Focused MV5 execution tests pass 42/42 expectations.
- New cell-projection sources are immutable, permit held-out constants only for
  cell/pseudobulk use, and are rejected by gene topology.
- Exact landscape rows are finite and complete; 750/750 are marked exact and
  all-active-level.
- Private raw counts, source PDFs, reviewer material, private job bundles, and
  distance bundles remain ignored. `example_run.r` remains untracked.
- No push, public-default change, new data, fusion, or endpoint evaluation was
  performed.

## Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | Approve cell; revise/narrow gene feasibility |
| Correctness demonstrated? | Technical pass for completed strata |
| Computation feasible? | One tissue yes; full run no under current architecture |
| Biological interpretation permitted? | Prohibited until endpoint gate |
| Next action | MV5-C2/MV5-D resource-safe execution specification |
