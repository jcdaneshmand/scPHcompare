# MV5-M post-retrieval benchmark-gap and gate audit specification v1

## Document control

| Field | Value |
|---|---|
| Contract ID | `mv05m_benchmark_gap_gate_v1` |
| Date | 2026-08-10 |
| Status | Executed and independently validated |
| Source revision | `a3a98895d47b99f8a12691b10a2320b82d195433` |
| Decision type | no-outcome benchmark-axis prioritization |
| Selected next sprint | MV5-N clustering contract and complete-matrix resource gate only |

## 1. Purpose and boundary

MV5-M determines what should follow the completed MV5-E/K/L biological
retrieval sequence. It inventories unresolved scientific axes, inspects actual
code and artifact readiness, applies a fixed weighted gate, and authorizes one
bounded next sprint.

MV5-M does not compute a biological or technical outcome. It does not inspect
MV5-L tissue-specific performance to select a tissue or method. It does not run
clustering, technical-mixing, robustness, integration, gene, fusion, new-data,
or optimization work. Confidential reviewer wording and the local PDFs remain
unpublished.

## 2. Evidence basis

The audit uses:

- the dissertation aims and future directions already mapped in
  `docs/PROJECT_EVIDENCE.md` (notably dissertation PDF pages 18-20, 46-65,
  82-84, and 127-133);
- the preprint's integration and sample-clustering aims already summarized in
  the same ledger;
- generalized reviewer workstreams from the Git-ignored response matrix,
  without quoting or identifying confidential report text;
- the frozen MV-05 statistical benchmark hierarchy;
- accepted MV5-E, MV5-K, and MV5-L audits;
- current production code and pair-scope flags; and
- measured distance-engine and resource evidence.

The external methodological framing remains consistent with scIB's separation
of biological conservation and batch removal
([Luecken et al., 2022](https://doi.org/10.1038/s41592-021-01336-8)) and with
evidence that integration can remove biological variation
([Zhang et al., 2025](https://doi.org/10.1038/s41587-024-02463-1)). These
references motivate separate axes; they do not determine a result.

## 3. Audited axes

Nine axes are assessed independently:

1. biological-conservation retrieval;
2. technical/batch mixing;
3. label-free sample clustering;
4. robustness/sensitivity;
5. integration-method expansion;
6. gene topology;
7. cell/gene fusion;
8. external validation/new data; and
9. optimization/Rust.

For each axis the audit records its estimand, primary unit, current evidence,
dominant validity risk, missing computation, and disposition. A completed axis
is not rescored as a candidate next sprint.

## 4. Frozen selection criteria

Every candidate receives an integer score from 0 to 4 on six criteria. Higher
is better. The weights are frozen before scoring:

| Criterion | Weight |
|---|---:|
| Scientific value | 3 |
| Reviewer relevance | 2 |
| Identifiability/validity | 3 |
| Existing-artifact readiness | 2 |
| Resource feasibility | 1 |
| Protection from outcome-driven selection | 2 |

Candidates failing the validity gate cannot be selected regardless of total.
Eligible candidates rank by weighted total, then reviewer relevance, artifact
readiness, and canonical candidate ID.

## 5. Axis dispositions

### Biological conservation retrieval

Complete for the frozen existing-data endpoint. MV5-E and MV5-K separately
showed no H0/H1 advantage over matched energy; MV5-L showed no favorable
representation difference-in-differences. This axis remains complete with its
null/negative result and cannot be retuned.

### Technical mixing

Scientifically essential but currently blocked at the sample-level primary
design. Every eligible study is nested in one tissue, so study and tissue are
strongly confounded. The accepted LOSO query-to-training pairs exclude the
held-out study by definition and therefore contain no same-study comparator
for the held-out query. Training-only or cell-neighborhood diagnostics may be
possible, but they require a new identifiability contract and must remain
separate from biological conservation.

### Label-free sample clustering

High scientific and reviewer relevance and central to the dissertation-era
project. The current production landscape artifacts explicitly support only
query-to-training retrieval, not training distance matrices. Nevertheless all
90 per-fold cell views and persistence diagrams exist, the stable-k helper is
implemented, and query-to-training distances needed for held-out assignment
already exist. A leakage-safe training-cluster/held-out-prediction design is
therefore plausible but must be frozen and resource-gated first.

### Robustness

High validity and the second-ranked next axis. The current five seeds assess
cell-subsample repetition, but cell count, feature count, coordinate dimension,
and distance sensitivity are not a frozen full-data grid. Because primary
outcomes are known, future robustness is sensitivity evidence only and cannot
rescue or replace the null result.

### Remaining axes

Additional integration methods are deferred until the single integration
baseline's full benchmark structure is complete. Gene topology is blocked by
fold-specific eligibility and the MV5-C integrated-gene failure. Fusion is
blocked until both views pass independently. External validation is important
but remains deferred under the existing-data-first rule. Optimization/Rust is
deferred because measured current stages pass resource gates and the existing
Rust decision is negative.

## 6. Deterministic decision

| Candidate | Weighted score | Validity gate | Rank/disposition |
|---|---:|---|---|
| Label-free clustering contract/resource gate | 45 | Pass | Selected |
| Retrieval robustness/sensitivity | 43 | Pass | 2 |
| External validation/new data | 40 | Pass | 3; deferred existing-data-first |
| Technical mixing evaluation | 38 | **Fail** | blocked by current identifiability |
| Integration-method expansion | 37 | Pass | 4 |
| Optimization/Rust | 35 | Pass | 5; no bottleneck trigger |
| Gene topology | 30 | Pass | 6; incomplete view readiness |
| Cell/gene fusion | 18 | **Fail** | blocked components |

MV5-N is selected. The score is a prioritization device, not a performance
metric or biological result.

## 7. Exact clustering pair scope

An inductive clustering design does not need distances among held-out queries.
It needs training-training distances to fit clusters and the already accepted
query-to-training distances to assign held-out samples.

Across 15 folds and five seeds:

- existing query-to-training pairs: 35,350 per representation;
- missing training-training pairs: 262,675 per representation; and
- missing H0/H1 component rows: 525,350 per representation.

The exact counts were independently reconstructed from the accepted MV5-E fold
geometry. Linear projection from accepted component-row throughput gives
landscape-only lower bounds of approximately 8.655 worker-hours for SCT and
4.280 worker-hours for integrated. These are not execution forecasts: they
exclude baseline construction, validation, I/O, repeat, clustering, and failed
or slow-tail groups.

## 8. MV5-N authorization

MV5-N may only:

1. freeze training-only PAM/k-medoids as the primary clustering method;
2. freeze candidate `k=2:10` and select the smallest k within one Monte Carlo
   SE of maximum five-seed adjusted-Rand stability, without labels;
3. define held-out assignment to frozen training medoids from immutable
   query-to-training distances;
4. specify a fair average-linkage sensitivity and keep current spectral code
   ineligible pending its own implementation gate;
5. inventory exact training-pair requests and cache reuse;
6. profile bounded label-closed representative/minimum/maximum folds under
   time, memory, and storage caps; and
7. validate synthetic clustering/assignment behavior and deterministic resume.

MV5-N may not run full training-matrix production, open tissue/study/approach
labels for clustering outcomes, select methods/tissues, report ARI/NMI, promote
spectral clustering, or begin gene/fusion/new-data/optimization work.

## 9. Acceptance

MV5-M passes when the scoring scale, scores, dispositions, fold pair counts,
confidential-source boundary, no-outcome counters, and selected next sprint are
independently reproduced; all seven production artifacts repeat byte for byte;
focused/full tests pass; package check is clean; and the public decision record
is committed locally without private artifacts.
