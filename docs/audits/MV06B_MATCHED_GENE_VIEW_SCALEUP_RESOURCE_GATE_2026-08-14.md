# MV6-B matched gene-view scale-up and resource gate

| Field | Value |
|---|---|
| Date | 2026-08-14 |
| Prefreeze commit | `2e903f4` |
| Inventory implementation | `c7f48e9` |
| Pre-output schema correction | `2b2531b` |
| Contract | `mv06b_matched_gene_view_scaleup_resource_gate_v1` |
| Accepted private records | 75 MV5-D1 + 75 MV5-G |
| Outcome-label state | Closed |
| Decision | `stop_contract_revision_required` |

## Outcome

The accepted 90-sample blocked benchmark cannot currently support the strict
matched gene-view axis required for MV6 blocked fusion. This is a scientific
input-contract boundary, not a landscape failure and not evidence against gene
topology itself.

All 150 accepted resource-ledger and private cache identities verified. The
inventory then found:

| Representation | View instances | Exact 500-gene panel | Incomplete panel | Affected groups | Missing feature instances | Maximum missing in one view | Additional boundary |
|---|---:|---:|---:|---:|---:|---:|---|
| SCT fold | 6,750 | 6,679 | 71 | 31 / 75 | 111 | 4 | Variance remains unresolved for the other 379 held-out views |
| Inductive integration | 6,750 | 6,679 | 71 | 31 / 75 | 111 | 4 | Accepted artifacts define cell coordinates, not corrected gene expression, for all 6,750 views |

The 71 incomplete instances are held-out sample-seed views. They cannot use the
full training-selected panel. Mapping an absent feature to the training mean is
valid for the accepted cell-coordinate route because the coordinate contributes
zero to within-sample cell distances. It is invalid for strict gene topology
because the resulting gene point is constant and Pearson correlation-chord
geometry is undefined.

The other 379 held-out SCT instances have all 500 panel features, but their
per-gene within-sample variance was not certified by the retained D1 cell-view
payload. This variance question could be resolved from accepted D0 matrices;
doing so cannot repair the 71 already incomplete instances, so the gate did not
open or scan an additional 450 matrix caches.

## Why integrated gene topology is not already present

MV5-G stores 384-cell by 30-PC transfer-projection coordinates. Those are valid
`cell_topology_v1` inputs. They are not a gene-by-cell matrix and must not be
transposed or relabeled as gene topology.

MV5-G retains query active-feature identities for provenance, but it does not
retain a fold-safe corrected gene-expression assay for training and held-out
samples. Creating one would require a new definition and rerunning at least
part of the integration calculation. SCT expression also cannot be silently
called an integrated gene representation.

## Work that did not run

The strict candidate would require, per representation:

| Unit | Count |
|---|---:|
| Candidate gene-view instances | 6,750 |
| H0/H1 diagram components | 13,500 |
| Directed held-out-query/training-reference pairs | 35,350 |
| Dimension-specific landscape distances | 70,700 |
| Training-fitted cell/gene H0/H1 scales | 300 |
| Five-weight fusion pair rows | 176,750 |

Runtime is deliberately not extrapolated from an invalid input contract. Zero
gene views, PH jobs, landscapes, fusion rows, clusters, or biological outcomes
were computed.

## Validation

The independent validator does not call the production inventory finalizer. It
rehashes all production outputs and both accepted resource ledgers, then
reconstructs:

1. all 75 fold-seed axes and 6,300/450 training/query counts;
2. every SCT group missingness count from the accepted D1 ledger;
3. integrated group maxima from the accepted MV5-G ledger and the shared raw
   feature-availability provenance;
4. both representation totals, workload axes, and the prospective stop rule;
5. the zero-execution and public-safe schema boundaries.

All 12 categories pass. A clean repeat reverified all 150 accepted private
records and reproduced all six production outputs with identical bytes and
SHA-256 hashes. The temporary repeat directory was removed after comparison;
the accepted evidence remains in `docs/audits/mv06b-scaleup-evidence/`.

The focused MV6-B test context passes 10/10 expectations, and the complete
source-loaded suite passes with only the two existing intentional optional-Rust
and build-boundary skips. A local installed-package `--as-cran` check built and
installed successfully and passed the new MV6-B tests, but it is not recorded
as a package-check pass: 16 existing PH-subprocess/toy-baseline assertions fail
because the installed child process exits before writing a diagram. The same
functional tests pass source-loaded, and MV6-B changes neither launcher nor toy
baseline. This environment-specific installed-test issue remains a separate
publication-readiness follow-up.

## Contract choices requiring project-owner input

The gate intentionally does not select among these materially different
scientific estimands:

1. **Matched global-core dual view.** Freeze a label-closed gene panel that is
   present and nonconstant across the full existing-data axis, then rerun both
   cell and gene topology on that same panel. This preserves the cleanest
   cell-versus-gene comparison but changes the accepted cell benchmark,
   introduces explicitly transductive technical harmonization, and requires a
   new bounded resource profile before any large run.
2. **Complete-case fold-specific gene view.** Keep each training-derived panel
   and evaluate only held-out instances that support it. This changes the
   estimand, removes 71/450 held-out instances across 31 groups, and risks
   technology/study-dependent selection; it is unsuitable as the default
   without a confounding analysis and explicit missing-view estimand.
3. **Sample-specific active gene sets.** Drop unavailable/constant genes per
   sample. This preserves more samples but makes point identities and counts
   differ, confounding topology with feature availability. It is not
   recommended for primary fusion under the current framework.
4. **New inductive integrated-gene representation.** Define and retain a
   fold-safe corrected expression space before gene PH. This is substantial new
   methods work, separate from the accepted transfer-projection cell route.
5. **Close full fusion as not estimable for now.** Retain MV6-A as an exploratory
   pilot, report cell and pilot gene views separately, and proceed to an MV-07
   cell-focused robustness/confounding synthesis. This preserves the accepted
   blocked benchmark but does not support a full gene/fusion claim.

The recommended path depends on the intended paper claim. If full dual-view
fusion is central, option 1 is the most interpretable next research program and
must begin with a panel/estimand prefreeze plus bounded profiling. If preserving
the completed blocked benchmark and reaching a defensible manuscript sooner is
the priority, option 5 is the safer disposition. Options 2 and 3 should not be
primary; option 4 is best treated as a later methods extension.

## Gate disposition

G-MV6 remains open. MV6-05 blocked fusion evaluation and MV6-06 advanced fusion
remain closed. MV-07 cannot freeze a dual-view confirmatory specification until
the project owner selects the gene/fusion estimand or explicitly accepts a
cell-focused scope.

The dissertation-aligned landscape definition is unchanged: H0 and H1 remain
separate; finite positive-persistence intervals retain every active consecutive
level; integration remains exact or error-controlled; and no universal grid or
level cap is introduced.
