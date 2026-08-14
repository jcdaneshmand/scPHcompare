# MV-02 orientation-safe constructors and analytical fixtures — 2026-08-05

## Outcome

G-MV2 passes as a technical implementation gate. The package now has a non-default typed route that makes the scientific object, observation axis, coordinate axis, metric, identifiers, transformations, and cache identity explicit before persistent homology can run.

This sprint did not process biological data, generate eligible pilot diagrams, activate corrected behavior in the historical pipeline, change the landscape definition, or authorize manuscript claims.

## Implemented API

`R/dual_view_topology.R` adds:

| Function | Contract role |
|---|---|
| `select_matched_cells()` | Deterministic RNG-isolated sampling after canonical cell-ID sorting |
| `new_dual_view_source()` | Typed standardized `genes_by_cells` source with immutable axis, sample, fit-scope, representation, seed, and digest identity |
| `fit_cell_topology_pca()` | Deterministic shared PCA in conventional cells-by-genes orientation, canonicalized by sample ID |
| `construct_cell_topology_view()` | Named cells-by-shared-PC point cloud under `cell_topology_v1` |
| `construct_gene_topology_view()` | Named explicit `dist` object under Pearson correlation-chord geometry for `gene_topology_v1` |
| `validate_topology_view()` | Class, contract, axes, identifiers, payload, provenance, and cache validation |
| `run_topology_view_ph()` | Corrected PH dispatch accepting only typed cell/gene view objects |
| `run_legacy_matrix_ph()` | Explicitly acknowledged, warned, scientifically ineligible historical rows-as-points compatibility route |

All functions are documented and exported, but none is wired into `process_datasets_PH()` as a default.

## Contract profiles

The `scientific` profile admits only:

- 500 named genes;
- 384 named cells;
- a 30-PC cell-view contract;
- `scientific_eligible = TRUE` at the source/view boundary.

Reduced tests must explicitly request `analytical_fixture`, provide their expected dimensions, and are permanently stamped `scientific_eligible = FALSE`. A reduced object therefore cannot be mistaken for an eligible MV-03 artifact even if its numerical calculation succeeds.

`new_dual_view_source()` accepts an already standardized matrix rather than silently extracting or transforming a Seurat assay. Representation-specific SCT/integrated extraction, feature-panel fitting, and fit-scope standardization remain explicit MV-03 inputs. The constructor rejects anonymous, transposed, missing-axis, duplicated-ID, nonfinite, wrong-shape, and zero-variance-gene inputs.

## Typed view identities

Every corrected view records:

- `view_id` and contract version/profile;
- scientific-eligibility flag;
- source cache key;
- sample, cohort, representation, fit scope, and subsample seed;
- input, point, and coordinate axis roles;
- ordered point and coordinate IDs;
- point-metric ID and transformations;
- payload SHA-256 and a contract-prefixed cache key.

The cache identity is recomputed during validation. Mutating the matrix, PCA basis, point metric, axis identity, identifiers, transformations, or payload without constructing a new object fails as stale rather than silently reusing a key.

Corrected result keys begin with `corrected_topology_result_v1:`. Historical compatibility keys begin with `legacy_topology_result_v0:`. The two artifact families cannot collide.

## Cell-view implementation

Compatible typed sources are sorted by sample ID before pooling, so caller list order does not alter the fitted model. They must share cohort, representation, fit scope, subsample seed, ordered genes, and contract. Because every source has the same cell count, each fit sample contributes equally.

`stats::prcomp()` is run over pooled cells by genes with `center = FALSE` and `scale. = FALSE` because fitted standardization is an upstream source contract. The retained variance-weighted scores are projected with the same named loadings for every sample. Fewer than the requested usable PCs, incompatible genes, stale PCA identity, nonfinite scores, or zero-norm cell scores are hard failures.

Numerically coincident cells are retained because they can occur in a real point cloud, but their count is recorded in diagnostics rather than hidden.

## Gene-view implementation

Each gene is centered over its matched cells and divided by its L2 norm. The explicit distance object is ordinary Euclidean distance between those unit vectors:

`||u_g - u_h|| = sqrt(2 * (1 - correlation(g, h)))`.

The fixture verifies symmetry, non-negativity, zero diagonal, the `[0, 2]` bound, and every triangle inequality. Joint cell-column permutation leaves the distance matrix unchanged. Because ordered cell identity remains provenance, the permuted source receives a different cache key even though its invariant gene geometry agrees numerically.

Numerically duplicated nonconstant genes become coincident points and are diagnosed. Constant or effectively zero-variance genes fail at the typed source boundary.

## Corrected and legacy PH boundary

`run_topology_view_ph()` refuses bare matrices, bare distance objects, unrecognized classes, stale payloads, and altered identities. It dispatches:

- the cell view as its explicit named cells-by-PC point cloud;
- the gene view as its explicit named `dist` object.

The v1 entry point enforces `max_dim = 1` and complete filtration with `threshold = -1`, then records the ripserr version, coefficient field, point count, axes, metric, view key, and diagram digest. Because ripserr 0.3.0 omits the essential H0 class from its returned matrix, the corrected route adds exactly one `(dimension = 0, birth = 0, death = Inf)` interval when absent and records whether it did so. Provenance also records finite, infinite, zero-persistence, and invalid interval counts. Landscape evaluation remains responsible for excluding the infinite interval as specified.

The existing `process_datasets_PH()` calculation remains the historical implementation. It now warns unless the caller explicitly acknowledges the legacy orientation, and its returned provenance is stamped `topology_contract_id = legacy_gene_view_v0` and `scientific_eligible = FALSE`. This warning/stamp does not repair or reclassify historical artifacts.

## Analytical fixtures

The reduced 4-gene by 5-cell fixture is intentionally asymmetric in its axes:

| View | Points | Coordinates/input | Observed finite H0 deaths |
|---|---:|---|---:|
| `cell_topology_v1` | 5 cells | 2 shared PCs | 4 |
| `gene_topology_v1` | 4 genes | explicit correlation-chord `dist` | 3 |

Thus a transpose or axis substitution cannot pass while yielding an apparently equivalent typed result. Both fixture results are scientifically ineligible.

The test file also covers deterministic and unequal cell subsampling, transposition, anonymous/duplicated IDs, missing/reordered genes, constant genes, duplicated gene/cell values, mismatched seeds and PCA models, cell-permutation invariance, metric properties, stable repeats, cache invalidation, bare-input rejection, full-shape enforcement, and legacy/corrected separation.

## Validation

- Focused MV-02 test file: 12 test blocks passed.
- Complete package suite: all 11 test files passed.
- Source-package build/check under locked Ubuntu R 4.4.1: `Status: OK`, 0 errors, 0 warnings, 0 notes.
- Locked ripserr version used by the fixture: 0.3.0.

An initial outer diagnostic guard returned a nonzero shell status after `R CMD check` had completed because it treated the returned check object fields as simple vectors. The authoritative check log itself completed successfully. A subsequent clean wrapper invocation also completed with exit status 0 and `Status: OK`; no package-check error, warning, or note occurred.

## G-MV2 disposition

| Criterion | Result |
|---|---|
| Bare assay matrices excluded from corrected PH | Pass |
| Axis/transposition/identifier failures explicit | Pass |
| Cell and gene analytical diagrams distinct | Pass |
| Mathematical gene metric verified | Pass |
| Stable hashes and mutation invalidation | Pass |
| Legacy/corrected cache families separated | Pass |
| Full tests and source-package check | Pass |
| Biological eligibility or conclusions | Not evaluated; prohibited in MV-02 |

MV-03 is next. It must implement representation-native extraction and fit-scope standardization, reconcile fresh QC against the frozen pilot manifest, construct scientifically eligible full-shape sources, add monitored/resource-capped typed PH execution, and stop at the prespecified feasibility boundaries.
