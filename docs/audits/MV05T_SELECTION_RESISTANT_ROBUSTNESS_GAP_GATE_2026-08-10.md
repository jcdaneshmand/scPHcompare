# MV5-T selection-resistant robustness/gap gate

Date: 2026-08-10

Accepted outcome base: `64d84d9`

Prospective gate commit: `7f6784d`

New robustness outcomes: zero

## 1. Decision

MV5-T approves a bounded resource admission for three outcome-independent
robustness families and rejects immediate full execution. The admitted families
are nested cell count, 20-PC truncation, and cosine-chord point geometry. They
produce four one-factor-at-a-time configurations. The next queue contains 24
resource-only groups: three fixed folds by two representations by four
configurations at seed `20260805`.

No MV5-S value chose these families. They were already named by the dual-view
contract and MV5-M robustness plan, target distinct alternative explanations,
and can reuse accepted coordinates without refitting or leakage. Labels remain
closed; no new diagram, distance, retrieval endpoint, cluster, or outcome was
computed in this gate.

## 2. Candidate audit

| Candidate family | Weighted score | Gate disposition |
|---|---:|---|
| Nested 192/256 cell counts | 55 | admit bounded resource pilot |
| 20-PC truncation | 53 | admit bounded resource pilot |
| Cosine-chord geometry | 53 | admit bounded resource pilot |
| Existing H0/H1/distance panel | 48 | already complete; no new execution |
| 250/1,000-feature panels | 42 | defer; requires training-only refit |
| 50-PC refit | 40 | blocked; accepted caches contain 30 PCs |
| Additional cell seeds | 36 | reject; five seeds already complete |
| Clustering/k expansion | 34 | reject/defer post-outcome search; current roles complete |
| Small-study eligibility rerun | 30 | defer; requires complete fold rerun |
| Integration-method expansion | 40 | defer; requires native inductive gate |

Scores use frozen weights for scientific value, identifiability, outcome-
selection safety, artifact readiness, resources, and generalized reviewer
relevance. Confidential reviewer wording is neither reproduced nor tracked.

## 3. Admitted configurations

| Configuration | Cells | Coordinates | Point geometry |
|---|---:|---:|---|
| `nested_cells_192_pc30_euclidean_v1` | 192 | 30 | Euclidean |
| `nested_cells_256_pc30_euclidean_v1` | 256 | 30 | Euclidean |
| `cells384_pc20_euclidean_v1` | 384 | 20 | Euclidean |
| `cells384_pc30_cosine_chord_v1` | 384 | 30 | row-unit-normalized Euclidean chord |

Every row differs from the accepted 384-cell/30-PC/Euclidean baseline by one
factor only. Nested cell order hashes sample ID, seed, and cell ID and is shared
across representations. The 20-PC condition truncates the existing training-fit
decomposition; it does not refit PCA. Zero-norm cells are a hard failure for
cosine chord.

## 4. Source and coordinate readiness

The source freeze contains 164 identities: 14 committed contracts,
implementations, manifests, audits, and resource records plus SHA-256/size
identities for all 75 accepted SCT and 75 accepted integrated private coordinate
files. Private contents and absolute paths are not published or tracked.

The three deterministic admission folds are:

- minimum training: `large_loso_v1:SRA779509`, n=65;
- median training with lexical tie break: `large_loso_v1:SRA703206`, n=86;
- maximum training with lexical tie break: `large_loso_v1:SRA713577`, n=89.

At seed `20260805`, all 270 selected sample pairs pass exact 384-by-30 shape,
sample-axis, cell-ID, finite-value, first-20-coordinate, and nested 192-in-256
checks. Across both representations, zero selected cells have zero coordinate
norm.

## 5. Resource gate

Accepted measurements give 7,720.432 worker-seconds of PH and 6,267.045
worker-seconds of landscape-distance work per complete two-representation
setting, or 13,987.477 seconds (3.885 hours). Four full settings project to
15.542 worker-hours before assembly/validation. Conservative repeated staging
projects to 10.18 GB.

Full execution is therefore not authorized. The 24-group admission is limited
to one worker, two worker-hours, 600 seconds and 4 GiB RSS per group, and 2 GiB
new private storage. It must measure whether streaming/shared staging can make a
later full retrieval-sensitivity queue safe.

## 6. Validation and reproducibility

Focused prospective tests pass 16/16; the pre-gate repository suite passes 847
checks with zero failures/errors/warnings and two expected CRAN-only skips. The
builder freezes and validates all source axes, candidate/configuration roles,
the exact 24-unit unopened queue, 12 required validation categories, and ten
abort rules.

A clean gate regeneration reproduces all 10 public CSV artifacts byte-for-byte.
Every queue row reports `labels_opened=FALSE`, `outcomes_computed=FALSE`, and
`admission_executed=FALSE`.

## 7. Decision table

| Question | Disposition |
|---|---|
| Scientific contract coherent? | approve |
| Correctness demonstrated? | pass for gate construction and source readiness |
| Computation feasible? | bounded admission only; full grid not yet authorized |
| Biological interpretation permitted? | prohibited in MV5-T |
| Next action | execute only the 24-group MV5-U resource admission |

MV5-T stops before robustness production or outcomes. The next sprint may
construct and execute only the admitted coordinate/PH/landscape/energy resource
pilot with labels closed. Full robustness, method ranking, spectral promotion,
gene topology, fusion, new data, optimization/Rust, package-default changes,
manuscript claims, PDF/reviewer/private-cache tracking, `example_run.r`, and
pushing remain closed.
