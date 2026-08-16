# MV8-D HCA structural/QC closure

## Decision

MV8-D is complete with disposition
`blocked_exact500_annotation_incompatibility`. This is a correct prospective
gate result, not a computation failure and not a biological result.

The exact eight HCA files are structurally valid, every biological unit has far
more than the required 384 post-QC cells, and all five accepted 124-sample
reference transforms remain cryptographically and internally valid. The frozen
500-gene external-projection contract nevertheless fails because only 475 panel
genes have an exact Ensembl-stable-ID match in the HCA feature reference. The
same 25 genes are absent from all eight files. No PCA coordinate, PH diagram,
landscape, topology distance, clustering result, biological endpoint, or
manuscript claim was calculated.

## What was acquired

The project owner authorized exactly the MV8-C primary manifest. Eight files
were downloaded sequentially into the private cache, verified before opening,
and atomically published. Their combined observed size is exactly 202,770,089
bytes, and all eight observed SHA-256 values equal the frozen values. There are
no substituted donors, extra files, or partial downloads.

The raw H5 files remain outside the Git worktree. `.gitignore` now also excludes
H5, H5AD, and H5Seurat payloads. Public evidence contains no expression value,
cell barcode, selected-cell list, local absolute path, or transient signed URL.

## Structural and QC result

Every file exposes a valid feature-by-barcode CSC matrix with:

- 33,538 Gene Expression features;
- 737,280 raw barcode columns;
- nonnegative integer counts;
- unique barcode and Ensembl axes;
- `GRCh38` genome metadata;
- `Single Cell 3' v2` chemistry metadata; and
- the same feature-identifier axis across all eight donors.

The prospective legacy-comparable QC was applied exactly as frozen: 200-feature
entry, three-cell feature entry, 500--9,000 detected genes, mitochondrial
percentage at most 20, ribosomal percentage greater than 5, and final detection
in more than three retained cells. Post-QC cell counts range from 3,403 to
4,707, with median 4,165. All eight donors therefore pass the 384-cell gate.

This entry rule is explicitly a historical deterministic heuristic, not a claim
of modern statistical droplet calling. A later modern-QC sensitivity would need
its own prospective contract.

## Why the panel gate blocks

The exact stable-ID intersection is 475/500, or 95%. Every one of those 475
genes survives final QC in every donor. The remaining 25 stable IDs are absent
from all eight H5 feature references, so this is a shared annotation/reference
mismatch rather than stochastic donor-level expression loss. Only one missing
gene has an exact symbol-only candidate, and its Ensembl ID differs; no symbol
rescue was applied.

The H5 metadata identify `GRCh38` but not the precise Ensembl/GENCODE annotation
release. Official 10x documentation records that Cell Ranger reference releases
change gene IDs and gene names and remove some annotations through filtering.
That makes annotation-release incompatibility the evidence-supported
explanation, although the exact HCA build provenance is not encoded in these
files. See the [10x reference release notes](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes)
and [10x reference build guidance](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references).

Zero filling, symbol-only substitution, or silently dropping 25 input rows from
the existing 500-gene transform would change the scientific object after
inspection. All are rejected.

## Reference-transform result

All five accepted seed-specific reference bundles pass repeat validation:

- 124 fit samples per seed;
- the exact ordered 500-gene panel;
- 500-element center and scale axes;
- 500 x 30 PCA rotations;
- source-file, record, panel, and PCA-model cache identities unchanged; and
- no HCA expression access and no reference refit during validation.

The reference itself is therefore sound. It simply cannot be applied as the
exact frozen 500-gene transform to a query with 25 unavailable genes without an
unapproved imputation or object change.

## Recovery recommendation requiring owner approval

The recommended recovery is `common475_reference_only_refit`:

1. Freeze the exact ordered 475 stable-ID intersection already established by
   this audit; do not rank or select genes using outcomes or HCA expression.
2. Rebuild five seed-specific transforms using only the accepted 124-sample
   reference SCT caches. HCA values must not contribute to reference centering,
   scaling, PCA, feature ranking, or tuning.
3. Normalize each HCA donor/seed under the already accepted SCTransform policy,
   apply the frozen common panel, and project into the corresponding new
   reference-only 30-PC model.
4. Reconstruct both cell and gene topology for reference and query samples on
   the same 475-gene object. The accepted 500-gene PH artifacts cannot be mixed
   with the new object.
5. Recompute 1,240 reference PH jobs plus 80 HCA PH jobs. Prior MV7-H evidence
   makes this workload feasible but does not eliminate the need for a new
   resource prefreeze and exact fallback receipts.
6. Compute only query-to-reference distances for the external endpoint. H0 and
   H1 remain separate; landscapes use every consecutive active level, exclude
   essential H0, and use exact or error-controlled squared-L2 with no primary
   fixed grid or level cap.
7. Keep labels closed until normalization, coordinates, PH, landscapes,
   distances, and prediction locks validate.

This result would be an external *harmonized-panel replication*, not a literal
replication of the accepted 500-gene analysis. The 5% panel reduction may affect
gene H1 more strongly than cell topology, so both view-specific stability and a
500-versus-475 reference sensitivity must be reported rather than assumed.

The scientifically valid alternative is to preserve the exact-500 contract and
stop HCA topology work. Switching datasets is possible later, but any candidate
must pass the same header/identifier gate before a large expression download.

## Reproducibility

Production and clean-repeat outputs are byte-identical for all seven core
ledgers. The deterministic dossier assembly repeats byte-identically, and the
independent validator passes 17/17 checks. Its SHA-256 is
`b69a606a843686598e558c49083ba6703731605cdc2113f27c791af0798b4011`.
The complete repository suite passes 2,062 expectations with zero failures or
warnings and four established environment-dependent skips.

The public evidence is in
`docs/audits/mv08d-hca-structural-qc-evidence/`. The primary decision is
`mv08d-decision.csv`; the complete 25-gene discrepancy is
`mv08d-missing-panel-summary.csv`; and the four explicit recovery dispositions
are in `mv08d-recovery-options.csv`.

Dr. Eric C. Rouchka and Dr. Akshitkumar Mistry remain credited in project
metadata. MV8-D makes no final authorship or CRediT-role decision.
