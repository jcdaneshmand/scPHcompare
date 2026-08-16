# MV8-A corrected publication figures

Date: 2026-08-16

## Decision

MV8-A is complete. Seven corrected publication figures were generated as
editable SVG and 300-dpi PNG pairs, their eleven public data tables were
published, a clean second render reproduced all 29 output files byte-for-byte,
and independent validation passed 11/11 checks. Visual review of all seven PNG
files found no clipping, overlap, illegible labels, or unresolved legacy
landscape language.

This authorizes author review of the figures. It does not authorize manuscript
submission, claim promotion, external-data download, or new persistent
homology.

## Scientific contract retained

The figures preserve the accepted dissertation-aligned definition:

- biological samples are the comparison and clustering units;
- cells and genes are distinct within-sample observation views;
- finite positive intervals are retained and the essential H0 class is
  excluded;
- every consecutive active landscape level enters the scientific distance;
- H0 and H1 are calculated and reported separately;
- squared-L2 integration is exact or error-controlled;
- no fixed scientific grid or universal landscape-level cap is used; and
- the unweighted within-view H0/H1 composite remains secondary and
  descriptive.

Figure 3 evaluates six landscape levels on a 600-point grid only for display.
Its labels explicitly distinguish that display grid from the all-level,
exact-or-error-controlled scientific calculation.

## Figure family

1. Corrected sample-comparison architecture, including the separate cell and
   gene views and secondary composite branch.
2. Cohort accountability and confounding, including the 127/124/90 flow and
   perfect approach nesting.
3. A fixed hash-validated finite persistence diagram and consecutive H0/H1
   landscape illustration.
4. All-pair H1 contribution distributions and all H0/composite partition
   concordances, with the two medians and one non-identical cell-average seed
   directly labeled.
5. All 54 label-free stability rows with jackknife ribbons, one-SE thresholds,
   and the selected k=2 marked for all six representations.
6. All 120 descriptive outcome units with the non-estimable primary-90
   approach endpoint stated explicitly.
7. All 30 PAM-versus-average-linkage partition comparisons without favorable
   algorithm selection.

Figure 8 remains deliberately deferred to a cross-stage synthesis so that the
90-sample blocked estimand is not visually conflated with the 124-sample
descriptive estimand.

## Reproducibility and privacy

The accepted prospective implementation head is
`912df83acb3a655becf61dd1dc86e7969d5f063f`. The renderer independently checked
the accepted MV7-H hash for the fixed PH artifact and the MV7-I hash for the
private all-pair H1 table. Public tables contain neither the fixed sample
identity nor confidential reviewer material.

Independent validation confirmed:

- seven SVG and seven PNG files with manifest-matching hashes and dimensions;
- live vector text with no embedded raster or legacy first-level/100-grid
  claim;
- complete row families of 12/12, 16/5, 734/7,200, 600/20, 54, 120, and 30;
- explicit corrected-landscape language;
- no sample identifiers, confidential review text, manuscript PDF, or
  dissertation PDF in the public bundle;
- no new data, PH, p-values, ranking, or claim authorization; and
- all 29 output files reproduced byte-for-byte in a clean repeat.

The complete source-loaded package test suite also passed. Its only four skips
were the established optional-Rust, build-context, and CRAN subprocess/toy
guards; there were no failures or warnings.

## Visual review

All seven PNGs were opened and inspected after structural validation. The
review corrected and then rechecked deterministic jitter, vector-safe color
encoding, log-percentage labels, direct Figure 4 annotations, explicit Figure
5 one-SE thresholds, selected-k markers, and concise numeric legends. The
final rerender passed the same 11/11 independent checks.

## Durable evidence

- Prefreeze specification:
  `docs/specifications/MV08A_CORRECTED_PUBLICATION_FIGURE_PREFREEZE_V1.md`
- Renderer:
  `scripts/render_mv08a_corrected_publication_figures.R`
- Independent validator:
  `scripts/validate_mv08a_corrected_publication_figures.R`
- Figures and public data:
  `docs/figures/mv08a-corrected-publication/`
- Validation evidence:
  `docs/audits/mv08a-figure-evidence/`

## Next authorized action

Proceed to the read-only external-dataset admission audit. The audit must
include the separate 25-sample GSE120221 bone-marrow technical-validation
cohort as a named first candidate, distinct from the 31 bone-marrow samples in
the heterogeneous cohort and outside the 127/124/90 flow. It should assess
observation independence, accession overlap, metadata provenance, raw-input
availability, cell-view compatibility, fixed-gene-panel overlap, H0/H1
estimands, and resource requirements before any new calculation.

No external download or new PH is authorized by this closure.
