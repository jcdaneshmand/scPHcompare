# MV8-B read-only external-dataset admission audit prefreeze v1

Date: 2026-08-16  
Status: prospective metadata/source audit; no external expression download and no new PH

## Purpose

MV8-B determines which existing or external datasets could answer a named
generalization question without selecting data from favorable topological
outcomes. It begins with the separately stored 25-library GSE120221
bone-marrow analysis used in the dissertation-era workflow, then audits a
bounded registry of independent bone-marrow and multi-tissue candidates.

This sprint may inspect local source axes, public accession metadata, data-use
terms, file formats, donor/sample structure, and published study design. It may
not download external expression matrices, calculate PCA or PH, calculate a
landscape or distance, open an external biological outcome, or promote a
manuscript claim.

## Questions and estimands

The audit separates four prospective roles:

1. `technical_reprocessing_validation`: the same public biological libraries
   processed through a distinct historical source/preprocessing path;
2. `same_tissue_external_replication`: independent bone-marrow donors queried
   against the frozen existing-corpus reference;
3. `composition_shift_sensitivity`: independent but deliberately sorted or
   lineage-restricted bone-marrow cells, reported as a sensitivity rather than
   a direct whole-marrow replication; and
4. `multi_tissue_external_generalization`: independent donor-tissue samples
   from a prospectively bounded multi-organ atlas.

No role can be upgraded after PH or retrieval outcomes are seen. A technical
reprocessing cohort cannot be called external biological validation.

## Incumbent GSE120221 rule

The 25 locally stored validation inputs are audited first. The audit must map
every validation GSM, BioSample, SRS, SRR, SRA submission, SRA study, and
BioProject using official NCBI metadata and compare those identifiers with the
31 bone-marrow samples inside the retained 124-sample corpus.

If any validation biological library appears in the 124-sample corpus, that
library is not independent. If all 25 map to the 25 `SRA779509` corpus rows,
the cohort is retained as a separately processed technical/reproducibility
population only. Its separate files, filtering, and legacy results remain
historically meaningful, but they cannot establish unseen-dataset
generalization.

The audit may compare local matrix dimensions, feature axes, and normalized
cell-barcode sets. It may not export expression values or private file paths.

## Hard admission gates

A candidate is eligible for a future primary external calculation only when
all applicable gates are resolved before expression download:

1. **Independence:** zero known accession, BioSample, library, donor, or source
   cohort overlap with the frozen 124-sample corpus. Unknown donor overlap is
   a blocker, not presumed independence.
2. **Observation unit:** individual biological donor/sample units can be
   reconstructed. Pooled libraries require a public, auditable demultiplexing
   map; technical replicates must remain linked to their donor.
3. **Assay and specimen:** human single-cell transcriptomic counts are
   available. Whole, unsorted bone-marrow mononuclear cells are required for a
   direct bone-marrow replication; selected HSPCs or immune subpopulations are
   sensitivity-only.
4. **Depth:** each admitted biological sample is expected to retain at least
   384 cells after a frozen QC rule and must never bypass the historical
   250-cell exclusion rule. Samples below 384 are non-estimable for the fixed
   matched-depth analysis.
5. **Counts and identifiers:** raw or minimally processed nonnegative count
   matrices, stable gene identifiers, cell barcodes, and sample/donor metadata
   are publicly obtainable under explicit access terms.
6. **Fixed gene panel:** the frozen 500-gene MV7-FP panel is tested without
   outcome-informed replacement. The primary external gene view requires all
   500 genes after a prospectively specified identifier normalization. Missing
   genes make the gene view non-estimable unless a separate, pre-outcome
   sensitivity contract is approved. The cell view is evaluated separately.
7. **Label firewall:** tissue, study, donor, disease, treatment, age, and
   approach labels may establish eligibility and blocking, but may not fit the
   transform, select cells/genes/levels, tune methods, or choose a candidate
   after results.
8. **Resource feasibility:** download size, source conversion, five 384-cell
   seeds, two views, H0/H1 PH, all-active-level landscapes, query-to-reference
   distances, storage, and peak memory are estimated before authorization.
9. **Use rights:** public access and reuse/citation requirements are recorded.
10. **Auditability:** official record URLs, access date, fact provenance,
    unresolved fields, and an explicit disposition are published without
    confidential reviewer text or private local paths.

## View-specific continuation

Cell and gene topology are separate views and receive separate admission
states. The cell view uses the frozen 124-sample feature/center/scale/PCA
reference and a prospective held-out projection route; external samples do not
refit the reference. The gene view uses the unchanged 500-gene panel and the
same 384 selected cells. Failure of one view does not authorize changing the
other.

Any future topology must preserve the accepted landscape contract:

- finite positive-persistence intervals only;
- essential H0 excluded;
- every consecutive active landscape level;
- H0 and H1 calculated and reported separately;
- exact or error-controlled squared-L2 integration on dimension-specific
  support;
- no universal grid and no universal landscape-level cap; and
- streamed or chunked computation.

The unweighted H0/H1 composite remains secondary and descriptive.

## Frozen candidate registry

The bounded registry is
`docs/specifications/mv08b-dataset-candidate-registry-v1.csv`. It contains:

- the existing GSE120221 25-library cohort;
- the HCA HematopoieticImmuneCellAtlas as the direct independent
  bone-marrow candidate;
- GSE185381 healthy-control records as a conditional pooled/demultiplexing
  candidate;
- GSE180298 healthy HSPCs as a composition-shift sensitivity;
- GSE212092 untreated baseline bone marrow as a small external feasibility
  candidate; and
- Tabula Sapiens/GSE201333 as the bounded multi-tissue candidate.

This registry closes candidate discovery for MV8-B. A later dataset requires a
new prospective amendment explaining why none of these candidates can answer
the named estimand.

## Allowed dispositions

- `preferred_pending_download_authorization`
- `reserve_pending_metadata_resolution`
- `technical_reprocessing_only`
- `composition_shift_sensitivity_only`
- `small_feasibility_only`
- `multi_tissue_later_stage`
- `reject_overlap`
- `reject_unit_not_reconstructable`
- `reject_incompatible_input`

Every disposition must name the estimand it can and cannot support. Unknown
facts remain explicit; they are not converted to passes.

## Resource policy

Planning uses observed MV7-H production measurements, not a generic benchmark.
The initial projection assumes five seeds, two views, one cell PH and one gene
PH per sample-seed, plus source preparation and streamed query-to-reference
landscape distances. The primary PH engine remains Ripserr. Exact GUDHI is
only the already approved resource fallback under the 12-GiB cap. A candidate
that cannot fit the cap is staged or rejected; the landscape definition is not
truncated to make it fit.

## Privacy and publication boundary

Public evidence may contain public accessions, aggregate donor/sample counts,
file formats, public URLs, source hashes, and aggregate resource projections.
It must not contain dissertation/manuscript PDFs, confidential reviewer text,
private paths, expression values, cell-level labels, or protected donor
attributes.

Dr. Eric C. Rouchka and Dr. Akshitkumar Mistry remain in the project credit
registry. Dataset admission does not decide author order or CRediT roles.

## Acceptance criteria

- all six frozen candidates have official-source facts and explicit unknowns;
- all 25 GSE120221 validation libraries are accession-mapped against the main
  corpus and local source compatibility is summarized;
- biological independence and separate storage/processing are not conflated;
- direct replication, composition sensitivity, technical reprocessing, and
  multi-tissue roles remain distinct;
- cell and gene admission states are separate;
- the fixed 500-gene panel and corrected all-active-level H0/H1 landscape
  contract are visible in every continuation decision;
- observed MV7 resource evidence is used for projections;
- a production/repeat audit bundle is byte-identical and independently
  validated;
- no external expression data, new PH, external outcome, or confidential
  material is introduced; and
- the final decision stops before download and calculation authorization.

## Gate

Passing MV8-B authorizes only an owner/author-team decision on one precisely
specified next calculation. External expression download, a new PH run, a
change to the fixed panel, a new data-use obligation, and a material manuscript
claim each remain separately gated.

