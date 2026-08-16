# MV8-B read-only external-dataset admission audit

Date: 2026-08-16

## Decision

MV8-B is complete. The bounded six-candidate audit passed 10/10 independent
checks, and a clean repeat reproduced all eight public evidence files
byte-for-byte. No external expression object was downloaded, no new PCA,
persistent homology, landscape, distance, clustering, or biological outcome
was calculated, and no manuscript or author-role decision was made.

The separately stored GSE120221 bone-marrow analysis remains in the plan, but
its role is now narrower and more accurate: it is a same-library technical
reprocessing population, not an independent external biological cohort. The
HCA HematopoieticImmuneCellAtlas adult bone-marrow cohort is the preferred
independent candidate, pending a separate metadata/download/compute decision.

## GSE120221 accession and source finding

Official NCBI run metadata closes all 25 historical validation libraries to
the 25 `SRA779509`/`SRP162214`/`PRJNA492227` rows already present in the
retained 124-sample corpus. Every GSM, BioSample, SRS, and SRR maps exactly;
all 25 pairs also share normalized local cell barcodes. Separate storage and
preprocessing therefore do not establish biological independence.

The cohort represents 25 libraries from 20 healthy donors. All 25 historical
validation inputs and all 25 paired main-corpus inputs pass the fixed 384-cell
depth gate. The historical validation metadata range is 789--6,990 post-QC
cells, and the paired main-corpus range is 874--7,261.

The two topology views receive different dispositions:

- **Cell view:** admissible for a same-library technical/reprocessing check,
  subject to a prospective donor/replicate aggregation rule before any new
  calculation.
- **Gene view:** non-estimable under the primary frozen 500-gene contract.
  Every historical validation matrix retains 467 panel genes, whereas every
  paired main-corpus matrix retains all 500. The panel is not reduced or
  replaced. A 467-gene analysis would require a separate pre-outcome
  sensitivity contract and could not be reported as the primary gene view.

This finding corrects earlier repository wording that treated the separately
processed 25-library path as an additional validation sample set outside the
127/124/90 flow. The processing path is separate; the biological libraries
are not.

## Candidate dispositions

| Candidate | Prospective role | MV8-B disposition |
|---|---|---|
| GSE120221, 25 libraries | Same-library technical reproducibility | Retain for cell-view technical reprocessing only; primary gene view non-estimable at 467/500 frozen genes |
| HCA HematopoieticImmuneCellAtlas adult bone marrow | Independent same-tissue replication | Preferred pending exact manifest, donor/library map, per-donor depth, identifier and 500-gene closure, byte/checksum inventory, and download authorization |
| GSE185381 healthy controls | Conditional same-tissue replication | Reserve pending disease-safe donor reconstruction and pooled/hashed CITE-seq demultiplexing |
| GSE180298 healthy HSPCs | Restricted-composition sensitivity | Retain as sensitivity only; sorted HSPCs cannot establish whole-marrow replication |
| GSE212092 day-zero bone marrow | Small external feasibility | Retain for feasibility only; three baseline donors do not meet the frozen minimum of four |
| Tabula Sapiens/GSE201333 | Multi-tissue generalization | Defer to a later stage requiring prospective donor-tissue unit reconstruction |

The HCA project is preferred because the associated adult bone-marrow analysis
reports more than 100,000 cells across eight healthy donors, which directly
matches the independent whole-marrow estimand better than the pooled,
restricted-composition, small, or multi-tissue alternatives. It is not yet
admitted: all hard gates remain unresolved until exact file-level metadata can
be inspected under a separately authorized retrieval.

## Resource envelope

For eight HCA adult bone-marrow donor units, the existing five-seed/two-view
design projects a maximum of 80 PH jobs. Using observed MV7-H timings gives
about 0.280 median worker-hours or 0.303 conservative p95 worker-hours for PH
itself. The measured PH peak was 6,654,734,336 bytes and the 124-sample source
fit peak was 806,625,280 bytes. These estimates exclude download, conversion,
reference projection, landscape distance streaming, and storage, so they do
not yet authorize execution. Ripserr remains primary and the exact GUDHI
resource fallback remains bounded by the approved 12-GiB cap.

## Landscape and view contract retained

No dataset may change the accepted dissertation-aligned landscape definition:

- biological sample is the comparison and clustering unit;
- cells and the fixed 500 genes are separate within-sample observation views;
- only finite positive-persistence intervals are retained;
- the essential H0 class is excluded;
- every consecutive active landscape level is included;
- H0 and H1 are calculated and reported separately;
- squared-L2 integration is exact or error-controlled;
- no universal grid or landscape-level cap is introduced; and
- computation is streamed or chunked, with the unweighted H0/H1 composite
  remaining secondary and descriptive.

Failure of the GSE120221 gene view does not authorize shrinking the panel,
changing the cell view, truncating landscapes, or combining H0 and H1.

## Reproducibility, privacy, and credit

The accepted prospective implementation head is
`7229f08032e4dddaaf8106c61778b736aa273d59`. The audit hashes eight frozen
input sources, including private metadata by hash only, and publishes no
private path, expression value, cell-level label, dissertation/manuscript PDF,
or confidential peer-review text. Seven generated data artifacts are
independently rehashed by their manifest; the manifest and those seven files
are identical in the production and repeat bundles.

The complete source-loaded package test suite also passed. Its only four skips
were the established optional-Rust, build-context, and CRAN subprocess/toy
guards; there were no failures or warnings.

Dr. Eric C. Rouchka and Dr. Akshitkumar Mistry remain in the project credit
registry. Dataset admission does not decide author order or final CRediT roles.

## Official sources

- [GSE120221 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221)
- [GSE120221 source study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6328018/)
- [HCA HematopoieticImmuneCellAtlas project](https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79)
- [HCA adult bone-marrow analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC6296228/)
- [GSE185381 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381)
- [GSE180298 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180298)
- [GSE212092 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212092)
- [Tabula Sapiens/GSE201333 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201333)

## Durable evidence

- Prefreeze specification:
  `docs/specifications/MV08B_READ_ONLY_DATASET_ADMISSION_AUDIT_PREFREEZE_V1.md`
- Candidate registry:
  `docs/specifications/mv08b-dataset-candidate-registry-v1.csv`
- Official-source fact table:
  `docs/audits/mv08b-dataset-source-facts-v1.csv`
- Evidence bundle:
  `docs/audits/mv08b-dataset-admission-evidence/`
- Independent validation:
  `docs/audits/mv08b-dataset-admission-validation/mv08b-independent-validation.csv`
- Builder and validator:
  `scripts/build_mv08b_read_only_dataset_admission_audit.R` and
  `scripts/validate_mv08b_read_only_dataset_admission_audit.R`

## Next authorized action

Prepare an owner-reviewable HCA metadata/download/compute dossier. It must name
the exact adult-bone-marrow files, donor and technical-replicate mapping,
per-donor expected depth, count and feature identifier state, frozen-panel
closure rule, bytes and checksums, data-use obligations, reference-only cell
projection, maximum job queue, storage, and stop/resume criteria.

This closure does **not** authorize external expression download, new PH, a
467-gene sensitivity, manuscript claim promotion, or material author-team
decisions.
