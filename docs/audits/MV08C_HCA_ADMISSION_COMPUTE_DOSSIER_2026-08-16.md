# MV8-C HCA admission and compute dossier closure

Date: 2026-08-16  
Status: complete at the metadata-closed download-candidate boundary  
Prospective cohort prefreeze: `cc281cd`  
Accepted implementation and evidence head: `6c6f600c9a3844ad08f36e0ba8b0e97fcab4c0b3`

## Decision

The HCA project `cc95ff89-2e68-4a08-a234-480eca21ce79` is accepted as a
`metadata_closed_download_candidate` for independent adult whole-bone-marrow
validation. The exact primary retrieval set is eight one-donor H5 count files,
one for each of eight unique donor/specimen/suspension/sequencing-process
chains. This decision does not admit either scientific view and does not
authorize an expression-file download, PCA, PH, landscape, distance,
clustering, tissue-retrieval endpoint, or manuscript claim.

The next owner decision is whether to authorize the exact eight-file
structural/QC download sprint. No substitution or outcome-informed file
selection is allowed.

## Official-source closure

The prospective query contract was committed before the compact manifest was
retrieved. Two independent metadata acquisitions returned byte-identical
official inputs:

- HCA Data Browser API/OpenAPI version `19.0`, catalog `dcp60`;
- one project entity;
- 24 whole-bone-marrow mononuclear sample entities;
- 63 selected bone-marrow hematopoietic-cell entities; and
- eight exact one-donor H5 manifest rows.

The [HCA project page](https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79)
reports access granted, project label `HematopoieticImmuneCellAtlas`, INSDC
project accession `ERP122984`, BioStudies accession `S-SUBS12`, 1,453,784
project cells, 28 project donors, and an update date of 2025-02-14. These are
project-wide values, not the primary eight-file post-QC depth.

The [HCA quick-start guide](https://data.humancellatlas.org/guides/quick-start-guide)
defines the compact TSV as metadata rather than the actual data payload. Raw
API JSON and the compact manifest remain ignored inputs; their hashes and
sanitized conclusions are public.

## Cohort accountability

The marrow inventory closes without topology, labels, or outcomes:

| Component | Count | Primary status | Reason |
|---|---:|---|---|
| One-donor whole-marrow mononuclear samples | 8 | metadata candidate | one biological donor per file |
| Eight-donor whole-marrow pools | 16 | excluded | donor-specific topology is not reconstructable |
| Selected marrow hematopoietic-cell samples | 63 | excluded from primary | composition-restricted rather than whole marrow |
| Aggregate marrow H5AD | 1 | excluded | multiple samples/pools combined |

The eight admitted metadata rows contain eight unique file UUIDs, versions,
SHA-256 checksums, donor documents, specimen documents, cell suspensions, and
sequencing processes. All use `10x 3' v2`; the linked H5 summaries are
accessible, non-intermediate count matrices. The compact chain represents the
library through a unique sequencing process and library-preparation protocol;
it does not expose a separate library document identifier.

The exact H5 payload totals `202,770,089` bytes (`0.1888443613` GiB), safely
below the frozen 2-GiB retrieval cap. The official `7,000` cells per sample is
a metadata loading estimate, not a post-QC cell count and not evidence that the
384-cell gate passes.

## Independence conclusion and limits

Exact and punctuation-insensitive comparisons against every frozen main-corpus
metadata field found zero overlaps for the HCA project UUID/accessions, donor
codes/documents, specimen documents, sequencing processes, file UUIDs, and
file names. The HCA project accession is also distinct from the main corpus's
SRA study identifiers.

The defensible conclusion is **no known metadata overlap**. It is not proof
that pseudonymized records cannot represent the same real person. The compact
manifest does not expose BioSample accessions, so direct BioSample comparison
is explicitly not testable from this frozen metadata format. Those limits are
retained in the public independence ledger and in any later manuscript text.

## Access, privacy, and attribution

The project page states that downloaded/exported data is governed by the
[HCA Data Release Policy](https://www.humancellatlas.org/data-release-policy/)
and licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The
[HCA Data Use Agreement](https://data.humancellatlas.org/about/data-use-agreement)
also identifies CC BY 4.0 for this release. Later use must attribute the
project, link the license, identify modifications, and preserve the applicable
release-policy obligations.

Public artifacts contain stable file locators, UUIDs, versions, byte counts,
checksums, aggregate counts, and one-way hashes needed for auditability. They
contain no donor ages/sex, contact information, signed URLs, access tokens,
private paths, expression values, cell-level metadata, dissertation/paper
PDFs, or confidential reviewer text.

Dr. Eric C. Rouchka and Dr. Akshitkumar Mistry remain in the project credit
registry. MV8-C makes no author-order or CRediT-role decision.

## View-specific admission boundary

Both proposed views use the same frozen ordered 500-gene source matrix and the
same deterministic 384 cells per donor-seed. Metadata closes the eight-unit
cohort and file identities, but the following hard gates remain unresolved
until a separately authorized exact download:

- H5 object structure, count type, and features-by-cells orientation;
- stable feature identifiers and unique mappings;
- presence of all 500 ordered MV7-FP genes without padding or replacement;
- at least 384 usable post-QC cells in every donor; and
- exact reuse of the accepted reference center, scale, and 30-PC rotation
  without external-data refitting.

Failure of the feature or depth gate blocks both primary views because both
share the same source matrix. Cell topology would then use cells as
observations in frozen 30-PC Euclidean space; gene topology would use the 500
genes as observations under Pearson chord distance. Neither coordinate system
has been calculated for HCA in this sprint.

## Landscape and compute contract

The accepted dissertation-aligned landscape contract remains unchanged:

- finite positive-persistence intervals only;
- essential H0 excluded;
- all consecutive active levels;
- H0 and H1 computed and reported separately;
- exact or error-controlled squared-L2 integration;
- no universal grid or landscape-level cap; and
- streamed or chunked execution.

If all eight files later pass structural/QC admission, the maximum primary
queue is 8 donors x 5 deterministic seeds x 2 views = **80 view-level PH
jobs**, followed by 8 x 124 x 5 x 2 views x 2 homology dimensions = **19,840
query-to-reference component distances**. No query-query distances, reference
refits, or reference recomputations are planned. Ripserr remains primary;
exact GUDHI is available only under the approved fallback policy and 12-GiB
cap.

The proposed later endpoint is same-tissue retrieval against the immutable
124-sample reference: 31 bone-marrow references versus 93 other references.
Cell-H0, cell-H1, gene-H0, and gene-H1 remain separate primary reporting
components; within-view H0/H1 composites are secondary and cross-view fusion
remains excluded. Labels stay closed until coordinates, PH, landscapes,
distances, and prediction receipts pass validation.

## Reproducibility and validation

The production and clean-repeat dossier assemblies are byte-identical for all
15 generated files. The independent validator passes all 15 categories:
source identity, project identity, component accountability, exact file
manifest, unit linkage, independence limits, access terms, view gates,
resource queue, execution/stop-resume policy, landscape contract, decision
boundary, artifact hashes, byte-repeat identity, and privacy/confidentiality.
The validation ledger itself repeats byte-identically with SHA-256
`43f9f12ed8add38c4c477ab250d51774799a75a003a9790c4777f65ff5d8df6f`.

Focused tests pass after publication. The complete repository suite passes
`2,023` checks with zero failures/warnings and two established skips (the
optional MV5-BC Rust library is absent, and public audit documents are
intentionally excluded from R builds).

## Exact next sprint

MV8-D should be a bounded exact-file structural/QC admission sprint only:

1. obtain owner authorization for the eight frozen H5 files;
2. download them sequentially using the published stable locators;
3. verify every byte count and SHA-256 before opening a file;
4. inspect H5 structure, axes, feature identifiers, and count type;
5. apply the frozen QC definition to establish the real usable-cell count;
6. test exact ordered 500-gene coverage and mapping;
7. validate the immutable external-projection interface; and
8. publish a pass/block ledger before authorizing any PCA or PH.

The sprint must stop on any manifest, checksum, structure, mapping, panel,
depth, transform, privacy, or reproducibility failure. Such a stop is a
technical gate, not a scientific exclusion and not a reason to substitute a
different donor after seeing results.
