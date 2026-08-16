# MV8-C HCA admission and compute dossier prefreeze v1

Date: 2026-08-16  
Status: prospective metadata and execution-plan audit; no expression-file
download and no scientific calculation

## Purpose and inherited decision

MV8-C converts the HCA candidate selected by MV8-B into an owner-reviewable,
file-exact admission and compute dossier. Candidate discovery remains closed.
No topology, outcome, or favorable file content may influence the cohort rule.

The fixed project is HCA
`cc95ff89-2e68-4a08-a234-480eca21ce79`, project label
`HematopoieticImmuneCellAtlas`, catalog `dcp60`. The official index reports a
mixture of cord blood, peripheral blood, selected marrow cells, individual
whole-marrow controls, and genetic pools. Metadata inspection before this
prefreeze established two marrow components:

- 24 whole-bone-marrow mononuclear sample entities: eight one-donor controls
  and sixteen eight-donor pooled entities; and
- 63 individually linked, selected `bone marrow hematopoietic cell` entities.

Those are eligibility facts, not topological outcomes. This prefreeze fixes the
primary subset before any count matrix, feature axis, barcode axis, PCA, PH,
landscape, distance, clustering, or tissue-retrieval result is opened.

## Exact primary cohort rule

An HCA entity enters the primary admission manifest only if every predicate is
true:

1. project ID is `cc95ff89-2e68-4a08-a234-480eca21ce79`;
2. catalog is `dcp60`;
3. organism is `Homo sapiens`;
4. organ is `bone marrow`;
5. selected cell type is `mononuclear cell of bone marrow`;
6. sample entity is a specimen and is publicly accessible;
7. donor count is exactly one;
8. library construction is `10x 3' v2`;
9. the file is a non-intermediate H5 `Count matrix`; and
10. exactly one distinct donor, specimen, cell suspension, library, and H5
    count file can be linked without a pooled-donor edge.

The expected metadata result is eight files representing eight biological
donors. Technical replicates, if discovered, remain linked to their donor and
cannot become additional biological units.

## Prespecified exclusions

- All sixteen whole-marrow entities with `donorCount = 8` are excluded from
  the primary cohort because genetic pooling does not provide one matrix per
  independent donor.
- All 63 `bone marrow hematopoietic cell` entities under organ part
  `bone marrow` are excluded from direct whole-marrow replication because they
  are composition-restricted. They may be considered only in a separately
  prefrozen sensitivity sprint.
- Project-level or component-level H5AD matrices are excluded because they
  aggregate multiple samples and cannot be treated as one biological unit.
- FASTQ, BAM, loom, CSV, ADT, cell-hashing, and intermediate files are not in
  the minimal primary retrieval set.
- No file may be substituted after its matrix contents or topology are seen.

## Current sprint boundary

MV8-C may retrieve and hash:

- the official OpenAPI schema;
- official project and sample-index JSON;
- a compact metadata manifest for the exact eight-file subset;
- public access terms; and
- HTTP metadata needed to name size, version, and checksum.

It may not issue a GET or Range request for an H5, H5AD, FASTQ, BAM, loom, or
other expression/sequence payload. It may not inspect an H5 header or feature
axis, read expression values, fit or project PCA, select cells, calculate PH or
landscapes, compute distances or clusters, open tissue-retrieval outcomes, or
promote a manuscript claim.

## Metadata admission gates

The dossier can recommend a later structural download only when all of the
following metadata gates pass:

1. exactly eight manifest rows and eight unique H5 file UUIDs;
2. exactly eight unique one-donor units with no pooled-donor membership;
3. one specimen, cell suspension, and library chain per admitted file, or an
   explicit donor-linked technical-replicate rule;
4. all files are accessible, non-intermediate H5 count matrices;
5. every file has a stable version, byte count, and SHA-256 checksum;
6. the summed H5 payload is within the 2-GiB metadata-stage download cap;
7. no project, accession, library, BioSample, or donor-code overlap is found
   with the frozen 124-sample corpus;
8. the HCA Data Release Policy, Data Use Agreement, CC BY 4.0 obligation, and
   project citation are recorded; and
9. the exact manifest and every official input are hash-addressed without
   publishing expiring credentials or protected donor attributes.

Public pseudonyms may establish uniqueness but are not published when a hash
or aggregate count is sufficient. Absence of known overlap is not described as
proof of real-world identity separation beyond the available pseudonymized
metadata.

## Structural gates after a separately authorized download

Passing this sprint does not establish view eligibility. A later, separately
authorized structural-download gate must verify, before reading expression
values or running PCA/PH:

- downloaded bytes and SHA-256 match the frozen manifest;
- the H5 object is a nonnegative integer count matrix with an unambiguous
  features-by-cells orientation;
- stable symbols and/or Ensembl identifiers map uniquely;
- each donor retains at least 384 usable cells under the frozen QC rule;
- all 500 ordered MV7-FP genes are present without zero-padding or
  outcome-informed replacement; and
- the exact reference center, scale, 30-PC rotation, gene order, and identity
  can be reused without refitting.

Because both the cell PCA and gene topology use the same frozen 500-gene
source matrix, incomplete panel closure blocks both primary views. A separate
sensitivity contract may be proposed before outcomes, but cannot replace the
primary definition.

## Frozen future compute queue

If and only if all eight files pass the later structural gate, the maximum
primary workload is:

- eight biological donors;
- five deterministic 384-cell seeds per donor;
- two separate views per donor-seed;
- one H0/H1 PH record per view, for 80 PH jobs total;
- 8 x 124 x 5 x 2 x 2 = 19,840 query-to-reference component landscape
  distances; and
- no query-query distance, external refit, or reference recomputation in the
  primary validation.

The existing 124-sample reference artifacts remain immutable. Every external
query is projected through the accepted center, scale, and 30-PC rotation.
Gene topology uses the unchanged 500-gene panel and the same 384 cells.

The prespecified primary endpoint is per-query same-tissue retrieval AUROC
when the 31 retained bone-marrow references are ranked against the other 93
reference samples. Top-1 and top-5 bone-marrow retrieval, and the median
bone-marrow-versus-other distance margin, are secondary. Endpoints are
reported for cell-H0, cell-H1, gene-H0, and gene-H1 separately. The unweighted
within-view H0/H1 composite remains secondary and descriptive; cross-view
fusion remains outside this sprint.

Labels are used only after all query coordinates, PH, landscapes, distances,
and immutable receipts pass validation. No endpoint selects the panel, seed,
view, homology dimension, distance, or candidate.

## Landscape contract

Any later admitted calculation must retain:

- finite positive-persistence intervals only;
- essential H0 excluded;
- every consecutive active landscape level;
- H0 and H1 calculated and reported separately;
- exact or error-controlled squared-L2 integration;
- no universal fixed grid or landscape-level cap; and
- streamed or chunked execution.

## Resource and stop policy

Planning reuses measured MV7-H resource evidence. Ripserr remains primary;
exact GUDHI is permitted only under the approved resource-fallback contract
and 12-GiB cap. The later execution stops without scientific exclusion if any
file, donor, feature, depth, transform, checksum, memory, timeout, or
reproducibility gate fails.

Current-sprint outcomes are limited to:

- `metadata_closed_download_candidate`;
- `blocked_manifest_or_donor_mapping`;
- `blocked_access_or_checksum`;
- `blocked_independence`;
- `blocked_resource_plan`; or
- `not_primary_whole_marrow`.

No current-sprint outcome authorizes an expression download or calculation.

## Privacy, publication, and credit

Public evidence may contain official URLs, UUIDs for project/files/sample
entities, aggregate counts, bytes, hashes, and pseudonym-free linkage states.
It must not contain expiring signed URLs, access tokens, private paths,
expression values, cell-level metadata, protected donor attributes,
dissertation/manuscript PDFs, or confidential peer-review text.

Dr. Eric C. Rouchka and Dr. Akshitkumar Mistry remain in the project credit
registry. This dossier does not decide author order or final CRediT roles.

## Acceptance criteria

- the four official metadata queries are frozen before the compact manifest
  is retrieved;
- the eight primary files and their donor/specimen/library relations are
  resolved without opening a matrix payload;
- pooled and selected-cell exclusions are explicit and reproducible;
- exact bytes, versions, checksums, use terms, and unresolved structural gates
  are documented;
- cell and gene admission states remain separate but both honor the frozen
  500-gene source contract;
- the 80-job and 19,840-distance maximum queues are reproduced from inputs;
- production and repeat dossier bundles are byte-identical and independently
  validated;
- no expression payload, PCA, PH, landscape, distance, clustering, outcome,
  claim, or author decision is introduced; and
- the final decision stops at an owner review of the proposed structural
  download.

