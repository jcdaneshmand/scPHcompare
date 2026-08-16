# MV8-E reference reconciliation and panel-sensitivity prefreeze

## Decision

MV8-E identifies the HCA quantification reference exactly and closes the
identifier-crosswalk question. The eight HCA feature axes are the official
Cell Ranger 3.0.0 GRCh38 reference obtained by applying the documented 10x
gene-biotype filter to Ensembl release 93. Filtered Ensembl IDs, names, order,
count, and newline-delimited ID-axis SHA-256 match the HCA H5 reference.

The 25 genes absent from every HCA count matrix are not retired or incorrectly
crosswalked IDs. All occur in the unfiltered Ensembl-93 annotation and all
still resolve through current Ensembl lookup, but their Ensembl-93 biotypes
were outside the official Cell Ranger 3.0.0 allow-list. They were therefore
not quantified in these HCA matrices. Symbol substitution and stable-ID
crosswalking cannot reconstruct their counts.

The exact ordered 475-gene intersection is now frozen. The next authorized
calculation is the reference-only 500-versus-475 sensitivity specified in
`MV08E_REFERENCE_RECONCILIATION_AND_PANEL_SENSITIVITY_PREFREEZE_V1.md`.
HCA FASTQ download and raw-read reprocessing remain closed until that
sensitivity has passed its independent gates and a separate custom-reference
prefreeze is accepted.

## Exact reference fingerprint

The official Ensembl release 93 GRCh38 GTF has SHA-256
`810e3bb63bb24bd5a005b14397d69280dda34c41a612fb86a18f3c4836fce57d`
and contains 58,395 gene rows. Applying the documented 17-biotype Cell Ranger
3.0.0 allow-list yields 33,538 genes, exactly the HCA feature count. Both the
filtered GTF axis and HCA feature axis hash to
`a6a914db1d218e2c3ae4c6680ac8fc546ba52ee4cfaedf423eb635c826c085be`;
the gene-name axes also match in exact order.

All 500 frozen panel stable IDs occur in unfiltered Ensembl 93. Exactly 475
occur in the filtered reference. The ordered `common475_v1` feature axis hashes
to `b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba`.
The original 500-gene panel remains the existing-data reference object; this
475-gene object is its prespecified harmonization sensitivity and external-data
candidate.

The 25 filtered annotations comprise:

- 17 processed pseudogenes;
- two unprocessed pseudogenes;
- one polymorphic pseudogene;
- one transcribed processed pseudogene; and
- four processed transcripts.

The supporting authoritative sources are the
[10x Cell Ranger reference-build steps](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps),
the [Ensembl release 93 GTF](https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz),
the [HCA bone-marrow publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC6296228/),
and the [Ensembl stable-ID lookup API](https://rest.ensembl.org/documentation/info/lookup).

## Why the sensitivity precedes raw reads

The reference-only sensitivity measures the actual consequence of removing
the 25 filtered features while holding the 124 samples, five 384-cell draws,
normalization, two topology orientations, PH definition, landscape definition,
and label firewall fixed. It will refit the 475-gene cell PCA and reconstruct
both cell and gene PH; it will not derive a 475 result by deleting rows from an
already-computed 500-gene persistence object.

The accepted 500-gene source, PH, and all 20 landscape groups remain intact.
All 170 added-sample SCT caches remain intact. The 450 older primary-sample SCT
caches are no longer present, but all 90 retained individual source files are
available. Those caches must be reproduced through the accepted MV5-D0 v2
pipeline and match their public logical and byte identities before the
sensitivity may begin. A mismatch is an environment/provenance stop, not a
reason to relax equality.

If all four separate primary components—cell-H0, cell-H1, gene-H0, and
gene-H1—show high harmonization stability, a 475-gene HCA run is a defensible
external harmonized-panel replication. Material sensitivity recommends
exact-500 raw-read reprocessing. Mixed results require component-specific
reporting and an owner decision. The project owner's preference for raw-read
reprocessing is recorded and remains eligible even if the 475 result is
stable.

## Raw-read feasibility boundary

The exact eight one-donor metadata records expose 48 FASTQ files totaling
85,034,239,918 bytes, or 79.194307 GiB. None has been downloaded. A later raw
path must prospectively freeze the custom GTF and filter, quantifier/version,
10x 3' v2 read handling, multimapping and overlap policy, compute/storage caps,
and comparison wording before acquisition or counting. Reprocessing HCA alone
also retains a quantification-pipeline asymmetry with the original 124 samples;
that limitation must be measured and disclosed rather than implied away.

## Landscape and label boundaries

The landscape object is unchanged: finite positive intervals, essential H0
excluded, every consecutive active level, H0/H1 separate, exact or
error-controlled squared-L2, no primary fixed grid or level cap, and streamed
or chunked execution. HCA expression and biological metadata remain unused.
No PCA, PH, landscape, distance, clustering, biological endpoint, or manuscript
claim was calculated in this reconciliation stage.

## Reproducibility and public/private boundary

The deterministic reconciliation repeats byte-for-byte for all nine production
artifacts. Independent validation passes 12/12 categories; its SHA-256 is
`ce171149043ea540f9a1050713a940309f4f3d0848fd300e407746cd1bb60374`.
Public evidence is under
`docs/audits/mv08e-reference-reconciliation-evidence/` and contains no
expression value, cell barcode, local absolute path, raw GTF, HCA H5 payload,
FASTQ payload, signed URL, biological outcome, or private reviewer material.

Dr. Eric C. Rouchka and Dr. Akshitkumar Mistry remain in the project credit
registry. MV8-E makes no final authorship-order or CRediT-role decision.
