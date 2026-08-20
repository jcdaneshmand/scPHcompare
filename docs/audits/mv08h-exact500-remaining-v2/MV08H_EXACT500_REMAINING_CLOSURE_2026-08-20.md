# MV8-H exact-500 remaining-unit closure

Date: 2026-08-20

The seven remaining adult bone-marrow HCA units completed the frozen Cell Ranger 8.0.1 exact-500 raw-read contract serially under the four-core/32-GiB policy. The corrected independent validator passed **68/68 checks** across all seven units.

## Verified gates

- Input binding: seven units, six manifest FASTQ files per unit; manifest SHA-256 `b315e097ef9b2f8cf1cad31cda08c4a440919641f9fff1674cf09e86928b2136`.
- Reference binding: 19 files, 20,765,871,518 bytes, tree SHA-256 `5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c`.
- Runtime: `cellranger cellranger-8.0.1`.
- Panel identities: exact500 `48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e`; common475 `b7b802ca862a63d7a4bbcaeab5af1192577663992a5ebde831371b6efafbc0ba`.
- Every unit: 33,563 feature IDs, 500/500 exact-panel IDs, 475/475 common-panel IDs, and at least 384 frozen-QC-eligible cells.
- Every unit: normal Cell Ranger success marker, empty stderr, required output structure, and passing resource evidence.

## Firewalls

This closure publishes only aggregate validation, input-binding identities, execution/resource evidence, and artifact hashes. Matrices, expression values, barcodes, donor metadata, labels, outcomes, persistence diagrams, landscapes, clustering results, fusion, manuscript files, and deletion remain closed.

Validator source SHA-256: `13968894aab073355a64c183716d7824b136a227c975b5125d0b42577a91df84`.
