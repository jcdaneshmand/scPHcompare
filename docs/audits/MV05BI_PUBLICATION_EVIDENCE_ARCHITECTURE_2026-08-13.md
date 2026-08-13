# MV5-BI publication evidence architecture

## Outcome

MV5-BI accepts a publication-standard split between a lean GitHub software
repository and a separate DOI-bearing derived-evidence dataset. It implements
and validates the local predeposit package but performs no upload, DOI
reservation, push, PR, release, or history replacement.

The design follows current computational-biology and repository guidance:
source code remains visible and reviewable in GitHub; exact manuscript
software is frozen from a tagged release; generated evidence is preserved in
a versioned archival record with a persistent identifier; and both objects are
linked and cited by their exact version DOIs.

## Frozen source and classification

- Publication base: `3910b15ce6d3f88197847a8aba94b184e7d9c034`.
- Evidence source: `53c5f2e0864aff288d34871fb8272da15ac660ca`.
- Generated audit tables: 995 files and 261,595,043 bytes.
- Derived R result objects: 141 files and 5,958,215 bytes.
- Total: 1,136 files and 267,553,258 uncompressed bytes.

The classifier includes changed `docs/audits/**/*.csv`,
`docs/audits/**/*.csv.gz`, and `results/**` paths. It keeps source, tests,
locks, human-readable Markdown reports, specifications, and compact
publication manifests in Git. It operates on Git blobs at the frozen source
revision and therefore cannot ingest `example_run.r`, local PDFs, reviewer
correspondence, private caches, or mutable working-tree content.

## Deterministic bundle

Two independent builds produced a 161,963,473-byte ZIP with identical SHA-256:

`0cd4fa24c407cfe3611c1e424889f38c15bb66d1784039c3d1291325155afa62`

The verifier independently reopened the ZIP and matched all 1,136 evidence
members to the Git-resident manifest by path, size, and SHA-256. The archive
also contains its exact manifest, a README, and an explicitly provisional
Zenodo metadata draft. The draft follows the public preprint order of Jonah
Daneshmand, Julia H. Chariker, Akshitkumar Mistry, and Eric C. Rouchka while
marking creator/contributor classification, order, ORCIDs, affiliations,
CRediT roles, and license as author-team decisions.

## Lean software proof

A source tree was materialized from the same Git revision with every evidence
member omitted. It contained no `results` directory and no generated audit
CSV/CSV.GZ. The R source package built successfully, installed, loaded,
executed its examples, and ran its complete test suite. `R CMD check
--no-manual --no-build-vignettes` completed with `Status: OK` using the
previously accepted locked library. The final check used R 4.4.1 on Ubuntu
22.04.4 and completed in 497.9 seconds.

This proves the generated evidence is supplementary to the distributable
software rather than a hidden runtime dependency.

## Reviewable Git publication stack

After the evidence paths are excluded, an exact dependency-ordered planner
produces eight source slices. Each is bounded below conservative targets of
150 files, 15,000 changed lines, and 1,000,000 raw diff bytes. Observed maxima
are 138 files, 14,966 changed lines, and 854,989 raw bytes. The full boundary
ledger is `docs/publication/scphcompare-publication-stack-v1.csv`.

These are source-history boundaries, not branches yet. The later reconstruction
must reproduce the final Git-resident tree exactly and validate each stacked
branch independently before any push.

## Citation and credit correction

The root citation metadata no longer uses a placeholder DOI. It now identifies
the existing bioRxiv preprint DOI `10.1101/2025.07.24.666637` and its
authoritative author order: Jonah Daneshmand, Julia H. Chariker, Akshitkumar
Mistry, and Eric C. Rouchka. Package metadata credits the latter three as
contributors while retaining Jonah Daneshmand as package author and
maintainer. Final revised-manuscript authorship and CRediT roles remain an
author-team decision.

## Decision and next gate

The reasonable next sprint is a **non-destructive stacked-branch
reconstruction rehearsal** in temporary worktrees. It should reconstruct the
eight lean slices without the archived evidence, prove final Git-resident tree
equivalence, run slice-appropriate checks, and retain the original 169-commit
branch as a local safety reference. It must stop before replacing the working
branch, authenticating GitHub, pushing, creating PRs, reserving a DOI, or
uploading the verified ZIP.
