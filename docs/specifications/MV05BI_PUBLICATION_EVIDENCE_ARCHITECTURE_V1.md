# MV5-BI publication evidence architecture v1

## Decision

Use a lean GitHub repository plus two linked Zenodo objects:

- a versioned **software** DOI created from a tagged GitHub release; and
- a versioned **derived-evidence dataset** DOI containing deterministic audit
  tables and derived R result objects.

This follows common computational-biology practice: GitHub remains the living,
reviewable source repository, while a DOI-bearing archive freezes the exact
software and evidence used by a manuscript version.

## Git-resident classes

- R, Python, shell, PowerShell, C, and Rust source;
- package metadata, dependency lock, public ABI header, and workflows;
- package/unit/native tests and intentionally small analytical fixtures;
- specifications, architecture decisions, plans, and human-readable reports;
- a compact member-level evidence manifest with path, byte count, Git blob,
  SHA-256, source revision, category, and media type;
- archive build/verification and evidence-retrieval instructions.

## DOI-archive classes

- every changed `docs/audits/*.csv` and `docs/audits/*.csv.gz` generated table;
- every changed `results/**` derived R result object;
- an archive-local README, exact manifest, and metadata draft.

The archive classifier is path- and revision-bound. It does not include source
code, Markdown reports, private caches, reviewer correspondence, PDFs, raw
biological data, credentials, or untracked files.

## Artifact rules

1. Build from Git blobs at an immutable source commit.
2. Sort members bytewise by repository-relative path.
3. Use fixed ZIP metadata and deterministic compression.
4. Record SHA-256 and Git blob identity for every evidence member.
5. Verify archive membership, sizes, and SHA-256 before deposit.
6. Repeat the build and require identical archive SHA-256.
7. Reserve a DOI only after the author team reviews creator order, roles,
   ORCIDs, affiliations, licensing, and manuscript relationship.
8. Cite the exact version DOI in the manuscript; retain the concept DOI only
   for discovery across versions.

## Credit boundary

The evidence metadata draft provisionally follows the public preprint order:
Jonah Daneshmand, Julia H. Chariker, Akshitkumar Mistry, and Eric C. Rouchka.
This preserves project-owner direction that Dr. Rouchka and Dr. Mistry receive
credit without silently omitting a public preprint author. The author team must
approve creator/contributor classification, ordering, affiliations, ORCIDs, and
CRediT roles before any deposit or software release.

## Git publication rules

The accumulated modernization history is reconstructed into dependency-ordered
publication slices after evidence paths are omitted. Each target slice should
stay below 250 changed files, 15,000 changed lines, and 750,000 raw diff bytes,
leaving margin below GitHub's display limits. Every slice must pass its relevant
tests independently. The final source tree must match the accepted local source
revision for every Git-resident path.

## Closed actions

MV5-BI does not authorize a Zenodo draft, DOI reservation, upload, public
release, GitHub push, pull request, history replacement, R-default change,
Rust adoption, new biological calculation, partition, or manuscript claim.
