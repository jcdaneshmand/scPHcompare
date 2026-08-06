# Phase 9 Subplan — Release and Archive

## Objective

Publish an installable, reproducible, citable artifact aligned with the accepted submission package while protecting confidential and restricted material.

## Task ledger

| ID | Task | Required evidence/output | Acceptance test | Status |
|---|---|---|---|---|
| P9-01 | Define release contents | Inclusion/exclusion manifest and license/data-rights review | Every artifact has distribution authority and retention decision | `not_started` |
| P9-02 | Run privacy/confidentiality audit | Git-history/tree/release-bundle scan and signed checklist | No private reviewer material, secrets, identifiers, or restricted data | `not_started` |
| P9-03 | Run clean-machine reproduction | Fresh supported systems, commands, timing/hardware, outputs/checksums | Installation and documented example pass; full-analysis limits are explicit | `not_started` |
| P9-04 | Finalize public documentation | README, tutorials, configuration, restart, data access, hardware/runtime, troubleshooting | New user can identify quick and full reproduction routes | `not_started` |
| P9-05 | Finalize metadata | Version, changelog, CITATION, authors/CRediT, license, software/data availability | Metadata matches approved authorship and manuscript | `not_started` |
| P9-06 | Create release candidate | Signed/tagged candidate and immutable artifact manifest | Candidate matches tested commit and checksum bundle | `not_started` |
| P9-07 | Archive and mint DOI | Approved archive record and DOI | Archive resolves, metadata is correct, artifacts match release | `not_started` |
| P9-08 | Publish GitHub release | Semantic version, notes, archive/DOI links | Public release matches candidate exactly | `not_started` |
| P9-09 | Preserve reproducibility record | Final environment, manifests, logs, known limitations, archival locations | Independent audit can reconstruct what was released and why | `not_started` |
| P9-10 | Evaluate G9 | Final release sign-off | Installable, reproducible, citable, private-safe, manuscript-aligned | `not_started` |

## Gate G9 checklist

- [ ] Release contents and rights are approved.
- [ ] Confidentiality and secret scans pass, including Git history and bundle.
- [ ] Clean-machine reproduction passes supported routes.
- [ ] Version, citation, authorship, archive, and manuscript metadata agree.
- [ ] DOI archive and GitHub release match by checksum.
- [ ] Final limitations and support expectations are published.
