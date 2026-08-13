# MV5-BH pre-push branch reconciliation

## Outcome

The branch is Git-reconciled and privacy-clean, but it is not accepted for a
single push/PR. A fresh fetch found `origin/main` unchanged at `3910b15`; that
revision is exactly the merge base of the audited local head `af81621`. There
are zero remote-only commits and 166 local commits. No merge or rebase is
needed.

Publication remains closed because the proposed comparison contains 1,869
files, 530,811 insertions, and 260.69 MiB of added or modified blobs. GitHub's
current documented diff display is limited to 300 files, 20,000 loadable
lines, and 1 MiB of raw diff data. The branch exceeds each review boundary by
a wide margin. Its largest object is 25.61 MiB, so no regular-Git 100 MiB hard
limit is breached, but nine objects exceed 10 MiB and much of the branch is
programmatically generated scientific evidence.

## Git and authentication state

- `origin` remains `https://github.com/jcdaneshmand/scPHcompare.git` for fetch
  and push.
- Fetch completed without authentication because the repository is public.
- GitHub CLI is not authenticated. This is not a local Git failure; it becomes
  relevant only after publication is otherwise accepted.
- No push, upstream assignment, PR, hosted run, artifact upload, or release
  occurred.

## Privacy and history audit

All 166 proposed commit trees were scanned, rather than only the final tree.
No credential signature was detected. No PDF, `docs/private`, `tmp`,
`provenance`, or `example_run.r` path appears in the proposed history. The
paper and dissertation occur only as explicit local/ignored source references
and recorded evidence hashes; their contents are absent. Reviewer
correspondence is absent.

The 141 tracked RDS files are confined to the intentionally public
`results/mv03` and `results/mv04` evidence and total 5.68 MiB. The two files
whose names contain `private` are public ledgers describing repeat status or
opaque cache identity; they do not contain the private RDS payloads.

## Portability remediation

The audit found 23 new R analysis/validation scripts that pinned the local WSL
checkout with `setwd("/mnt/e/Repositories/Jonah/PH Pipeline Repo/scPHcompare")`.
They now derive the repository root from their own `Rscript --file` argument,
including R's `~+~` encoding for spaces. All 300 repository R scripts parse,
and `build_mv05au_corrected_consumer_prefreeze.R` succeeds when invoked from
`/tmp`, outside the repository. A repository test now rejects absolute
string-literal `setwd()` calls and verifies that bootstrap users decode paths
with spaces.

The old absolute paths in
`inst/extdata/inputs/metadata_MultiTissueAnalysis.csv` predate this branch and
remain a separate public-fixture cleanup item; they are not part of the
proposed diff and were not silently changed here.

## Logical slice feasibility

Six dependency-ordered scientific/engineering slices were measured. Four fit
under the 300-file display ceiling, but the two central existing-data and
robustness slices contain 521 and 632 changed files and approximately 122.54
and 119.34 MiB of blobs. Merely opening six stacked PRs would therefore not
solve reviewability.

The preferred publication design is:

1. keep source, tests, specifications, compact ledgers, manifests, and hashes
   in Git;
2. place bulky generated evidence in an immutable scientific archive or
   versioned release assets, with checksums and retrieval instructions kept in
   Git;
3. construct smaller dependency-ordered publication branches whose source
   diffs stay under GitHub's review boundaries;
4. run the four-host Rust candidate workflow only after those branches and
   their evidence identities are frozen.

Git LFS remains possible, especially for genuine binary inputs, but it would
not by itself make a 1,869-file PR reviewable and would add a collaborator and
bandwidth dependency. It is therefore not selected automatically.

## Decision and next gate

MV5-BH accepts the fast-forward reconciliation, privacy boundary, and script
portability remediation. It rejects a single-PR publication plan. The next
gate requires the owner's choice to externalize bulky generated evidence (the
recommended route) or deliberately adopt Git LFS/full in-repository evidence.
No history rewrite, evidence removal, authentication, push, or PR should occur
before that choice.
