# MV5-BH pre-push branch reconciliation contract v1

## Purpose

MV5-BH decides whether the accumulated local modernization branch is safe and
reviewable enough to publish for hosted candidate certification. It is a
local-only publication gate. It does not authorize a push, pull request,
release, Rust adoption, R-default change, added-seed calculation, partition,
or scientific claim.

## Frozen source revisions

- Remote base after a fresh fetch: `3910b15ce6d3f88197847a8aba94b184e7d9c034`
  (`origin/main`).
- Audited local source head before MV5-BH remediation:
  `af81621fe51aac5443c2db44cee34ffc8d7e8781`.
- The remote base must be an ancestor of the audited source head. Any future
  remote divergence requires a new reconciliation audit.

## Required checks

1. Refresh `origin` without pushing and record base, head, merge base, and
   ahead/behind counts.
2. Inspect every proposed commit for forbidden paths, private source material,
   reviewer correspondence, credential signatures, and local checkout paths.
3. Inventory changed files, lines, blob bytes, binary evidence, and largest
   objects against current GitHub repository and diff guidance.
4. Test logical publication slices at stable scientific/engineering
   boundaries; do not claim that a slice is reviewable merely because it is
   below the hard push limit.
5. Remove newly introduced absolute checkout roots and verify scripts from a
   working directory outside the repository.
6. Preserve `example_run.r` as untracked and preserve all local PDFs and
   confidential correspondence outside Git history.
7. Record a continuation decision before authentication, push, PR creation,
   or hosted CI.

## Abort rules

Stop before publication when any of the following occurs:

- `origin/main` is not the recorded merge base;
- a forbidden or confidential path appears anywhere in proposed history;
- a credential signature appears anywhere in proposed history;
- a proposed regular-Git object exceeds GitHub's enforced object limit;
- the complete change cannot be meaningfully reviewed under GitHub's diff
  limits without an explicit artifact and PR-partition policy;
- publication would require silently dropping auditable evidence;
- GitHub authentication, release publication, or a production-default change
  would be needed to complete a local check.

## Acceptance rule

MV5-BH may accept local reconciliation while withholding publication. A clean
fast-forward relationship and privacy pass are necessary but not sufficient:
reviewability and generated-evidence storage must also have an owner-approved
policy. Until then the branch remains local.
