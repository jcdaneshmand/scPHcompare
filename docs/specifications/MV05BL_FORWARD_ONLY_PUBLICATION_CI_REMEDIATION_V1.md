# MV5-BL forward-only publication CI remediation v1

## Purpose

MV5-BL repairs two GitHub Actions expression-context errors discovered only
after the accepted eight-slice publication stack was pushed. The correction is
CI infrastructure only. It preserves every original scientific candidate
commit in ancestry, leaves the stacked review boundaries unchanged, and uses
normal fast-forward pushes without force.

## Frozen pre-remediation state

- publication base: `3910b15ce6d3f88197847a8aba94b184e7d9c034`;
- draft PRs: `#111` through `#118` in P01-to-P08 order;
- original candidate commits: the eight `rehearsal_commit` values in
  `docs/publication/scphcompare-publication-candidates-v1.csv`;
- local candidate branches and their `origin/*` counterparts must both equal
  those original candidate commits before remediation;
- all eight PRs must remain draft and retain their frozen base/head pairs;
- the primary development worktree may contain only the previously untracked
  `example_run.r`, which must never be staged or enter a candidate.

## Root cause and exact correction

GitHub rejected both workflows before creating jobs because the `runner`
context is unavailable in job-level `env`. The official `actionlint` v1.7.12
diagnostic identified:

- `.github/workflows/r-package-ci.yml:21`;
- `.github/workflows/rust-accelerator-certification.yml:61`.

The only accepted source changes are:

1. `${{ runner.tool_cache }}/scphcompare-renv-library` becomes
   `${{ github.workspace }}/../scphcompare-renv-library` in the R workflow;
2. `${{ runner.tool_cache }}/scphcompare-renv-${{ matrix.target }}` becomes
   `${{ github.workspace }}/../scphcompare-renv-${{ matrix.target }}` in the
   P08 Rust workflow.

The replacement remains outside the checked-out package directory, is
cross-platform, and uses a context permitted at job-level `env`.

## Forward-only stack construction

1. Create retention refs for all eight original candidate heads.
2. Add one ordinary commit to P01 containing only the R workflow correction.
3. Merge each updated parent branch into P02 through P08 with a merge commit.
   The original candidate commit remains the first-parent ancestor and the
   corrected parent becomes the second parent.
4. Add one ordinary commit to P08 containing only the Rust workflow correction.
5. Never rebase, reset, amend, cherry-pick, replace a published ref, or
   force-push.
6. Push P01 through P08 normally only after local validation. Each remote old
   head must be an ancestor of the new local head.

## Validation contract

- Run checksum-verified official `actionlint` v1.7.12 on the corrected R and
  Rust workflow files and require zero diagnostics.
- Require each original candidate commit to remain an ancestor of its updated
  branch.
- Require each updated parent candidate to be an ancestor of its child.
- Require the P02-through-P07 PR tree deltas to remain identical to the
  original scientific slice deltas.
- Require P01 to differ from its original accepted tree only by the one-line R
  workflow correction.
- Require P08 to differ from its original accepted slice only by the one-line
  Rust workflow correction after accounting for the inherited R correction.
- Require no `results/**`, generated evidence CSV/CSV.GZ, PDF, local example,
  or other untracked artifact to enter any candidate.
- After normal pushes, independently verify all PR numbers, titles, draft
  status, base/head names, current head SHAs, and review-boundary bodies.
- Monitor the R package check on all eight PRs and the no-release guard plus
  four-platform Rust candidate matrix on P08.

## Documentation contract

The original promotion ledger remains immutable historical evidence. A new
MV5-BL remediation ledger records original and current commit/tree identities,
ancestry, exact source corrections, actionlint results, push/run/PR URLs, and
hosted check outcomes. PR bodies must distinguish the original deterministic
scientific boundary from the current CI-remediated head.

## Abort and closed actions

Abort on remote drift, unexpected PR topology, a non-fast-forward push, an
additional source delta, failed local validation, or an unrelated dirty file.
MV5-BL does not authorize merge, branch deletion, force-push, tag, release,
binary publication, DOI/Zenodo action, default change, Rust adoption, new
biological calculation, or manuscript claim.
