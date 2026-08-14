# Safe Synchronization Procedure

## Purpose

Bring local `main` up to date without losing untracked planning/source material or masking divergence. This is a procedure, not authorization to execute synchronization.

## Current cached assessment

Local `main` is at `64ec4fd`; cached `origin/main` is two commits ahead at `3910b15`. The cached incoming content changes only `R/PH_Functions.R`, and there are no current untracked-path collisions. The cached ref is stale (last updated 2026-01-05), so these findings must be refreshed.

## Fresh verification update

The upstream state was refreshed on 2026-08-03 at 21:06:47 America/New_York. `origin/main` remained at `3910b15`; local `HEAD` is its direct ancestor, and exact, case-insensitive, and path-prefix collision checks all passed. See `FRESH_UPSTREAM_VERIFICATION_2026-08-03.md`. Steps 3–6 below are complete for this point in time but must still be repeated if material time or local state changes before integration.

## Preconditions

- [ ] Owner reviews P0-01 state capture and artifact inventory.
- [ ] Public untracked files have an explicit keep/track/relocate decision.
- [x] The two source PDFs are local-only and explicitly Git-ignored by owner decision.
- [ ] Confidential materials remain ignored and appear only in ignored-file checks.
- [ ] Protected artifacts verify against the manifests.
- [ ] No unexpected tracked modifications or deletions exist.
- [ ] A recoverable copy or intentional commit exists for work that must survive.

## Procedure

1. **Repeat the local snapshot.** Record branch, `HEAD`, short status including ignored files, and public untracked paths.
2. **Verify protected hashes.** Stop if an irreplaceable artifact changed without an expected explanation.
3. **Fetch only after approval.** Fetch/prune `origin`; record the new `origin/main` SHA and fetch evidence. Do not pull yet. Completed for the 2026-08-03 verification sprint.
4. **Recompute divergence.** Record `rev-list --left-right --count HEAD...origin/main` and the left/right commit list.
5. **Inspect every incoming commit.** Record author/date/subject/stat and review the full diff, with special attention to scientific behavior, package metadata, ignore rules, and paths shared with untracked work.
6. **Re-run collision analysis.** Intersect incoming changed paths with untracked paths. Also inspect rename/delete operations and case-only filename changes.
7. **Choose the integration mode.** If local `main` still has zero unique commits, tracked files are clean, and collisions are absent, use a fast-forward-only update. Otherwise stop and design a branch-specific merge/rebase plan.
8. **Update with fast-forward protection.** Use an operation that refuses non-fast-forward history. Do not auto-stash confidential or untracked material.
9. **Verify post-update state.** Record new `HEAD`, status, changed-file summary, protected hashes, ignore behavior, and whether untracked work remains intact.
10. **Run proportional checks.** For the freshly verified logging-only change, run the relevant tests plus package parse/check as feasible. Expand testing if a later fetch contains other changes.
11. **Record the outcome.** Append task/revision/evidence/results to `docs/plans/WORK_LOG.md`; update P0-02 and the decision log.

## Stop conditions

Stop without pulling or merging if any of the following occurs:

- current upstream contains an unexpected force update or local branch has unique commits;
- incoming paths collide with untracked files;
- protected hashes change unexpectedly;
- confidential files cease to be ignored;
- an incoming change alters scientific definitions or reference outputs without an approved audit;
- distribution status of source PDFs remains unresolved and the preservation method would expose them.

## Verified conditional outcome

The fresh fetch revealed no commits beyond `3910b15`. The source PDFs are now explicitly excluded from Git, and `example_run.r` is preserved untracked on hold pending a provenance/usefulness check. The repository is technically eligible for a fast-forward-only update after the remaining public planning work is intentionally preserved. That conclusion is conditional and is not authorization to integrate.
