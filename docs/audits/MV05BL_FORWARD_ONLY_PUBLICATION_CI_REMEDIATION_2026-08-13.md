# MV5-BL forward-only publication CI remediation

## Outcome

The accepted eight-slice publication stack is live as draft PRs #111 through
#118. GitHub initially rejected both new workflow files before job creation
because `runner.tool_cache` is not available in job-level `env`. MV5-BL applied
the two prospectively frozen one-line context corrections, preserved all eight
original scientific candidate commits in ancestry, and atomically advanced all
eight remote branches by ordinary fast-forward pushes without force.

## Root cause and remediation

Checksum-verified official `actionlint` v1.7.12 identified the invalid context
at `.github/workflows/r-package-ci.yml:21` and
`.github/workflows/rust-accelerator-certification.yml:61`. The R library paths
now use `${{ github.workspace }}/../...`, which is legal at job level, remains
outside the package checkout, and is cross-platform. Both corrected workflow
files pass `actionlint` with zero diagnostics.

P01 received one ordinary R-workflow commit. P02 through P08 received the
corrected parent through ordinary merge commits. P08 then received a separate
one-line Rust-workflow commit. Retention refs under
`refs/codex/remediation/mv05bl/pre/P01..P08` preserve the eight pre-remediation
heads in addition to the existing two rehearsal ref sets.

The first executing Rust matrix exposed a separate PowerShell portability
defect after all macOS ARM64 Rust and fixture gates had passed: bare `R CMD`
resolved to PowerShell's `R`/`Invoke-History` alias instead of the installed R
application. The prospectively frozen MV5-BM correction makes that P08-only
step resolve the external R application explicitly. It was added and pushed as
another ordinary forward commit; no earlier commit or branch was rewritten.

The first executing Ubuntu R checks then consistently exposed missing hosted
system dependencies, with `igraph` unable to load `libglpk.so.40`. MV5-BN added
the five prerequisites that `renv` itself reported for the frozen lockfile to
the Ubuntu R workflow and the Linux Rust candidate. The R-workflow correction
was propagated through the stack with ordinary merge commits, and the P08 Rust
correction remained a separate forward commit. A second atomic non-force push
advanced all eight branches.

The replacement macOS ARM64 Rust candidate then passed restoration, all Rust
gates, the native harness, and all five fixtures before its nested `R CMD check`
could not see the verified renv library. MV5-BO bound `R_LIBS_USER` to that
already-restored library for the nested package checks. This was a one-line
P08-only forward commit and normal non-force push; P01 through P07 continued
their hosted runs without interruption.

P07 then proved the general R workflow had the same post-restore visibility
boundary: restoration and package build passed, but `R CMD check` retained the
setup action's temporary library path. MV5-BP publishes the verified
`RENV_PATHS_LIBRARY` through `GITHUB_ENV` after successful restoration. The
one-line P01 correction was propagated through all children with ordinary
merge commits, scientifically revalidated, and atomically pushed without force.

The next P08 Windows candidate passed locked restoration, Rust identity,
formatting, all four unit tests, strict Clippy, and its first optimized build,
but two clean MSVC DLL builds were not byte-identical. MV5-BQ preserved the
strict raw SHA-256 equality requirement and made the Windows linker inputs
reproducible with `/Brepro` plus `/PDBALTPATH:%_PDB%`. This workflow-only
change advanced P08 by another ordinary fast-forward commit; it did not change
the Rust kernel, R source, fixtures, tolerances, dependency locks, or any
scientific slice.

P02 then proved that `R_LIBS_USER` alone did not stop the repository
`.Rprofile` from reactivating renv around later commands. Its restore report
was exact at 265/265 records, but `R CMD check` searched the platform-suffixed
project path computed by `renv/activate.R` rather than the explicit directory
that had just been restored and verified. MV5-BR exports that directory as
both `R_LIBS` and `R_LIBS_USER` and disables redundant autoloading only after
successful verification. The correction was added at P01, propagated through
P08 by ordinary merge commits, added separately to P08's Rust workflow, and
atomically pushed without force.

The corrected Windows Rust candidate then passed exact DLL byte identity,
native C ABI execution, dependency normalization, and all direct analytical,
corruption, and fallback fixtures. Its accelerator-absent full package check
exposed a pre-existing cross-platform launcher defect: Windows temporary paths
were inserted into `Rscript -e` source, so `\U` was parsed as an R escape before
the real-PH child reached `ripserr`. MV5-BS passes the two paths and two numeric
parameters as process arguments read through `commandArgs()` while preserving
the exact PH call. The WSL real-child regression passed 27/27 before this
separate P08 forward commit was pushed normally without force.

## Scientific-boundary validation

- Every original candidate commit is an ancestor of its current branch head.
- Every updated parent branch is an ancestor of its updated child.
- P02 through P07 retain patch-identical accepted scientific slices.
- P08 retains its patch-identical accepted scientific slice before the separate
  Rust workflow correction.
- P01 differs from its original tree only in the one-line R workflow path.
- P08's inherited merge differs from the original P08 tree only in that R line,
  and its later commits differ from that merge only in the audited Rust
  workflow path and portable R-application invocation.
- No landscape, topology, data, analysis, result, PDF, generated evidence, or
  local example file changed or entered a candidate.

## Remote publication verification

GitHub account `jcdaneshmand` had `ADMIN` access and workflow-enabled repository
write scope. Immediately before push, all eight remote heads still equaled the
frozen candidates and `main` still equaled
`3910b15ce6d3f88197847a8aba94b184e7d9c034`. One atomic, non-force push advanced
all eight candidate branches. Independent post-push inspection proved:

- PRs #111 through #118 remain drafts;
- titles, base branches, head branches, and current head SHAs are exact;
- all eight remote bodies are byte-equal to the committed local templates;
- each body distinguishes the original deterministic scientific boundary from
  its current CI-remediated head and tree.

The exact per-slice identities, PR URLs, and hosted run URLs are recorded in
`docs/publication/mv05bl-remote-publication-ledger-v1.csv`.

## Hosted checks

GitHub now parses both workflows and creates real jobs. Superseded runs caused
by the atomic parent/head updates were cancelled by the workflow concurrency
policy; those cancellations are not candidate failures. The final retained
general R-package runs succeeded for P01 through P08. Each restored and
verified the exact 265-record lock, passed `R CMD check`, and completed the
realistic H0 and H1 fixture step.

P08's final Rust candidate run `31750611952` succeeded overall. Its
`no-release-guard`, `linux-x86-64`, `windows-x86-64`, `macos-arm64`, and
`macos-x86-64` jobs all succeeded. Every native row passed pinned toolchain
identity, format/unit/strict-Clippy checks, two byte-identical clean release
builds, native ABI execution, dependency normalization, all five direct R
analytical/corruption/fallback fixtures, accelerator-absent and
accelerator-present full R package checks, manifest generation, and short-lived
candidate-evidence upload. The separate final P08 general R run `31750612007`
also succeeded. Exact run URLs and conclusions are frozen in the companion
ledger.

Earlier failed runs remain auditable diagnostics rather than retained
candidate evidence: workflow context parsing, the PowerShell `R` alias, missing
Ubuntu system libraries, post-restore library visibility, nondeterministic
default MSVC link inputs, and the Windows path-in-source PH child launcher were
each isolated before a prospective forward-only remediation. No scientific
parameter or landscape definition was changed in those remediations.

## Closed actions

All PRs remain draft. No merge, retarget, branch deletion, force-push, tag,
release, binary publication, Rust default/adoption, DOI/Zenodo action, default
branch change, new biological calculation, or manuscript claim occurred.
