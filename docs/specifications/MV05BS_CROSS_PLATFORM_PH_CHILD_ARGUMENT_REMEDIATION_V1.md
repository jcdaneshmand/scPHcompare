# MV5-BS cross-platform PH child-argument remediation contract v1

Date frozen: 2026-08-13
Status at freeze: prospective P08 forward fix; draft and unmerged

## Triggering evidence

P08 Rust candidate run `31747769480` proved that the MV5-BQ Windows linker
correction worked: the two clean optimized DLLs were byte-identical, the native
C ABI harness passed, and all five direct R analytical/corruption/fallback
fixtures passed. The accelerator-absent `R CMD check` then passed 1,345 test
assertions but failed nine assertions in the pre-existing real-PH subprocess
test because the nominal success child returned `ph_child_error` before writing
a diagram.

The accepted child launcher constructs an `Rscript -e` expression by inserting
`dataset_file` and `pd_file` between R string quotes. On Windows, ordinary paths
such as `C:\Users\...` therefore place `\U` and other backslash escapes in R
source text. Parsing fails before the expression's `tryCatch`, so the valid
numeric fixture never reaches `ripserr`. The failure is independent of the Rust
accelerator, landscape calculation, dependency library, and PH fixture values.

## Prospective correction

Keep the child expression static and pass exactly four trailing `Rscript`
arguments through `processx`:

1. dataset RDS path;
2. output persistence-diagram RDS path;
3. maximum homology dimension; and
4. filtration threshold.

The child reads the unmodified values with `commandArgs(trailingOnly = TRUE)`,
converts only the two numeric parameters, reads the same dataset, calls the same
`ripserr::vietoris_rips()` function with the same arguments, and writes the
same diagram destination. R's official `Rscript` interface defines trailing
arguments for both script files and expressions supplied with `-e`; `processx`
passes each vector element without source-code interpolation.

## Placement and boundaries

- Add this as a separate P08 forward commit because hosted cross-platform
  certification discovered it after the exact accepted P08 slice was already
  published. Preserve that accepted slice and every earlier remediation in
  ancestry; do not rewrite or force-push.
- Change only the child command construction in `R/PH_Functions.R` plus the
  audit/publication records required to trace it.
- Do not change input orientation, topology view, persistence dimension,
  threshold selection, coefficient field, `ripserr` version, landscapes,
  distances, clustering, tolerances, retry/timeout behavior, provenance
  schema, dependencies, fixtures, APIs, or defaults.
- Do not skip or weaken the Windows real-PH test or any Rust/ABI/R/package gate.
- Do not merge, tag, release, publish a binary, change the default branch,
  perform a DOI/Zenodo action, or make a biological/manuscript claim.

## Acceptance gate

1. R parses the modified source and the package's complete existing tests pass
   on the final P08 hosted R check.
2. P08's final Rust candidate run passes no-release guard plus Linux x86-64,
   Windows x86-64, macOS ARM64, and macOS x86-64 in full.
3. Windows specifically passes exact clean-DLL identity, native C ABI, all five
   direct fixtures, the real-PH subprocess test, and both accelerator-absent
   and accelerator-present `R CMD check` modes.
4. P08's accepted original scientific candidate remains an ancestor and its
   pre-remediation slice remains available at the immutable retention ref.
5. The final P08 head/tree, run IDs, outcomes, body, and unchanged `main` are
   recorded and independently rechecked before publication-goal closure.

## Nonclaims

This correction transports existing PH child inputs safely across operating
systems. It does not alter a persistence diagram or landscape definition,
constitute new scientific evidence, certify a public binary release, or
authorize production Rust adoption.
