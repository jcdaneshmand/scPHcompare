# MV5-BQ Windows MSVC reproducible-link remediation contract v1

Date frozen: 2026-08-13
Status at freeze: prospective, forward-only CI correction; P08 remains draft and unmerged

## Triggering evidence

Hosted P08 Rust candidate run `31745846731`, job `94600203305`, passed the
exact locked R restore, Rust 1.97.1 host-identity check, formatting, all four
Rust unit tests, strict Clippy, and its first optimized Windows build. The
second clean `x86_64-pc-windows-msvc` build produced a different raw DLL hash,
so the workflow stopped before the ABI, R-fixture, and package-check gates.

The failure is confined to the Windows MSVC link-product identity assertion.
It does not authorize a weaker normalized-binary comparison, removal of the
second clean build, a numerical tolerance change, or any scientific/source
change.

## Prospective correction

For the Windows matrix row only, set `RUSTFLAGS` before both clean builds to:

```text
-C link-arg=/Brepro -C link-arg=/PDBALTPATH:%_PDB%
```

`/Brepro` requests a reproducible MSVC link product. `/PDBALTPATH:%_PDB%`
prevents the absolute hosted-runner PDB path from being embedded in the DLL,
while retaining the actual PDB filename. Non-Windows rows receive no new link
arguments.

The acceptance criterion remains exact SHA-256 equality of two independently
cleaned release DLL builds on the same pinned Windows runner/toolchain. The
second build must still pass before the candidate is copied or any downstream
ABI, analytical, fallback, package-check, manifest, or upload step runs.

## Frozen boundaries

- Change only the build-excluded Rust candidate-certification workflow and
  audit/publication records needed to trace this correction.
- Do not change Rust or R source, landscape definitions, inputs, fixtures,
  tolerances, dependency locks, package contents, APIs, defaults, or evidence.
- Preserve the exact four-host matrix, Rust 1.97.1, R 4.4.1, locked Cargo and R
  dependencies, read-only permissions, no-release guard, seven-day candidate
  artifacts, and all downstream certification gates.
- Advance P08 by ordinary fast-forward commit and push only. Do not force-push,
  rewrite, merge, tag, release, publish a binary, change the default branch, or
  make a manuscript claim.

## Acceptance gate

1. The workflow parses and passes the pinned local `actionlint` validator.
2. Static inspection proves Windows alone receives both exact linker flags and
   that the raw two-build SHA-256 equality assertion remains unchanged.
3. P08's scientific slice remains patch-identical to its retained pre-publication
   ref when workflow files are excluded.
4. A new hosted P08 Rust run passes no-release guard plus all four native jobs,
   including exact Windows DLL byte identity and every downstream gate.
5. P08's general R-package check also passes for the final head.
6. The final P08 head/tree, run IDs, outcomes, PR metadata, and unchanged `main`
   are recorded in the publication ledger before this goal is closed.

## Authority and nonclaims

This correction operationalizes the already accepted two-clean-build contract;
it does not certify a release, satisfy the separate glibc 2.17 release baseline,
authorize production Rust adoption, or alter the canonical R implementation.
Only a completely green final hosted matrix can close this remediation.
