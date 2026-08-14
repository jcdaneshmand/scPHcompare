# MV6-F stage-1 validation and admission audit

| Field | Result |
|---|---|
| Date | 2026-08-14 |
| Disposition | `pass_stage2_label_closed_only` |
| Frozen group | `mv06f_group_v1:fbed6ad04f8243313ed439ecb5f29ddd43326a478d9b60fb21ff84be70b6ebf1` |
| Fold / seed | `large_loso_v1:SRA779509` / `20260807` |
| Sample roles | 25 held-out / 65 training |
| Queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Implementation root | `599074b3cd078cf27eb4a85148eb1df2ce3f84a5bdfd3160617b80a78f78c05e` |
| Rust SHA-256 | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Outcome state | Closed; no fusion, clustering, or biological-outcome jobs |

## Scope and scientific contract

Stage 1 executed exactly the prospectively selected maximum MV6-F group. It
created 180 typed cell/gene PH records, 360 H0/H1 diagram records, 1,625
biological query-to-training pairs, and 6,500 separate view/dimension landscape
rows. The calculation retained the approved dissertation-aligned landscape
definition:

- H0 and H1 are calculated and stored separately;
- only finite positive-persistence intervals enter a landscape;
- every active consecutive landscape level is included, with zero padding;
- squared L2 distance is exact or independently error-controlled;
- no universal level cap or uniform evaluation grid is applied; and
- grouped/streamed execution changes resource use, not the mathematical object.

The 12 independent Persim rows exercised 25 through 1,856 active levels and up
to 3,418,015 critical points. Every row reported `all_active_levels = TRUE`,
`exact = TRUE`, and `level_cap_applied = FALSE`.

## Resource gate

| Run | Elapsed | Peak process-tree RSS | Private root | Result |
|---|---:|---:|---:|---|
| Primary | 344.799 s | 3,045,322,752 B | 10,346,588 B | Pass |
| Clean repeat | 356.902 s | 3,106,680,832 B | 10,346,587 B | Pass |
| Frozen cap | 1,200 s | 6,442,450,944 B | 10,737,418,240 B | — |

Both runs completed with exit status zero and a complete atomic group
directory. The external monitor measured the complete process tree, not only
the parent R process. No cap was approached closely enough to require a revised
envelope.

## Determinism and resume

The three scientific artifacts were byte-identical across a clean repeat:

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `diagrams.rds` | 7,484,632 | `d9066b21809e0e647b69c50ef3d5f74aa9d8639ee5f37f2b2762b299be2a59ec` |
| `diagram-manifest.csv` | 63,070 | `89d09fd92c9384b0c00e504366b0df5c482590fced8de3f10ec0787b27353921` |
| `distances.csv` | 2,797,261 | `776bf94392cb6c37a0fdc9243032f1f007ea580d497252bc4cb172ece1dd6a4c` |

Run-local timing/status files are intentionally not required to match between
clean executions. A subsequent resume reused the validated primary directory
and left all five primary files unchanged in SHA-256, byte count, and
modification time; the machine-readable identity ledger is
`mv06f-stage1-resume.csv`.

## Independent landscape validation

The oracle selection was fixed without outcomes at three interval-depth strata
for each of the four cell/gene by H0/H1 components: four minimum, four median,
and four maximum rows. It is balanced at six cell/six gene and six H0/six H1.

- R passed 12/12 rows. Ten used the exact breakpoint stream and two deep rows
  used partitioned adaptive QUADPACK. Maximum absolute Rust-versus-R error was
  `5.00222085975111e-12`, below the largest certified tolerance
  `3.87846524909276e-07`.
- Grouped Persim 0.3.8 passed 12/12 rows against both Rust and R. Maximum error
  was `8.185452315956354e-12` versus Rust and `3.183231456205249e-12` versus R.
- The clean repeat passed 3/3 scientific artifact comparisons.

The validation environment was Ubuntu/WSL with R 4.4.1 and a dedicated Python
3.12.11 environment containing Persim 0.3.8, NumPy 2.2.6, and SciPy 1.15.2.
The R and Persim gates were rerun successfully after the post-merge P1
integrity remediation was applied locally.

An initial Persim staging attempt exposed duplicate interval rows when one
diagram appeared in more than one selected pair. The R staging script was
corrected to de-duplicate by immutable diagram identity and interval values;
the corrected R gate and independent Persim gate then passed completely. This
was an oracle-staging defect, not a production-distance discrepancy.

## Post-merge review and CI gate

The four late P1 review findings from merged PR 112 are remediated by draft PR
120. The changes add actual frozen-file hashing, live process-cap enforcement,
profile-derived scientific eligibility, and shared-PCA standardization
identity. The complete local R suite passed in 164.9 seconds with only two
established optional skips. GitHub Actions run 31818846372 passed in 34m48s,
Codex found no further major issue on the remediation commit, and all four
original review threads are resolved. These integrity changes do not alter the
landscape definition or the frozen MV6-F numerical implementation root.

## Decision and next boundary

Stage 1 satisfies the prospective resource, atomicity, repeat, R-oracle,
Persim-oracle, and immutable-resume gates. The remaining 74 frozen MV6-F groups
may therefore proceed as label-closed stage 2 under the existing queue,
implementation, Rust, resource, and abort contracts.

This decision does **not** authorize fusion evaluation, clustering, outcome
access, public Rust adoption, a release, or manuscript claims. Stage 2 must be
validated and summarized before blocked fusion begins. Draft PR 120 remains
unmerged pending the project owner's public merge decision.

## Public evidence

- `docs/audits/mv06f-stage1-evidence/mv06f-stage1-resource.csv`
- `docs/audits/mv06f-stage1-evidence/mv06f-stage1-repeat-resource.csv`
- `docs/audits/mv06f-stage1-evidence/mv06f-stage1-scientific-repeat.csv`
- `docs/audits/mv06f-stage1-evidence/mv06f-stage1-resume.csv`
- `docs/audits/mv06f-stage1-evidence/mv06f-stage1-oracle-selection.csv`
- `docs/audits/mv06f-stage1-evidence/mv06f-stage1-r-oracles.csv`
- `docs/audits/mv06f-stage1-evidence/mv06f-stage1-persim-oracles.csv`
- `scripts/validate_mv06f_stage1_r_oracles.R`
- `scripts/validate_mv06f_stage1_persim.py`

Private diagrams, interval staging, logs, caches, the Rust binary, and source
PDFs remain excluded from Git.
