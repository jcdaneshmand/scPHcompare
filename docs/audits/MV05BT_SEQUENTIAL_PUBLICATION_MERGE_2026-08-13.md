# MV5-BT sequential publication merge

Date: 2026-08-13 (America/New_York)
Status: accepted; eight-slice merge and final-main validation complete

## Authorization

The project owner approved D1-D6 in
`docs/publication/MERGE_READINESS_DECISION_DOSSIER_2026-08-13.md` and
explicitly authorized the sequential merge procedure in D7:

> I approve D1-D6 and authorize the sequential merge procedure in D7.

That approval authorized the dependency-ordered merge of draft PRs #111 through
#118, a separate P09 documentation PR, and no release, tag, DOI, binary,
default-engine, calculation, or manuscript action.

## Outcome

PRs #111 through #118 were merged into `main` in order using GitHub merge
commits. Each child was retargeted to `main` only after its parent merge
completed. Before every merge:

- the head SHA was required to equal the frozen current candidate;
- the new `main` commit and tree were required to equal the preceding accepted
  slice;
- every current-base-to-head path/status row was compared with the 824-row
  merge-readiness inventory;
- the retained successful R run was bound to the unchanged head;
- P08 additionally required the retained successful Rust run and each of its
  five named successful jobs;
- the candidate branch was verified after merge and deliberately preserved.

All eight gates passed. No conflict, unexpected path, head movement, tree
difference, or failing retained check occurred.

## Merge ledger

| Slice | PR | Merge commit | Accepted resulting tree |
|---|---:|---|---|
| P01 | #111 | `7a6c7ef03440664ad9e9b84360b97a0d284263d5` | `bc581aa0ee1c4a6ddff50a6bcd726f322dbb729f` |
| P02 | #112 | `6941ae386296beeb184178e73e743cdadc5eda99` | `ed60fc692e5680498eea23c7b93e69809a353f12` |
| P03 | #113 | `f136eeff30109a09b8b39fa788f75923e704dc21` | `e2e3bc85c623bb5cb478359695b10d44d9c5dd64` |
| P04 | #114 | `adf2aeb59aeeb603d0bb6eae941ac0d580d8e7fb` | `b015eff511e57991d233ec77872714a86313df26` |
| P05 | #115 | `e0d36e4240843d25dd7c90506df10eb370f4b46c` | `7d8682b03e5d47cf9609d49bf5140a1e31303ab0` |
| P06 | #116 | `3d8b1b44985ba20b3c6c67fea3eea089e70e7e6a` | `675f226e6abd3e0bbd0e936044db408b94015921` |
| P07 | #117 | `4c5e0530ddda588729c4ff8fd6222276bc602626` | `7a7ec8cdae06394f3411ade51fa2a5bc8d3ae357` |
| P08 | #118 | `5e5ecb7d955579561c53246d89250b1978c67d0d` | `1fbf3350a4f15fe0209cc0a9eb1c70e8ef3fd7f6` |

The exact per-slice identities, timestamps, retained runs, and gate flags are in
`docs/publication/mv05bt-sequential-merge-ledger-v1.csv`.

## Scientific and landscape assurance

P07 received an additional check that all 38 landscape-key paths in its
accepted inventory remained present and exact after retargeting. The resulting
P07 and final P08 trees preserve the accepted landscape contract:

- every finite positive-persistence interval;
- every consecutive active level with zero padding;
- separate H0 and H1 primary results;
- exact or independently validated error-controlled squared-L2 integration;
- no universal uniform grid or level cap;
- secondary H0/H1 combination only with its components;
- R as the canonical/default engine.

P08 retained the audited Windows PH-child launcher commit
`98b3cc1f6b835d60bf92fdcd7cb273841a48855e`, R run
`31750612007`, and Rust run `31750611952`. The Rust run's
`no-release-guard`, `linux-x86-64`, `windows-x86-64`, `macos-arm64`,
and `macos-x86-64` jobs were each completed successfully before merge.
Rust remains optional and default-off; no candidate binary was published.

## Post-merge main CI

The R workflow runs on every push to `main`. Its concurrency policy cancelled
the seven intermediate merge-push runs as each newer merge arrived. Those
cancellations are expected and are not scientific failures.

The final retained main run is:

- run: `31764726529`
- head: `5e5ecb7d955579561c53246d89250b1978c67d0d`
- URL: <https://github.com/jcdaneshmand/scPHcompare/actions/runs/31764726529>
- initial status at audit write: `in_progress`
- final status: `completed`
- final conclusion: `SUCCESS`
- duration: `26m34s`

P09 subsequently passed its retained R-package check at run `31766356585`, was
merged as PR #119 at `d0192d35a4ab52006aa83b0ad3b0ad6a19f066cb`,
and produced exact tree `188e24a7b0c588b1926a2d99f43ba9e48e938053`.
The exact merged-main run `31768108363` completed successfully in 22m19s,
including the privacy guard, exact dependency restore, package check, realistic
H0/H1 fixtures, and synthetic audit-evidence upload. All nine publication
branches remained preserved; no release, tag, DOI, binary publication, or
default change followed.

## Remote state and closed actions

Immediately after P08:

- live `main` is `5e5ecb7d955579561c53246d89250b1978c67d0d`;
- its tree is exactly
  `1fbf3350a4f15fe0209cc0a9eb1c70e8ef3fd7f6`;
- all eight PRs report `MERGED`, `base=main`, and their frozen head SHA;
- all eight candidate branches still exist at their frozen current heads;
- the default branch remains `main`;
- the repository remains public;
- the sole existing tag/release remains pre-existing `v0.9.9`, pointing to
  commit `82f75a8610c37f5706a8ac8ef5fa38d4dd999cdf` and published in 2025.

No candidate branch deletion, force push, squash/rebase merge, new tag, new
release, Rust binary publication, Rust adoption/default, DOI/Zenodo action,
default-branch change, biological calculation, data change, or manuscript claim
occurred.

## P09 boundary

P09 is documentation-only. It may publish:

- the merge-readiness owner decision dossier and exact file inventory;
- the remote publication/CI remediation audit and run ledger;
- the remediation specifications needed to interpret the final candidate heads;
- this sequential-merge audit and ledger;
- a discoverability update to `docs/publication/README.md`.

P09 must not modify R, Rust, workflows, dependency locks, scientific
specifications, results, generated evidence, PDFs, reviewer correspondence, or
`example_run.r`. It remains a draft until its documentation diff and CI are
verified.
