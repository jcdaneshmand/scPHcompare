# MV6-F stage-two serial 12-GiB policy prefreeze

| Field | Frozen value |
|---|---|
| Date | 2026-08-14 |
| Parent prefreeze commit | `eaa3001` |
| Scope | The unchanged 74-row stage-two queue, executed serially |
| Queue root | `f5471633e21d229eeabecadf12989dece2a3a7ab5b5d09f4584b0c3b6410bb5d` |
| Scientific implementation root | `5a1258e87eb30c367648daa4bb02d1aec1f4b40dd3799a91fbb7067e9558d292` |
| Scientific runner SHA-256 | `1128523f948f5adf29f1165300cd71f53cf93eb1335d0778cdf5160e9b404ca3` |
| Rust library SHA-256 | `51d3fca4808c6620c46c7b1100ddfad4753af02d37f87ec3eb278008846b160d` |
| Policy SHA-256 | `7a9751d32ae12e62412c2e8b147e940e2b54012f9204383a05af4ccd8201ace4` |
| Driver / monitor SHA-256 | `b39320a70f40f24fb54f92fb18b36382a86993b97940a9a0ccdf5bd91373cc9e` / `80e65e1db00c2c73f6ebf18629b39ba0e7ff45188b5c3f6c2dcef2cb6cdc0007` |
| Resource policy | One worker; 12 GiB RSS/group; 1,800 s/group; 10 GiB private storage; no automatic retry |
| Scientific firewall | Labels closed; fusion, clustering, outcomes, and claims remain zero/closed |

## Evidence for the amendment

The original serial monitor stopped execution-order group 2 at 8,747,204,608
bytes, just above its prospective 8-GiB per-group ceiling. It published no
group directory, launched no later group, and preserved the failure metric.

The separately prefrozen diagnosis then ran that exact group with only the
process-tree ceiling raised to 12 GiB. It completed in 423.049569 seconds at a
9,575,215,104-byte peak, returned exit status zero, and published a complete
20,628,449-byte group directory that passes the frozen atomic validator. This
locates the first stop in the execution resource ceiling rather than in the
scientific calculation or numerical result.

## Frozen execution policy

The policy contains the exact 74 stage-two queue rows in execution order. Each
row binds the parent failure and successful diagnostic, queue root,
implementation root, scientific runner, Rust library, exception monitor, and
serial driver. The driver validates pre-existing completed groups, writes a
one-row authorization for each missing group, invokes the already frozen
exception monitor, and stops immediately on any nonzero exit, resource breach,
missing metric, or invalid atomic directory. It performs no automatic retry.

This is an execution-policy amendment, not a scientific-root revision. The
cell and gene views, H0 and H1 diagrams, exact all-active-level landscapes,
pair inventory, folds, seeds, panel, and inputs are unchanged. The successful
diagnostic group is reused only after full frozen-directory validation.

## Admission boundary

Serial stage-two production may proceed under this policy. If any group exceeds
12 GiB, exceeds 1,800 seconds, breaches the private-storage ceiling, fails
numerically, or cannot publish and validate its atomic directory, execution
must stop before the next group. Complete-production validation and a full
immutable-resume check are mandatory before any fusion, clustering, label
opening, result synthesis, or claim.

Machine-readable evidence is in
`docs/audits/mv06f-stage2-production-evidence/`. Private outputs and logs remain
excluded from Git.
