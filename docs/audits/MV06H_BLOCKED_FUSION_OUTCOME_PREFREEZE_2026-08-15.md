# MV6-H blocked fusion outcome prefreeze audit

## Disposition

`prediction_lock_pass_outcomes_not_yet_executed`

The complete corrected-root MV6-G prediction corpus is now bound by a durable,
label-closed MV6-H prediction manifest. Independent validation passes 13/13
categories and a clean rebuild reproduces all nine public lock artifacts
byte-for-byte. Metadata values, biological outcomes, advanced fusion,
clustering, and claims were not opened or computed in this sprint.

## Frozen prediction boundary

- implementation commit: `f2b8da4c1f8fc5d5f2d62a8088c9a24b260bf74a`;
- prediction-manifest root:
  `c752408f064800d4539b475455e5bca7e9fcb010019e90a1fc261824c865f4fd`;
- group manifest:
  `21f8371e9c9d40072cde8b81cad2b36eca5109d0da38b6c5fc8d32ef4e319dc2`;
- public prediction-lock file:
  `21e445b4167ec7b870d34da127301956afe2dcaacdcb2dd461b2a2a822566605`;
- authoritative metadata expected SHA-256:
  `e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`.

The metadata file path was not supplied to the builder or validator. The
external source appears only as its expected hash with
`opened_as_outcome_source = FALSE`.

## Bound corpus and analysis

- 75 fold-seed groups, 15 folds, five technical seeds;
- 375 group artifacts independently rehashed;
- 318,150 ranking rows and 4,050 query-method rank sequences;
- nine fixed methods with `fusion_gene_weight_050` as the sole primary fusion;
- required comparators `cell_composite` and `gene_composite`;
- cross-study tissue MRR primary and fixed 1-NN balanced accuracy supportive;
- 2,000 tissue-stratified study-block bootstraps, seed `20260815`;
- 9,999 paired held-out-study sign flips, seed `20260816`;
- Holm adjustment across exactly two primary MRR contrasts;
- seeds averaged within biological samples rather than treated as independent.

All method, endpoint, contrast, inference, source, implementation, and firewall
registries are hash-bound by the public lock.

## Validation

The independent validator did not call the production builder. It reconstructed
the complete group axis, all file hashes, every rank sequence under ascending
distance plus canonical training-sample ID, the nine-method panel, two-endpoint
registry, two-contrast F1 family, inference seeds, and ten closed firewall
rules. It passed 13/13 categories. A second clean build reproduced 9/9 public
lock artifacts byte-for-byte. Focused tests passed 21/21; the complete package
suite passed 1,580 tests with four established skips and no failure or warning.

## Bounded pre-lock correction

The first builder invocation stopped before lock publication because Ubuntu
Git could not resolve the linked worktree's absolute Windows `.git` pointer.
A process-local `GIT_DIR`/`GIT_WORK_TREE` binding resolved that environment
boundary without changing global configuration. The next attempt traversed all
groups but stopped before lock publication because `file.path()` had discarded
the logical names of the eight-file hash vector. The partial label-closed
manifests were quarantined under ignored `tmp/`. Name preservation was corrected
in the builder, runtime verifier, and independent validator at `f2b8da4`; the
focused suite passed again before the clean build. Neither failed attempt read
metadata or computed an outcome.

## Durable execution boundary

The outcome runner requires:

1. the current Git HEAD to equal the explicitly supplied durable lock commit;
2. all lock artifacts to be tracked and clean;
3. all 75 prediction groups to hash-match before label access;
4. an atomic lock receipt to be written before metadata hashing or parsing;
5. metadata and all prediction hashes to remain unchanged afterward.

## Next gate

After this audit and the lock evidence are committed, run the fixed MV6-H
blocked fusion outcome evaluation at that exact commit. Independently validate
the label join, endpoints, seed/sample/tissue aggregation, bootstrap and sign
matrices, two contrasts, Holm adjustment, and post-label source hashes before
deciding G-MV6 or beginning MV-07 synthesis.
