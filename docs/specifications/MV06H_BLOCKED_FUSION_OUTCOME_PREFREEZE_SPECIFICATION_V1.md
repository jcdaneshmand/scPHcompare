# MV6-H blocked fusion outcome prefreeze specification v1

## Document control

| Field | Frozen value |
|---|---|
| Contract | `mv06h_blocked_fusion_outcome_prefreeze_v1` |
| Date | 2026-08-15 |
| Accepted predictions | Complete corrected-root MV6-G corpus |
| Samples / studies / tissues | 90 / 15 / 5 |
| Folds / seeds | 15 LOSO folds / 5 fixed technical seeds |
| Methods | Nine nonselective MV6-G methods |
| Primary fusion | `fusion_gene_weight_050` |
| Required comparators | `cell_composite`, `gene_composite` |
| Outcome state during this sprint | Closed |

## 1. Purpose

MV6-H creates the durable prediction-before-label boundary for the complete
dual-view corpus and freezes the exact blocked outcome implementation. It does
not inspect metadata values or calculate outcomes during prefreeze. Once the
manifest, independent validation, and repeat evidence are committed, the
unchanged runner may evaluate the prespecified tissue endpoints at that exact
commit.

## 2. Prediction manifest

The manifest binds all 75 group identities and the SHA-256 and byte count of
every training-distance, scale, ranking, metric, and status artifact. It also
binds the accepted complete queue, group inventory, independent validation,
445-file immutable-resume ledger, canonical resource metrics, method/endpoint/
contrast/inference contracts, and every MV6-H implementation source.

The expected authoritative metadata SHA-256 is
`e2b6f1e7d7dc05ba8e80627aad4ea1677c3b718a844397acd06596d7e22a18d0`.
Only that hash is recorded during prefreeze: no metadata path is accepted and
no metadata file is opened.

## 3. Fixed predictions and endpoints

The 318,150 already accepted ranking rows are immutable predictions. No scale,
distance, rank, tie rule, weight, method, endpoint, tissue, or comparator may
change after labels open. The primary endpoint is cross-study tissue MRR; fixed
1-NN balanced accuracy is supportive. Seeds are averaged within each biological
sample, samples within tissue, and the five tissue means equally.

The F1 family contains exactly `fusion_gene_weight_050 - cell_composite` and
`fusion_gene_weight_050 - gene_composite` for MRR. Both point estimates must be
positive for any fusion-benefit interpretation. All nine methods receive point
estimates and intervals; descriptive methods cannot replace the primary.

## 4. Frozen inference

- 2,000 tissue-stratified held-out-study bootstrap replicates;
- bootstrap seed `20260815`;
- 9,999 paired held-out-study sign flips;
- randomization seed `20260816`;
- two-sided Monte Carlo p-values with plus-one correction;
- Holm adjustment across exactly the two F1 MRR contrasts;
- held-out biological sample blocked by study as the inferential unit;
- technical seeds never increase the independent sample count.

## 5. Durable label boundary

Outcome execution requires the current Git HEAD to equal an explicitly supplied
40-character lock commit. Every lock artifact must be tracked and clean. The
runner verifies and loads all predictions, then atomically writes a lock receipt
before hashing or parsing metadata. It rehashes all ranking artifacts and the
metadata after evaluation. Any drift, partial output, untracked lock, dirty lock,
or label-before-receipt state aborts.

## 6. Stop boundary

Prefreeze authorizes only construction, independent reconstruction, byte-repeat
testing, and durable commit of the label-closed prediction lock. Outcome
execution requires a later sprint at that exact commit. Advanced or learned
fusion, clustering, package defaults, release, external data, and manuscript
claims remain closed regardless of the blocked outcome result.
