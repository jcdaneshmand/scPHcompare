# MV5-Q label-closed clustering-artifact production

Date: 2026-08-10  
Accepted distance-production base: `eb5c812`  
Final prospective prefreeze: `b4bd052`  
Accepted private production root: `tmp/mv05q/production-v2`  
Outcome-label state: closed

## Outcome

MV5-Q completed the label-closed clustering-artifact layer frozen by MV5-N and
made prospective in the MV5-Q specification. All 150 analysis groups completed:
15 leave-one-study-out folds, two cell representations, five frozen distances,
and five seed-specific training matrices per analysis group.

This is a methodological artifact milestone, not a biological result. The
pipeline did not open tissue, approach, class, biological, or technical labels
and did not calculate agreement against any label.

## Landscape and distance definitions retained

MV5-Q did not alter the accepted persistence-landscape definition. H0 and H1
remain separate exact L2 distances over all consecutive active landscape levels,
with no fixed level cap and no fixed uniform-grid approximation. The five
clustering inputs per representation are:

1. exact all-active-level H0 landscape distance;
2. exact all-active-level H1 landscape distance;
3. descriptive raw `sqrt(H0^2 + H1^2)` without component rescaling;
4. representation-matched cell-distribution energy distance;
5. shared training-standardized pseudobulk Euclidean distance.

The explicit alias audit proves that every training/query method pair resolves
to exactly 35,350 accepted query rows across all 15 folds and all five seeds.
The SCT query representation is `sct_fold`, while its corresponding training
representation is `sct_whole`. Both analysis representations reuse the exactly
verified shared pseudobulk distance source.

## Label-closed clustering execution

For each analysis group, PAM was fit at every integer k from 2 through 10 in
each of the five frozen seed matrices: 6,750 candidate fits total. The selected
k is the smallest value within one jackknife SE of the maximum mean agreement
among the ten seed-partition pairs. No held-out sample or label participates in
fitting or selection.

At selected k, seed-specific PAM partitions are primary and average linkage is
retained only as sensitivity. Cluster IDs are canonicalized by lexicographically
sorting each cluster's sorted member-ID signature. Held-out PAM assignment uses
the nearest frozen training medoid with lexicographic medoid-ID ties. Held-out
average-linkage assignment uses minimum mean dissimilarity to the frozen
training cluster with canonical member-signature ties.

## Accepted artifact counts

| Artifact | Observed |
|---|---:|
| Analysis groups | 150 |
| Analysis-matrix contexts | 750 |
| Source training-matrix validations | 675 |
| Candidate PAM fits | 6,750 |
| Private candidate-partition rows | 567,000 |
| Stability rows | 1,350 |
| Selected PAM/average training-partition rows | 126,000 |
| Held-out assignment rows | 9,000 |
| Selected k range | 2 to 10 |

The selected-k range is reported only as an execution-axis check. MV5-Q does
not rank representations, distance components, methods, folds, or tissues.

## Independent validation

The independent validator re-read every candidate, stability, selected-
partition, and held-out artifact and reproduced:

- all 150 complete k grids and one-SE selections;
- all canonical sorted-member cluster IDs;
- exactly one frozen medoid per PAM cluster;
- all PAM and average-linkage held-out assignments;
- every public row and source/queue/completion hash boundary;
- closed-label and no-biological-outcome counters.

The canonical maximum-size fold, `large_loso_v1:SRA713577` with 89 training
samples, was repeated cleanly for all ten representation-distance paths. All 40
candidate, stability, selected-partition, and held-out artifacts were
byte-identical. A full explicit resume reused all 150 groups. Before/after
snapshots showed all 753 artifact paths, SHA-256 hashes, byte sizes, and
modification timestamps unchanged, with zero rebuild.

## Resource evidence

The 150 clustering groups used `134.125` aggregate elapsed seconds. Maximum
single-group time was `1.838` seconds and peak process RSS was `1219973120`
bytes, below the frozen 900-second and 4-GiB group limits. Source revalidation
is intentionally I/O-heavy because it re-hashes all 4,565 accepted MV5-P units;
that validation time is not included in clustering-group elapsed time.

## Corrective execution history

Two prospective checks stopped unsafe execution without altering accepted
upstream artifacts:

1. The first one-group pilot exposed a wrapper handoff that projected away
   identity columns after validation. No clustering artifact had been written;
   the handoff was corrected and the prefreeze was regenerated.
2. The first full attempt stopped after nine immutable groups when the SCT
   pseudobulk alias incorrectly used the integrated bundle's representation
   name. That root is not accepted. The alias was corrected, a ten-alias
   Cartesian prefreeze audit and ten-path pilot were added, and production was
   restarted in the accepted `production-v2` root.

During independent validation, embedded carriage-return separators in average-
linkage member signatures were observed to normalize to line feeds on CSV read.
Validation therefore normalizes platform line endings before comparing the same
ordered member signature; cluster, distance, tie choice, and stored hashes are
unchanged.

## Public-safety boundary

Public MV5-Q tables contain only frozen identifiers, label-free seed stability,
cluster/medoid/assignment artifacts, provenance, validation, and resources.
They contain no outcome labels or label-agreement endpoints. The private
candidate cache and stopped/pilot roots remain Git-ignored. `example_run.r`,
PDFs, reviewer material, and private source caches remain untracked.

## Code and package verification

- Parent MV5-N focused tests: 29/29 passed.
- MV5-Q focused tests: 15/15 passed.
- Complete package test suite: 831/831 passed with zero failures, warnings, or
  skips.
- Six MV5-Q operational R scripts: 6/6 parse.
- Clean tracked-source `R CMD check`: `Status: OK`.
- Source tarball SHA-256:
  `934a89a3ab4f0c08e952b42cdb8f90c6a7fb22d84d1f41c478cebcc3b24cbe02`.
- Check-log SHA-256:
  `a5048ffaa5cb849655c206e586cb918e7e49243edd2853b33e56950b31360964`.

## Gate and next action

MV5-Q passes its technical completion gate. The next sprint should prospectively
freeze MV5-R prediction-locked clustering-outcome evaluation before opening any
accepted label source. That specification must bind label identities and hashes,
define training-alignment and held-out-generalization endpoints, freeze
missingness and multiplicity handling, preserve PAM as primary and average
linkage as sensitivity, and prohibit tuning after outcomes are visible.

Gene topology, cell/gene fusion, spectral promotion, robustness expansion, new
data, Rust/optimization, package-default changes, and manuscript claims remain
outside that next sprint.
