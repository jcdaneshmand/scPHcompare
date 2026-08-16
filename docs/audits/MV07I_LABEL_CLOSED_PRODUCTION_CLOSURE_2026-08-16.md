# MV7-I label-closed production closure

Date: 2026-08-16

## Decision

MV7-I label-closed matrix and clustering production is complete. Independent
validation passes 13/13 and authorizes only a separate MV7-I outcome
prefreeze. Metadata outcomes, manuscript claims, new distances, cross-view
fusion, external data, and public release of the private result artifacts
remain unauthorized.

## Executed estimand

The run preserved the frozen landscape estimand: finite positive intervals,
essential H0 excluded, all active consecutive landscape levels, H0 and H1
reported separately, exact squared-L2 integration, no uniform grid, no level
cap, and streamed groups. The label-closed stage reconstructed the 124-sample
matrices from all 20 validated MV7-H distance files after checking every source
SHA-256 against the completed MV7-H inventory.

Six representations were produced: cell H0, cell H1, gene H0, gene H1, and the
two prespecified within-view unweighted H0/H1 composites. Across five seeds,
this yielded 30 matrices, 45,756 pair summaries, 15,252 H1-contribution
summaries, 270 candidate PAM fits, 54 stability rows, and 7,440 selected PAM
and average-linkage assignment rows. No tissue, study, approach, clinical, or
other biological metadata were loaded.

## Resource result

The production run completed in 46.622 seconds with peak process-tree RSS of
147,791,872 bytes. The clean repeat completed in 49.395 seconds with peak RSS
of 147,828,736 bytes. Both were far below the frozen 3,600-second and
4,294,967,296-byte caps. Each private artifact bundle occupies 25,952,976
bytes.

## Independent validation

The validator independently:

- checked all 20 source hashes and all 7,626 ordered pairs per group;
- reconstructed all 30 seed-specific matrices and six median matrices;
- reproduced every pair and H1-contribution summary;
- refit and reproduced all 270 PAM candidates;
- recomputed all ten pairwise seed ARIs per candidate k, five-seed
  delete-one-seed jackknife standard errors, and all six one-SE selections;
- reproduced the selected PAM and average-linkage partitions;
- verified all completion hashes;
- confirmed seven byte-identical scientific artifacts across independent
  private roots; and
- proved immutable resume by preserving the hash, size, and modification time
  of every production file.

Public evidence contains only checks, hashes, counts, resource measurements,
and the authorization decision. The scientific result files remain under the
ignored private `tmp/` roots.

## Failure and correction record

The first validation invocation stopped before computation because a
three-component `file.path()` call dropped vector names required by a named
artifact lookup. Commit `98e10ad` restored explicit names. The next complete
recomputation reached only the final firewall check; its broad substring
matcher incorrectly treated `linkage` as an age field. Commit `2497833`
changed the screen to delimiter-aware field matching and added boundary tests.
Neither attempt published a validation directory or modified the immutable
production and repeat artifacts.

## Durable evidence

- Prefreeze specification:
  `docs/specifications/MV07I_DESCRIPTIVE_TOPOLOGY_CLUSTERING_PREFREEZE_SPECIFICATION_V1.md`
- Prefreeze audit: `docs/audits/MV07I_DESCRIPTIVE_PREFREEZE_2026-08-16.md`
- Public validation evidence: `docs/audits/mv07i-label-closed-validation/`
- Production implementation: commit `7fea36f`
- Independent validator and resume hardening: commit `f4d0fa0`
- Validator corrections: commits `98e10ad` and `2497833`

## Next authorized action

Build and independently validate a prospective outcome prefreeze that fixes
the metadata endpoints, population strata, reporting order, association
statistics, multiplicity policy, uncertainty summaries, confounding warnings,
and interpretation limits before joining metadata to the immutable partitions.
The outcome prefreeze itself may inspect metadata structure and counts but may
not inspect association results or select favorable representations, seeds,
algorithms, tissues, studies, or approaches.
