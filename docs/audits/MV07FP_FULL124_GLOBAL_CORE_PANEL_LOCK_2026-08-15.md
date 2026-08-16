# MV7-FP full-124 global-core panel lock audit

Date: 2026-08-15

Status: complete; MV7-G prefreeze authorized

## Scope and boundary

MV7-FP applied the unchanged MV6-C availability and within-cache SCT-variance
selection rule to the complete label-closed 124-sample corpus: 620 cache
records (124 samples by five seeds), each containing the frozen 384-cell
selection. It did not fit PCA, construct typed views, calculate PH or
landscapes, open labels, or calculate biological outcomes.

The authoritative input freeze is
`docs/audits/mv07fp-prefreeze-evidence-v4/`, committed at `805bd57`. It binds
the exact 620 cache hashes and sizes, the 124-by-five axis, the accepted
90-derived comparison panel, the selector/runtime sources, one streamed panel
inventory, a 2,700-second elapsed cap, an 8-GiB process-RSS cap, and the
panel-only label firewall. Its independent prefreeze validation passed 7/7.

## Scientific result

The accepted run verified 620/620 cache identities and measured:

- 34,970 features in the union and 2,388 in the exact intersection;
- 2,262 retained-category intersection features and 126 excluded technical
  features;
- zero features with nonfinite variance in any cache;
- 58 features with variance not strictly above machine epsilon in every
  cache;
- 2,205 eligible unique canonical genes for 500 requested slots; and
- an eligibility margin of 1,705 genes.

The exact ordered 500-gene panel has canonical SHA-256
`48e881ee753893bfaecd7101fc16fbcb552dd30a2a791fedfc7204d7b322a32e`.
That digest is reproducible from a row-name-free payload containing only
integer `panel_order`, character `feature_id`, and character `gene`.

Compared with the accepted 90-sample panel, 465 of 500 features overlap
(Jaccard 0.8691588785). The order is not identical, as expected after adding
34 samples to a transductive technical-harmonization fit.

Across the five prespecified seeds, the seed-specific top-500 panels overlap
the global panel by 492 to 495 genes. Jaccard values range from 0.9685039370
to 0.9801980198, and Spearman correlations between global and seed-specific
candidate ranks range from 0.9992371607 to 0.9993309225.

## Resource result

The accepted production inventory took 986.152 elapsed seconds, peaked at
463,958,016 bytes process RSS, and read 4,163,901,516 bytes of frozen cache
artifacts. Both prospective resource gates passed. The clean exact repeat took
approximately 907 seconds and reproduced the same panel digest.

## Independent and deterministic validation

The standalone production validator independently reopened all 620 caches and
reconstructed the exact feature intersection, eligibility rule, ordered panel,
canonical panel digest, comparison, decision/resource state, and artifact
manifest. It passed 9/9 categories.

The exact same production runner was then executed from the same `805bd57`
commit into ignored temporary storage. Seven deterministic public artifacts
matched byte-for-byte and SHA-for-SHA; only the runtime-dependent resource file
was excluded prospectively.

The original 9/9 validator checked the five seed axes but did not independently
reconstruct every seed metric. A supplemental independent validator therefore
reopened all 620 caches, reconstructed the exact ordered global panel and all
five seed-specific summaries, and passed 5/5 categories. Counts and overlaps
matched exactly; Jaccard and Spearman values matched at `1e-14` tolerance.

The complete repository R test suite passed across 360 test result groups with
four established environment-dependent skips and no failure or warning.

## Retained corrections and failed attempts

Four pre-acceptance attempts are retained in history:

1. The first standalone runner did not source the inherited variance helper
   and stopped before cache inventory or public output. The source binding was
   corrected at `e23413c`.
2. The next runner reached the first sparse SCT matrix and stopped because the
   inherited helper called base `rowMeans()`/`rowSums()` on an unattached S4
   Matrix class. Commit `ff53e73` added a sparse-safe Matrix path and a
   dense-versus-sparse numerical regression test without changing the sample
   variance formula.
3. A complete provisional run selected 500 genes, but independent validation
   rejected its digest because the serialized data frame included hidden row
   names that were absent from the public CSV. That provisional evidence was
   not committed and was removed before regeneration.
4. Commit `027a128` canonicalized the digest payload and added row-name
   invariance tests. The authoritative v4 freeze and accepted calculation were
   rebuilt from this corrected contract.

These corrections changed runtime binding and content-address
reproducibility, not the feature-availability rule, variance formula, ranked
gene order, cache contents, labels, or biological outcomes. Earlier v1-v3
prefreeze directories remain public audit history; v4 is authoritative.

## Public and private separation

Public evidence contains only the ordered panel, aggregate eligibility and
stability summaries, resource totals, hashes, Boolean gates, and validation
ledgers. It contains no expression matrices, per-cache variances or ranks,
cell identifiers, private source paths, labels, outcomes, PDFs, or generated
binaries. Private caches and the clean repeat remain under ignored `tmp/`
state.

## Decision

MV7-FP passes. The next authorized stage is a separate MV7-G prefreeze for the
six-sample, five-seed, dual-view typed-view/PH sentinel. Full PH, pairwise
landscape calculation, labels, outcomes, clustering claims, and manuscript
claim promotion remain closed until their later gates.

The revised landscape definition remains unchanged: finite positive intervals,
essential H0 excluded, all consecutive active levels, H0 and H1 separate,
exact or error-controlled squared-L2 distance, no universal grid or level cap,
and streamed or chunked computation.
