# MV8-H Cell Ranger 8.0.1 custom-reference closure

## Outcome

The prospectively frozen four-core/32-GiB `cellranger mkref` run completed
successfully in 13,603 seconds. Independent validation passes
15/15 checks. The final reference contains 19 regular files,
20,765,871,518 regular-file bytes, zero partials, and complete required STAR
index components. Its deterministic tree SHA-256 is `5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c`.

## Exact scientific inputs

- The embedded primary-assembly FASTA is byte-identical to the frozen
  Ensembl-93 input: 3,151,425,857 bytes and SHA-256
  `78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c`.
- The embedded gzip-compressed GTF decompresses to the exact frozen
  target-complete annotation: 1,099,737,654 bytes and SHA-256
  `e28e4c4faf0dd76884d5e94c481fce2db43ad303968067c1276092a234727182`.
- The annotation contains 33,563 unique genes and
  2,565,751 feature records. All 500 exact-panel stable IDs and all
  475 common-panel stable IDs are present.
- `reference.json` binds Cell Ranger 8.0.1, four threads, 32 GiB, the exact
  genome name, input filenames, and Cell Ranger's matching SHA-1 identities.

## Resource closure

The private monitor recorded 450 samples. Peak subprocess-tree RSS was
30.494 GiB, below the selected 32-GiB allocation. Peak run-tree
size was 27,819,306,242 bytes, below 50 GiB. Minimum free space was
1,987,831,930,880 bytes, above 1 TiB. No resource breach, automatic kill,
or launcher deletion occurred. Cell Ranger's own normal temporary-file cleanup
is reflected in the changing run-tree samples.

The prospectively committed launcher sampled RSS correctly but did not include
RSS in its composite `resource_breach_detected` Boolean. This independent
closure therefore applies the intended 32-GiB memory gate directly to all 450
samples; the observed 30.494-GiB peak passes. The launcher is corrected
forward-only so any future RSS overage also marks a formal breach. This does not
change the executed command, reference, or scientific inputs.

## Scope and next gate

This closure does **not** execute or authorize `cellranger count`. It opens only
a prospective one-complete-unit count-sentinel prefreeze. The remaining seven
units, QC, matrices, PCA, PH, persistence landscapes, clustering, labels,
outcomes, and deletion remain closed.

The corrected topology contract is unchanged: cell and gene observations are
separate typed views; H0 and H1 remain separate; essential H0 is excluded; and
landscapes use every consecutive active level with exact or error-controlled
squared-L2 integration, no fixed grid, and no universal level cap.

The executed launcher was prospectively committed as `42fdc91` before the run.
Public evidence contains only identities and aggregates; private paths, licensed
runtime contents, expression, barcodes, donor attributes, labels, and outcomes
are not published.
