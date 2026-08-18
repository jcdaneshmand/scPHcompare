# MV8-H Cell Ranger 8.0.1 count-sentinel prefreeze

## Outcome

One complete unit is prospectively frozen without running `cellranger count`.
The deterministic, label-closed selection rule chooses the complete unit with
the greatest total compressed FASTQ bytes, breaking any tie by original unit
order. The unique result is `HCA_BM_002`: six exact files totaling
11,249,623,632 bytes. All six live cache files pass independent size
and SHA-256 checks. Selection did not inspect expression, QC, donor attributes,
labels, or outcomes.

## Exact execution contract

```text
cellranger count --id=mv08h_count_sentinel_hca_bm_002 --transcriptome=<verified_custom_reference> --fastqs=<verified_hca_bm_002_fastq_directory> --sample=MantonBM2_HiSeq_9 --chemistry=SC3Pv2 --expect-cells=7000 --include-introns=false --create-bam=false --nosecondary --localcores=4 --localmem=32 --disable-ui
```

The live 19-file custom reference independently retains its exact
20,765,871,518-byte tree and SHA-256
`5e2aff9e7154e6b02f98552a4419bd48edce66e617e579ae562e714f79199f1c`. Cell Ranger 8.0.1 remains bound by runtime
tree SHA-256 `aafd39e293e0ba9d14dba3896a6aeda077304531a2702d26bda0c62c4688fdf3` and launcher SHA-256
`4ee3a1670b4f14c826004fe8e17b4759e1edc701b15ff2e9623753bf1b34d4d6`.

The fixed contract is **4 cores** and **32 GiB** (four-core/32-GiB), a
downward-only resource amendment from the earlier 16-core/64-GiB maximum.
The existing 80-GiB process-tree RSS ceiling, 200-GiB workspace ceiling, and
1-TiB free-space floor remain. The elapsed observation ceiling is 96 hours so the
deliberately smaller core allocation may take longer. Monitoring is
non-destructive: a breach invalidates scientific admission and preserves artifacts
for owner review; it does not kill or delete.

## Scientific and privacy stop

This prefreeze does **not** execute or authorize count. It does not open a
matrix, evaluate QC, select cells, normalize data, run PCA, compute PH or
persistence landscapes, cluster, access labels/outcomes, process the remaining
seven units, or delete anything. Public closure may contain only file/runtime/
reference identities, matrix shapes and feature-axis identities, aggregate
resource measurements, and run status—not private paths, expression values,
barcodes, donor attributes, labels, or outcomes.

The corrected downstream landscape definition is unchanged: cells and genes
are separate typed observation views; H0 and H1 remain separate; essential H0
is excluded; and landscapes use every consecutive active level with exact or
error-controlled squared-L2 integration, no fixed grid, and no universal level cap.
A successful future sentinel can open only a separate structural/QC review
and remaining-unit decision.
