# MV5-AZ label-closed stability/resampling prefreeze

Date: 2026-08-12
Status: accepted; acceleration benchmark authorized; partitions and added seeds closed

## Decision

The existing resampling artifacts do not identify a primary partition for the
accepted MV5-AY matrices. MV03 provides both views and five seeds but only two
samples per stratum. MV05C provides six samples and five seeds, but its fold-
specific axes and incomplete gene/representation coverage do not match the
current target. MV5-AY provides the intended complete axes and corrected
landscape contract but only one seed. These assets remain useful for technical
reproducibility, method reference, acceleration equivalence, and baseline-axis
freezing, respectively.

The first primary stability target is prospectively frozen as four independent
large ten-sample panels: SCT-cell, SCT-gene, integration-cell, and integration-
gene. Five identical sample axes are required. A future launch would retain the
accepted seed 20260805 and add four seeds: 160 diagrams and 720 unordered pairs,
including 360 high-depth gene-H1 pairs under the R baseline.

For each representation/view/H0-or-H1 stratum, PAM partitions are compared
across five seeds for `k=2:9`. Stability is the mean of the ten pairwise adjusted
Rand indices. Uncertainty is the delete-one-seed jackknife standard error of
that statistic. The smallest `k` within one jackknife SE of the maximum mean is
selected. A deterministic medoid/tie policy must still be frozen before launch.
H0 and H1 remain separate; combined distance, fusion, labels, outcomes, and
oracle `k` cannot select the partition.

## Speed decision

MV5-BA is authorized to benchmark the corrected Persim critical-pair engine
against accepted MV5-AY results. Persim's built-in norm remains prohibited.
Adoption requires analytical, exact-corpus, adaptive-certificate, determinism,
memory, dependency-lock, and at least threefold median high-depth speed gates.
If equivalence passes without the speed threshold, the accepted R path remains.
Rust starts only if the mature path misses throughput/packaging or the projected
full run breaches its 40-worker-hour, 12-wall-hour, two-process, 2-GiB/process,
or 20-GiB-retained caps—and it must pass the identical equivalence corpus.

## Validation

The evidence builder verifies the 80-row MV03, 150-row MV05C, and accepted
56-diagram/204-pair MV5-AY inventories before writing any decision. Twelve
independent prefreeze categories pass twice. All six public evidence files are
byte-identical across clean builds.

No partition, extra seed, biological outcome, source object, workflow default,
or legacy artifact changed in this sprint.
