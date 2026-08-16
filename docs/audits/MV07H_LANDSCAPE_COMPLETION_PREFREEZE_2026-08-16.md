# MV7-H landscape-completion validation prefreeze

Date: 2026-08-16

Status: accepted; three completion repeats authorized

## Completed private corpus

All 20 frozen landscape groups are present and hash-verified: five seeds for
each of cell-H0, cell-H1, gene-H0, and gene-H1. They contain 152,520 total
component rows. The 19-group stage-two run completed in 2,988.372 charged
seconds; combined with the primary stress group, landscape production used
3,652.144 of the frozen 23,232.367-second budget. Peak monitored process-tree
RSS was 204,963,840 bytes, and the private landscape root occupied 94,770,517
bytes at completion.

## Prospective validation selections

Commit `90f9798` freezes the completion-validation inputs before any new repeat
or full-corpus numerical oracle is run. The prefreeze passes 8/8 checks and
selects without biological labels or distance values.

The existing byte-identical stress repeat covers gene-H1. Three additional
component repeats are selected by maximum MV7-G sentinel interval burden, then
queue order:

- seed 20260805 gene-H0;
- seed 20260808 cell-H1; and
- seed 20260805 cell-H0.

One numerical oracle per component is selected by maximum combined finite
interval count, then pair ID. This yields combined depths of 766 for cell-H0,
1,025 for cell-H1, 998 for gene-H0, and 5,865 for gene-H1. The gene-H1 oracle
therefore independently challenges the heavy fallback seed rather than merely
reusing the earlier stress-seed oracle.

The repeat runner enforces the original 3,600-second and 12-GiB group limits,
uses one worker and zero retries, and accepts a repeat only when its entire
distance CSV is byte-identical to the frozen production input. Dimension
combination, clustering, labels, outcomes, and claims remain closed pending
complete-corpus independent validation.
