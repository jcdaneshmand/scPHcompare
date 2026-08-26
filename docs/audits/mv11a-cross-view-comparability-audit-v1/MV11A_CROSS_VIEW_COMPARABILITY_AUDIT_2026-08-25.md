# MV11-A cross-view comparability audit

## Decision

A direct comparison of the current MV10 gene outcome with the historical MV7
cell outcome would confound topology view with clustering scope: MV7 used only
PAM/average at K=2, whereas MV10 uses five methods and selected K=2 for H0 and
K=3 for H1. The newer all-QC SCT and Pearson-residual production contains
internal gene topology only, so it cannot yet supply a matched newer cell view.

The historical selected-fit matrix bundle is nevertheless a valid bounded
bridge. It contains the same 124 samples, five seeds, and separate cell/gene
H0/H1 exact-landscape matrices. MV11 should therefore first expand only the
ten historical cell matrices to the unchanged MV10 five-method, K=2:10,
five-seed clustering contract. After independent closure, a separate
prefreeze may compare historical cell and gene partitions at the prespecified
common K=2 and K=3 sensitivities. View-native selected K may be reported
separately but may not be mistaken for a controlled cell-versus-gene contrast.

## Scope boundaries

- H0 and H1 remain separate; the stored composites are excluded.
- PAM remains primary; average, complete, and DIANA are sensitivities; single
  linkage remains a diagnostic.
- No label, outcome, method choice, view ranking, fusion, biological claim, or
  manuscript claim enters the cell clustering stage.
- New all-QC cell coordinates, PH, and landscapes remain deferred. They become
  eligible only under a separate source/PCA/topology contract if the matched
  historical benchmark establishes a concrete need.
- The matrix bundle remains private; only hashes, schemas, aggregate clustering
  diagnostics, and closure evidence may be public.

## Next action

Implement and prospectively freeze MV11-B/C/D for the ten immutable historical
cell H0/H1 matrices: 450 fits, 55,800 private assignments, 450 quality rows,
90 stability rows, two PAM-K rows, and 900 method-agreement rows. Stop before
labels or cross-view comparison until the cell benchmark independently closes.
