# Publication and reproducibility assets

scPHcompare uses two linked, versioned research objects:

1. **Software:** the lean GitHub repository, archived from a tagged GitHub
   release to Zenodo and assigned a version-specific software DOI.
2. **Derived evidence:** a separate Zenodo dataset record containing generated
   audit tables and derived R result objects, with its own version-specific DOI.

The Git repository retains all source code, tests, dependency locks,
specifications, human-readable audit reports, compact manifests, checksums,
and retrieval tooling. Bulky generated evidence is reconstructed into one
deterministic ZIP because Zenodo recommends ZIP packaging for deposits with
more than 20 files.

No DOI shown in this directory is valid until the owner creates and reviews a
Zenodo draft. Creator order, ORCIDs, affiliations, CRediT roles, licenses, and
the relationship to a revised manuscript require author-team approval. The
project must preserve credit for Dr. Eric Rouchka and Dr. Akshitkumar Mistry.

## Local build

The evidence builder reads immutable Git blobs, not mutable working-tree
files. A typical predeposit build is:

```text
python scripts/build_publication_evidence_bundle.py build \
  --base <published-base> \
  --revision <accepted-source-commit> \
  --output-dir tmp/publication/evidence \
  --git-manifest docs/publication/scphcompare-derived-evidence-manifest-v1.csv
```

Verify the resulting archive independently before upload:

```text
python scripts/build_publication_evidence_bundle.py verify \
  --archive tmp/publication/evidence/scphcompare-derived-evidence-<sha>.zip \
  --manifest docs/publication/scphcompare-derived-evidence-manifest-v1.csv
```

The archive, DOI, and release remain unpublished until the author team reviews
metadata and the GitHub publication stack has passed CI.

The current local predeposit identities are recorded in:

- `scphcompare-derived-evidence-manifest-v1.csv` — all 1,136 archived members;
- `scphcompare-derived-evidence-bundle-v1.sha256` — deterministic ZIP identity;
- `ZENODO_EVIDENCE_METADATA_DRAFT.json` — explicitly provisional metadata;
- `scphcompare-publication-stack-v1.csv` — eight bounded Git publication slices;
- `DATA_AND_CODE_AVAILABILITY_TEMPLATE.md` — manuscript-ready placeholder text.
