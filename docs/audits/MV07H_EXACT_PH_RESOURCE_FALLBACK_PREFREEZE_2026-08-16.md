# MV7-H exact PH resource-fallback prefreeze

Date: 2026-08-16

Status: accepted; PH-only resume authorized

## Authorization and purpose

The project owner approved Ripserr-primary/TDA-GUDHI exact resource fallback
with a 12 GiB fallback cap. This amendment permits completion of the frozen
124-sample PH corpus without reducing or omitting gene H1 and without changing
the subsequent persistence-landscape definition.

The accepted implementation is commit
`936acf4bd2e3aa3a2268c2dbeb8e2559e85302e0`, with implementation root
`8f8fbefe782f76bbf42095651576b27167944d1727936949092cb3bcd1136ec8`.
The exact base-prefreeze root is
`8a1e1a81925ebe4081417f00c4358d9329b664099aa483fff0a8b3667c8d7f6d`,
and the resource-gate root is
`c2e7e7971cf5b65125ace533f8aa20605717c31ff8fc5f0bea470d1e2f987843`.

## Exact policy

Ripserr remains primary for every queued cell and gene view. One GUDHI attempt
is eligible only for `gene_topology_v1` after the corresponding primary attempt
has disposition `rss_cap_exceeded` and has published no output. Other failure
types remain fatal. The fallback uses one worker, an 1,800-second elapsed cap,
and a 12,884,901,888-byte process-tree RSS cap.

GUDHI receives the exact frozen 500-gene distance matrix and complete
Vietoris--Rips filtration through H1 over field 2. Only its engine-specific
capped essential H0 representation is normalized to the contract's single
infinite essential H0 interval. Every selected fallback must pass the existing
typed-view, MST, finite-positive-H1, hash, and atomic-output validators and must
be repeated byte-for-byte.

The primary and fallback attempt ledgers are separate. Public PH metrics record
the selected engine and version; public attempt evidence retains both the
failed primary and successful fallback aggregate resource receipts.

## Resume identity

The original private ledger SHA-256 is
`f651410707790e8a72d3a24c53814dd66a7bd39721af527e1f912c14d2739899`.
It contains 21 attempts: five completed source bundles, eight completed cell PH
records, seven completed gene PH records, and one gene-side Ripserr RSS-cap
failure. The 20 completed source/PH artifacts independently match their receipt
hashes and byte sizes. Their aggregate manifest root is
`19596999f168130c6be1364eecc11bc2d3fce0d15d2b67a5f8e85f09f85215c6`.
The failed output remains absent.

The amendment therefore reuses five source and fifteen PH artifacts and leaves
1,225 PH jobs to complete. It authorizes no landscape, clustering, label,
outcome, combination, or claim job.

## Validation

- Focused fallback suite: 40 expectations, zero failures or warnings.
- Complete repository suite: 1,704 expectations in 370 result groups, zero
  failures or warnings, four established skips.
- Real hard-view fallback: valid 93,031-byte MV7-H PH record with SHA-256
  `ac78ae74f263cd7c66f67e2d42e4a012071904060407669878063870ec95ec36`.
- Real hard-view repeat: byte- and SHA-identical.
- Independent prefreeze validation: 12/12 categories pass for both builds.
- Prefreeze repeat: 33/33 artifacts byte-identical.

The next authorized action is exact-head resumption of full PH only. The one
prospectively selected landscape stress group remains closed until the complete
1,240-record corpus, all fallback repeats, and the independent full-PH
validator pass.
