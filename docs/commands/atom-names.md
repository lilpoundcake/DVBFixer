# `dvbfixer atom-names`

Converts residue and atom names in a legacy PDB for an AMBER or CHARMM
GROMACS consumer. Unlike `rename`, this command preserves explicit protonation
variants and applies target-specific atom naming. Unlike `convert`, it does not
convert glycan nomenclature or infer chemical connectivity.

The source is never modified. V1 accepts `.pdb` and `.ent` input only, requires
an explicit output, and rejects multi-model structures, naming collisions, and
ambiguous residue identities before publishing an output.

```bash
dvbfixer atom-names input.pdb -o amber.pdb \
  --target-ff amber --profile gromacs
```

Use `--dry-run` to validate and summarize a conversion without writing the PDB.
Use `--report-json report.json` for a machine-readable report on either success
or a validation failure.

Variant overrides are supplied as a JSON array. Residue identity includes the
case-sensitive chain ID, residue number, and insertion code:

```json
[
  {"chainId": "H", "residueNumber": "82", "insertionCode": "A", "variant": "HIE"}
]
```

```bash
dvbfixer atom-names input.pdb -o amber.pdb --target-ff amber \
  --variant-overrides variants.json --report-json naming-report.json
```

The report uses `schemaVersion: 1` and includes the request, success state,
stable error details when applicable, output publication state, diagnostics,
summary counts, and per-coordinate changes. Exit code `2` denotes validation or
conversion failure; publication and other I/O failures use exit code `1`.

See the generated [CLI reference](../reference/atom-names.md) for all options.
