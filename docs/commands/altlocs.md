# `dvbfixer altlocs`

[User guide](../README.md) · [Command index](index.md) · [CLI reference](../reference/altlocs.md)

Analyze alternate locations and select exactly one labeled state per residue.
Common atoms without altloc are retained. Selection does not repair chemistry,
add atoms, change occupancy, or optimize coordinates. Input is never overwritten.

## Quick start

```bash
dvbfixer altlocs input.pdb --policy occupancy -o selected.pdb --report-json selection.json
```

PDB, ENT, CIF and mmCIF input are accepted; output is PDB. CIF uses the shared
normalization boundary and reports any chain mapping. Output and report paths
must be new and different from the source. Inspect the result before modeling.

## Policies

| Policy | Meaning |
|---|---|
| `occupancy` | Highest mean occupancy of labeled heavy atoms; no atomwise mixing |
| `first` | First encountered labeled state, not necessarily A |
| `prefer` | Prefer one label, with explicit missing-label fallback |
| `manual` | Every alternate residue needs an explicit selection |

Every automatic policy accepts manual overrides. Occupancy ties choose the first
encountered state and warn. Missing occupancies remain unresolved. Different atom
inventories are reported; no claim of conformer completeness is made. Matching
labels across different residues are not assumed to describe correlated states.

Prefer B, using occupancy only where B is absent:

```bash
dvbfixer altlocs input.pdb --policy prefer --prefer B --fallback occupancy -o selected.pdb
```

Without `--fallback`, absent preferred states require manual decisions.

## Manual overrides

```bash
dvbfixer altlocs input.pdb --policy first --select C:82A:B --select c:15:A -o selected.pdb
```

Selectors are case-sensitive and preserve insertion codes. `_` selects a blank
chain. Use an exact JSON plan when multiple MODELs or repeated occurrences make
the shorthand ambiguous. Duplicate selections are rejected.

## Analysis and JSON plans

Omit output to analyze only:

```bash
dvbfixer altlocs input.pdb --policy manual --report-json analysis.json
```

The report lists states, labeled/shared serials, occupancies, warnings and
unresolved positions. Copy each exact residue `key` and choose one label in a
plan bound to the report's `sourceSha256`. Schematic plan (replace the digest):

```text
{
  "schemaVersion": 1,
  "sourceSha256": "<SHA-256 from analysis>",
  "selections": {"[1,\"C\",\"82\",\"A\",1]": "B"}
}
```

The key contains model ordinal, chain, residue number, insertion code and
occurrence. The GUI creates this plan automatically.

```bash
dvbfixer altlocs input.pdb --policy manual --selection-plan plan.json -o selected.pdb
```

`--dry-run` never publishes coordinates. Apply is refused until all alternate
residues have exactly one state. A stale source digest requires fresh analysis.

## Records and limitations

Chosen coordinates, occupancy and atom serials survive unchanged. Altloc columns
are cleared. ANISOU/SIGATM/SIGUIJ follow their atoms. CONECT retains only surviving
endpoints; stale records are dropped rather than repointed. LINK/SSBOND referring
to removed state chemistry are dropped with a warning. Incompatible global
connectivity across MODELs and ambiguous repeated link identities are refused.
Headers are retained except stale MASTER counts. Conflicting shared/labeled atom
identities are errors, not silently merged. Source-linked repair sidecars are not
rewritten; prepare the selected structure to regenerate matching metadata.

## GUI

Open **Altlocs** in the DVBfixer panel, choose a Workspace structure and a policy,
then **Analyze**. Select exactly one mutually exclusive checkbox per residue.
Mol* shows the focused residue's alternatives and common atoms in All, Selected,
or Compare mode. Preview does not modify the source or the two main viewers.
**Apply selection** runs this CLI and registers a new structure with a retained
selection report. Unresolved states disable Apply.

## Batch mode

Directory batch input is not supported in this version. Run once per explicit
structure; plans are source-specific. See [batch matrix](../batch-mode.md#support-by-tool).
