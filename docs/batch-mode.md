# Batch mode

Batch mode runs one supported single-structure command independently for every
`.pdb` or `.ent` file in a directory (`split` also accepts `.gro`). Outputs are
written under a separate directory; `--recursive` preserves relative paths and
processing continues after individual failures unless `--fail-fast` is used.

```bash
dvbfixer zbs --input-dir structures --output-dir fixed --recursive
```

The shared keys are `--input-dir DIR`, `--output-dir DIR`, `--recursive`, and
`--fail-fast`.

## Support by tool

| Tool | Batch input | Notes |
|---|---:|---|
| `split` | Yes | Accepts PDB, ENT, and GRO files. |
| `renumber` | Yes | Renumbers each structure independently. |
| `model` | Yes | Models each structure independently. |
| `pull` | Yes | Pulls each structure independently. |
| `prepare` | Yes | Prepares each structure independently. |
| `minimize` | Yes | Minimizes each structure independently. |
| `protonate` | Yes | Protonates each structure independently. |
| `rename` | Yes | Renames each structure independently. |
| `convert` | Yes | Converts each structure independently. |
| `conect` | Yes | Rebuilds connectivity for each structure independently. |
| `puppet` | Yes | Applies the requested operation to each structure independently. |
| `diagnose` | Yes | Writes one diagnostic report per structure. |
| `zbs` | Yes | Runs the selected ZBS stages independently per structure. |
| `top` | No | Requires one explicit input per invocation. |
| `transplant` | No | Its source/target relationship is not a directory batch operation. |
| `cluster` | No | Operates on an explicitly selected trajectory or structure set. |
| `parametrize` | No | Requires explicit ligand/parameterization inputs. |
| `homology` | No | Requires explicit query, template, and alignment relationships. |
| `msa` | No | Operates on an explicitly selected sequence set. |
| `salign` | No | Operates on explicitly selected structures/sequences. |
| `doctor` | No | Checks the installation once; it has no structure input. |

Every tool repeats its status in `dvbfixer TOOL --help` and in its command
guide. Passing `--input-dir` to an unsupported tool fails with a message that
lists the supported commands.
