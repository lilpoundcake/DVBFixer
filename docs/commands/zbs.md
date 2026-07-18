# dvbfixer zbs — Full Pipeline

[← command index](index.md) · [← README](../../README.md)

Runs the complete preparation workflow in one command: **renumber → model → prepare → minimize → protonate → minimize**.

Six steps. The pipeline was simplified from seven — the old `protonate --no-hydrogens` rename-only step (and its second-protonate cleanup) has been dropped because `--no-hydrogens` leaves existing H atoms in wrong positions relative to the new HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN residue names. The full `protonate` (PROPKA + MolProbity Reduce for HIS tautomers + ASN/GLN flip detection + variant-aware OpenMM addHydrogens) now runs once between the two minimize passes; the second minimize just refines the H positions.

Intermediate files are cleaned up by default — use `--keep-interim` to preserve them. The `.dat` file flows from model (gap atoms) through prepare (merged with PDBFixer additions) to minimize (selective restraints). Water is removed by default. Each step can be skipped individually.

## Pipeline

| # | Step | What it does |
|---|---|---|
| 1 | `renumber` | SEQRES-based renumbering; removes insertion codes |
| 2 | `model` | Rebuild missing loops via Modeller (`--num-output` saves top-N candidates) |
| 3 | `prepare` | Fix missing atoms + preliminary hydrogens via PDBFixer |
| 4 | `minimize` (pass 1) | Relax heavy atoms using prepare's approximate hydrogens |
| 5 | `protonate` (**full**) | PROPKA pKa + MolProbity Reduce (`--protassign` default ON) + variant-aware addHydrogens. Places final H correctly with respect to the ASN/GLN flips and HIS tautomer choices. |
| 6 | `minimize` (pass 2) | Refine the freshly-placed hydrogens; heavy atoms stay restrained. Preserves AMBER variant names via minimize's `_input_variants` capture. Optional post-minimize refinement via `--refine {xtb, obminimize}`. |

## Usage

```bash
# Full pipeline
dvbfixer zbs input.pdb -v

# Skip terminal modeling and use vacuum minimization
dvbfixer zbs input.pdb --no-terminal --no-solvent -v

# Skip model and minimize steps
dvbfixer zbs input.pdb --skip-model --skip-minimize -v

# Custom output and pH
dvbfixer zbs input.pdb -o output.pdb --ph 6.5 -v

# With point mutations
dvbfixer zbs input.pdb --mutate A:39:ALA --mutate B:100:GLY -v

# Save top 3 candidate loop models from Modeller (best is used downstream)
dvbfixer zbs input.pdb --num-loops 4 --num-output 3 -v

# xtb refinement pass after minimize pass 2 (organic ligands)
dvbfixer zbs input.pdb --refine xtb --refine-heterogens-only -v

# Disable MolProbity Reduce in protonate (fall back to pure PROPKA + --his-default)
dvbfixer zbs input.pdb --no-protassign -v
```

## Options

Organised by which pipeline step each flag flows into.

### General
| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_zbs.pdb` | Final output PDB file |
| `--ph` | 7.0 | pH for protonation and hydrogen addition (used by prepare, minimize, protonate) |
| `--ff` | `auto` | Force field forwarded to prepare / minimize / protonate. Accepts a short name (`auto`, `amber`, `amber+glycam`, `charmm`, …) or explicit OpenMM XML paths. See [force-fields.md](../force-fields.md). |
| `--keep-water` | off | Keep water molecules (removed by default) |
| `--no-infer-conect` | off | Skip auto CONECT inference in prepare/minimize/protonate |
| `--keep-interim` | off | Keep all intermediate files (default: only final output) |
| `--align-to-input` / `--no-align-to-input` | **on** | After every pipeline step, Kabsch-align the output back to the ORIGINAL input on protein backbone atoms. Prevents accumulated rigid-body drift so residue-by-residue comparisons in a viewer line up. Pass `--no-align-to-input` for the legacy behaviour (each step's output in its own frame). |
| `-v`, `--verbose` | off | Print detailed progress for all steps |

### `renumber` (step 1)
| Flag | Default | Description |
|------|---------|-------------|
| `--skip-renumber` | off | Skip this step |

### `model` (step 2)
| Flag | Default | Description |
|------|---------|-------------|
| `--skip-model` | off | Skip this step |
| `--no-terminal` | off | Do not model N/C terminal residues |
| `--num-loops` | 2 | Number of loop refinement models per initial model |
| `--md-level` | fast | Modeller MD refinement level (`none`, `fast`, `slow`, `very_slow`, `slow_large`) |
| `--fasta` | none | FASTA file with complete sequence(s) |
| `--num-output` | 1 | Save top-N candidate models; zbs uses the best (`_1`) downstream, other candidates remain on disk |
| `--pin-input` / `--no-pin-input` | **on** | Freeze input residues during Modeller's MD refinement (only gap residues move). Prevents flanking drift. Pass `--no-pin-input` for legacy LoopModel behaviour (gap ±~3 residue flank mobile). |

### `prepare` (step 3)
| Flag | Default | Description |
|------|---------|-------------|
| `--skip-prepare` | off | Skip this step |
| `--strip-heterogens` | off | Strip heterogens (protein-only pipeline); default keeps them |
| `--no-heterogen-h` | off | Skip H addition on heterogens in prepare (default: add H BioLuminate-style) |
| `--mutate CHAIN:RESNUM:NEW_AA` | none | Mutate a residue; repeatable; use `del` for deletion |
| `--rename` | off | Canonicalise non-standard residue names before prepare/minimize |

### `minimize` (steps 4 & 6)
| Flag | Default | Description |
|------|---------|-------------|
| `--skip-minimize` | off | Skip both minimize steps |
| `--no-solvent` | off | Minimize in vacuum |
| `--rebuild-h` | off | Force `--rebuild-h` on both passes (default off; pass 2 already has correct H from protonate) |
| `--restraint-k` | 100.0 | Restraint force constant for original atoms |
| `--max-iter` | 1000 | Max minimization iterations per phase |
| `--platform` | auto | OpenMM platform (`CPU`, `CUDA`, `OpenCL`, `Reference`) |
| `--refine` | `none` | Post-minimize refinement pass in step 6: `xtb` (GFN-FF) or `obminimize` (UFF). Off by default. |
| `--refine-heterogens-only` | off | Restrict `--refine` pass to heterogen residues (protein frozen). Only meaningful with `--refine != none`. |
| `--parametrize-ligands` | off | Forward `--parametrize-ligands` to both minimize passes: GAFF2 + AM1-BCC for each heterogen residue with no template in the resolved `--ff`. See [force-fields.md](../force-fields.md). |

### `protonate` (step 5)
| Flag | Default | Description |
|------|---------|-------------|
| `--skip-protonate` | off | Skip this step (also disables pass 2 minimize since it needs post-protonate input) |
| `--no-protassign` | off | Disable MolProbity Reduce (falls back to pure PROPKA + `--his-default`). Default: run Reduce. |

## See also

- [`model`](model.md), [`prepare`](prepare.md), [`minimize`](minimize.md), [`protonate`](protonate.md) — the individual steps
- [Pipelines](../pipelines.md) — manual step-by-step recipes and alternatives

## How it works
Full pipeline in 6 steps: renumber → model → prepare → minimize (pass 1) → protonate (full) → minimize (pass 2). Interim files are deleted by default (`--keep-interim` to preserve). The old 7-step design used `protonate --no-hydrogens` twice (rename-only + name-restore) — that pattern was dropped because `--no-hydrogens` leaves existing H in wrong positions after variant renames. The new design runs one FULL `protonate` call between the two minimize passes (PROPKA + MolProbity Reduce for HIS tautomers and ASN/GLN flip detection via `--protassign` default ON + variant-aware `addHydrogens`). Pass-2 minimize no longer force-adds `--rebuild-h` because H are already correct from the protonate step; it just refines them. `minimize.py`'s `_input_variants` capture-restore path preserves AMBER variant names on write, so the old step-7 name-restore protonate is no longer needed. Water removed by default. `.dat` flows from model → prepare (merged) → minimize.

**New pass-through flags** (added after tools evolved past zbs): `--num-output N` (model — save top-N candidates; zbs picks `_1` best for downstream), `--pin-input` / `--no-pin-input` (model — freeze input residues during Modeller MD; default ON, forwards `--no-pin-input` when the user opts out), `--no-heterogen-h` (prepare — skip H on heterogens), `--rename` (prepare + minimize — canonicalise non-standard names first), `--refine {none, xtb, obminimize}` (minimize pass 2 — post-minimize refinement pass), `--refine-heterogens-only` (restrict refinement to heterogens), `--no-protassign` (protonate — disable MolProbity Reduce), `--no-infer-conect` (prepare + minimize + protonate — skip auto-CONECT inference), `--parametrize-ligands` (minimize — real GAFF2 templates for unknown ligands). Existing flags: `--ph`, `--ff`, `--num-loops`, `--md-level`, `--fasta`, `--no-terminal`, `--strip-heterogens`, `--mutate`, `--no-solvent`, `--rebuild-h`, `--restraint-k`, `--max-iter`, `--platform`, `--keep-water`, `--skip-{renumber,model,prepare,minimize,protonate}`, `--keep-interim`, `-v`.

**`--align-to-input` / `--no-align-to-input`** (default ON): after every pipeline step, `_maybe_align(out)` runs `dvbfixer.align.kabsch_align_pdb(out, original_input, out, selection='backbone')` in-place. Kabsch rotation + translation is computed on backbone atoms (N/CA/C/O of standard AAs) matched by `(chain, resseq, icode, atomname)` and applied to EVERY atom in the file — protein + heterogens stay in a consistent relative frame with the user's original input, no accumulated rigid-body drift from successive OpenMM minimizations. Alignment is pure numpy (~30 lines in `align.py`), no extra deps. Legacy behaviour (each step's output in its own frame) via `--no-align-to-input`. No standalone `dvbfixer align` subcommand — kept as a pipeline-internal helper per user preference.
