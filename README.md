# dvbfixer

A suite of Python CLI tools for preparing PDB (Protein Data Bank) structural biology files. Handles common issues with PDB files from MD simulations and the PDB database: missing chain IDs, antibody insertion codes, missing loops/residues, loop rebuilding with Modeller, multi-template homology modeling, energy minimization with selective restraints, protonation state assignment, GROMACS topology generation, GLYCAM glycoprotein transplanting, small molecule parametrization (GAFF2), and glycan conformational clustering from MD trajectories.

This README is the root of a manual-style documentation tree. Each subcommand has its own page under [`docs/commands/`](docs/commands/index.md); the [pipelines](docs/pipelines.md) page collects end-to-end recipes. For design notes see [`ARCHITECTURE.md`](ARCHITECTURE.md); for opinionated recipes and gotchas see [`BEST_PRACTICES.md`](BEST_PRACTICES.md).

## Quick start

```bash
micromamba create -f environment.yml
micromamba activate dvbfixer
pip install -e .
dvbfixer --help
```

Full install instructions (including the Modeller license step) are in [`docs/installation.md`](docs/installation.md).

## Commands

| Command | What it does |
|---------|--------------|
| [`split`](docs/commands/split.md) | Empirical chain splitting for GRO/PDB files without chain IDs (multi-MODEL aware) |
| [`renumber`](docs/commands/renumber.md) | SEQRES-based renumbering OR antibody schemes (Kabat/Chothia/IMGT/Martin/Aho/EU) |
| [`model`](docs/commands/model.md) | Rebuild missing loops/gaps with Modeller's LoopModel |
| [`prepare`](docs/commands/prepare.md) | PDBFixer-based missing-atom/H repair, plus substitution and deletion mutations |
| [`pull`](docs/commands/pull.md) | OpenMM partial minimization to form SS / glycosidic bonds |
| [`minimize`](docs/commands/minimize.md) | Energy minimization with selective restraints, optional xtb/obminimize refinement |
| [`protonate`](docs/commands/protonate.md) | PROPKA3 pKa prediction + AMBER residue renaming + H repair |
| [`rename`](docs/commands/rename.md) | Canonicalize residue names (AMBER/CHARMM/MSE → standard PDB) |
| [`top`](docs/commands/top.md) | GROMACS topology from PDB/GRO (AMBER, CHARMM, or ACPYPE pipeline) |
| [`transplant`](docs/commands/transplant.md) | Transplant molecules between PDB structures (GLYCAM glycoprotein workflow) |
| [`convert`](docs/commands/convert.md) | Convert between PDB/AMBER/GLYCAM and CHARMM naming (sugars + protonation variants) |
| [`cluster`](docs/commands/cluster.md) | Glycan conformational clustering from MD trajectories (GFDB-style) |
| [`homology`](docs/commands/homology.md) | Multi-template homology modeling with Modeller (antibody-aware) |
| [`parametrize`](docs/commands/parametrize.md) | GAFF2 + AM1-BCC/RESP small-molecule parametrization (GROMACS-ready) |
| [`puppet`](docs/commands/puppet.md) | Strip a PDB to backbone-only polyglycine (template / visualization) |
| [`zbs`](docs/commands/zbs.md) | Full pipeline: renumber → model → prepare → minimize → protonate → minimize |

The [command index](docs/commands/index.md) groups the same commands by workflow stage.

## Pipelines

End-to-end recipes — quick `zbs` one-liner, manual step-by-step, GLYCAM glycoprotein, CHARMM-GUI alternative, GROMACS topology export, glycan clustering, antibody homology, small-molecule parametrization — live in [`docs/pipelines.md`](docs/pipelines.md). For the all-in-one path, see [`zbs`](docs/commands/zbs.md).

## Known issues

See [`docs/known-issues.md`](docs/known-issues.md) for FF-template gaps (terminal ASH/GLH in AMBER14), mixed 1-4 scaling (ACPYPE workaround), Modeller alignment quirks, and other gotchas.

## More

- [`BEST_PRACTICES.md`](BEST_PRACTICES.md) — opinionated recipes and gotchas
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — module structure, key abstractions, design decisions
- [`CLAUDE.md`](CLAUDE.md) — internal reference covering every subcommand in detail (also used by Claude Code)
