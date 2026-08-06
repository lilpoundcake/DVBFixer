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

## Testing

The PDB and FASTA inputs consumed by pytest are tracked under
[`tests/fixtures/`](tests/fixtures/README.md). Verify their integrity with
`sha256sum -c tests/fixtures/MANIFEST.sha256`, then run the fast suite with
`pytest -m 'not slow'`. Tests marked `slow` exercise full preparation and ZBS
pipelines and may require Modeller, AmberTools, OpenBabel, and other optional
chemistry backends.

## Commands

| Command | What it does |
|---------|--------------|
| [`split`](docs/commands/split.md) | Empirical chain splitting or REMARK 350/BIOMT biological-assembly extraction |
| [`renumber`](docs/commands/renumber.md) | FASTA/SEQRES renumbering OR antibody schemes (Kabat/Chothia/IMGT/Martin/Aho/EU) |
| [`model`](docs/commands/model.md) | Rebuild missing loops/gaps with Modeller's LoopModel |
| [`prepare`](docs/commands/prepare.md) | PDBFixer-based missing-atom/H repair, optional SMILES-guided ligand chemistry, plus substitution and deletion mutations |
| [`pull`](docs/commands/pull.md) | OpenMM partial minimization to form SS / glycosidic bonds |
| [`minimize`](docs/commands/minimize.md) | Energy minimization with selective restraints, optional xtb/obminimize refinement |
| [`protonate`](docs/commands/protonate.md) | PROPKA3 pKa prediction + AMBER residue renaming + H repair |
| [`rename`](docs/commands/rename.md) | Canonicalize residue names (AMBER/CHARMM/MSE → standard PDB) |
| [`top`](docs/commands/top.md) | GROMACS topology from PDB/GRO (AMBER, CHARMM, or ACPYPE pipeline) |
| [`transplant`](docs/commands/transplant.md) | Transplant molecules between PDB structures (GLYCAM glycoprotein workflow) |
| [`convert`](docs/commands/convert.md) | Convert between PDB/AMBER/GLYCAM and CHARMM naming (sugars + protonation variants) |
| [`conect`](docs/commands/conect.md) | Infer missing CONECT records (SS, glycosidic, glycosylation) — runs automatically inside prepare/top/minimize/transplant/convert |
| [`cluster`](docs/commands/cluster.md) | Glycan conformational clustering from MD trajectories (GFDB-style) |
| [`homology`](docs/commands/homology.md) | Multi-template homology modeling with Modeller (antibody-aware) |
| [`parametrize`](docs/commands/parametrize.md) | GAFF2 + AM1-BCC/RESP small-molecule parametrization (GROMACS-ready) |
| [`puppet`](docs/commands/puppet.md) | Strip a PDB to backbone-only polyglycine (template / visualization) |
| [`diagnose`](docs/commands/diagnose.md) | Report structure-quality issues (missing atoms, coincident atoms, valence, clashes, chirality) — report-only |
| [`doctor`](docs/commands/doctor.md) | Report installed backends, external executables, and OpenMM platforms |
| [`zbs`](docs/commands/zbs.md) | Full pipeline: renumber → model → prepare → minimize (PROPKA + Reduce run inside prepare) |

The [command index](docs/commands/index.md) groups the same commands by workflow stage.

## Batch mode: folder input

Batch mode runs the selected single-structure command independently on every
supported structure in a directory. One broken input does not prevent the
remaining structures from being processed. Outputs go into a separate
directory, and recursive input preserves the relative directory layout.

```bash
dvbfixer zbs --input-dir structures --output-dir fixed \
  --recursive --no-solvent
```

Supported commands are `split`, `renumber`, `model`, `pull`, `prepare`,
`minimize`, `protonate`, `rename`, `convert`, `conect`, `puppet`, `diagnose`,
and `zbs`. Commands with multiple semantic inputs or multi-file topology
outputs (`transplant`, `cluster`, `parametrize`, and `top`) remain explicit
single-run workflows.

Batch processing continues after individual failures by default and prints a
success/failure summary. Add `--fail-fast` to stop at the first failure.
For `diagnose`, exit status 1 means the analysis completed but found at least
one ERROR-severity structural issue; batch output labels these as `FINDINGS`
rather than execution failures and points to the per-structure report.

## Pipelines

End-to-end recipes — quick `zbs` one-liner, manual step-by-step, GLYCAM glycoprotein, CHARMM-GUI alternative, GROMACS topology export, glycan clustering, antibody homology, small-molecule parametrization — live in [`docs/pipelines.md`](docs/pipelines.md). For the all-in-one path, see [`zbs`](docs/commands/zbs.md).

## Known issues

See [`docs/known-issues.md`](docs/known-issues.md) for FF-template gaps (terminal ASH/GLH in AMBER14), mixed 1-4 scaling (ACPYPE workaround), Modeller alignment quirks, and other gotchas.

## More

- [`BEST_PRACTICES.md`](BEST_PRACTICES.md) — opinionated recipes and gotchas
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — module structure, key abstractions, design decisions
- [`CLAUDE.md`](CLAUDE.md) — internal reference covering every subcommand in detail (also used by Claude Code)
