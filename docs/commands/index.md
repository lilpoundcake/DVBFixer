# Command index

[← README](../../README.md)

Every `dvbfixer` subcommand has its own page. The first column is alphabetical; the workflow group hints at where each command typically sits in a pipeline.

| Command | What it does | Workflow group |
|---------|--------------|----------------|
| [`cluster`](cluster.md) | Glycan conformational clustering from MD trajectories (GFDB-style) | Analysis |
| [`glycam`](glycam.md) | Convert between PDB/CHARMM and GLYCAM nomenclature (bidirectional) | Glycoprotein prep |
| [`homology`](homology.md) | Multi-template homology modeling with Modeller (antibody-aware) | Modeling |
| [`minimize`](minimize.md) | Energy minimization with OpenMM + selective restraints, optional xtb/obminimize refinement | Refinement |
| [`model`](model.md) | Rebuild missing loops/gaps with Modeller's LoopModel | Structure prep |
| [`parametrize`](parametrize.md) | GAFF2 + AM1-BCC/RESP small-molecule parametrization (GROMACS-ready) | Topology |
| [`prepare`](prepare.md) | PDBFixer-based missing-atom/H repair with BioLuminate-style heterogen H | Structure prep |
| [`protonate`](protonate.md) | PROPKA3 pKa prediction + AMBER residue renaming + H repair | Refinement |
| [`pull`](pull.md) | OpenMM partial minimization to form bonds (SS, glycosidic) | Refinement |
| [`puppet`](puppet.md) | Strip a PDB to backbone-only polyglycine (template / visualization) | Utilities |
| [`rename`](rename.md) | Canonicalize residue names (AMBER/CHARMM/MSE → standard PDB) | Utilities |
| [`renumber`](renumber.md) | SEQRES-based residue renumbering, removes insertion codes | Structure prep |
| [`split`](split.md) | Empirical chain splitting for GRO/PDB files without chain IDs | Structure prep |
| [`top`](top.md) | GROMACS topology from PDB/GRO (AMBER, CHARMM, or ACPYPE pipeline) | Topology |
| [`transplant`](transplant.md) | Transplant molecules between PDB structures (GLYCAM glycoprotein workflow) | Glycoprotein prep |
| [`zbs`](zbs.md) | Full pipeline: renumber → model → prepare → minimize → protonate → minimize | Pipeline |

## See also

- [Pipelines](../pipelines.md) — end-to-end recipes
- [Known issues](../known-issues.md)
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — recipe collection + gotchas
- [ARCHITECTURE.md](../../ARCHITECTURE.md) — module structure + key abstractions
