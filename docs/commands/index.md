# Command index

[← README](../../README.md)

Every `dvbfixer` subcommand has its own page. The first column is alphabetical; the workflow group hints at where each command typically sits in a pipeline.

Commands that accept directory input document it on their own page; see the
shared [Batch mode guide](../batch-mode.md) for supported tools and failure
behavior.

| Command | What it does | Workflow group |
|---------|--------------|----------------|
| [`cluster`](cluster.md) | Glycan conformational clustering from MD trajectories (GFDB-style) | Analysis |
| [`conect`](conect.md) | Infer missing CONECT records (SS, glycosidic, glycosylation) into a PDB | Utilities |
| [`convert`](convert.md) | Convert between PDB/AMBER/GLYCAM and CHARMM naming (sugars + protonation variants); bidirectional | Glycoprotein prep |
| [`diagnose`](diagnose.md) | Report structure-quality issues without modifying the input | Analysis |
| [`doctor`](doctor.md) | Report installed backends, executables, and OpenMM platforms | Utilities |
| [`homology`](homology.md) | Multi-template homology modeling with Modeller (antibody-aware) | Modeling |
| [`msa`](msa.md) | Multiple protein-sequence alignment using MAFFT, MUSCLE 5, or Clustal Omega | Modeling |
| [`minimize`](minimize.md) | Energy minimization with OpenMM + selective restraints, optional xtb/obminimize refinement | Refinement |
| [`model`](model.md) | Rebuild missing loops/gaps with Modeller's LoopModel | Structure prep |
| [`parametrize`](parametrize.md) | GAFF2 + AM1-BCC/RESP small-molecule parametrization (GROMACS-ready) | Topology |
| [`prepare`](prepare.md) | PDBFixer missing-atom/H repair; optional SMILES-guided ligands; substitution and deletion mutations | Structure prep |
| [`protonate`](protonate.md) | PROPKA3 pKa prediction + AMBER residue renaming + H repair | Refinement |
| [`pull`](pull.md) | OpenMM partial minimization to form bonds (SS, glycosidic) | Refinement |
| [`puppet`](puppet.md) | Strip a PDB to backbone-only polyglycine (template / visualization) | Utilities |
| [`rename`](rename.md) | Canonicalize residue names (AMBER/CHARMM/MSE → standard PDB) | Utilities |
| [`renumber`](renumber.md) | FASTA/SEQRES renumbering OR antibody schemes (Kabat/Chothia/IMGT/Martin/Aho/EU) | Structure prep |
| [`salign`](salign.md) | Structure-based multiple alignment and superposition using Modeller SALIGN | Modeling |
| [`split`](split.md) | Empirical chain splitting for GRO/PDB files without chain IDs (multi-MODEL aware) | Structure prep |
| [`top`](top.md) | GROMACS topology from PDB/GRO (AMBER, CHARMM, or ACPYPE pipeline) | Topology |
| [`transplant`](transplant.md) | Transplant molecules between PDB structures (GLYCAM glycoprotein workflow) | Glycoprotein prep |
| [`zbs`](zbs.md) | Full pipeline: renumber → model → prepare → minimize (PROPKA + Reduce run inside prepare) | Pipeline |

## See also

- [Batch mode](../../README.md#batch-mode-folder-input) — folder processing for single-structure commands

- [Force fields](../force-fields.md) — short-name aliases (`--ff amber`, `--ff charmm`, `--ff amber+glycam`, …), auto-detection rules, the two `--ff` namespaces
- [Pipelines](../pipelines.md) — end-to-end recipes
- [Known issues](../known-issues.md)
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — recipe collection + gotchas
- [ARCHITECTURE.md](../../ARCHITECTURE.md) — module structure + key abstractions
