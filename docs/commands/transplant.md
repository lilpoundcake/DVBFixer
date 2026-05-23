# dvbfixer transplant — Molecule Transplanting

[← command index](index.md) · [← README](../../README.md)

Transplants molecules from a graft PDB into an acceptor PDB. Designed for the GLYCAM glycoprotein workflow: extract glycosylation site residues from your protein, submit to GLYCAM-Web, then transplant the GLYCAM output (renamed protein residues + glycan trees) back into the full structure. Also works with CHARMM-GUI output — use simple transplant mode (`--donor` + `--select`) to copy glycan chains or other molecules from CHARMM-GUI PDB into your structure.

## Usage

```bash
# Graft workflow: replace donor residues in acceptor with GLYCAM output
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb

# With Kabsch superposition (if structures are not pre-aligned)
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb --superpose

# With OpenMM relaxation (AMBER+GLYCAM, 4-stage minimization)
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb --relax

# Export GROMACS topology via ACPYPE
dvbfixer transplant acceptor.pdb --donor donor.pdb --graft glycam_output.pdb --gromacs gmx_output/

# Simple transplant: copy selected residues from donor to acceptor
dvbfixer transplant acceptor.pdb --donor donor.pdb --select A:NAG

# CHARMM-GUI: transplant glycan chains from CHARMM-GUI output
dvbfixer transplant protein.pdb --donor charmm_gui.pdb --select G,H --superpose
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--donor` | (required) | Donor PDB: original residues from acceptor (identifies replacement sites) |
| `--graft` | none | Graft PDB: modified donor + added molecules (e.g. GLYCAM output) |
| `--select` | none | Selection for simple transplant: chain IDs, residue names, or ranges |
| `--align` | none | Chain mapping for superposition: DONOR:ACCEPTOR (e.g. H:H, repeatable) |
| `--superpose` | off | Enable Kabsch superposition (auto-detect chain mapping) |
| `--relax` | off | Run OpenMM minimization with AMBER+GLYCAM after transplant |
| `--relax-stages` | `1000:5000,...` | Relaxation stages: k1:iter1,k2:iter2,... (k in kJ/mol/nm2) |
| `--gromacs` | none | Export GROMACS topology via ACPYPE to specified directory |
| `-o`, `--output` | `<acceptor>_transplant.pdb` | Output PDB |
| `-v`, `--verbose` | off | Print detailed progress |

## Graft Workflow

1. **`--donor`**: Original protein residues extracted from acceptor (e.g. ASN307, SER308). Used to identify which acceptor residues to replace and provides CA atoms for alignment.
2. **`--graft`**: GLYCAM output with renamed residues (NLN/OLS/OLT) + glycan trees (4YB/VMB/UYB etc.). Preserves GLYCAM atom and residue names.
3. Optional `--superpose`: Kabsch alignment of donor to acceptor, same transform applied to graft.
4. Donor residues removed from acceptor, graft protein residues inserted at correct positions, glycans appended.
5. Non-protein residues with resseq backward jumps are split into separate chains (prevents duplicate residue numbers when multiple glycan trees share a graft chain).
6. CONECT records remapped via atom identity (chain, resseq, atomname).

## Relaxation (`--relax`)

4-stage energy minimization with AMBER14 + GLYCAM_06j-1:
- Protein heavy atoms restrained; glycans move freely
- Stages: k=1000 -> 100 -> 10 -> 0 kJ/mol/nm2
- Preprocessing: CYS->CYX for disulfides, bond addition for GLYCAM residues, hydrogen re-addition

## GROMACS Export (`--gromacs DIR`)

Same ACPYPE pipeline as `dvbfixer top --acpype`: OpenMM parametrization -> ParmEd -> ACPYPE with `[ pairs_nb ]` for mixed 1-4 scaling. Outputs `topol.top`, `.gro`, and `posre_*.itp` to the specified directory. The `.top` includes position restraints and water/ion moleculetypes.

## See also

- [`glycam`](glycam.md) — convert PDB sugar names to GLYCAM convention before `transplant`
- [`top`](top.md) — `--acpype` provides the same GROMACS export pipeline for stand-alone structures
- [`prepare`](prepare.md) — run before `transplant` to clean up the acceptor
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — full GLYCAM-Web glycoprotein workflow recipe
