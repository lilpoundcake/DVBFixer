# dvbfixer top — GROMACS Topology Generation

[← command index](index.md) · [← README](../../README.md)

Generates GROMACS topology files directly from PDB or GRO files by parsing force field RTP/ARN/R2B/TDB files in Python — no `pdb2gmx` or GROMACS installation required. GRO files are auto-converted via MDAnalysis. Supports AMBER99SB-ILDN and CHARMM36 force fields (bundled). Output is a modular set of `.itp` files: `ffparams.itp` (all FF parameters), `{chain}.itp` (each chain moleculetype), `water.itp`, `ions.itp`, `posre_*.itp` (position restraints), `interchain_ss.itp` (inter-chain SS bonds + protein-glycan bonds), and a compact `topol.top` with only `#include` directives — no external FF directory needed. Handles proteins, carbohydrates (CHARMM glycan topology with glycosidic bond detection and protein-glycan bonds), glycolipids (ceramide + sugar tree as single moleculetype from CHARMM-GUI output), small CGenFF molecules (ACET, ACEH — auto-detected via distance splitting), lipids, nucleic acids, and all other CHARMM molecule types (~2400+ residues). Water (SOL/HOH/WAT), ions, and buffer particles (BUF) in the input are auto-detected and added to `[ molecules ]`. Automatically splits chains with overlapping residue numbers (e.g. duplicate glycan trees from `transplant`).

## Usage

```bash
# Basic usage with AMBER (default) — writes topol.top, posre_*.itp, conf.pdb
dvbfixer top input.pdb

# GRO file input (auto-converted via MDAnalysis)
dvbfixer top system.gro --ff charmm

# Buffer system (ACET/ACEH auto-detected as small molecules)
dvbfixer top buffer.gro --ff charmm

# CHARMM36 force field
dvbfixer top input.pdb --ff charmm

# Explicit disulfide bonds and HIS protonation
dvbfixer top input.pdb --ss A:22:A:96 --his A:64:HID

# Protonate all ASP/GLU/HIS (CHARMM: ASPP/GLUP/HSP)
dvbfixer top input.pdb --ff charmm --protonate all

# Protonate specific residues
dvbfixer top input.pdb --ff charmm --protonate A:66,A:46:GLUP

# Glycolipid from CHARMM-GUI (ceramide + sugar tree as one moleculetype)
dvbfixer top glycolipid_charmm.pdb --ff charmm -o glycolipid_top/

# ACPYPE mode: full AMBER14+GLYCAM parametrization via OpenMM
# Handles mixed 1-4 scaling for glycoproteins via [ pairs_nb ]
dvbfixer top input.pdb --acpype

# Custom output
dvbfixer top input.pdb -o my_topology.top --pdb my_conf.pdb
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `topol.top` | Output .top file |
| `--ff` | `amber` | Force field: `amber` or `charmm` |
| `--ff-dir` | (bundled) | Custom force field directory |
| `--water` | `tip3p` | Water model: tip3p, spc, spce, tip4p |
| `--ignh` | off | Ignore hydrogens in input PDB |
| `--ss` | auto | Disulfide bond: CHAIN1:NUM1:CHAIN2:NUM2 (repeatable) |
| `--his` | auto | HIS protonation: CHAIN:NUM:STATE (HIE/HID/HIP, repeatable) |
| `--protonate` | off | Protonate residues: `all` for every ASP/GLU/HIS, or `CHAIN:NUM[:STATE],...` for specific. H placed via OpenMM Modeller with CHARMM/AMBER FF |
| `--merge` | off | Merge all chains into single moleculetype |
| `--pdb` | `conf.pdb` | Output PDB with topology-matched atom names |
| `--acpype` | off | Use ACPYPE pipeline (AMBER14+GLYCAM -> ParmEd -> GROMACS) with per-pair 1-4 scaling |
| `-v`, `--verbose` | off | Print detailed progress |

## RTP-based mode (default)

Parses force field RTP files to build topology from bond graph: resolves inter-residue bonds (`-C`, `+N`), enumerates angles/dihedrals/pairs algorithmically, copies impropers and CMAP from templates, applies terminal patches. Intra-chain disulfide bonds (SG-SG between CYS2 residues) are added explicitly to the chain topology with all derived angles, dihedrals, and 1-4 pairs; inter-chain SS bonds go to `interchain_ss.itp`. For CHARMM36, loads all molecule-type RTP files (aminoacids, carb, lipid, na, cgenff, ethers, metals, silicates, solvent — ~2400+ residue types). Non-protein chains are built without terminal patches. Glycan trees are detected from C1-O distances and built with proper glycosidic bond handling (HO removal, charge redistribution, atom type changes). Output is modular: `ffparams.itp` (all FF parameters including atomtypes, bonded params, water model params), chain `.itp` files, `water.itp`, `ions.itp`, `posre_*.itp` (position restraints, `#ifdef POSRES`), and `interchain_ss.itp` (if needed). `topol.top` contains only `#include` directives, `[ system ]`, and `[ molecules ]`. Ions and buffer particles (BUF) in the input PDB are auto-detected and added to the `[ molecules ]` section.

## ACPYPE mode (`--acpype`)

Uses OpenMM to parametrize with AMBER14 + GLYCAM_06j-1, converts via ParmEd to AMBER prmtop/inpcrd, then ACPYPE generates GROMACS `topol.top`/`.gro` with `[ pairs_nb ]` directive for per-pair 1-4 parameters. This solves the mixed 1-4 scaling problem (AMBER fudgeLJ=0.5 vs GLYCAM fudgeLJ=1.0) that GROMACS cannot express globally. Output includes position restraints (`#ifdef POSRES` / `#include "posre_{stem}.itp"` / `#endif`) and water/ion moleculetypes ready for `gmx solvate`/`genion`. Best for glycoprotein systems. Ignores `--ff`/`--water`/`--merge` flags.

## Glycolipid support (CHARMM36)

Automatically detects glycolipids — ceramide residues (CER1, CER160, CER180, etc.) covalently bonded to a sugar tree — and builds them as a single GROMACS moleculetype. Designed for CHARMM-GUI PDB output which uses 4-character residue names (CER1, BGLC, BGAL, ANE5, AFUC, etc.).

**How it works:**
- Detects sugar O atom within 2.0 A of ceramide C1S (CHARMM-GUI removes ceramide O1/HO1 at the bond site)
- Builds ceramide from `lipid.rtp` + sugar tree from `carb.rtp` into one moleculetype
- Handles linkage charge redistribution: ceramide O1+HO1 charge goes to C1S, sugar HO1 charge goes to O1
- Sugar O1 at ceramide junction gets type OC301 (linear ether); sugar-sugar junctions get OC3C61 (cyclic ether)
- Sialic acid (Neu5Ac/ANE5AC) handled correctly: links via C2, O2 removed at linkage
- Auto-detects BGAL with N-acetyl atoms (N, HN, CT) and remaps to BGALNA

**Supported ceramides:** CER160 (d18:1/16:0), CER180, CER181, CER2, CER200, CER220, CER240, CER241, CER3E.

**Example** (Fucosyl-GM1 ganglioside from CHARMM-GUI):
```bash
dvbfixer top FucGM1_Charmm.pdb --ff charmm -o gmx_top/ -v
# Output: Glycolipid_ _1.itp (251 atoms, 256 bonds, charge -1.0)
# Residues: CER160 + BGLC + BGAL + BGALNA + AFUC + ANE5AC
```

**Note:** Glycolipid support is CHARMM36 only. AMBER/GLYCAM+Lipid21 cannot be combined for glycolipids (incompatible atom type namespaces, no cross-FF parameters).

## See also

- [`parametrize`](parametrize.md) — GAFF2 small molecules to `#include` from `topol.top`
- [`transplant`](transplant.md) — `--gromacs` uses the same ACPYPE export pipeline
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — CHARMM-GUI alternative path
