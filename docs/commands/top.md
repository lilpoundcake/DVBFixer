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
| `--water` | `tip3p` | Water model: `tip3p`, `spc`, `spce`, `tip4p`, `tip4pew`, `opc`. With `--ff charmm` only `tip3p`/`spc`/`spce` are accepted (the CHARMM-tuned variants); `opc`/`tip4p`/`tip4pew` with CHARMM is rejected because CHARMM ions are fitted to CHARMM-TIP3P |
| `--ion-set` | `auto` | Ion LJ parameter set: `auto` picks the set matched to the water model (`jc-tip3p` for TIP3P, `jc-spce` for SPC/SPCE, `jc-tip4pew` for TIP4P/TIP4P-Ew, `lm-hfe-opc` for OPC). Override choices: `jc-tip3p`, `jc-spce`, `jc-tip4pew`, `lm-hfe-opc`, `lm-iod-opc`, `dang-legacy` (bundled Aqvist Na⁺/Dang Cl⁻ — pre-2008 behaviour, for backwards-compat). AMBER only — ignored with `--ff charmm`. Covers Na⁺/K⁺/Cl⁻/Ca²⁺/Mg²⁺/Zn²⁺ |
| `--ignh` | off | Ignore hydrogens in input PDB (strip all H and let the FF templates rebuild them) |
| `--keep-all-hydrogens` | off | Preserve every input H atom (default OFF: `HO1`/`HO2`/`HO3`/`HO4`/`HO6` at glycosidic linkage sites are stripped and their charge is redistributed onto the linked O — matches the CHARMM RTP template for a formed glycosidic bond). See warning below |
| `--ss` | auto | Disulfide bond: CHAIN1:NUM1:CHAIN2:NUM2 (repeatable) |
| `--his` | auto | HIS protonation: CHAIN:NUM:STATE (HIE/HID/HIP, repeatable) |
| `--protonate` | off | Protonate residues: `all` for every ASP/GLU/HIS, or `CHAIN:NUM[:STATE],...` for specific. H placed via OpenMM Modeller with CHARMM/AMBER FF |
| `--merge` | off | Merge all chains into single moleculetype |
| `--pdb` | `conf.pdb` | Output PDB with topology-matched atom names |
| `--acpype` | off | Use ACPYPE pipeline (AMBER14+GLYCAM -> ParmEd -> GROMACS) with per-pair 1-4 scaling |
| `-v`, `--verbose` | off | Print detailed progress |

## Preserving all input hydrogens (`--keep-all-hydrogens`)

By default, `dvbfixer top` implements the CHARMM glycosidic-bond
chemistry when it detects a sugar-sugar (or ceramide-sugar) linkage:
it strips the reducing-end `HO1`/`HO2`/`HO3`/`HO4`/`HO6` at the linked
position and redistributes that H's charge onto the linked O (which
also gets its atom type changed from `OC311` hydroxyl to `OC3C61`
ether). Similarly, ceramide `HO1` is dropped at the ceramide-sugar
junction. This matches the RTP template for a formed glycosidic bond
and gives an integer residue net charge.

`--keep-all-hydrogens` skips both the H strip AND the charge/type
redistribution — the input `HO*` atoms pass through verbatim and the
O keeps its `OC311` hydroxyl type + full RTP charge. Use this when:

- The input is a free reducing-end sugar (no glycosidic link at that
  position) whose H is REAL, not vestigial.
- You need charges to round-trip untouched to a downstream tool that
  expects the un-redistributed sum.

```bash
# Default — HO atoms stripped at detected linkages
dvbfixer top glycan.pdb --ff charmm -o gmx/

# Preserve every input H
dvbfixer top glycan.pdb --ff charmm --keep-all-hydrogens -o gmx/
```

**Warning**: at a REAL glycosidic linkage this produces an over-valent O
(the H is bonded to it AND the neighbouring sugar's anomeric C bonds to
the same O). `gmx grompp` may warn and the resulting MD energy will be
chemically wrong. Use only when you know the linkage detection is a
false positive or when your downstream tool expects the un-modified
charges. The default (`--keep-all-hydrogens` off) is the right choice
for normal glycoprotein / glycolipid MD.

Also honoured by `--acpype` mode: skips the `_GLYCAM_KEEP_H` allowlist
filter in `acpype_export.prepare_for_openmm` so any GLYCAM-side H that
was going to be dropped (e.g. `HD22` on `NLN`) passes through.

## RTP-based mode (default)

Parses force field RTP files to build topology from bond graph: resolves inter-residue bonds (`-C`, `+N`), enumerates angles/dihedrals/pairs algorithmically, copies impropers and CMAP from templates, applies terminal patches. Intra-chain disulfide bonds (SG-SG between CYS2 residues) are added explicitly to the chain topology with all derived angles, dihedrals, and 1-4 pairs; inter-chain SS bonds go to `interchain_ss.itp`. For CHARMM36, loads all molecule-type RTP files (aminoacids, carb, lipid, na, cgenff, ethers, metals, silicates, solvent — ~2400+ residue types). Non-protein chains are built without terminal patches. Glycan trees are detected from C1-O distances and built with proper glycosidic bond handling (HO removal, charge redistribution, atom type changes). Output is modular: `ffparams.itp` (all FF parameters including atomtypes, bonded params, water model params), chain `.itp` files, `water.itp`, `ions.itp`, `posre_*.itp` (position restraints, `#ifdef POSRES`), and `interchain_ss.itp` (if needed). `topol.top` contains only `#include` directives, `[ system ]`, and `[ molecules ]`. Ions and buffer particles (BUF) in the input PDB are auto-detected and added to the `[ molecules ]` section.

## Water-matched ions (AMBER)

By default, `--water` only swapped the water moleculetype; ion atom types were
loaded unchanged from the bundled `ffnonbonded.itp` (Aqvist Na⁺, Dang Cl⁻ — fit
for TIP3P). This caused stability problems when users combined newer water models
(OPC, TIP4P-Ew) with the legacy ions.

`dvbfixer top` now swaps the ion Lennard-Jones parameters to match the water
model. Set `--ion-set auto` (default) and the right set is chosen automatically:

| `--water` | Default ion set | Source |
|-----------|-----------------|--------|
| `tip3p`   | `jc-tip3p`      | Joung & Cheatham 2008 monovalents + Li-Merz 2013 12-6 HFE divalents |
| `spce`    | `jc-spce`       | same papers, SPC/E fit |
| `spc`     | `jc-spce` (alias, warns) | plain SPC was not parametrized by JC — using SPC/E |
| `tip4pew` | `jc-tip4pew`    | same papers, TIP4P-Ew fit |
| `tip4p`   | `jc-tip4pew` (alias, warns) | original TIP4P was not parametrized by JC — using TIP4P-Ew |
| `opc`     | `lm-hfe-opc`    | Sengupta-Li-Merz 2021 monovalents + Li-Song-Merz 2020 12-6 HFE divalents (OPC) |

`--ion-set` can override the default. `lm-iod-opc` swaps to Li-Merz's
ion-oxygen-distance-fit set (better first-shell RDF for OPC, slightly worse
hydration free energy). `dang-legacy` keeps the pre-2008 Aqvist/Dang values
for backwards compatibility with older runs.

Covered ions: **Na⁺, K⁺, Cl⁻, Ca²⁺, Mg²⁺, Zn²⁺**.

**CHARMM** uses its own ions (SOD/CLA/POT/CAL/MGA in the bundled `ions.itp`),
fitted to CHARMM-TIP3P. `--ion-set` is ignored with `--ff charmm`. Using
`--water opc/tip4p/tip4pew` with `--ff charmm` is rejected at the CLI level —
combine OPC with `--ff amber` instead.

OPC water itself is a new water choice. The bundled `FF/amber99sb-ildn-lipid21.ff/opc.itp`
implements the 4-site OPC water (Izadi-Anandakrishnan-Onufriev 2014) with the
M-site as a `virtual_sites3` particle.

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
