# dvbfixer glycam — Convert between PDB/CHARMM and GLYCAM Nomenclature

[← command index](index.md) · [← README](../../README.md)

Bidirectional converter between standard PDB/CHARMM sugar naming and GLYCAM force field nomenclature.

**Forward (default)**: PDB → GLYCAM. Renames sugars from standard PDB codes (NAG, BMA, MAN, GAL, FUC, SIA, …) to GLYCAM 3-char codes `[linkage][sugar][anomer]` (UYB, VMB, 0MA, 6LB, 0fA, 0SA, …). Detects glycosidic bonds from CONECT records (or distance-based fallback), determines linkage patterns, renames atoms to GLYCAM convention (hydroxyl `HO3→H3O`; N-acetyl `C7→C2N`, `O7→O2N`, `C8→CME`, `H2N→HN2`; methyl `HT1→H1M` etc.). Sialic acid (SIA→0SA): full stereo-specific rename including methylene `H3/H32→H3A/H3E`, `H9/H92→H9R/H9S`, amide `HN5→H5N`, and methyl `H11/H112/H113→H1M/H2M/H3M`; the spurious `HO1B` (carboxylate H added by PDBFixer) is dropped because the GLYCAM template has no slot for it. Adds ROH cap at the reducing end unless `--no-roh`. Detects protein-linked glycans and renames `ASN→NLN`, `SER→OLS`, `THR→OLT`.

**Reverse (`--to-charmm`)**: GLYCAM → standard PDB / CHARMM-compatible. Strips GLYCAM linkage characters, inverts atom-name renames, drops ROH/OME caps, reverts `NLN/OLS/OLT → ASN/SER/THR`. Output uses standard 3-char PDB sugar codes (NAG/NDG/BMA/MAN/GAL/FUL/SIA/…) accepted by both CHARMM-GUI and `dvbfixer top --ff charmm` (the latter maps PDB → CHARMM RTP names via `PDB_TO_CARB`). Linkage information is preserved in CONECT records, not in the residue name.

Text-based — no OpenMM dependency. Handles input from PDB, CHARMM-GUI, or `dvbfixer prepare`.

**AMBER protonation-variant cleanup**: when an input PDB labels a residue with an AMBER variant name (`LYN`/`CYX`/`CYM`/`HID`/`HIE`) but still carries H atoms that the variant template doesn't have — common when a user manually renamed `LYS→LYN` to mark deprotonation but forgot to drop `HZ1` — glycam drops the extra H atoms during its rename pass so the output matches the AMBER template directly. Per-variant drops:

| Variant | Atoms dropped | Why |
|---------|---------------|-----|
| `LYN` | `HZ1`, `1HZ` | Deprotonated NZ has only HZ2 + HZ3 |
| `CYX` | `HG`, `HG1` | Disulfide-bonded SG has no H |
| `CYM` | `HG`, `HG1` | Deprotonated SG has no H |
| `HID` | `HE2` | HD1-only tautomer |
| `HIE` | `HD1` | HE2-only tautomer |

`HIP`, `ASH`, `GLH` are not in the table — they ADD an H rather than miss one. Vanilla `LYS`/`CYS`/`HIS` are untouched.

## Usage

```bash
# Forward: PDB → GLYCAM (writes input_glycam.pdb)
dvbfixer glycam glycan.pdb -v

# Reverse: GLYCAM → standard PDB / CHARMM-compatible (writes input_charmm.pdb)
dvbfixer glycam glycam_in.pdb --to-charmm -v

# Without ROH cap at reducing end
dvbfixer glycam glycan.pdb --no-roh

# Custom output
dvbfixer glycam glycan.pdb -o glycam_output.pdb -v
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_glycam.pdb` (or `_charmm.pdb` with `--to-charmm`) | Output file path |
| `--no-roh` | off | Do not add ROH cap at the reducing end (forward only) |
| `--to-charmm` | off | Reverse direction: GLYCAM → standard PDB / CHARMM-compatible naming |
| `-v`, `--verbose` | off | Print conversion details |

## GLYCAM Naming Convention

Each sugar residue gets a 3-character name: `[linkage][sugar][anomer]`

**Linkage code** (1st character): `0`=terminal, `2`-`9`=single position, `V`=O3+O6, `W`=O3+O4, `U`=O4+O6, `Z`=O2+O3, `X`=O2+O6, `Y`=O2+O4 (multi-linkage for branching sugars).

**Sugar code** (2nd character): `G`=glucose, `L`=galactose, `M`=mannose, `Y`=GlcNAc, `V`=GalNAc, `f`=fucose (lowercase=L-sugar), `S`=Neu5Ac, `X`=xylose, `R`=ribose, etc.

**Anomer code** (3rd character): `A`=alpha, `B`=beta.

## Example

```
BGC(res1, child at O4)     -> 4GB
GAL(res2, children at O3+O4) -> WLB
SIA(res6, terminal)        -> 0SA
NGA(res3, child at O3)     -> 3VB
GAL(res4, child at O2)     -> 2LB
FUC(res5, terminal)        -> 0fA
+ ROH cap at reducing end
```

## See also

- [`transplant`](transplant.md) — pair with `glycam` for the GLYCAM-Web glycoprotein workflow
- [`top`](top.md) — `--acpype` consumes GLYCAM-named input; `--ff charmm` consumes the `--to-charmm` output
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — full GLYCAM glycoprotein recipe
