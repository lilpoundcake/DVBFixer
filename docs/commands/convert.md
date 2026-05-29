# dvbfixer convert — Convert between PDB/AMBER/GLYCAM and CHARMM naming

[← command index](index.md) · [← README](../../README.md)

> Previously called `dvbfixer glycam` — the old name still works and emits a deprecation notice.

Bidirectional converter for sugar nomenclature AND protein protonation variants. Both directions are idempotent — running on already-correct input is a no-op.

**Default direction (`--to-amber`)** — converts PDB / CHARMM-named input to AMBER-friendly naming:
- Sugars: PDB codes (NAG, BMA, MAN, GAL, FUC, SIA, …) → GLYCAM 3-char codes `[linkage][sugar][anomer]` (UYB, VMB, 0MA, 6LB, 0fA, 0SA, …)
- Glycoprotein residues: `ASN→NLN`, `SER→OLS`, `THR→OLT` when sidechain-bonded to a sugar
- Protonation variants: CHARMM `HSD→HID`, `HSE→HIE`, `HSP→HIP`, `ASPP→ASH`, `GLUP→GLH`, `LSN→LYN` (and LSN's `HZ1→HZ2`, `HZ2→HZ3` to match the AMBER LYN template)
- Atom names: hydroxyl `HO3→H3O`, N-acetyl `C7→C2N`, `O7→O2N`, `C8→CME`, methyl `HT1→H1M`, etc.
- Sialic acid (`SIA→0SA`): stereo-specific atom renames including `H3/H32→H3A/H3E`, `H9/H92→H9R/H9S`, `HN5→H5N`, methyl `H11/H112/H113→H1M/H2M/H3M`; the spurious `HO1B` (PDBFixer-added) is dropped
- ROH cap added at the reducing end (suppress with `--no-roh`)

Output is consumable by `prepare`, `minimize`, `protonate`, `top --acpype`.

**Reverse direction (`--to-charmm`)** — converts AMBER/GLYCAM-named input to CHARMM-friendly naming:
- Sugars: GLYCAM codes → standard 3-char PDB / CHARMM-GUI 4-char codes (NAG/NDG/BMA/MAN/GAL/FUL/SIA/… or BGLCNA/BMAN/AMAN/ANE5AC/… where `top --ff charmm` maps via `PDB_TO_CARB`)
- Glycoprotein residues: `NLN/OLS/OLT → ASN/SER/THR`; ROH/OME caps dropped
- Protonation variants: AMBER `HID→HSD`, `HIE→HSE`, `HIP→HSP`, `ASH→ASPP`, `GLH→GLUP`, `LYN→LSN` (and LYN's `HZ2→HZ1`, `HZ3→HZ2` to match the CHARMM LSN convention); `CYX→CYS` (CHARMM uses DISU patch via SSBOND)
- CYM stays as CYM (CHARMM36 has a `[ CYM ]` residue)
- Linkage information preserved via CONECT records, not in the residue name

Output is consumable by `top --ff charmm`.

**Accidental wrong-FF naming** is handled automatically: if the input has the wrong FF's convention (e.g. AMBER LYN in a structure being converted with `--to-charmm`), it's renamed to the target convention. If the input is already in the target convention (e.g. CHARMM LSN + `--to-charmm`), nothing changes.

Text-based — no OpenMM dependency. Handles input from PDB, CHARMM-GUI, or `dvbfixer prepare`.

**Stale-H cleanup**: when an input labels a residue with an AMBER variant name but still carries H atoms the AMBER template doesn't have (e.g. a user manually renamed `LYS→LYN` without dropping `HZ1`), `convert` drops the extra H atoms during its rename pass:

| Variant | Atoms dropped | Why |
|---------|---------------|-----|
| `LYN` | `HZ1`, `1HZ` | AMBER LYN keeps HZ2 + HZ3 only |
| `CYX` | `HG`, `HG1` | Disulfide-bonded SG has no H |
| `CYM` | `HG`, `HG1` | Deprotonated SG has no H |
| `HID` | `HE2` | HD1-only tautomer |
| `HIE` | `HD1` | HE2-only tautomer |

`HIP`, `ASH`, `GLH` are not in the table — they ADD an H rather than miss one. Vanilla `LYS`/`CYS`/`HIS` are untouched.

## Usage

```bash
# Default direction: PDB/CHARMM → GLYCAM + AMBER (writes input_amber.pdb)
dvbfixer convert input.pdb -v

# Explicit (same as default)
dvbfixer convert input.pdb --to-amber -v

# Reverse: GLYCAM/AMBER → CHARMM (writes input_charmm.pdb)
dvbfixer convert input.pdb --to-charmm -v

# Without ROH cap at reducing end (default direction only)
dvbfixer convert glycan.pdb --no-roh

# Custom output
dvbfixer convert input.pdb -o output.pdb -v
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_amber.pdb` (or `_charmm.pdb` with `--to-charmm`) | Output file path |
| `--to-amber` | (default) | PDB/CHARMM → GLYCAM (sugars) + AMBER (protonation variants). Mutually exclusive with `--to-charmm`. |
| `--to-charmm` | off | Reverse direction: GLYCAM/AMBER → CHARMM. |
| `--no-roh` | off | Do not add ROH cap at the reducing end (default direction only) |
| `-v`, `--verbose` | off | Print conversion details |

## GLYCAM Naming Convention

Each sugar residue gets a 3-character name: `[linkage][sugar][anomer]`

**Linkage code** (1st character): `0`=terminal, `2`-`9`=single position, `V`=O3+O6, `W`=O3+O4, `U`=O4+O6, `Z`=O2+O3, `X`=O2+O6, `Y`=O2+O4 (multi-linkage for branching sugars).

**Sugar code** (2nd character): `G`=glucose, `L`=galactose, `M`=mannose, `Y`=GlcNAc, `V`=GalNAc, `f`=fucose (lowercase=L-sugar), `S`=Neu5Ac, `X`=xylose, `R`=ribose, etc.

**Anomer code** (3rd character): `A`=alpha, `B`=beta.

## Protonation-variant mapping table

| State | AMBER (`--to-amber`) | CHARMM (`--to-charmm`) | Glycam-stage atom changes |
|---|---|---|---|
| Protonated K (cationic) | LYS | LYS | none |
| Neutral K | LYN | LSN | `--to-charmm`: HZ2→HZ1, HZ3→HZ2 (and reverse) |
| HIS, ND1 only | HID | HSD | none |
| HIS, NE2 only | HIE | HSE | none |
| HIS, doubly protonated | HIP | HSP | none |
| Protonated D | ASH | ASPP | none |
| Protonated E | GLH | GLUP | none |
| Disulfide C | CYX | CYS (+ DISU patch) | none |
| Deprotonated C | CYM | CYM | none |

Methylene H renames (HB2/HB3 ↔ HB1/HB2 etc.) and backbone H ↔ HN are NOT applied here — `top.py`'s ARN + methylene-shift passes handle them during topology generation.

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

- [`transplant`](transplant.md) — pair with `convert` for the GLYCAM-Web glycoprotein workflow
- [`top`](top.md) — `--acpype` consumes default `convert` output; `--ff charmm` consumes the `--to-charmm` output
- [`prepare`](prepare.md) / [`minimize`](minimize.md) / [`protonate`](protonate.md) — all consume the default (AMBER) direction
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — full GLYCAM glycoprotein recipe
