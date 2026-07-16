# Force fields in dvbfixer

[← README](../README.md) · [← command index](commands/index.md)

`dvbfixer` uses OpenMM's bundled force fields for every energy-based step (`prepare`, `minimize`, `protonate`, `pull`, `zbs`). To spare you from typing raw XML paths, all of those tools accept a short-name alias via `--ff <name>`, and by default they auto-detect the right FF from the input's residue names.

`top` (GROMACS topology generation) uses a **different** short-name namespace — see [Two `--ff` namespaces](#two---ff-namespaces) below.

## `--ff` for OpenMM-using tools

Applies to: `prepare`, `minimize`, `protonate`, `pull`, `zbs`.

### Short-name aliases

| Short name       | Expands to (OpenMM XMLs)                                                                              | When to use                                                                                          |
|------------------|-------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| `auto`           | (chosen by auto-detection — see below)                                                                | **Default** for every tool.                                                                          |
| `amber`          | `amber19/protein.ff19SB.xml`, `amber19/tip3p.xml`                                                     | Plain protein, ff19SB. Same as `amber19`.                                                            |
| `amber19`        | `amber19/protein.ff19SB.xml`, `amber19/tip3p.xml`                                                     | Explicit ff19SB.                                                                                     |
| `amber14`        | `amber14/protein.ff14SB.xml`, `amber14/tip3p.xml`                                                     | Older ff14SB. Needed as base for GLYCAM.                                                             |
| `amber+glycam`   | `amber14-all.xml`, `amber14/GLYCAM_06j-1.xml`, `amber14/tip3pfb.xml`                                  | Glycoproteins with GLYCAM-named residues (NLN/OLS/OLT, 3-char sugar codes like UYB/4YB/VMB).         |
| `amber14+glycam` | same as `amber+glycam`                                                                                 | Alias.                                                                                               |
| `amber+lipid`    | `amber14-all.xml`, `amber14/lipid17.xml`, `amber14/tip3p.xml`                                         | Membrane proteins (Lipid17 phospholipids).                                                           |
| `amber+nucleic`  | `amber14-all.xml`, `amber14/DNA.OL15.xml`, `amber14/RNA.OL3.xml`, `amber14/tip3p.xml`                 | Protein-nucleic acid complexes (OL15 DNA, OL3 RNA).                                                  |
| `charmm`         | `charmm36.xml`, `charmm36/water.xml`                                                                  | CHARMM36 (covers protein + carbohydrates + lipids + nucleic acids natively). Same as `charmm36`.     |
| `charmm36`       | `charmm36.xml`, `charmm36/water.xml`                                                                  | Explicit CHARMM36.                                                                                   |
| `charmm2024`     | `charmm36_2024.xml`, `charmm36/water.xml`                                                             | Newer CHARMM36 release (2024).                                                                       |

You can also pass **explicit XML paths** for backward compatibility:

```bash
dvbfixer minimize input.pdb --ff amber14-all.xml amber14/tip3p.xml -v
```

Any argument containing `.xml` or a path separator is treated as an explicit list and passed straight to `openmm.app.ForceField(*paths)` unchanged.

### Auto-detection rules

When you say `--ff auto` (the default) or don't pass `--ff`, dvbfixer scans the input PDB's residue names once and picks the FF from these markers:

**CHARMM markers** (any hit → `charmm`, unambiguous FF-prep signal):

- **Protonation variants**: `HSD`, `HSE`, `HSP`, `ASPP`, `GLUP`, `LSN`
- **CHARMM-GUI 4-char sugars**: `BGLC`, `AGLC`, `BMAN`, `AMAN`, `BGAL`, `AGAL`, `BFUC`, `AFUC`, `BGLCNA`, `AGLCNA`, `BGALNA`, `AGALNA`, `ANE5`, `BNE5`, `ANE5AC`, `BNE5AC`, `AIDO`, `BIDO`
- **Ceramides**: `CER1`, `CER160`, `CER180`, `CER181`, `CER2`, `CER200`, `CER220`, `CER240`, `CER241`, `CER3E`

**GLYCAM markers** (any hit → `amber+glycam`):

- **Glycoproteins**: `NLN`, `OLS`, `OLT`
- **Reducing-end caps**: `ROH`, `OME`, `TBT`, `CMET`
- **Sugar residues**: any 3-char resname matching the GLYCAM `[linkage][sugar][anomer]` pattern — e.g. `4YB`, `UYB`, `VMB`, `0YA`, `0SA`, `2MA`

**Ambiguous — PDB Chemical Component Dictionary sugar names**: `NAG`, `NDG`, `BMA`, `MAN`, `FUC`, `FUL`, `GAL`, `BGC`, `GLC`, `SIA` (etc.). These appear in raw crystal PDBs, in GLYCAM output *before* renaming, and in CHARMM-GUI output *before* CHARMM renaming — they tell dvbfixer nothing about which FF the rest of the file was prepared for. Neither `amber14+GLYCAM` nor `charmm36.xml` has templates for these bare PDB names. **No auto-selection is made from PDB sugar names alone** — you get the plain `amber` default with a warning to convert first:

```
FF: amber  (PDB-standard sugar name(s) detected (FUC, GAL, MAN) with no
  FF-specific markers — cannot auto-select. Run `dvbfixer convert
  --to-amber` (→ GLYCAM UYB/VMB/...) or `dvbfixer convert --to-charmm`
  (→ BGLCNA/BMAN/...) before this step, or pass --ff explicitly)
  → amber19/protein.ff19SB.xml amber19/tip3p.xml
```

**Precedence when both CHARMM and GLYCAM markers appear** (rare — hand-assembled structures): CHARMM wins, because CHARMM protonation names are unambiguous full-file-prep signals.

**Upgrade behaviour**: if you explicitly asked for `amber` (or `amber14` / `amber19`) but the input clearly has CHARMM or GLYCAM markers, dvbfixer *upgrades* your choice and prints why:

```
FF: amber+glycam  (upgraded from 'amber' → 'amber+glycam' because GLYCAM
   residue(s) detected (4YB, NLN, VMB))
  → amber14-all.xml amber14/GLYCAM_06j-1.xml amber14/tip3pfb.xml
```

Downgrading (e.g. `--ff charmm` on a plain-protein input) is respected — no auto-override.

### The startup banner

Every OpenMM-using tool prints two lines at startup so you always know what FF is in play:

```
FF: <alias>  (<reason if auto-selected or upgraded>)
  → <expanded XML list>
```

### Overriding

- `--ff auto` — explicit auto-detection (also the default).
- `--ff <short-name>` — pick a specific alias.
- `--ff <xml> [<xml> …]` — pass raw OpenMM XML paths verbatim.

### Adding a new short name

Edit `FF_ALIASES` in `src/dvbfixer/ffutils.py`:

```python
FF_ALIASES = {
    ...
    'my-ff': ['some/openmm.xml', 'some/water.xml'],
}
```

If the new FF has residue names that unambiguously identify it, extend the marker sets (`_CHARMM_MARKERS`, GLYCAM residue lists, etc.) and add a branch to `detect_ff_from_pdb`.

## `--ff` for `top` (GROMACS topology)

Applies to: `top` **only**. Different namespace — it doesn't load OpenMM XML files. Instead, it parses **bundled GROMACS FF directories** at `FF/amber99sb-ildn-lipid21.ff/` and `FF/charmm36_ljpme-jul2022.ff/` via the RTP parser in `rtp_parser.py`.

| `top --ff` value | What it loads                                          |
|------------------|--------------------------------------------------------|
| `amber`          | `FF/amber99sb-ildn-lipid21.ff/` — AMBER99SB-ILDN + Lipid21 |
| `charmm`         | `FF/charmm36_ljpme-jul2022.ff/` — CHARMM36 with LJ-PME    |

`top` also has `--acpype` mode, which uses OpenMM (AMBER14 + GLYCAM) → ParmEd → ACPYPE and ignores `--ff` entirely.

## Two `--ff` namespaces (side-by-side)

| Aspect                 | OpenMM tools (`prepare`, `minimize`, `protonate`, `pull`, `zbs`) | `top`                                          |
|------------------------|------------------------------------------------------------------|------------------------------------------------|
| Backend                | OpenMM `ForceField(*xmls)`                                       | GROMACS RTP parser (bundled `FF/*.ff/` dirs)   |
| Short names            | `auto`, `amber`, `amber+glycam`, `charmm`, `charmm2024`, …       | `amber`, `charmm`                              |
| Explicit path          | `--ff a.xml b.xml …`                                             | `--ff-dir /path/to/custom.ff/`                 |
| Auto-detection         | Yes (this doc)                                                   | No — user picks                                |
| `charmm` maps to       | `charmm36.xml` + water XML                                       | Full bundled CHARMM36 RTP directory            |

The two are separate because they consume completely different file formats: OpenMM parses XML; GROMACS parses `.rtp` / `.atp` / `.itp`. The bundled GROMACS FF dirs let dvbfixer emit topologies that don't need any external FF installation on the target machine.

## Water models

The OpenMM aliases pick a default water XML that matches the FF (e.g. `amber` uses `tip3p`, `amber+glycam` uses `tip3pfb`). If you need a different water model, pass the water XML explicitly:

```bash
dvbfixer minimize input.pdb --ff amber19/protein.ff19SB.xml amber19/opc.xml
```

For `top`, water is a separate `--water` argument (`tip3p|spc|spce|tip4p|tip4pew|opc`) and ions come from `--ion-set` — see [`top`](commands/top.md).

## Examples

```bash
# Default: auto-detect
dvbfixer minimize glycoprotein.pdb -v
# → FF: amber+glycam  (GLYCAM residue(s) detected (4YB, VMB))
#   → amber14-all.xml amber14/GLYCAM_06j-1.xml amber14/tip3pfb.xml

dvbfixer minimize plain_protein.pdb -v
# → FF: amber  (no non-standard residues detected)
#   → amber19/protein.ff19SB.xml amber19/tip3p.xml

# Force CHARMM
dvbfixer prepare input.pdb --ff charmm -v

# Force a newer CHARMM release
dvbfixer minimize input.pdb --ff charmm2024 -v

# Backward-compat: explicit XML paths still work
dvbfixer minimize input.pdb --ff amber14-all.xml amber14/tip3p.xml -v

# Whole pipeline with a specific FF
dvbfixer zbs input.pdb --ff amber+glycam -v
```
