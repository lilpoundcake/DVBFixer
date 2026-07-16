# dvbfixer protonate — PROPKA3 Protonation Assignment

[← command index](index.md) · [← README](../../README.md)

Runs PROPKA3 to predict per-residue pKa values, then renames titratable residues to their correct protonation state at the target pH and adds the corresponding hydrogen atoms. Uses AMBER force field naming conventions (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN). Existing hydrogens are stripped and re-added by OpenMM based on the renamed residue templates.

**GLYCAM support**: PROPKA3 doesn't recognize GLYCAM glycoprotein residues (NLN/OLS/OLT), so they are renamed to ASN/SER/THR in a temp PDB (heterogens stripped) before PROPKA runs; pKa results are mapped back to the original input. When GLYCAM residues are detected, the FF auto-switches to AMBER14 + GLYCAM_06j-1 + tip3pfb (ff19SB has no GLYCAM templates → crash); `add_glycam_bonds(positions=...)` populates the missing bonds and `glycam-hydrogens.xml` provides the H definitions. PROPKA renames are filtered out for NLN/OLS/OLT positions (sugar-bonded sidechains, different chemistry). AMBER variant names (HID/HIE/HIP/CYX/etc.) already present in the input are preserved by scanning the raw PDB text before OpenMM normalizes them. After writing, `fix_atom_hetatm_records` rewrites HETATM lines for protein residues back to ATOM.

## Protonation Logic

| Residue | Condition | Renamed to | Description |
|---------|-----------|------------|-------------|
| HIS | pKa > pH | HIP | Doubly protonated (charged) |
| HIS | pKa < pH | HIE (default) | Neutral, Ne2 protonated |
| HIS | pKa < pH | HID (option) | Neutral, Nd1 protonated |
| ASP | pKa > pH | ASH | Protonated (neutral) |
| GLU | pKa > pH | GLH | Protonated (neutral) |
| CYS | pKa >= 90 | CYX | Disulfide bridge |
| CYS | pKa < pH | CYM | Deprotonated thiolate |
| LYS | pKa < pH | LYN | Neutral |

## Usage

```bash
# Basic usage — writes input_prot.pdb
dvbfixer protonate input.pdb

# Custom pH
dvbfixer protonate input.pdb --ph 6.5

# Full pKa summary table
dvbfixer protonate input.pdb --summary

# Show non-standard protonation changes
dvbfixer protonate input.pdb -v

# Use HID as default neutral histidine tautomer
dvbfixer protonate input.pdb --his-default HID
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_prot.pdb` | Output file path |
| `--ph` | 7.0 | Target pH |
| `--his-default` | HIE | Default neutral HIS tautomer (HIE or HID) |
| `--cys-disulfide-pka` | 90.0 | pKa threshold for CYS -> CYX assignment |
| `--protassign` / `--no-protassign` | **ON** | Run MolProbity Reduce to optimise HIS tautomers and detect ASN/GLN side-chain flips (see below). Pass `--no-protassign` to disable |
| `--protassign-binary` | auto | Override path to the `reduce` binary |
| `--no-hydrogens` | off | Only rename residues, do not add/fix hydrogen atoms |
| `--ff` | `auto` | Force field for hydrogen addition. Accepts a short name (`auto`, `amber`, `amber+glycam`, `charmm`, …) or explicit OpenMM XML paths. See [force-fields.md](../force-fields.md). |
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
| `--summary` | off | Print full pKa table |
| `-v`, `--verbose` | off | Print non-standard protonation changes |

## ProtAssign-style optimisation (`--protassign`, default ON)

PROPKA only predicts pKa — it doesn't pick between the two neutral HIS
tautomers (HID vs HIE) or detect ASN/GLN side-chain flips. dvbfixer
wraps **MolProbity Reduce** (Word, Lovell & Richardson 1999) to make
those decisions from local H-bond geometry and van der Waals clash
scoring — analogous to Schrödinger's ProtAssign preprocessor.

**This runs by default on every `protonate` invocation.** Pass
`--no-protassign` to disable (PROPKA-only mode, the legacy behaviour).

What it does:

- **HIS tautomer choice** — for each HIS, Reduce decides whether HID
  (Nδ1-protonated), HIE (Nε2-protonated), or HIP (both, charged) gives
  the best local H-bond network. dvbfixer overlays Reduce's choice on
  top of PROPKA's pKa-driven rename. **PROPKA still wins HIP** (pKa-driven
  charged-state decision is more reliable than Reduce's local count);
  Reduce wins HID vs HIE picks for neutral tautomers.
- **ASN flip** — for each ASN, Reduce evaluates whether swapping OD1
  and ND2 atom positions improves the H-bond network. ~20% of PDB
  entries have ASN/GLN flipped (Popowicz 2007, NQ-Flipper); the
  carbonyl O and amide N have nearly indistinguishable electron density
  at typical crystallographic resolution. If Reduce flips it, dvbfixer
  swaps the OD1/ND2 coordinates in the input before adding hydrogens
  so the HD21/HD22 protons end up on the correct nitrogen.
- **GLN flip** — same as ASN, on OE1/NE2.

```bash
# Default: ProtAssign runs (HIS tautomer + ASN/GLN flip detection)
dvbfixer protonate input.pdb -v
# -v shows which HIS got re-tautomerised + how many ASN/GLN got flipped

# Disable ProtAssign (PROPKA-only mode — pH-driven HIS tautomers,
# no flip detection, no Reduce binary required)
dvbfixer protonate input.pdb --no-protassign -v
```

Example output with `-v`:

```
Running MolProbity Reduce for HIS tautomers + ASN/GLN flip optimisation...
  --protassign: 7 HIS tautomer override(s) from Reduce
  --protassign: 3 ASN flip(s), 6 GLN flip(s)
    ASN flip at C:387
    ASN flip at C:392
    ASN flip at L:210
    GLN flip at A:38
    ...
```

**Requires** the `reduce` binary. It's bundled with AmberTools (already
in the dvbfixer env). If missing, install via
`conda install -c conda-forge ambertools`, pass `--protassign-binary
PATH` to point at a local build, or pass `--no-protassign` to skip
the optimisation (PROPKA-only mode).

**When to use**:
- Crystal structures (X-ray, EM) where ASN/GLN orientations may be
  mis-assigned at the deposit stage.
- Any time the H-bond network of the starting structure matters
  (MD initialization, docking, free-energy calculations).

**When NOT to use**:
- If your input has been through Schrödinger PrepWizard or
  MolProbity Reduce already (no benefit — same algorithm).
- Cyclic peptides / non-standard residues (Reduce ignores them; no harm).

**Reference**: Word JM, Lovell SC, Richardson JS, Richardson DC (1999).
"Asparagine and glutamine: using hydrogen atom contacts in the choice
of side-chain amide orientation." *J Mol Biol* 285, 1735.
[PubMed 9917408](https://pubmed.ncbi.nlm.nih.gov/9917408/).

## See also

- [`prepare`](prepare.md) — run before `protonate` for structure repair
- [`minimize`](minimize.md) — re-runs after `protonate` to relax hydrogens
- [`rename`](rename.md) — reverse AMBER names back to canonical PDB names
