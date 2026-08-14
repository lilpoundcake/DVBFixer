# dvbfixer protonate — PROPKA3 Protonation Assignment

[← command index](index.md) · [← README](../../README.md)

Runs PROPKA3 to predict per-residue pKa values, then renames titratable residues to their correct protonation state at the target pH and adds the corresponding hydrogen atoms. Existing hydrogens are stripped and re-added by OpenMM based on the renamed residue templates.

**FF-aware output naming** (since v0.3): with `--ff amber` (default) the output uses AMBER variant names (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN). With `--ff charmm` the output uses CHARMM36 equivalents (HSD/HSE/HSP, CYM, CYS) — AMBER-only variants that CHARMM36 doesn't ship a template for (ASH/GLH/LYN) are folded back to their standard parent (ASP/GLU/LYS) with a WARNING that the PROPKA-requested charge state can't be expressed with the shipped CHARMM XML. Internally, `addHydrogens` always runs with AMBER-style names because that's the only vocabulary OpenMM's `hydrogens.xml` speaks; the CHARMM output rewrite happens after H addition finishes.

**Text-level variant pre-rename** (fixes 2BNQ-class crashes): OpenMM's `PDBFile` parser only recognises standard AA residue names when inferring peptide bonds. A mid-chain LYN/HIE/HSE (etc.) silently blocks the peptide bond from the previous residue's C to its N, and downstream `addHydrogens` then fails on the ADJACENT residue with a misleading `No template found for residue X (ASN). ... missing 1 C atom externally bonded` error. Protonate now rewrites all AMBER + CHARMM variant names → standard parent (HIS/ASP/GLU/CYS/LYS) in the raw PDB TEXT before OpenMM parses. Variant names are restored on the output topology after `addHydrogens` completes so the final PDB carries HID/HIE/HSE/LYN etc. correctly.

**Input sanitization**: `ffutils.sanitize_protein_hetatm` runs first — rewrites protein-residue HETATM → ATOM lines and drops spurious mid-chain TER records. Same helper `prepare` uses. No-op on clean inputs.

**Unsupported PROPKA-titratable groups**: PROPKA also reports pKas for TYR (deprotonated tyrosinate) and ARG (neutral arginine). Neither AMBER14/19 nor CHARMM36 ships a TYD/TYN/ARN variant template — there's no way to add those H atoms via the stock `addHydrogens` flow. In `--verbose` mode you'll see them logged as "Unsupported PROPKA-titratable groups at pH X"; no rename is attempted. Terminal N+/C- pKas are handled automatically by OpenMM's NXXX/CXXX terminal patches.

**GLYCAM support**: PROPKA3 doesn't recognize GLYCAM glycoprotein residues (NLN/OLS/OLT), so they are renamed to ASN/SER/THR in a temp PDB (heterogens stripped) before PROPKA runs; pKa results are mapped back to the original input. When GLYCAM residues are detected, the FF auto-switches to AMBER14 + GLYCAM_06j-1 + tip3pfb (ff19SB has no GLYCAM templates → crash); `add_glycam_bonds(positions=...)` populates the missing bonds and `glycam-hydrogens.xml` provides the H definitions. PROPKA renames are filtered out for NLN/OLS/OLT positions (sugar-bonded sidechains, different chemistry). After writing, `fix_atom_hetatm_records` rewrites HETATM lines for protein residues back to ATOM.

**Diagnostic on `addHydrogens` failure**: wrapped with `explain_template_error` so failures surface chain/resseq/icode + atom-set diff + prev/next residue instead of OpenMM's bare topology-index message.

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

### Cysteine: `CYS`, `CYM`, and `CYX`

`--cys-disulfide-pka` does **not** set or override a cysteine's predicted
pKa. It only recognizes PROPKA's unusually high disulfide sentinel. With the
default value `99.99`, each cysteine is classified independently as follows:

- predicted pKa below `--ph` → `CYM`, a deprotonated free thiolate with no SG
  hydrogen;
- predicted pKa at or above `--cys-disulfide-pka` → `CYX`, the AMBER
  disulfide state;
- otherwise → ordinary protonated `CYS` with an SG hydrogen.

Explicit/inferred SG–SG connectivity takes precedence and assigns both bonded
cysteines as `CYX`. To disable all pKa-driven states, including `CYM`, use
`--no-propka`; changing `--cys-disulfide-pka` does not disable thiolate
assignment.

OpenMM calls every cysteine hydrogen state without SG–H `CYX`, even when the
chemical state is a free thiolate. Dvbfixer therefore passes `CYM` internally
as `CYX` only while placing hydrogens, then restores `CYM` before selecting the
AMBER force-field template and writing the output. Thus a final `CYM` name is
intentional and retains the thiolate charge model; it is not evidence that the
disulfide threshold was ignored.

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
| `--cys-disulfide-pka` | 99.99 | Disulfide-sentinel cutoff for CYS → CYX; does not override predicted pKa or disable pKa < pH → CYM. |
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

## How it works
Runs PROPKA3 for pKa prediction, renames residues to AMBER protonation names (HID/HIE/HIP, ASH, GLH, CYM, CYX, LYN) based on target pH. `--no-hydrogens` for text-based rename only (used in pipeline). Default mode strips H and re-adds via OpenMM with `variants` parameter. **Full GLYCAM support**: PROPKA3 doesn't know GLYCAM glycoprotein residues (NLN/OLS/OLT) — `_sanitize_for_propka` writes a temp PDB with NLN→ASN / OLS→SER / OLT→THR rename and heterogens stripped, runs PROPKA on that, maps results back to the original input keyed by (chain, resnum, icode). `_scan_glycam_residues` in `main()` auto-switches `args.ff` from ff19SB to `['amber14-all.xml', 'amber14/GLYCAM_06j-1.xml', 'amber14/tip3pfb.xml']` when GLYCAM residues are detected and the user didn't override `--ff` (ff19SB has no GLYCAM templates → crash). PROPKA protonation renames for NLN/OLS/OLT positions are FILTERED OUT (different chemistry — sidechain N/O is sugar-bonded). AMBER variant names (HID/HIE/HIP/CYX/etc.) already present in input are preserved by scanning raw PDB text before `PDBFile` normalizes them. `_add_hydrogens_to_output` GLYCAM mode: keep heterogens in topology, build FF via `create_forcefield_with_openff`, load `glycam-hydrogens.xml`, run `add_glycam_bonds(positions=...)` for intra-residue + protein-glycan + sugar-sugar bonds, then `Modeller.addHydrogens(forcefield, pH, variants=...)`. Calls `fix_atom_hetatm_records` on final output so HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/NLN/OLS/OLT are written as ATOM records.

**ProtAssign-style optimisation (`--protassign` / `--no-protassign`, default ON)**: wraps MolProbity Reduce (Word, Lovell & Richardson 1999) to pick HIS tautomers and detect ASN/GLN side-chain flips from local H-bond network + clash scoring. **Runs by default on every protonate invocation** (since June 2026). `argparse.BooleanOptionalAction` generates both `--protassign` (explicit ON) and `--no-protassign` (disable, PROPKA-only mode). Reduce binary at `/tmp/mamba/envs/dvbfixer/bin/reduce` (v4.10) — bundled with AmberTools. `_find_reduce_binary` searches override→PATH→env bin dir. `_run_reduce` invokes `reduce -build -flip -quiet -noheterogens` and captures stdout. `_parse_reduce_decisions` determines HIS tautomer by which H atoms Reduce placed (HD1+HE2→HIP, HD1→HID, HE2→HIE) and detects ASN/GLN flips by coordinate diff (OD1_reduce within 0.1 Å of ND2_input → flipped; same for OE1/NE2). `_apply_flips_to_pdb_text` swaps OD1↔ND2 / OE1↔NE2 / HD21↔HD22 / HE21↔HE22 coordinate fields in raw PDB text before the file is passed to `_add_hydrogens_to_output`. **Merge logic**: Reduce's HIS picks overlay PROPKA's renames, but PROPKA's HIP wins (pKa-driven charged-state decision more reliable than Reduce's local count). GLYCAM positions (NLN/OLS/OLT) are filtered from Reduce's picks/flips. Missing reduce binary → hard error mentioning `conda install ambertools` AND `--no-protassign` (escape hatch). The default flip from OFF→ON is a deliberate behaviour change; scripts that depend on the old byte-output behaviour must add `--no-protassign` explicitly.

**Input sanitization + text-level variant pre-rename** (fixes crash on 2BNQ-style inputs where a mid-chain LYN blocks OpenMM's peptide-bond inference): `main()` calls `ffutils.sanitize_protein_hetatm(input_path)` first — rewrites protein-residue HETATM → ATOM lines (no-op on clean input, shared with `prepare`). Then `_add_hydrogens_to_output` calls `_text_rename_variants_to_parent(input_path)` which rewrites LYN/HID/HIE/HIP/ASH/GLH/CYX/CYM (+ CHARMM's HSD/HSE/HSP/ASPP/GLUP/LSN) → standard parent (LYS/HIS/…) in the raw PDB TEXT before OpenMM parses. Reason: OpenMM's `PDBFile` parser only recognises standard AA names when inferring peptide bonds — a mid-chain LYN blocks the peptide bond FROM the previous residue's C TO the LYN's N, and downstream `addHydrogens` then fails on the adjacent residue (typically the ASN/whatever BEFORE the LYN) with a confusing "residue X missing 1 C atom externally bonded" error. The rename returns a saved map `(chain, resid, icode) → original_variant` which is merged into `_saved` before addHydrogens; `_restore_variants_post_addhydrogens` restores the variant names on the OUTPUT topology so LYN/HIE/etc. survive to the final PDB. `addHydrogens` also now wrapped with `explain_template_error` for chain/resseq/atom-set diagnostics on failure (instead of OpenMM's bare topology-index message).

**FF-aware output naming (`--ff charmm`)**: `_remap_amber_variants_to_charmm_in_pdb` rewrites the AMBER-style variant names in the final output PDB to CHARMM36 equivalents: HID→HSD, HIE→HSE, HIP→HSP, CYX→CYS (SS via SSBOND), CYM→CYM. ASH/GLH/LYN have NO OpenMM-`charmm36.xml` template — those get folded back to standard ASP/GLU/LYS with a clear WARNING that the PROPKA-requested charge state can't be expressed with the shipped CHARMM XML. `--ff amber` (default) leaves AMBER names as-is. All rewrites use `f"{name:<3s}"` (fixed 3-char left-align) so column alignment is preserved. Selection driven by `_is_charmm_ff(args.ff)` which scans the resolved XML list for `charmm`.

**Unsupported PROPKA titratable groups**: PROPKA also reports pKas for TYR (deprotonated tyrosinate) and ARG (neutral arginine). Neither AMBER14/19 nor CHARMM36 ships a TYD/TYN/ARN variant template — no way to add those H atoms via the stock `addHydrogens` flow. `decide_protonation` skips them; `main()` logs them as "Unsupported PROPKA-titratable groups at pH X" in `--verbose` mode so users know. Terminal N+/C- pKas are handled automatically by OpenMM's NXXX/CXXX terminal patches.

## Batch mode

`protonate` supports directory input with a common pH and force field:
`dvbfixer protonate --input-dir structures --output-dir protonated --ph 7.4`.
See [Batch mode](../batch-mode.md) for shared keys.
