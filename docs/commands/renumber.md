# dvbfixer renumber — SEQRES- or Antibody-Scheme Residue Renumbering

[← command index](index.md) · [← README](../../README.md)

Renumbers residues with one of two strategies:

- **Default (`--scheme seqres`)** — aligns ATOM records to the SEQRES section via subsequence matching. Removes insertion codes, makes resseq sequential per chain, preserves gap positions.
- **Antibody (`--scheme kabat|chothia|imgt|martin|aho|eu`)** — applies an antibody numbering convention via ANARCI for V-domains and bundled EU reference sequences for C-domains. Works on incomplete chains (Fv-only, Fc-only, partial domains). Falls back to SEQRES for chains that aren't antibodies. Fully local — no web service.

Updates **all** PDB sections that reference residue numbers:
- ATOM, HETATM, TER
- HELIX, SHEET, SSBOND, LINK, CISPEP
- HET, DBREF, SEQADV
- CONECT (atom serial remapping)
- REMARK 465, 500, 610

## Antibody schemes

| Scheme | V-domain | C-domain | Notes |
|--------|----------|----------|-------|
| `kabat`   | Kabat (1991) via ANARCI | EU positions (Edelman 1969) | CDR insertions at 27/52/82/100 in heavy. |
| `chothia` | Chothia via ANARCI | EU | Same as Kabat with CDR loop boundaries adjusted. |
| `imgt`    | IMGT (1-128) via ANARCI | EU shifted to avoid V/C collision | CDR1=27-38, CDR2=56-65, CDR3=105-117. |
| `martin`  | Martin / Honnegger via ANARCI | EU | Structural alignment basis. |
| `aho`     | Aho (1-149) via ANARCI | EU shifted | Same as Martin in practice. |
| `eu`      | Kabat (V-domain EU = Kabat) | EU | Canonical IgG numbering for whole chain. |

C-domain numbering is always EU — Kabat/Chothia/Martin/Aho don't define C-domain positions, so EU is used regardless of the V-scheme. When the V-scheme extends past EU position 117 (IMGT/Martin/Aho), the EU C-domain numbers are shifted forward to avoid collisions and a warning is printed.

ANARCI is an optional dependency: only required when `--scheme` is non-default. The default seqres path has no extra requirements.

## Usage

```bash
# Basic usage (SEQRES-based) — writes input_renum.pdb
dvbfixer renumber input.pdb

# Antibody Kabat numbering (CDR insertions at H100A/B/C etc.)
dvbfixer renumber antibody.pdb --scheme kabat -v

# EU numbering across the full IgG (V 1-113, CH1 118-215, etc.)
dvbfixer renumber igg.pdb --scheme eu -v

# Mixed per-chain schemes
dvbfixer renumber bispecific.pdb --scheme seqres \
    --chain-scheme H:kabat --chain-scheme L:kabat -v

# IMGT (V 1-128); EU C-domain auto-shifted
dvbfixer renumber antibody.pdb --scheme imgt -v

# Fc fragment (CH2+CH3 only) — partial chain handled automatically
dvbfixer renumber fc.pdb --scheme eu -v   # starts CH2 at 231

# Verbose output showing alignment details and gaps
dvbfixer renumber input.pdb -v

# Custom output
dvbfixer renumber input.pdb -o renumbered.pdb
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_renum.pdb` | Output file path |
| `--scheme` | `seqres` | Numbering scheme: `seqres`, `kabat`, `chothia`, `imgt`, `martin`, `aho`, `eu`. Antibody schemes use ANARCI for V-domains + bundled EU references for C-domains. |
| `--chain-scheme` | none | Per-chain override (e.g. `H:kabat`). Repeatable. Wins over `--scheme`. |
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
| `--rename` | off | Rename non-canonical residues (AMBER/CHARMM) to standard names before processing |
| `-v`, `--verbose` | off | Print alignment details and gap positions |

## How It Works

### Default SEQRES path

1. Parses SEQRES records to get the full sequence per chain
2. Extracts unique (resSeq, iCode, resname) tuples from ATOM records
3. Aligns ATOM residues to SEQRES as a subsequence — each ATOM residue name is matched to the next occurrence in SEQRES
4. Assigns new sequential numbering based on SEQRES position (position 1 = first SEQRES residue)
5. Non-SEQRES residues (waters, ligands) are numbered sequentially after the last SEQRES position
6. Chains without SEQRES entries are renumbered sequentially from 1
7. All PDB sections are updated with the new numbering

### Antibody path (`--scheme` ≠ `seqres`)

1. For each chain, build the 1-letter AA sequence from ATOM records.
2. Run ANARCI on the sequence with the chosen scheme (Kabat / Chothia / IMGT / Martin / Aho; EU uses Kabat). If a V-domain is found, place its residues at the ANARCI-assigned (resseq, iCode) positions — preserving CDR insertion codes.
3. Take the residues NOT placed by ANARCI (typically the post-V-domain residues) and align them against three bundled human reference sequences (IgG1 heavy constant, Cκ, Cλ) via semi-global Needleman-Wunsch. The best-scoring reference wins; placed residues get their EU positions.
4. If V and C numbering would collide (IMGT-V ends at 128, EU CH1 starts at 118), shift the EU numbering forward by `(max_V_resseq + 5 - first_C_EU)` and emit a warning.
5. Chains where neither ANARCI nor the C-domain alignment placed anything fall back to the SEQRES path.

## See also

- [`split`](split.md) — apply before `renumber` when the input has no chain IDs
- [`model`](model.md) — rebuild missing loops (consumes renumbered output)
- [`rename`](rename.md) — text-only rename of non-canonical residues
- [`homology`](homology.md) — antibody-aware homology modeling (also uses ANARCI)
