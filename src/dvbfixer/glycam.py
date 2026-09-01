"""Convert PDB glycan nomenclature to GLYCAM force field naming.

GLYCAM uses 3-character residue codes encoding [linkage][sugar][anomer]:
  - Linkage: 0=terminal, 2-9=single position, V/W/U/Z/X/Y=multi-position
  - Sugar: G=Glc, L=Gal, M=Man, Y=GlcNAc, V=GalNAc, f=Fuc, S=Neu5Ac, ...
  - Anomer: A=alpha, B=beta (lowercase sugar code = L-sugar)

Detects glycosidic bonds from CONECT records (or distance-based fallback),
determines linkage patterns, and renames residues and atoms accordingly.
Optionally adds ROH cap at the reducing end.
"""

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

# PDB residue name -> (GLYCAM sugar code, GLYCAM anomer code)
# Lowercase sugar code = L-sugar
PDB_TO_GLYCAM = {
    'BGC': ('G', 'B'),   # beta-D-glucopyranose
    'GLC': ('G', 'A'),   # alpha-D-glucopyranose
    'GAL': ('L', 'B'),   # beta-D-galactopyranose
    'BGA': ('L', 'A'),   # alpha-D-galactopyranose
    'MAN': ('M', 'A'),   # alpha-D-mannopyranose
    'BMA': ('M', 'B'),   # beta-D-mannopyranose
    'AMA': ('M', 'A'),   # alpha-D-mannopyranose (alt PDB code)
    'NAG': ('Y', 'B'),   # N-acetyl-beta-D-glucosamine
    'NDG': ('Y', 'A'),   # N-acetyl-alpha-D-glucosamine
    'BGL': ('Y', 'B'),   # beta-GlcNAc (alt PDB code)
    'NGA': ('V', 'B'),   # N-acetyl-beta-D-galactosamine
    'A2G': ('V', 'A'),   # N-acetyl-alpha-D-galactosamine
    'FUC': ('f', 'A'),   # alpha-L-fucopyranose
    'FUL': ('f', 'B'),   # beta-L-fucopyranose
    'AFU': ('f', 'A'),   # alpha-L-fucopyranose (alt PDB code)
    'SIA': ('S', 'A'),   # alpha-Neu5Ac (sialic acid)
    'XYS': ('X', 'B'),   # beta-D-xylopyranose
    'XYP': ('X', 'A'),   # alpha-D-xylopyranose
    'RIB': ('R', 'B'),   # beta-D-ribose
    'GCU': ('Z', 'B'),   # beta-D-glucuronic acid
    'IDS': ('U', 'A'),   # alpha-L-iduronic acid
    'RAM': ('h', 'A'),   # alpha-L-rhamnose
}

# Multi-linkage position sets -> GLYCAM linkage code
MULTI_LINKAGE = {
    frozenset({2, 3}): 'Z',
    frozenset({2, 4}): 'Y',
    frozenset({2, 6}): 'X',
    frozenset({3, 4}): 'W',
    frozenset({3, 6}): 'V',
    frozenset({4, 6}): 'U',
    frozenset({2, 3, 4}): 'T',
    frozenset({2, 3, 6}): 'S',
    frozenset({2, 4, 6}): 'R',
    frozenset({3, 4, 6}): 'Q',
    frozenset({2, 3, 4, 6}): 'P',
}

# PDB atom name -> GLYCAM atom name (per PDB residue type)
# Universal hydroxyl H rename: PDB HOx → GLYCAM HxO (applies to ALL sugars)
_HYDROXYL_H_RENAME = {
    'HO1': 'H1O', 'HO2': 'H2O', 'HO3': 'H3O', 'HO4': 'H4O',
    'HO6': 'H6O', 'HO7': 'H7O', 'HO8': 'H8O', 'HO9': 'H9O',
}

# N-acetyl group rename: PDB standard → GLYCAM
_NACETYL_RENAME_PDB = {
    'C7': 'C2N', 'O7': 'O2N', 'C8': 'CME',           # PDB standard (NAG/NGA)
}
# N-acetyl rename for CHARMM-GUI style names
_NACETYL_RENAME_CHARMM = {
    'N': 'N2', 'HN': 'H2N',                            # amide N + H
    'C': 'C2N', 'O': 'O2N',                            # carbonyl C + O
    'CT': 'CME', 'HT1': 'H1M', 'HT2': 'H2M', 'HT3': 'H3M',  # methyl
}

# N-acetyl sugar PDB names (these get both PDB and CHARMM N-acetyl renames)
_NACETYL_SUGARS = {'NAG', 'NDG', 'BGL', 'NGA', 'A2G'}

# Per-residue specific renames (overrides). Applied AFTER the universal
# hydroxyl H rename and BEFORE the N-acetyl maps.
GLYCAM_ATOM_MAP = {
    # Sialic acid (Neu5Ac → 0SA / *SA in GLYCAM). The GLYCAM template
    # uses stereo-specific naming for the two methylene H atoms (H3A/H3E
    # on C3, H9R/H9S on C9), the amide H (H5N on N5), and the N-acetyl
    # methyl H atoms (H1M/H2M/H3M on CME, plus C10→C5N / C11→CME / O10→O5N).
    # We map both PDB-style (H3, H32) and CHARMM-style (H31, H32) input
    # H-names; the actual stereo assignment is arbitrary because OpenMM
    # only matches by name, then minimization refines positions.
    'SIA': {
        'C10': 'C5N', 'C11': 'CME', 'O10': 'O5N',
        # Methylene H on C3 (PDB: H3+H32 or H31+H32)
        'H3':  'H3A', 'H31': 'H3A', 'H32': 'H3E',
        # Methylene H on C9 (PDB: H9+H92 or H91+H92)
        'H9':  'H9R', 'H91': 'H9R', 'H92': 'H9S',
        # Amide H on N5 (HN5 → H5N; H5 on C5 stays as-is)
        'HN5': 'H5N',
        # N-acetyl methyl H atoms on CME (PDB H11/H112/H113)
        'H11':  'H1M', 'H112': 'H2M', 'H113': 'H3M',
    },
}

# Atom names to DROP from specific residues (e.g. PDBFixer-added H on a
# carboxylate). Applied during the rename pass.
GLYCAM_ATOM_DROP = {
    # Carboxylate -COO- has no H; PDBFixer sometimes adds one on O1B.
    'SIA': {'HO1B', 'HO1A'},
}

# Atom names to DROP from AMBER protonation-variant protein residues so the
# residue's atom set matches the AMBER template. Applied during glycam's
# rename pass — common when a user manually renamed e.g. LYS → LYN to mark
# deprotonation but didn't strip the HZ1 atom that LYN doesn't carry.
#
# Atom names verified from the AMBER ff14SB / ff19SB XML templates in
# OpenMM's data dir:
#   LYS NZ: HZ1, HZ2, HZ3.
#   LYN NZ: HZ2, HZ3 (HZ1 is the one stripped by deprotonation).
#
# Note: OpenMM's `hydrogens.xml` uses an inconsistent convention — it gates
# HZ3 by `variant="LYS"` (so addHydrogens with variant=LYN produces HZ1+HZ2,
# the OPPOSITE of the AMBER template). prepare and minimize patch this up
# with a post-addHydrogens HZ1→HZ3 rename. For glycam (text-only PDB
# rename) the AMBER template's atom set is the ground truth, so:
#   LYN: drop HZ1   (AMBER LYN keeps HZ2 + HZ3)
#   CYX/CYM: drop HG   (no SG-H)
#   HID: drop HE2   (HD1-only tautomer)
#   HIE: drop HD1   (HE2-only tautomer)
#   HIP / ASH / GLH have no drops.
#   NLN/OLS/OLT — handled by glycam's protein-link path (drops HD22/HG/HG1).
PROTEIN_VARIANT_ATOM_DROP = {
    'LYN': {'HZ1', '1HZ'},
    'CYX': {'HG', '1HG', 'HG1'},
    'CYM': {'HG', '1HG', 'HG1'},
    'HID': {'HE2'},
    'HIE': {'HD1'},
}

# Protein residues that can be glycosylated
PROTEIN_TO_GLYCAM = {
    'ASN': ('NLN', 'ND2'),   # N-linked glycosylation
    'SER': ('OLS', 'OG'),    # O-linked
    'THR': ('OLT', 'OG1'),   # O-linked
}

# Sialic acid residues (anomeric carbon is C2, not C1)
SIALIC_ACID_RESIDUES = {'SIA'}


# ---------------------------------------------------------------------------
# Reverse direction (GLYCAM → CHARMM)
# ---------------------------------------------------------------------------

# GLYCAM sugar letter → (alpha-anomer PDB name, beta-anomer PDB name).
# Output uses standard 3-char PDB codes (NAG/BMA/MAN/GAL/...). These are
# recognized by both CHARMM-GUI (as input) and `dvbfixer top --ff charmm`
# (via the PDB_TO_CARB mapping which translates PDB names to CHARMM RTP
# names like BGLCNA/BMAN/AMAN internally). The linkage character of the
# GLYCAM 3-char code is dropped — linkage info is preserved via CONECT
# records, exactly as in CHARMM-GUI PDB output.
_SUGAR_LETTER_TO_CHARMM = {
    'G': ('GLC', 'BGC'),   # glucose: alpha → GLC, beta → BGC
    'L': ('BGA', 'GAL'),   # galactose: alpha → BGA, beta → GAL
    'M': ('MAN', 'BMA'),   # mannose: alpha → MAN, beta → BMA
    'Y': ('NDG', 'NAG'),   # GlcNAc: alpha → NDG, beta → NAG
    'V': ('A2G', 'NGA'),   # GalNAc: alpha → A2G, beta → NGA
    'f': ('FUC', 'FUL'),   # L-fucose: alpha → FUC, beta → FUL
    'S': ('SIA', 'SIA'),   # Neu5Ac: only one PDB code (sialic always α)
    'X': ('XYP', 'XYS'),   # xylose: alpha → XYP, beta → XYS
    'R': ('RIB', 'RIB'),   # ribose (PDB has only one code)
    'Z': ('GCU', 'GCU'),   # glucuronic acid
    'U': ('IDS', 'IDS'),   # iduronic acid
    'h': ('RAM', 'RAM'),   # L-rhamnose
}

# GLYCAM glycoprotein residue → standard protein residue name.
GLYCAM_TO_STANDARD_PROTEIN = {'NLN': 'ASN', 'OLS': 'SER', 'OLT': 'THR'}

# AMBER protonation-variant residue name → CHARMM equivalent. Used by
# `glycam --to-charmm` so the output is directly consumable by
# `dvbfixer top --ff charmm` and downstream CHARMM-aware minimization.
# Source: AMBER ff14SB/ff19SB XML + CHARMM36 aminoacids.rtp (verified).
PROTONATION_AMBER_TO_CHARMM = {
    'HID': 'HSD', 'HIE': 'HSE', 'HIP': 'HSP',
    'ASH': 'ASPP', 'GLH': 'GLUP',
    'LYN': 'LSN',
    'CYX': 'CYS',   # CHARMM uses CYS + DISU patch, applied via SSBOND
    # CYM stays as CYM — CHARMM36 has a [ CYM ] residue.
}

# Reverse direction (default `glycam` forward path).
PROTONATION_CHARMM_TO_AMBER = {
    'HSD': 'HID', 'HSE': 'HIE', 'HSP': 'HIP',
    'ASPP': 'ASH', 'GLUP': 'GLH',
    'LSN': 'LYN',
}

# Per-residue atom renames at the glycam stage. Methylene H shifts
# (HB2/HB3 ↔ HB1/HB2 etc.) and backbone H ↔ HN are NOT in here — top.py
# already handles them during topology generation (methylene_shift +
# aminoacids.arn). The only AMBER↔CHARMM rename that's asymmetric on the
# atom set (so top.py can't infer it) is LYN/LSN's NH2-H pair:
#
#   AMBER LYN: HZ2 + HZ3 (HZ1 absent)
#   CHARMM LSN: HZ1 + HZ2 (HZ3 absent)
#
# So LYN → LSN needs HZ2→HZ1, HZ3→HZ2 (and the reverse for LSN → LYN).
# Apply in a single atomic pass to avoid the HZ2→HZ1 then HZ3→HZ2 collision.
PROTONATION_ATOM_RENAME_TO_CHARMM = {
    'LYN': {'HZ2': 'HZ1', 'HZ3': 'HZ2'},
}
PROTONATION_ATOM_RENAME_TO_AMBER = {
    'LSN': {'HZ1': 'HZ2', 'HZ2': 'HZ3'},
}

# PDB 3-char names that are N-acetyl sugars (have an N-acetyl group).
_NACETYL_PDB_NAMES = {'NAG', 'NDG', 'NGA', 'A2G'}

# Inverse of _HYDROXYL_H_RENAME (universal hydroxyl H rename, all sugars).
# GLYCAM H<n>O → PDB HO<n>
_REV_HYDROXYL_H = {v: k for k, v in _HYDROXYL_H_RENAME.items()}

# Inverse of _NACETYL_RENAME_PDB (N-acetyl group: GLYCAM → standard PDB)
# for N-acetyl sugars (NAG/NDG/NGA/A2G).
# C2N→C7, O2N→O7, CME→C8
_REV_NACETYL_PDB = {v: k for k, v in _NACETYL_RENAME_PDB.items()}
# Additional atom renames for standard PDB NAG-like residues:
_REV_NACETYL_PDB.update({
    'N2': 'N2',   # amide N — standard PDB keeps N2
    'H2N': 'HN2', # amide H — PDB uses HN2 (vs CHARMM HN)
    'H1M': 'H81', 'H2M': 'H82', 'H3M': 'H83',  # methyl H's: PDB H81/H82/H83
})

# Per-residue reverse atom rename. SIA has the N-acetyl group on C5,
# not the standard C7 position, so it has a different layout.
_REV_GLYCAM_ATOM_MAP = {
    'SIA': {
        # N-acetyl group on C5 (PDB SIA uses C10/CT11/O10 — but CHARMM-GUI
        # and most PDB readers actually accept C5N/CME/O5N for sialic. We
        # output standard PDB SIA naming: C10/C11/O10).
        'C5N': 'C10', 'CME': 'C11', 'O5N': 'O10',
        'H1M': 'H111', 'H2M': 'H112', 'H3M': 'H113',
        # Amide H on N5 (PDB SIA uses HN5)
        'H5N': 'HN5',
        # C3 methylene (PDB SIA: H31/H32)
        'H3A': 'H31', 'H3E': 'H32',
        # C9 methylene (PDB SIA: H91/H92)
        'H9R': 'H91', 'H9S': 'H92',
    },
}


def _glycam_to_charmm_resname(glycam_name):
    """Translate a GLYCAM 3-char sugar code to its CHARMM 4-char equivalent.
    Returns the CHARMM name, or None if not a recognized GLYCAM sugar code.
    """
    if len(glycam_name) != 3:
        return None
    sugar_letter = glycam_name[1]
    anomer = glycam_name[2]
    pair = _SUGAR_LETTER_TO_CHARMM.get(sugar_letter)
    if pair is None:
        return None
    if anomer == 'A':
        return pair[0]
    if anomer == 'B':
        return pair[1]
    return None


def _rename_atom_reverse(glycam_resname, charmm_resname, atom_name):
    """Rename a GLYCAM atom name back to standard PDB / CHARMM convention.

    Applies in order:
    1. Residue-specific reverse map (e.g. SIA's C5N→C10, H3A→H31)
    2. N-acetyl group reverse map — only for N-acetyl sugars
       (NAG/NDG/NGA/A2G in PDB; SIA uses its own map in step 1).
    3. Universal hydroxyl H rename reverse (H<n>O → HO<n>).
    """
    # 1. Residue-specific (sialic acid: keyed by destination PDB name SIA)
    if charmm_resname == 'SIA':
        specific = _REV_GLYCAM_ATOM_MAP.get('SIA', {})
        if atom_name in specific:
            return specific[atom_name]

    # 2. N-acetyl reverse — for NAG/NDG/NGA/A2G
    if charmm_resname in _NACETYL_PDB_NAMES:
        if atom_name in _REV_NACETYL_PDB:
            return _REV_NACETYL_PDB[atom_name]

    # 3. Universal hydroxyl H
    if atom_name in _REV_HYDROXYL_H:
        return _REV_HYDROXYL_H[atom_name]

    return atom_name


_PASSTHROUGH_EXCLUDED_RECORDS = (
    'ATOM  ', 'HETATM', 'CONECT', 'TER', 'END', 'MODEL', 'ENDMDL',
)


def _extract_passthrough_header_lines(input_path):
    """Return every line of `input_path` that isn't ATOM/HETATM/CONECT/
    TER/END/MODEL/ENDMDL — SEQRES, HELIX, SHEET, CRYST1, HEADER, TITLE,
    LINK, SSBOND, REMARK, etc. — preserving original order.

    Mirrors `align.py`'s `_apply_transform_preserving_headers`: PDB
    metadata that isn't atom-indexed or residue-name-dependent (SEQRES
    in particular — `model.py`'s gap-filling needs it) has no reason to
    be dropped just because atoms got renamed/renumbered. CONECT is
    excluded here because both `convert_to_glycam`/`convert_to_charmm`
    already regenerate it correctly with remapped serials — passing
    through the OLD one too would just duplicate/conflict. LINK is
    passed through as-is even though its residue names could go stale
    after a rename (e.g. NAG → a GLYCAM code) — CONECT (regenerated) is
    the bond source every downstream OpenMM-based tool actually reads;
    LINK is informational.
    """
    with open(input_path) as f:
        lines = f.readlines()
    return [line for line in lines
            if not line.startswith(_PASSTHROUGH_EXCLUDED_RECORDS)]


def convert_to_charmm(input_path, output_path, verbose=False):
    """Convert a GLYCAM-named PDB to CHARMM-compatible naming.

    Inverse of `convert_to_glycam`:
    - GLYCAM 3-char sugar codes (UYB/4YB/0SA/VMB/...) → CHARMM 4-char
      names (BGLCNA/ANE5AC/BMAN/...). Linkage info is lost from the
      residue name but preserved via CONECT records.
    - GLYCAM glycoprotein residues (NLN/OLS/OLT) → standard ASN/SER/THR.
    - Atom names reverted via the inverse of the forward rename maps.
    - ROH/OME caps stripped (CHARMM doesn't model reducing-end caps).

    Returns: number of sugar residues converted.
    """
    from dvbfixer.ffutils import is_glycam_sugar

    atoms, _residues, _bond_graph, _links = _parse_pdb(input_path)

    # Group atoms into residues, preserving order
    residues = {}  # (chain, resseq, icode) → list of atom indices
    res_order = []
    for i, atom in enumerate(atoms):
        key = (atom['chain'], atom['resseq'], atom['icode'])
        if key not in residues:
            residues[key] = []
            res_order.append(key)
        residues[key].append(i)

    # Classify each residue; collect rename decisions
    res_new_name = {}  # key → new resname (or None to skip residue)
    n_sugars = 0
    n_protein = 0
    n_caps_dropped = 0
    for key in res_order:
        first_atom = atoms[residues[key][0]]
        rn = first_atom['resname']
        if rn in {'ROH', 'OME', 'TBT', 'CMET'}:
            # GLYCAM reducing-end cap — drop in CHARMM output
            res_new_name[key] = None
            n_caps_dropped += 1
            if verbose:
                print(f"  Dropped {rn} cap at {key[0]}:{key[1]}{key[2].strip()} "
                      f"(CHARMM has no equivalent)")
            continue
        if rn in GLYCAM_TO_STANDARD_PROTEIN:
            res_new_name[key] = GLYCAM_TO_STANDARD_PROTEIN[rn]
            n_protein += 1
            if verbose:
                print(f"  {key[0]}:{rn}{key[1]} -> {res_new_name[key]}")
            continue
        if is_glycam_sugar(rn):
            charmm_name = _glycam_to_charmm_resname(rn)
            if charmm_name is None:
                if verbose:
                    print(f"  WARNING: GLYCAM sugar {rn} at {key[0]}:{key[1]} "
                          f"has no CHARMM equivalent; keeping original name")
                res_new_name[key] = rn
            else:
                res_new_name[key] = charmm_name
                n_sugars += 1
                if verbose:
                    print(f"  {key[0]}:{rn}{key[1]} -> {charmm_name}")
            continue
        # AMBER protonation variant → CHARMM equivalent
        if rn in PROTONATION_AMBER_TO_CHARMM:
            res_new_name[key] = PROTONATION_AMBER_TO_CHARMM[rn]
            if verbose:
                print(f"  {key[0]}:{rn}{key[1]} -> {res_new_name[key]} "
                      f"(AMBER→CHARMM protonation variant)")
            continue
        # Standard amino acid or unknown — keep as is
        res_new_name[key] = rn

    # Apply atom renames within sugar residues
    for key in res_order:
        new_rn = res_new_name.get(key)
        if new_rn is None:
            continue
        orig_rn = atoms[residues[key][0]]['resname']
        if not is_glycam_sugar(orig_rn):
            continue
        for idx in residues[key]:
            atoms[idx]['name'] = _rename_atom_reverse(
                orig_rn, new_rn, atoms[idx]['name'])

    # Apply protonation-variant atom renames (LYN→LSN swaps HZ atoms).
    # Single atomic pass per residue to avoid HZ2→HZ1 then HZ3→HZ2 collision.
    for key in res_order:
        orig_rn = atoms[residues[key][0]]['resname']
        atom_map = PROTONATION_ATOM_RENAME_TO_CHARMM.get(orig_rn)
        if not atom_map:
            continue
        for idx in residues[key]:
            old_name = atoms[idx]['name']
            if old_name in atom_map:
                atoms[idx]['name'] = atom_map[old_name]
                if verbose:
                    print(f"    Renamed {key[0]}:{orig_rn}{key[1]} "
                          f"{old_name} → {atom_map[old_name]}")

    # Write output: emit atom lines for residues we keep, remap serials,
    # then append CONECT records with remapped serials.
    serial = 0
    old_to_new_serial = {}
    out_lines = []
    for key in res_order:
        new_rn = res_new_name.get(key)
        if new_rn is None:
            continue
        for idx in residues[key]:
            atom = atoms[idx]
            serial += 1
            old_to_new_serial[atom['serial']] = serial
            line = _format_atom_line(atom, new_rn, serial)
            out_lines.append(line)

    # CONECT records — remap serials, drop bonds involving dropped atoms.
    with open(input_path) as f:
        for line in f:
            if not line.startswith('CONECT'):
                continue
            parts = []
            s = line[6:]
            while len(s) >= 5:
                chunk = s[:5].strip()
                if chunk:
                    try:
                        parts.append(int(chunk))
                    except ValueError:
                        pass
                s = s[5:]
            if len(parts) < 2:
                continue
            new_parts = [old_to_new_serial.get(p) for p in parts]
            if any(p is None for p in new_parts):
                continue  # bond involves a dropped atom
            line_out = 'CONECT' + ''.join(f'{p:5d}' for p in new_parts) + '\n'
            out_lines.append(line_out)
    out_lines.append('END\n')

    header_lines = _extract_passthrough_header_lines(input_path)
    with open(output_path, 'w') as f:
        f.writelines(header_lines + out_lines)

    print(f"Converted {n_sugars} GLYCAM sugar(s) and {n_protein} glycoprotein "
          f"residue(s) → CHARMM naming. Dropped {n_caps_dropped} cap(s).")
    return n_sugars


def _parse_pdb(path):
    """Parse PDB file into atoms, residues, and bond graph."""
    with open(path) as f:
        lines = f.readlines()

    atoms = []          # list of atom dicts (preserving order)
    serial_to_idx = {}  # serial -> index in atoms list
    bond_graph = defaultdict(set)  # serial -> set of bonded serials
    link_records = []   # parsed LINK records

    # CHARMM-GUI uses 4-character residue names (ASPP, GLUP, BGLC, ANE5AC, ...).
    # In PDB layout the 4th character extends into col 21 — the standard
    # chain-ID column. Detect this case via the trailing letter at col 20
    # (0-indexed): when col 17:21 spells a known 4-char resname, the chain
    # ID is at col 22 (shifted by 1). Otherwise the standard col-22 chain
    # applies. Known 4-char names accepted here:
    _CHARMM_4CHAR_RESNAMES = {
        'ASPP', 'GLUP',                           # protonated D / E
        'BGLC', 'AGLC', 'BMAN', 'AMAN', 'BGAL', 'AGAL', 'BFUC', 'AFUC',
        'BGLCNA', 'AGLCNA',                       # GlcNAc α/β
        'BGALNA', 'AGALNA',                       # GalNAc α/β
        'ANE5', 'ANE5AC', 'BNE5', 'BNE5AC',       # sialic acid variants
        'AIDO', 'BIDO',
        'CER1', 'CER160', 'CER180', 'CER181',
        'CER2', 'CER200', 'CER220', 'CER240', 'CER241', 'CER3E',
    }

    def _read_resname_chain(line):
        """Return (resname, chain) handling both 3-char and CHARMM 4-char layouts."""
        candidate4 = line[17:21].strip()
        if candidate4 in _CHARMM_4CHAR_RESNAMES:
            return candidate4, line[21]
        return line[17:20].strip(), line[21] if len(line) > 21 else ' '

    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            serial = int(line[6:11])
            name = line[12:16].strip()
            resname, chain = _read_resname_chain(line)
            resseq = int(line[22:26])
            icode = line[26] if len(line) > 26 else ' '
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) > 76 else name[0]

            atom = {
                'serial': serial, 'name': name, 'resname': resname,
                'chain': chain, 'resseq': resseq, 'icode': icode,
                'x': x, 'y': y, 'z': z, 'element': element,
                'line': line,
            }
            serial_to_idx[serial] = len(atoms)
            atoms.append(atom)

        elif line.startswith('CONECT'):
            serials = []
            s = line[6:]
            while len(s) >= 5:
                chunk = s[:5].strip()
                if chunk:
                    serials.append(int(chunk))
                s = s[5:]
            if len(serials) >= 2:
                src = serials[0]
                for dst in serials[1:]:
                    bond_graph[src].add(dst)
                    bond_graph[dst].add(src)

        elif line.startswith('LINK  '):
            # LINK records: atom1 (13-26) -- atom2 (43-56)
            try:
                name1 = line[12:16].strip()
                resname1 = line[17:20].strip()
                chain1 = line[21]
                resseq1 = int(line[22:26])
                name2 = line[42:46].strip()
                resname2 = line[47:50].strip()
                chain2 = line[51]
                resseq2 = int(line[52:56])
                link_records.append({
                    'name1': name1, 'resname1': resname1,
                    'chain1': chain1, 'resseq1': resseq1,
                    'name2': name2, 'resname2': resname2,
                    'chain2': chain2, 'resseq2': resseq2,
                })
            except (ValueError, IndexError):
                pass

    # Build residue dict: (chain, resseq) -> list of atom indices
    residues = defaultdict(list)
    for i, atom in enumerate(atoms):
        residues[(atom['chain'], atom['resseq'])].append(i)

    return atoms, residues, bond_graph, link_records


def _detect_glycosidic_bonds(atoms, residues, bond_graph):
    """Detect glycosidic bonds from CONECT records.

    Returns:
        glyco_bonds: list of (parent_reskey, child_reskey, linkage_position)
        protein_links: list of (protein_reskey, sugar_reskey)
    """
    glyco_bonds = []
    protein_links = []

    # Build serial -> atom lookup
    serial_to_atom = {a['serial']: a for a in atoms}

    seen_bonds = set()
    for src_serial, neighbors in bond_graph.items():
        if src_serial not in serial_to_atom:
            continue
        src = serial_to_atom[src_serial]
        src_reskey = (src['chain'], src['resseq'])

        for dst_serial in neighbors:
            if dst_serial not in serial_to_atom:
                continue
            dst = serial_to_atom[dst_serial]
            dst_reskey = (dst['chain'], dst['resseq'])

            if src_reskey == dst_reskey:
                continue  # intra-residue

            bond_key = tuple(sorted([src_serial, dst_serial]))
            if bond_key in seen_bonds:
                continue
            seen_bonds.add(bond_key)

            # Check for protein-sugar bond
            if src['resname'] in PROTEIN_TO_GLYCAM and dst['resname'] in PDB_TO_GLYCAM:
                glycam_name, link_atom = PROTEIN_TO_GLYCAM[src['resname']]
                if src['name'] == link_atom:
                    protein_links.append((src_reskey, dst_reskey))
                    continue
            if dst['resname'] in PROTEIN_TO_GLYCAM and src['resname'] in PDB_TO_GLYCAM:
                glycam_name, link_atom = PROTEIN_TO_GLYCAM[dst['resname']]
                if dst['name'] == link_atom:
                    protein_links.append((dst_reskey, src_reskey))
                    continue

            # Check for sugar-sugar glycosidic bond: C1-O (or C2-O for sialic acid)
            anomeric = None
            oxygen = None
            if _is_anomeric_carbon(src) and dst['name'].startswith('O'):
                anomeric, oxygen = src, dst
            elif _is_anomeric_carbon(dst) and src['name'].startswith('O'):
                anomeric, oxygen = dst, src

            if anomeric and oxygen:
                child_reskey = (anomeric['chain'], anomeric['resseq'])
                parent_reskey = (oxygen['chain'], oxygen['resseq'])
                # Extract linkage position from oxygen name (O4 -> 4)
                o_name = oxygen['name']
                try:
                    pos = int(o_name[1:])
                    glyco_bonds.append((parent_reskey, child_reskey, pos))
                except ValueError:
                    pass

    return glyco_bonds, protein_links


def _detect_glycosidic_bonds_by_distance(atoms, residues):
    """Fallback: detect glycosidic bonds by C1-O distance < 2.0 A."""
    glyco_bonds = []
    protein_links = []

    # Collect anomeric carbons and oxygens per residue
    anomeric_atoms = {}  # reskey -> atom dict for C1 (or C2 for SIA)
    oxygen_atoms = defaultdict(list)  # reskey -> list of O atoms

    for atom in atoms:
        reskey = (atom['chain'], atom['resseq'])
        if _is_anomeric_carbon(atom):
            anomeric_atoms[reskey] = atom
        if atom['name'].startswith('O') and atom['resname'] in PDB_TO_GLYCAM:
            oxygen_atoms[reskey].append(atom)

    # Check each anomeric C against all O atoms on other residues
    for child_key, c_atom in anomeric_atoms.items():
        for parent_key, o_list in oxygen_atoms.items():
            if child_key == parent_key:
                continue
            for o_atom in o_list:
                dx = c_atom['x'] - o_atom['x']
                dy = c_atom['y'] - o_atom['y']
                dz = c_atom['z'] - o_atom['z']
                dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                if dist < 2.0:
                    try:
                        pos = int(o_atom['name'][1:])
                        glyco_bonds.append((parent_key, child_key, pos))
                    except ValueError:
                        pass

    # Check protein-sugar links
    for atom in atoms:
        if atom['resname'] not in PROTEIN_TO_GLYCAM:
            continue
        _, link_atom = PROTEIN_TO_GLYCAM[atom['resname']]
        if atom['name'] != link_atom:
            continue
        prot_key = (atom['chain'], atom['resseq'])
        for sugar_key, c_atom in anomeric_atoms.items():
            dx = atom['x'] - c_atom['x']
            dy = atom['y'] - c_atom['y']
            dz = atom['z'] - c_atom['z']
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < 2.0:
                protein_links.append((prot_key, sugar_key))

    return glyco_bonds, protein_links


def _merge_glycosidic_bonds(conect_result, distance_result):
    """Merge CONECT-derived and distance-derived glycosidic bond
    detection.

    CONECT-derived bonds win where both agree, but CONECT/LINK
    annotation in real deposited PDBs is sometimes incomplete for a
    subset of genuine glycosylation sites — an annotation gap, not
    evidence the bond doesn't exist. Any bond the distance-based
    detector finds for a residue not already covered by a CONECT
    result is added rather than discarded; a residue with ANY CONECT
    bond record at all is otherwise treated as fully trustworthy and
    not second-guessed by geometry.
    """
    conect_glyco, conect_links = conect_result
    dist_glyco, dist_links = distance_result

    covered_children = {child for _parent, child, _pos in conect_glyco}
    merged_glyco = list(conect_glyco)
    for parent, child, pos in dist_glyco:
        if child not in covered_children:
            merged_glyco.append((parent, child, pos))
            covered_children.add(child)

    covered_sugars = {sugar for _prot, sugar in conect_links}
    merged_links = list(conect_links)
    for prot, sugar in dist_links:
        if sugar not in covered_sugars:
            merged_links.append((prot, sugar))
            covered_sugars.add(sugar)

    return merged_glyco, merged_links


def _is_anomeric_carbon(atom):
    """Check if atom is an anomeric carbon (C1, or C2 for sialic acid)."""
    if atom['resname'] in SIALIC_ACID_RESIDUES:
        return atom['name'] == 'C2'
    return atom['name'] == 'C1'


def _determine_linkage_code(positions):
    """Convert set of child linkage positions to GLYCAM linkage code."""
    if not positions:
        return '0'
    if len(positions) == 1:
        return str(next(iter(positions)))
    key = frozenset(positions)
    if key in MULTI_LINKAGE:
        return MULTI_LINKAGE[key]
    # Unknown multi-linkage — use lowest position as fallback
    print(f"  WARNING: Unknown multi-linkage combination {sorted(positions)}, "
          f"using '{min(positions)}'", file=sys.stderr)
    return str(min(positions))


def _rename_atom(pdb_resname, atom_name):
    """Rename a PDB atom name to GLYCAM atom name.

    Applies in order:
    1. Residue-specific overrides (GLYCAM_ATOM_MAP)
    2. N-acetyl renames (PDB standard + CHARMM-GUI style) for NAG/NGA/BGL/NDG/A2G
    3. Universal hydroxyl H rename (HOx → HxO) for all sugars
    """
    # 1. Residue-specific
    specific = GLYCAM_ATOM_MAP.get(pdb_resname, {})
    if atom_name in specific:
        return specific[atom_name]

    # 2. N-acetyl renames
    if pdb_resname in _NACETYL_SUGARS:
        if atom_name in _NACETYL_RENAME_PDB:
            return _NACETYL_RENAME_PDB[atom_name]
        if atom_name in _NACETYL_RENAME_CHARMM:
            return _NACETYL_RENAME_CHARMM[atom_name]

    # 3. Universal hydroxyl H rename
    if atom_name in _HYDROXYL_H_RENAME:
        return _HYDROXYL_H_RENAME[atom_name]

    return atom_name


def _format_atom_line(atom, new_resname, new_serial, new_resseq=None):
    """Format an ATOM/HETATM PDB line with updated names."""
    name = atom['name']
    # PDB atom name formatting: 4-char names start at col 13, shorter at col 14
    if len(name) < 4:
        name_field = f" {name:<3s}"
    else:
        name_field = f"{name:<4s}"

    resseq = new_resseq if new_resseq is not None else atom['resseq']
    icode = atom.get('icode', ' ')
    if icode.strip() == '':
        icode = ' '

    element = atom.get('element', name[0])

    # Standard amino acids → ATOM, everything else → HETATM
    _STD_AA = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HIE', 'HID', 'HIP', 'ASH', 'GLH', 'CYX', 'CYM', 'LYN',  # AMBER variants
        'HSD', 'HSE', 'HSP', 'ASPP', 'GLUP', 'LSN',  # CHARMM variants
        'NLN', 'OLS', 'OLT',  # GLYCAM protein residues
    }
    record = 'ATOM  ' if new_resname in _STD_AA else 'HETATM'

    # CHARMM 4-character residue names (ASPP, GLUP, LSN [3 char actually],
    # ANE5AC, ...) extend into the col-20 gap between resname and chain ID,
    # leaving altLoc (col 17) still as a space. Layout per PDB convention:
    #   cols 13-16: atom name (`name_field`)
    #   col 17:     altLoc (always space here)
    #   cols 18-20: resname (3 chars), or 18-21 if resname is 4 chars
    #   col 21:     space (3-char path) or chain ID (4-char path, immediately
    #               after the 4-char resname)
    #   col 22:     chain ID (3-char path) -- shifted left by 1 for 4-char
    if len(new_resname) >= 4:
        resname_block = f" {new_resname[:4]:<4s}"       # altLoc + 4-char + no gap
    else:
        resname_block = f" {new_resname:>3s} "          # altLoc + 3-char + gap
    return (
        f"{record}{new_serial:5d} {name_field}"
        f"{resname_block}{atom['chain']}{resseq:4d}{icode}"
        f"   {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
        f"                      {element:>2s}  \n"
    )


def convert_to_glycam(input_path, output_path, add_roh=True, verbose=False):
    """Convert PDB glycan nomenclature to GLYCAM naming."""
    atoms, residues, bond_graph, link_records = _parse_pdb(input_path)

    # Detect glycosidic bonds. Always run the distance-based detector
    # too, even when CONECT records exist — CONECT/LINK annotation in
    # real deposited PDBs is sometimes incomplete for a subset of
    # genuine glycosylation sites (an annotation gap, not evidence the
    # bond doesn't exist). Previously trusting "any CONECT present" as
    # an all-or-nothing gate silently dropped those under-annotated
    # sites entirely — their sugar tree ended up floating, unbonded to
    # the protein, in the final output.
    dist_glyco, dist_links = _detect_glycosidic_bonds_by_distance(atoms, residues)
    if bond_graph:
        conect_glyco, conect_links = _detect_glycosidic_bonds(
            atoms, residues, bond_graph)
        glyco_bonds, protein_links = _merge_glycosidic_bonds(
            (conect_glyco, conect_links), (dist_glyco, dist_links))
        if verbose:
            n_extra_glyco = len(glyco_bonds) - len(conect_glyco)
            n_extra_links = len(protein_links) - len(conect_links)
            extra = (f" (+{n_extra_glyco} more by distance, filling CONECT gaps)"
                     if n_extra_glyco > 0 else "")
            print(f"  Detected {len(conect_glyco)} glycosidic bonds from CONECT{extra}")
            if n_extra_links > 0:
                print(f"  +{n_extra_links} protein-sugar link(s) found by "
                      f"distance, not present in CONECT/LINK")
    else:
        glyco_bonds, protein_links = dist_glyco, dist_links
        if verbose:
            print(f"  Detected {len(glyco_bonds)} glycosidic bonds by distance")

    if verbose and protein_links:
        print(f"  Detected {len(protein_links)} protein-sugar links")

    # Build parent->children map and child->parent map
    children_positions = defaultdict(set)  # reskey -> set of linkage positions
    child_to_parent = {}                   # child_reskey -> parent_reskey
    for parent_key, child_key, pos in glyco_bonds:
        children_positions[parent_key].add(pos)
        child_to_parent[child_key] = parent_key

    # Identify sugar residues and protein-linked sugars
    sugar_reskeys = set()
    for reskey, atom_indices in residues.items():
        resname = atoms[atom_indices[0]]['resname']
        if resname in PDB_TO_GLYCAM:
            sugar_reskeys.add(reskey)

    protein_linked_sugars = {sugar_key for _, sugar_key in protein_links}

    # Find reducing-end sugars (roots): sugars with no parent sugar and not protein-linked
    root_sugars = set()
    for reskey in sugar_reskeys:
        if reskey not in child_to_parent and reskey not in protein_linked_sugars:
            root_sugars.add(reskey)

    # Build GLYCAM residue names
    glycam_names = {}  # reskey -> GLYCAM 3-char name
    for reskey in sugar_reskeys:
        atom_idx = residues[reskey][0]
        pdb_resname = atoms[atom_idx]['resname']

        if pdb_resname not in PDB_TO_GLYCAM:
            if verbose:
                print(f"  WARNING: Unknown sugar {pdb_resname} at "
                      f"{reskey[0]}:{reskey[1]}, skipping")
            continue

        sugar_code, anomer_code = PDB_TO_GLYCAM[pdb_resname]
        linkage_code = _determine_linkage_code(children_positions.get(reskey, set()))
        glycam_name = linkage_code + sugar_code + anomer_code
        glycam_names[reskey] = glycam_name

        if verbose:
            positions = children_positions.get(reskey, set())
            pos_str = ','.join(f'O{p}' for p in sorted(positions)) if positions else 'terminal'
            print(f"  {reskey[0]}:{pdb_resname}{reskey[1]} -> {glycam_name} ({pos_str})")

    # Rename atoms in sugar residues + drop atoms that GLYCAM templates
    # don't have (e.g. PDBFixer-added HO1B on sialic acid carboxylate).
    n_dropped = 0
    for reskey in sugar_reskeys:
        drop_set = GLYCAM_ATOM_DROP.get(atoms[residues[reskey][0]]['resname'], set())
        kept_idx = []
        for atom_idx in residues[reskey]:
            atom = atoms[atom_idx]
            if atom['name'] in drop_set:
                n_dropped += 1
                continue
            atom['name'] = _rename_atom(atom['resname'], atom['name'])
            kept_idx.append(atom_idx)
        residues[reskey] = kept_idx
    if verbose and n_dropped:
        print(f"  Dropped {n_dropped} atom(s) not in GLYCAM templates "
              f"(e.g. HO1B on sialic carboxylate)")

    # Handle protein-linked residues
    protein_renames = {}
    for prot_key, sugar_key in protein_links:
        prot_idx = residues[prot_key][0]
        prot_resname = atoms[prot_idx]['resname']
        if prot_resname in PROTEIN_TO_GLYCAM:
            glycam_prot_name, _ = PROTEIN_TO_GLYCAM[prot_resname]
            protein_renames[prot_key] = glycam_prot_name
            if verbose:
                print(f"  {prot_key[0]}:{prot_resname}{prot_key[1]} -> {glycam_prot_name}")

    # Write output PDB
    serial = 0
    out_lines = []
    roh_atoms = []  # (atom_dict, resseq) for ROH cap

    # Determine residue order from input
    reskey_order = []
    seen_reskeys = set()
    for atom in atoms:
        rk = (atom['chain'], atom['resseq'])
        if rk not in seen_reskeys:
            seen_reskeys.add(rk)
            reskey_order.append(rk)

    # Track old serial -> new serial for CONECT remapping
    serial_map = {}
    roh_serial_map = {}  # old serial -> new serial for ROH atoms

    for reskey in reskey_order:
        atom_indices = residues[reskey]
        first_atom = atoms[atom_indices[0]]
        pdb_resname = first_atom['resname']

        # Determine output residue name. CHARMM protonation variants in the
        # input (HSD/HSE/HSP/ASPP/GLUP/LSN) are renamed to their AMBER
        # equivalents (HID/HIE/HIP/ASH/GLH/LYN) so the output is consumable
        # by AMBER-aware downstream tools (prepare, minimize, top --acpype).
        if reskey in glycam_names:
            out_resname = glycam_names[reskey]
        elif reskey in protein_renames:
            out_resname = protein_renames[reskey]
        elif pdb_resname in PROTONATION_CHARMM_TO_AMBER:
            out_resname = PROTONATION_CHARMM_TO_AMBER[pdb_resname]
            if verbose:
                print(f"  {reskey[0]}:{pdb_resname}{reskey[1]} -> "
                      f"{out_resname} (CHARMM→AMBER protonation variant)")
        else:
            out_resname = pdb_resname

        # Apply protonation-variant atom renames (LSN → LYN swaps HZ atoms).
        # Build a per-atom-index rewrite dict from the source CHARMM resname's
        # entry in PROTONATION_ATOM_RENAME_TO_AMBER; apply atomically so
        # HZ1→HZ2 and HZ2→HZ3 don't collide.
        proton_atom_map = PROTONATION_ATOM_RENAME_TO_AMBER.get(pdb_resname, {})

        # Strip atoms that the output residue's AMBER template doesn't carry.
        # Applies to AMBER protonation variants (LYN/CYX/CYM/HID/HIE) in the
        # input — e.g. a user manually renamed LYS → LYN to mark deprotonation
        # but left HZ1 in place; LYN's AMBER template has only HZ2/HZ3.
        variant_drop_set = PROTEIN_VARIANT_ATOM_DROP.get(out_resname, set())

        for atom_idx in atom_indices:
            atom = atoms[atom_idx]

            # For reducing-end sugars: extract O1 for ROH cap
            if add_roh and reskey in root_sugars and atom['name'] == 'O1':
                roh_atoms.append((dict(atom), reskey[1]))
                # Don't write O1 in the sugar itself
                continue

            # Skip HD22 on NLN (glycosylated ASN)
            if reskey in protein_renames and protein_renames[reskey] == 'NLN':
                if atom['name'] == 'HD22':
                    continue

            # Apply protonation-variant atom rename (e.g. LSN HZ1 → LYN HZ2).
            if atom['name'] in proton_atom_map:
                old_name = atom['name']
                atom['name'] = proton_atom_map[old_name]
                if verbose:
                    print(f"    Renamed {reskey[0]}:{pdb_resname}{reskey[1]} "
                          f"{old_name} → {atom['name']}")

            # Drop atoms that don't belong in the AMBER variant template.
            if atom['name'] in variant_drop_set:
                if verbose:
                    print(f"  Dropped {reskey[0]}:{out_resname}{reskey[1]} "
                          f"{atom['name']} (not in AMBER {out_resname} template)")
                continue

            serial += 1
            serial_map[atom['serial']] = serial
            line = _format_atom_line(atom, out_resname, serial)
            out_lines.append(line)

    # Add ROH cap residues at end
    for roh_atom, orig_resseq in roh_atoms:
        serial += 1
        roh_serial_map[roh_atom['serial']] = serial
        serial_map[roh_atom['serial']] = serial
        roh_atom['name'] = 'O1'
        line = _format_atom_line(roh_atom, 'ROH', serial, new_resseq=orig_resseq)
        out_lines.append(line)
        if verbose:
            print(f"  Added ROH cap at {roh_atom['chain']}:{orig_resseq}")

    # Write TER
    out_lines.append('TER\n')

    # Remap and write CONECT records
    for src_serial, neighbors in sorted(bond_graph.items()):
        if src_serial not in serial_map:
            continue
        new_src = serial_map[src_serial]
        new_neighbors = []
        for dst in sorted(neighbors):
            if dst in serial_map:
                new_neighbors.append(serial_map[dst])
        if new_neighbors:
            conect = f"CONECT{new_src:5d}"
            for n in new_neighbors:
                conect += f"{n:5d}"
            out_lines.append(conect + '\n')

    out_lines.append('END\n')

    header_lines = _extract_passthrough_header_lines(input_path)
    with open(output_path, 'w') as f:
        f.writelines(header_lines + out_lines)

    return len(glycam_names), len(protein_renames)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer convert",
        description="Convert structure naming between PDB/AMBER/GLYCAM and "
        "CHARMM conventions. Default direction (`--to-amber`): convert PDB / "
        "CHARMM-named input to GLYCAM sugar codes + AMBER protonation variants "
        "(consumable by `prepare`, `minimize`, `protonate`, `top --acpype`). "
        "With `--to-charmm`: reverse direction, output CHARMM-compatible names "
        "(consumable by `top --ff charmm`). Both directions are idempotent — "
        "running the tool on already-correct input is a no-op.",
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB, PDBx/mmCIF, or crystallographic CIF file")
    io.add_argument("-o", "--output",
                    help="Output PDB file (default: <input>_amber.pdb in "
                         "the default direction, <input>_charmm.pdb with "
                         "--to-charmm)")

    direction_grp = p.add_argument_group("Conversion direction")
    direction = direction_grp.add_mutually_exclusive_group()
    direction.add_argument("--to-amber", action="store_true",
                           help="PDB/CHARMM → GLYCAM (sugars) + AMBER "
                                "(protonation variants HID/HIE/HIP/ASH/GLH/LYN/"
                                "CYX/CYM). This is the default if no direction "
                                "flag is given.")
    direction.add_argument("--to-charmm", action="store_true",
                           help="Reverse direction: convert GLYCAM/AMBER-named "
                                "input to CHARMM-compatible naming "
                                "(BGLCNA/BMAN/AMAN/ANE5AC/... sugars, "
                                "HSD/HSE/HSP/ASPP/GLUP/LSN protonation "
                                "variants). Glycoprotein residues "
                                "NLN/OLS/OLT revert to ASN/SER/THR; ROH/OME "
                                "caps are dropped. Linkage info preserved via "
                                "CONECT records.")

    content = p.add_argument_group("Content selection")
    content.add_argument("--no-roh", action="store_true",
                         help="Do not add ROH cap at the reducing end (forward / "
                              "--to-amber direction only)")
    content.add_argument("--no-infer-conect", dest="no_infer_conect",
                         action="store_true",
                         help="Skip automatic CONECT inference (default: infer "
                              "missing glycosidic / glycosylation bonds so linkage "
                              "detection works on CONECT-less inputs).")

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument("-v", "--verbose", action="store_true",
                      help="Print conversion details")

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p, batch=True)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    default_suffix = "_charmm" if args.to_charmm else "_amber"
    output_path = (Path(args.output) if args.output
                   else input_path.with_stem(input_path.stem + default_suffix))

    if not args.no_infer_conect:
        from dvbfixer.pdbutils import _materialise_inferred_pdb
        input_path = Path(_materialise_inferred_pdb(
            input_path, verbose=args.verbose))

    if args.to_charmm:
        convert_to_charmm(input_path, output_path, verbose=args.verbose)
        print(f"Wrote {output_path}")
        return

    n_sugars, n_protein = convert_to_glycam(
        input_path, output_path,
        add_roh=not args.no_roh,
        verbose=args.verbose,
    )

    print(f"Converted {n_sugars} sugar residue(s)" +
          (f", {n_protein} protein residue(s)" if n_protein else ""))
    print(f"Wrote {output_path}")
