"""dvbfixer top — Generate GROMACS .itp/.top topology files from PDB.

Parses GROMACS force field RTP/ARN/R2B/TDB files directly and builds
correct topology with proper atom types, charges, bonds, angles,
dihedrals, impropers, and CMAP (CHARMM).
"""

import argparse
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from dvbfixer.rtp_parser import (
    parse_arn,
    parse_atomtypes,
    parse_r2b,
    parse_rtp,
    parse_tdb,
)

# ---------------------------------------------------------------------------
# Force field directories bundled with dvbfixer
# ---------------------------------------------------------------------------
FF_DIR = Path(__file__).parent.parent.parent / 'FF'
FF_CHOICES = {
    'amber': 'amber99sb-ildn-lipid21.ff',
    'charmm': 'charmm36_ljpme-jul2022.ff',
}

# Map standard PDB residue names to the GMX names used in .r2b
# (these are the "GMX" column in aminoacids.r2b)
PDB_TO_GMX = {
    'HIS': 'HISE',    # default HIS → HIE/HSE
    'CYS': 'CYS',
    'ASP': 'ASP',
    'GLU': 'GLU',
    'LYS': 'LYS',
    # AMBER protonation names
    'HIE': 'HISE', 'HID': 'HISD', 'HIP': 'HISH',
    'ASH': 'ASPH', 'GLH': 'GLUH',
    'CYX': 'CYS2', 'CYM': 'CYS',
    'LYN': 'LYSN',
    # CHARMM protonation names
    'HSE': 'HISE', 'HSD': 'HISD', 'HSP': 'HISH',
    'ASPP': 'ASPH', 'GLUP': 'GLUH',
    'LSN': 'LYSN',
    # Non-canonical → canonical
    'MSE': 'MET',
}

# Standard amino acids that can appear in PDB
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL',
}

# PDB→CHARMM atom name mapping for sugars
# PDB uses different naming for NAG acetyl group and some other atoms
CARB_ATOM_MAP = {
    'NAG': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    'NDG': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    'NGA': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    'A2G': {
        'N2': 'N', 'C7': 'C', 'O7': 'O', 'C8': 'CT',
        'H81': 'HT1', 'H82': 'HT2', 'H83': 'HT3',
        'HN2': 'HN',
    },
    # Mannose: PDB may use H11/H12 for C1 hydrogens
    'MAN': {'H11': 'H1'},
    'BMA': {'H11': 'H1'},
    # Galactose/Glucose: similar
    'GAL': {'H11': 'H1'},
    'GLC': {'H11': 'H1'},
}

# PDB sugar names -> CHARMM carb.rtp names
PDB_TO_CARB = {
    'NAG': 'BGLCNA',   # N-acetylglucosamine (beta)
    'NDG': 'BGLCNA',   # 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
    'BMA': 'BMAN',     # beta-mannose
    'MAN': 'AMAN',     # alpha-mannose
    'GAL': 'BGAL',     # beta-galactose
    'GLC': 'BGLC',     # beta-glucose
    'FUC': 'AFUC',     # alpha-fucose
    'FUL': 'BFUC',     # beta-fucose
    'SIA': 'ANE5AC',   # sialic acid (Neu5Ac)
    'NGA': 'BGALNA',   # N-acetylgalactosamine (beta)
    'A2G': 'AGALNA',   # N-acetylgalactosamine (alpha)
    'BGC': 'BGLC',     # beta-glucose
    'XYS': 'BXYL',     # beta-xylose
    'AFU': 'AFUC',     # alpha-L-fucose (alternate PDB code)
    'AMA': 'AMAN',     # alpha-mannose (alternate PDB code)
    'BGA': 'BGAL',     # beta-galactose (alternate PDB code)
    'BGL': 'BGLCNA',   # beta-N-acetylglucosamine (alternate PDB code)
    # CHARMM-GUI native names (already CHARMM RTP names, pass through)
    'BGLC': 'BGLC',
    'BGAL': 'BGAL',
    'AFUC': 'AFUC',
    'AMAN': 'AMAN',
    'BMAN': 'BMAN',
    'BGLCNA': 'BGLCNA',
    'BGALNA': 'BGALNA',
    'AGALNA': 'AGALNA',
    'ANE5AC': 'ANE5AC',
    'ANE5': 'ANE5AC',
    # Non-standard PDB names (from transplant/GLYCAM workflows)
    'AGL': 'AGAL',     # alpha-galactose (_resolve_sugar_rtp auto-detects AGALNA if N-acetyl)  # CHARMM-GUI short name for sialic acid
}

# PDB lipid names -> CHARMM lipid.rtp names (ceramides)
PDB_TO_LIPID = {
    'CER1': 'CER160',   # CHARMM-GUI: CER1 = d18:1/16:0
    'CER2': 'CER2',     # generic ceramide
    'CER3': 'CER3E',    # ceramide variant
}

# CHARMM lipid.rtp ceramide residue names
CERAMIDE_RTP = {
    'CER160', 'CER180', 'CER181', 'CER2', 'CER200',
    'CER220', 'CER240', 'CER241', 'CER3E',
}

# 4-char resnames from CGenFF/solvent that GROMACS writes into PDB cols 17-20
# (standard PDB uses 3-char in cols 18-20, but GROMACS left-justifies from col 17)
_KNOWN_4CHAR_RESNAMES = {
    'ACET', 'ACEH', 'ACEM',  # acetate/acetic acid/acetamide (CGenFF)
    'TIP3', 'SPC', 'SPCE',   # water models
}

# Water residue names (for counting SOL molecules in PDB)
_WATER_RESNAMES = {'SOL', 'HOH', 'WAT', 'TIP3', 'SPC', 'SPCE', 'TIP4', 'TIP5'}

# Known atom name differences between AMBER RTP and PDB/IUPAC naming.
# Key = RTP name, value = PDB name.
# Applied per-residue after detecting which convention the PDB uses.
_EXPLICIT_RENAMES = {
    # ILE: AMBER uses CD, PDB/IUPAC uses CD1
    'CD': 'CD1',
    'HD1': 'HD11', 'HD2': 'HD12', 'HD3': 'HD13',
    # C-terminal: AMBER RTP uses OC1/OC2, PDB uses OXT/O
    'OC1': 'OXT', 'OC2': 'O',
}


def _match_atom_names(rtp_names, pdb_names, arn_rtp_to_pdb=None):
    """Match RTP atom names to PDB atom names, handling naming conventions.

    PDB/OpenMM uses IUPAC naming (HB2/HB3), AMBER RTP uses old naming
    (HB1/HB2). For prochiral methylene hydrogens the numbering is shifted:
      AMBER HB1 = IUPAC HB2,  AMBER HB2 = IUPAC HB3

    Also handles:
    - ILE CD/HD1-3 → CD1/HD11-13
    - C-terminal OC1/OC2 → OXT/O
    - N-terminal H1 ← H (when only H exists in PDB)
    - ARN-based renames (e.g. CHARMM HN → PDB H)

    arn_rtp_to_pdb: dict[rtp_name -> pdb_name] from ARN file (reverse mapping).

    Returns dict[rtp_name -> pdb_name] for matched atoms.
    """
    rtp_set = set(rtp_names)
    pdb_set = set(pdb_names)
    mapping = {}
    used_pdb = set()

    # Pass 0a: ARN-based renames (e.g. CHARMM HN→H)
    if arn_rtp_to_pdb:
        for rtp_name in rtp_names:
            if rtp_name in arn_rtp_to_pdb:
                pdb_name = arn_rtp_to_pdb[rtp_name]
                if pdb_name in pdb_set and rtp_name not in pdb_set:
                    mapping[rtp_name] = pdb_name
                    used_pdb.add(pdb_name)

    # Pass 0b: explicit renames (ILE CD→CD1, C-terminal OC1→OXT, etc.)
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if rtp_name in _EXPLICIT_RENAMES:
            pdb_name = _EXPLICIT_RENAMES[rtp_name]
            if pdb_name in pdb_set and rtp_name not in pdb_set:
                mapping[rtp_name] = pdb_name
                used_pdb.add(pdb_name)

    # Detect if shift mapping applies for H-atom prefixes:
    # RTP has HB1,HB2 and PDB has HB2,HB3 → shift applies for HB prefix
    # Only shift when: lowest RTP number is NOT in PDB, but lowest+1 IS
    shift_prefixes = set()
    # Group RTP H-atoms by prefix
    prefix_groups = defaultdict(list)
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit() and rtp_name[0] == 'H':
            prefix = rtp_name[:-1]
            num = int(rtp_name[-1])
            prefix_groups[prefix].append(num)

    for prefix, nums in prefix_groups.items():
        min_num = min(nums)
        max_num = max(nums)
        min_name = f"{prefix}{min_num}"
        max_shifted = f"{prefix}{max_num + 1}"
        # Shift applies if: lowest RTP name not in PDB AND highest+1 IS in PDB
        # This means the PDB numbering is offset by +1 from RTP for this prefix
        # e.g. RTP {HB1,HB2}, PDB has HB3 but not HB1 → shift
        if min_name not in pdb_set and max_shifted in pdb_set:
            shift_prefixes.add(prefix)

    # N-terminal special case: RTP H1 → PDB H (when PDB has H but not H1)
    # N-terminal residues have H1/H2/H3 in RTP but PDB may have H/H2/H3
    if 'H1' in rtp_set and 'H1' not in pdb_set and 'H' in pdb_set and 'H' not in rtp_set:
        mapping['H1'] = 'H'
        used_pdb.add('H')

    # Pass 1: shift matching for detected prefixes
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit() and rtp_name[0] == 'H':
            prefix = rtp_name[:-1]
            if prefix in shift_prefixes:
                num = int(rtp_name[-1])
                shifted = f"{prefix}{num + 1}"
                if shifted in pdb_set and shifted not in used_pdb:
                    mapping[rtp_name] = shifted
                    used_pdb.add(shifted)
                elif num == 1:
                    # H1 → H (N-terminal special case)
                    base = prefix
                    if base in pdb_set and base not in used_pdb:
                        mapping[rtp_name] = base
                        used_pdb.add(base)

    # Pass 2: exact matches for remaining
    for name in rtp_names:
        if name not in mapping and name in pdb_set and name not in used_pdb:
            mapping[name] = name
            used_pdb.add(name)

    # Pass 3: shift for any remaining unmatched
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit():
            prefix = rtp_name[:-1]
            num = int(rtp_name[-1])
            shifted = f"{prefix}{num + 1}"
            if shifted in pdb_set and shifted not in used_pdb:
                mapping[rtp_name] = shifted
                used_pdb.add(shifted)

    # Pass 4: singleton numbered atom → base name (e.g. CHARMM HG1 → PDB HG)
    for rtp_name in rtp_names:
        if rtp_name in mapping:
            continue
        if len(rtp_name) >= 2 and rtp_name[-1].isdigit():
            prefix = rtp_name[:-1]
            # Only if this is the sole RTP atom with this prefix
            count = sum(1 for n in rtp_names if n.startswith(prefix) and n != prefix
                        and len(n) > len(prefix) and n[len(prefix):].isdigit())
            if count == 1 and prefix in pdb_set and prefix not in used_pdb:
                mapping[rtp_name] = prefix
                used_pdb.add(prefix)

    return mapping


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class PDBResidue:
    chain_id: str
    resname: str
    resseq: int
    icode: str
    atoms: list = field(default_factory=list)  # [(name, x, y, z)]


@dataclass
class PDBChain:
    chain_id: str
    residues: list = field(default_factory=list)  # [PDBResidue]


@dataclass
class AtomEntry:
    index: int          # 1-based
    atom_type: str
    resnr: int          # 1-based residue number in chain
    resname: str
    atomname: str
    cgnr: int
    charge: float
    mass: float
    x: float = 0.0     # coordinates for PDB output
    y: float = 0.0
    z: float = 0.0
    chain_id: str = ' '
    orig_resseq: int = 0
    orig_resname: str = ''


@dataclass
class ChainTopology:
    name: str
    nrexcl: int
    atoms: list = field(default_factory=list)       # [AtomEntry]
    bonds: list = field(default_factory=list)        # [(i, j)]
    pairs: list = field(default_factory=list)        # [(i, j)]
    angles: list = field(default_factory=list)       # [(i, j, k)] or [(i, j, k, ftype)]
    dihedrals: list = field(default_factory=list)    # [(i,j,k,l)] or [(i,j,k,l,type_name)]
    impropers: list = field(default_factory=list)    # [(i, j, k, l)]
    cmap: list = field(default_factory=list)         # [(i, j, k, l, m)]


# ---------------------------------------------------------------------------
# PDB reader
# ---------------------------------------------------------------------------
def read_pdb_chains(path):
    """Read PDB file and extract chains with residues and atoms.

    Detects resseq backward jumps within a chain (e.g. two glycan trees
    with same chain ID and overlapping residue numbers) and splits them
    into separate sub-chains with generated chain IDs.
    """
    # First pass: collect lines per original chain ID, preserving order
    # Handle 4-char resnames (CHARMM-GUI style): if col 21 is not a space
    # and col 17-20 is not blank, the resname extends to col 21 and chain ID
    # is effectively blank.
    chain_lines = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            # Detect 4-char resnames: col 17-20 is the resname field,
            # col 21 is normally chain ID. If col 20 is not space and col 21
            # is not space either, it's likely a 4-char resname with no chain.
            resname_3 = line[17:20].strip()
            chain_id = line[21]
            if resname_3 and chain_id != ' ' and not chain_id.isalpha():
                # Could be 4-char resname (e.g. CER1, BGAL, ANE5, AGLC)
                resname_4 = line[17:21].strip()
                if (resname_4 in PDB_TO_LIPID or resname_4 in PDB_TO_CARB
                        or resname_4 in CERAMIDE_RTP
                        or (len(resname_4) == 4 and resname_3 not in STANDARD_AA)):
                    chain_id = ' '
            chain_lines[chain_id].append(line)

    # Second pass: split chains on resseq backward jumps
    chains = []
    used_ids = set(chain_lines.keys())
    # Pool of available chain IDs for sub-chains
    all_ids = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')

    for orig_id, lines in chain_lines.items():
        # Detect breaks (resseq goes backwards)
        segments = [[]]  # list of line groups
        prev_resseq = -999999
        for line in lines:
            resseq = int(line[22:26])
            if resseq < prev_resseq:
                segments.append([])
            segments[-1].append(line)
            prev_resseq = resseq

        for seg_idx, seg_lines in enumerate(segments):
            if seg_idx == 0:
                cid = orig_id
            else:
                # Assign a new chain ID
                cid = None
                for candidate in all_ids:
                    if candidate not in used_ids:
                        cid = candidate
                        break
                if cid is None:
                    cid = f'{orig_id}{seg_idx}'
                used_ids.add(cid)

            chain = PDBChain(chain_id=cid)
            for line in seg_lines:
                resname = line[17:20].strip()
                # Check for 4-char resname (CHARMM-GUI style or GROMACS output).
                # Use 4-char if it's in known sets, OR if the 3-char truncation
                # isn't a standard amino acid (catches all CHARMM carb/lipid names).
                resname_4 = line[17:21].strip()
                if len(resname_4) == 4:
                    if (resname_4 in PDB_TO_LIPID or resname_4 in PDB_TO_CARB
                            or resname_4 in CERAMIDE_RTP
                            or resname_4 in _KNOWN_4CHAR_RESNAMES
                            or resname not in STANDARD_AA):
                        resname = resname_4
                resseq = int(line[22:26])
                icode = line[26].strip()
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                key = (resseq, icode)
                if not chain.residues or (chain.residues[-1].resseq, chain.residues[-1].icode) != key:
                    chain.residues.append(PDBResidue(
                        chain_id=cid, resname=resname,
                        resseq=resseq, icode=icode,
                    ))
                chain.residues[-1].atoms.append((atom_name, x, y, z))

            chains.append(chain)

    return chains


def _parse_ion_names(ions_itp_path):
    """Parse moleculetype names from ions.itp file."""
    names = set()
    with open(ions_itp_path) as f:
        in_moltype = False
        for line in f:
            stripped = line.strip()
            if stripped == '[ moleculetype ]':
                in_moltype = True
                continue
            if in_moltype and stripped and not stripped.startswith(';'):
                names.add(stripped.split()[0])
                in_moltype = False
            if stripped.startswith('[') and 'moleculetype' not in stripped:
                in_moltype = False
    return names


def _count_molecules(pdb_path, mol_names):
    """Count ion/buffer molecules in PDB by residue name.

    Returns list of (name, count) in order of first appearance,
    preserving the PDB residue ordering for [ molecules ].
    """
    from collections import OrderedDict
    counts = OrderedDict()
    seen_residues = set()  # (chain, resseq, resname) to avoid double-counting atoms
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            resname = line[17:20].strip()
            # Also check 4-char resname (cols 17-20)
            resname4 = line[17:21].strip()
            name = None
            if resname in mol_names:
                name = resname
            elif resname4 in mol_names:
                name = resname4
            if name is None:
                continue
            chain = line[21]
            resseq = int(line[22:26])
            key = (chain, resseq, name)
            if key not in seen_residues:
                seen_residues.add(key)
                counts[name] = counts.get(name, 0) + 1
    return list(counts.items())


def _gro_to_pdb(gro_path):
    """Convert a GROMACS .gro file to a temporary PDB via MDAnalysis."""
    import tempfile
    import warnings
    import MDAnalysis as mda
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(gro_path))
        tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
        tmp.close()
        u.atoms.write(tmp.name)
    return Path(tmp.name)


def _split_chain_by_distance(chain, gap_cutoff=4.0):
    """Split a chain into sub-chains where consecutive residues are > gap_cutoff apart.

    Uses nearest-atom distance between consecutive residues (same approach as
    split_chains.py criterion 3). Molecules that are physically separate
    (e.g. ACET, ACEH in buffer) get split into individual chains.
    """
    import numpy as np

    if len(chain.residues) <= 1:
        return [chain]

    # Build coordinate arrays per residue
    res_coords = []
    for r in chain.residues:
        coords = np.array([(x, y, z) for _, x, y, z in r.atoms])
        res_coords.append(coords)

    # Find breaks: where nearest-atom distance > gap_cutoff
    breaks = [0]
    for i in range(1, len(chain.residues)):
        prev = res_coords[i - 1]
        cur = res_coords[i]
        # Nearest-atom distance
        diff = prev[:, None, :] - cur[None, :, :]
        min_dist = np.sqrt((diff ** 2).sum(axis=2)).min()
        if min_dist > gap_cutoff:
            breaks.append(i)

    if len(breaks) == 1:
        return [chain]  # no splits needed

    # Split into sub-chains
    result = []
    for bi in range(len(breaks)):
        start = breaks[bi]
        end = breaks[bi + 1] if bi + 1 < len(breaks) else len(chain.residues)
        sub = PDBChain(chain_id=chain.chain_id)
        sub.residues = chain.residues[start:end]
        result.append(sub)

    return result


def _count_water(pdb_path):
    """Count water molecules (SOL/HOH/WAT/TIP3) in PDB.

    Uses atom count / 3 (atoms per water) to handle resseq overflow
    in large systems where PDB wraps at 9999.
    """
    water_atoms = 0
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            resname = line[17:20].strip()
            resname4 = line[17:21].strip()
            if resname in _WATER_RESNAMES or resname4 in _WATER_RESNAMES:
                water_atoms += 1
    return water_atoms // 3


def _add_protonation_hydrogens(protein_chains, pdb_path, ff_type, verbose=False):
    """Add missing protonation H atoms using OpenMM Modeller with variants.

    Same approach as protonate.py: load PDB → strip H → strip non-protein →
    build variants list → PDBFixer fix missing atoms → Modeller.addHydrogens
    with CHARMM FF and variants → extract only needed protonation H coords.

    OpenMM variant names (ASH, GLH, HIP, HID, HIE) work with both
    charmm36.xml and amber14-all.xml force fields.
    """
    import tempfile
    import os

    # Map protonated names to OpenMM variant names
    _PROT_TO_VARIANT = {
        'ASPP': 'ASH', 'ASH': 'ASH', 'ASPH': 'ASH',
        'GLUP': 'GLH', 'GLH': 'GLH', 'GLUH': 'GLH',
        'HSP': 'HIP', 'HIP': 'HIP', 'HISH': 'HIP',
        'HSD': 'HID', 'HID': 'HID', 'HISD': 'HID',
        'HSE': 'HIE', 'HIE': 'HIE', 'HISE': 'HIE',
    }
    # Which H atoms we want for each protonated form
    _PROT_H_ATOMS = {
        'ASPP': {'HD2'}, 'ASH': {'HD2'}, 'ASPH': {'HD2'},
        'GLUP': {'HE2'}, 'GLH': {'HE2'}, 'GLUH': {'HE2'},
        'HSP': {'HD1', 'HE2'}, 'HIP': {'HD1', 'HE2'}, 'HISH': {'HD1', 'HE2'},
        'HSD': {'HD1'}, 'HID': {'HD1'}, 'HISD': {'HD1'},
        'HSE': {'HE2'}, 'HIE': {'HE2'}, 'HISE': {'HE2'},
    }

    # Collect residues that need protonation H
    need_h = {}  # (chain_id, resseq) -> (chain_ref, res_ref, prot_name, missing_h)
    for chain in protein_chains:
        for res in chain.residues:
            prot_name = res.resname.upper()
            if prot_name in _PROT_H_ATOMS:
                existing = {a[0] for a in res.atoms}
                missing = _PROT_H_ATOMS[prot_name] - existing
                if missing:
                    need_h[(chain.chain_id, res.resseq)] = (
                        chain, res, prot_name, missing)

    if not need_h:
        return

    if verbose:
        for (cid, rseq), (ch, res, pn, mh) in need_h.items():
            print(f"  Need H for {pn} {cid}:{rseq}: "
                  f"{', '.join(sorted(mh))}")

    from openmm.app import ForceField, Modeller, PDBFile
    from openmm import unit
    from pdbfixer import PDBFixer
    from dvbfixer.protonate import _strip_hydrogens
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    # Build variant lookup: (chain_id, resseq) -> OpenMM variant name
    variant_lookup = {}
    for (cid, rseq), (ch, res, pn, mh) in need_h.items():
        variant = _PROT_TO_VARIANT.get(pn)
        if variant:
            variant_lookup[(cid, rseq)] = variant

    # Load PDB with OpenMM (same approach as protonate.py)
    pdb = PDBFile(str(pdb_path))

    # Strip existing hydrogens
    topology, positions = _strip_hydrogens(pdb.topology, pdb.positions)

    # Strip non-protein residues (glycans, ligands, etc.)
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    to_delete = [res for res in topology.residues() if res.name not in known]
    if to_delete:
        modeller = Modeller(topology, positions)
        modeller.delete(to_delete)
        topology, positions = modeller.topology, modeller.positions

    # Build variants list
    def _build_variants(topo):
        vlist = []
        for res in topo.residues():
            cid = res.chain.id
            try:
                rseq = int(res.id)
            except ValueError:
                vlist.append(None)
                continue
            key = (cid, rseq)
            if key in variant_lookup:
                vlist.append(variant_lookup[key])
            else:
                vlist.append(None)
        return vlist

    variants = _build_variants(topology)

    # Use PDBFixer to fix any missing heavy atoms
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.pdb', delete=False) as f:
            PDBFile.writeFile(topology, positions, f, keepIds=True)
            tmp_path = f.name

        fixer = PDBFixer(filename=tmp_path)
        fixer.findMissingResidues()
        fixer.missingResidues = {}  # Only fix atoms, not residues
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        modeller = Modeller(fixer.topology, fixer.positions)

        # Rebuild variants for potentially reordered topology
        variants = _build_variants(modeller.topology)

        # Add H with proper geometry — use CHARMM or AMBER FF
        if ff_type == 'charmm':
            ff = ForceField('charmm36.xml', 'charmm36/water.xml')
        else:
            ff = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller.addHydrogens(ff, variants=variants)

        # Extract only the specific protonation H atoms we need
        added = 0
        for atom in modeller.topology.atoms():
            res = atom.residue
            cid = res.chain.id
            try:
                rseq = int(res.id)
            except ValueError:
                continue
            key = (cid, rseq)
            if key not in need_h:
                continue
            chain_ref, res_ref, prot_name, missing_h = need_h[key]
            if atom.name in missing_h:
                pos = modeller.positions[atom.index]
                xyz = pos.value_in_unit(unit.angstrom)
                res_ref.atoms.append((atom.name, xyz[0], xyz[1], xyz[2]))
                added += 1
                if verbose:
                    print(f"    Placed {atom.name} at {res_ref.resname} "
                          f"{cid}:{rseq} "
                          f"({xyz[0]:.3f}, {xyz[1]:.3f}, {xyz[2]:.3f})")

    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.unlink(tmp_path)

    if verbose:
        print(f"  Added {added} protonation H atom(s) via OpenMM Modeller")


def _extract_molecule_lines(pdb_path, mol_names):
    """Extract ATOM/HETATM lines for ion/buffer molecules from PDB."""
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            resname = line[17:20].strip()
            resname4 = line[17:21].strip()
            if resname in mol_names or resname4 in mol_names:
                lines.append(line)
    return lines


def read_ssbonds(path):
    """Read SS bonds from SSBOND records or auto-detect from CYS SG-SG distances.

    Returns [(chain1, resseq1, chain2, resseq2)].
    """
    ssbonds = []

    # Try SSBOND records first
    with open(path) as f:
        for line in f:
            if line.startswith('SSBOND'):
                ch1 = line[15]
                res1 = int(line[17:21])
                ch2 = line[29]
                res2 = int(line[31:35])
                ssbonds.append((ch1, res1, ch2, res2))

    if ssbonds:
        return ssbonds

    # Auto-detect from CYS SG atom distances (< 2.5 A)
    import math
    sg_atoms = []  # [(chain, resseq, x, y, z)]
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            resname = line[17:20].strip()
            atom_name = line[12:16].strip()
            if resname == 'CYS' and atom_name == 'SG':
                chain = line[21]
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                sg_atoms.append((chain, resseq, x, y, z))

    # Find SG pairs within 2.5 A
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            ch1, r1, x1, y1, z1 = sg_atoms[i]
            ch2, r2, x2, y2, z2 = sg_atoms[j]
            dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            if dist < 2.5:
                ssbonds.append((ch1, r1, ch2, r2))

    return ssbonds


def _is_ceramide(resname):
    """Check if a residue name is a ceramide (PDB or RTP name)."""
    return resname in PDB_TO_LIPID or resname in CERAMIDE_RTP


def _resolve_sugar_rtp(resname, atoms, residues_dict):
    """Resolve a sugar PDB/CHARMM name to its correct RTP entry.

    Handles special cases like BGAL with N-acetyl atoms -> BGALNA.
    """
    # First try direct RTP match (CHARMM-GUI native names)
    if resname in residues_dict:
        rtp_name = resname
    else:
        rtp_name = PDB_TO_CARB.get(resname)

    if rtp_name is None:
        return None

    # Auto-detect: BGAL/AGAL with N-acetyl atoms -> BGALNA/AGALNA
    if rtp_name in ('BGAL', 'AGAL'):
        pdb_atom_names = {a[0] for a in atoms} if isinstance(atoms, list) else atoms
        has_nacetyl = bool(pdb_atom_names & {'N', 'HN', 'CT'})
        if has_nacetyl:
            alt = 'BGALNA' if rtp_name == 'BGAL' else 'AGALNA'
            if alt in residues_dict:
                rtp_name = alt

    if rtp_name not in residues_dict:
        return None

    return rtp_name


def detect_glycan_links(chains):
    """Auto-detect glycosidic bonds from C1-O distances using parsed chain data.

    Detects:
    - Sugar-sugar links: C1 of one sugar within 2.0 A of Ox of another
    - Protein-sugar links: C1 of NAG within 2.0 A of ASN ND2
    - Sugar-ceramide links: C1 of sugar within 2.0 A of ceramide O1

    Returns list of GlycanLink tuples.
    """
    sugar_names = set(PDB_TO_CARB.keys())
    links = []

    # Also detect sugars already named with CHARMM RTP names
    charmm_sugar_names = set(PDB_TO_CARB.values())

    # Collect C1 atoms from sugars and ND2/O* from sugars+ASN
    c1_atoms = []
    target_atoms = []
    # Ceramide C1S atoms (for sugar O1-ceramide C1S detection)
    cer_c1s_atoms = []
    # Sugar O atoms (for reverse check against ceramide)
    sugar_o_atoms = []

    for chain in chains:
        ch = chain.chain_id
        for res in chain.residues:
            is_sugar = res.resname in sugar_names or res.resname in charmm_sugar_names
            # Sialic acid variants: anomeric carbon is C2 (not C1)
            is_sialic = res.resname in ('SIA', 'ANE5', 'ANE5AC', 'BNE5AC')
            for aname, x, y, z in res.atoms:
                if is_sugar and aname == 'C1' and not is_sialic:
                    c1_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))
                # Sialic acid links via C2
                if is_sugar and is_sialic and aname == 'C2':
                    c1_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

                if is_sugar and aname.startswith('O') and aname != 'O5':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))
                    sugar_o_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

                if res.resname == 'ASN' and aname == 'ND2':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

                # Collect ceramide C1S for sugar-ceramide bond detection
                if _is_ceramide(res.resname) and aname == 'C1S':
                    cer_c1s_atoms.append((ch, res.resname, res.resseq, x, y, z))

    # Find C1/C2-target pairs within 2.0 A (sugar-sugar and protein-sugar)
    for ch1, rn1, rs1, aname1, x1, y1, z1 in c1_atoms:
        for ch2, rn2, rs2, aname2, x2, y2, z2 in target_atoms:
            if ch1 == ch2 and rs1 == rs2:
                continue
            d = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            if d < 2.0:
                # Link: acceptor (has C1/C2) bonds to donor (has Ox/ND2)
                links.append((ch2, rs2, aname2, ch1, rs1, aname1))

    # Find sugar O - ceramide C1S bonds within 2.0 A
    # In glycolipids, sugar O1 bridges to ceramide C1S (ceramide's own O1/HO1
    # are already removed in CHARMM-GUI output)
    for ch1, rn1, rs1, aname1, x1, y1, z1 in sugar_o_atoms:
        for ch2, rn2, rs2, x2, y2, z2 in cer_c1s_atoms:
            d = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            if d < 2.0:
                # Link: ceramide (donor, has C1S) <- sugar (has Ox)
                links.append((ch2, rs2, 'C1S', ch1, rs1, aname1))

    return links


def build_glycan_trees(glycan_links, chains):
    """Build glycan trees from linkage list.

    Returns list of glycan trees, each a list of (chain_id, resseq, resname)
    in topological order (root first).
    Also detects ceramide-sugar links for glycolipids.
    """
    # Build adjacency: sugar -> list of child sugars
    sugar_residues = set()
    charmm_sugar_names = set(PDB_TO_CARB.values())
    for chain in chains:
        for res in chain.residues:
            if res.resname in PDB_TO_CARB or res.resname in charmm_sugar_names:
                sugar_residues.add((chain.chain_id, res.resseq))

    # Identify ceramide residues
    ceramide_residues = set()
    for chain in chains:
        for res in chain.residues:
            if _is_ceramide(res.resname):
                ceramide_residues.add((chain.chain_id, res.resseq))

    # Find root sugars (connected to protein ASN, not to another sugar)
    children = defaultdict(list)  # parent -> [child]
    parent_of = {}  # child -> parent
    protein_links = []  # [(protein_chain, protein_resseq, protein_atom, sugar_chain, sugar_resseq)]
    ceramide_links = []  # [(cer_chain, cer_resseq, cer_atom, sugar_chain, sugar_resseq, sugar_atom)]

    for don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom in glycan_links:
        acc_key = (acc_ch, acc_rs)
        don_key = (don_ch, don_rs)
        if don_key in sugar_residues and acc_key in sugar_residues:
            # sugar-sugar link: donor is parent, acceptor is child
            children[don_key].append((acc_key, don_atom, acc_atom))
            parent_of[acc_key] = (don_key, don_atom)
        elif don_key in ceramide_residues and acc_key in sugar_residues:
            # ceramide-sugar link — do NOT add to parent_of so the sugar
            # becomes a root of its glycan tree (ceramide handled separately)
            ceramide_links.append((don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom))
        elif acc_key in sugar_residues:
            # protein-sugar link
            protein_links.append((don_ch, don_rs, don_atom, acc_ch, acc_rs))

    # Find roots (sugars with no sugar parent)
    roots = []
    for ch, rs in sugar_residues:
        if (ch, rs) not in parent_of:
            roots.append((ch, rs))

    # BFS from each root to build tree
    trees = []
    for root in roots:
        tree = []
        # Find linkage info for this root
        link_atoms = {}  # (ch, rs) -> set of linked O positions
        queue = [root]
        visited = set()
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            tree.append(node)
            for child, donor_atom, child_atom in children.get(node, []):
                # Record which O on parent is linked
                if node not in link_atoms:
                    link_atoms[node] = set()
                link_atoms[node].add(donor_atom)
                # Child's anomeric O is linked: O1 for C1 linkage, O2 for C2 linkage
                if child_atom in ('C1', 'C2'):
                    linked_o = 'O1' if child_atom == 'C1' else 'O2'
                    if child not in link_atoms:
                        link_atoms[child] = set()
                    link_atoms[child].add(linked_o)
                queue.append(child)

        # Mark root's O1 as linked if bonded to ceramide
        for cl in ceramide_links:
            if (cl[3], cl[4]) == root:
                if root not in link_atoms:
                    link_atoms[root] = set()
                link_atoms[root].add('O1')

        # Also mark root's O1 as linked (bonds to protein or is reducing end)
        for pl in protein_links:
            if (pl[3], pl[4]) == root:
                if root not in link_atoms:
                    link_atoms[root] = set()
                link_atoms[root].add('O1')

        # Filter ceramide_links relevant to this tree
        tree_set = set(tree)
        tree_cer_links = [cl for cl in ceramide_links
                          if (cl[3], cl[4]) in tree_set]

        trees.append((tree, link_atoms, protein_links, tree_cer_links))

    return trees


# ---------------------------------------------------------------------------
# Topology builder
# ---------------------------------------------------------------------------
class TopologyBuilder:
    def __init__(self, ff_dir, ff_type='amber', verbose=False):
        self.ff_dir = Path(ff_dir)
        self.ff_type = ff_type
        self.verbose = verbose

        # Parse all FF files
        rtp_file = self.ff_dir / ('aminoacids.rtp' if ff_type == 'amber' else 'aminoacids.rtp')
        self.bonded_types, self.residues = parse_rtp(rtp_file)

        # For CHARMM, also parse merged.rtp if it exists (has more residues)
        merged_rtp = self.ff_dir / 'merged.rtp'
        if merged_rtp.exists():
            bt2, res2 = parse_rtp(merged_rtp)
            # merged.rtp is the main file for CHARMM, aminoacids.rtp may be subset
            self.residues.update(res2)
            self.bonded_types = bt2

        # Also load carb.rtp for sugar residues
        carb_rtp = self.ff_dir / 'carb.rtp'
        if carb_rtp.exists():
            carb_bt, carb_res = parse_rtp(carb_rtp)
            self.residues.update(carb_res)
            self.carb_bonded_types = carb_bt
        else:
            self.carb_bonded_types = None

        # Load all other molecule-type RTP files (CHARMM has lipids, NA, etc.)
        for rtp_name in ['lipid.rtp', 'na.rtp', 'cgenff.rtp', 'ethers.rtp',
                         'metals.rtp', 'silicates.rtp', 'solvent.rtp']:
            rtp_path = self.ff_dir / rtp_name
            if rtp_path.exists():
                _, extra_res = parse_rtp(rtp_path)
                self.residues.update(extra_res)

        self.r2b = parse_r2b(self.ff_dir / 'aminoacids.r2b')

        # Load all R2B files
        for r2b_name in ['carb.r2b', 'lipid.r2b', 'na.r2b', 'cgenff.r2b',
                         'ethers.r2b', 'metals.r2b', 'silicates.r2b', 'solvent.r2b']:
            r2b_path = self.ff_dir / r2b_name
            if r2b_path.exists():
                extra_r2b = parse_r2b(r2b_path)
                self.r2b.update(extra_r2b)

        self.arn = parse_arn(self.ff_dir / 'aminoacids.arn')
        self.atom_masses = parse_atomtypes(self.ff_dir / 'atomtypes.atp')

        # Terminal patches (CHARMM uses these, AMBER has empty TDB)
        n_tdb = self.ff_dir / 'aminoacids.n.tdb'
        c_tdb = self.ff_dir / 'aminoacids.c.tdb'
        self.n_patches = parse_tdb(n_tdb) if n_tdb.exists() else {}
        self.c_patches = parse_tdb(c_tdb) if c_tdb.exists() else {}

        # Build reverse ARN: (resname, ff_name) -> gromacs_name
        # For PDB->FF mapping we need (resname, pdb_name) -> ff_name
        # ARN file format is: resname gromacs_name ff_name
        # gromacs_name = what PDB uses, ff_name = what RTP uses
        self.arn_reverse = {}
        for (resname, gmx_name), ff_name in self.arn.items():
            self.arn_reverse[(resname, ff_name)] = gmx_name

    def _resolve_resname(self, pdb_resname, position, chain_ss_residues):
        """Map PDB residue name to RTP building block name.

        position: 'nter', 'cter', 'mid', 'twter' (both terminals = single residue chain)
        chain_ss_residues: set of resseq numbers involved in SS bonds
        """
        # Normalize non-canonical names
        gmx_name = PDB_TO_GMX.get(pdb_resname, pdb_resname)

        # Check if in r2b mapping
        if gmx_name not in self.r2b:
            # Try as-is (some residues use their PDB name directly)
            if pdb_resname in self.r2b:
                gmx_name = pdb_resname
            elif pdb_resname in self.residues:
                return pdb_resname  # Direct RTP match, no terminal variant needed
            else:
                return None

        main, nter, cter, twter = self.r2b[gmx_name]

        if position == 'twter' and twter != '-':
            return twter
        elif position == 'nter' and nter != '-':
            return nter
        elif position == 'cter' and cter != '-':
            return cter
        else:
            return main

    def _map_atom_name(self, rtp_resname, pdb_atom_name):
        """Map a PDB atom name to the name used in the RTP entry.

        ARN maps gromacs_name -> ff_name. PDB names match gromacs_names.
        We need to convert PDB name to RTP name (which is the ff_name).
        """
        # Check specific residue mapping first
        key = (rtp_resname, pdb_atom_name)
        if key in self.arn:
            return self.arn[key]

        # Check wildcard mapping
        key = ('*', pdb_atom_name)
        if key in self.arn:
            return self.arn[key]

        # No mapping needed — name is the same
        return pdb_atom_name

    def _get_pdb_name(self, rtp_resname, rtp_atom_name):
        """Get the PDB/GROMACS name for an RTP atom (reverse of _map_atom_name)."""
        key = (rtp_resname, rtp_atom_name)
        if key in self.arn_reverse:
            return self.arn_reverse[key]
        # Check wildcard
        for (resn, ff_name), gmx_name in self.arn_reverse.items():
            if resn == '*' and ff_name == rtp_atom_name:
                return gmx_name
        return rtp_atom_name

    def build_chain(self, chain, ss_residues=None, ss_pairs=None):
        """Build topology for a single chain.

        ss_residues: set of resseq involved in SS bonds in this chain.
        ss_pairs: list of (resseq1, resseq2) intra-chain SS bond pairs.
        """
        if ss_residues is None:
            ss_residues = set()
        if ss_pairs is None:
            ss_pairs = []

        n_res = len(chain.residues)
        if n_res == 0:
            return None

        # Step 1: Resolve RTP names for each residue
        rtp_names = []
        for i, res in enumerate(chain.residues):
            if n_res == 1:
                pos = 'twter'
            elif i == 0:
                pos = 'nter'
            elif i == n_res - 1:
                pos = 'cter'
            else:
                pos = 'mid'

            # Check SS bond
            pdb_name = res.resname
            if pdb_name == 'CYS' and res.resseq in ss_residues:
                pdb_name = 'CYX'

            rtp_name = self._resolve_resname(pdb_name, pos, ss_residues)
            if rtp_name is None:
                print(f"WARNING: Residue {res.resname} {res.chain_id}:{res.resseq} "
                      f"not found in force field", file=sys.stderr)
                return None

            if rtp_name not in self.residues:
                print(f"WARNING: RTP entry '{rtp_name}' not found for "
                      f"{res.resname} {res.chain_id}:{res.resseq}", file=sys.stderr)
                return None

            rtp_names.append(rtp_name)
            if self.verbose:
                label = f"  {res.chain_id}:{res.resname}{res.resseq}"
                if rtp_name != res.resname:
                    label += f" -> {rtp_name}"
                print(label)

        # Step 2: Build atom list from RTP entries, matched to PDB
        chain_top = ChainTopology(
            name=f"Protein_chain_{chain.chain_id.strip() or 'X'}",
            nrexcl=self.bonded_types.nrexcl,
        )

        # Map: (residue_index, rtp_atom_name) -> global atom index (1-based)
        atom_index_map = {}
        global_idx = 0
        cgnr_offset = 0

        for res_i, (res, rtp_name) in enumerate(zip(chain.residues, rtp_names)):
            rtp_res = self.residues[rtp_name]
            pdb_atom_names = {a[0] for a in res.atoms}
            rtp_atom_names = [a[0] for a in rtp_res.atoms]

            # Build per-residue ARN reverse mapping (rtp_name → pdb_name)
            arn_rtp_to_pdb = {}
            for rtp_aname in rtp_atom_names:
                # Check residue-specific ARN
                key = (rtp_name, rtp_aname)
                if key in self.arn_reverse:
                    arn_rtp_to_pdb[rtp_aname] = self.arn_reverse[key]
                # Check wildcard ARN
                for (resn, ff_name), gmx_name in self.arn_reverse.items():
                    if resn == '*' and ff_name == rtp_aname:
                        arn_rtp_to_pdb[rtp_aname] = gmx_name
                        break

            # Build RTP→PDB atom name mapping
            rtp_to_pdb = _match_atom_names(rtp_atom_names, pdb_atom_names,
                                           arn_rtp_to_pdb)

            # Build PDB atom coordinate lookup
            pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in res.atoms}

            # Protonation H atoms (HD2/HE2/HD1) are added to res.atoms
            # by _add_protonation_hydrogens() before build_chain is called.
            # They should already be in pdb_coords via the res.atoms loop above.

            for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
                pdb_name = rtp_to_pdb.get(atom_name)
                if pdb_name is None:
                    # Atom is in RTP but not in PDB — skip it
                    if self.verbose:
                        print(f"    Skipping {rtp_name}:{atom_name} "
                              f"(not in PDB {res.chain_id}:{res.resseq})")
                    continue

                global_idx += 1
                mass = self.atom_masses.get(atom_type, 0.0)
                x, y, z = pdb_coords.get(pdb_name, (0.0, 0.0, 0.0))
                chain_top.atoms.append(AtomEntry(
                    index=global_idx,
                    atom_type=atom_type,
                    resnr=res_i + 1,
                    resname=rtp_name,
                    atomname=pdb_name,
                    cgnr=cgnr + cgnr_offset,
                    charge=charge,
                    mass=mass,
                    x=x, y=y, z=z,
                    chain_id=res.chain_id,
                    orig_resseq=res.resseq,
                    orig_resname=res.resname,
                ))
                atom_index_map[(res_i, atom_name)] = global_idx

            cgnr_offset += max((cgnr for _, _, _, cgnr in rtp_res.atoms), default=0)

        # Step 3: Build bonds, resolving inter-residue references
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for a1, a2 in rtp_res.bonds:
                idx1 = self._resolve_atom_ref(a1, res_i, atom_index_map, n_res)
                idx2 = self._resolve_atom_ref(a2, res_i, atom_index_map, n_res)
                if idx1 is not None and idx2 is not None:
                    bond = (min(idx1, idx2), max(idx1, idx2))
                    chain_top.bonds.append(bond)

        # Deduplicate bonds
        chain_top.bonds = sorted(set(chain_top.bonds))

        # Step 3b: Add intra-chain SS bonds (SG-SG)
        # RTP CYS2 has CB-SG but not the inter-residue SG-SG bond
        resseq_to_resi = {res.resseq: i for i, res in enumerate(chain.residues)}
        for res1_seq, res2_seq in ss_pairs:
            resi1 = resseq_to_resi.get(res1_seq)
            resi2 = resseq_to_resi.get(res2_seq)
            if resi1 is None or resi2 is None:
                continue
            sg1 = atom_index_map.get((resi1, 'SG'))
            sg2 = atom_index_map.get((resi2, 'SG'))
            if sg1 is not None and sg2 is not None:
                bond = (min(sg1, sg2), max(sg1, sg2))
                chain_top.bonds.append(bond)
                if self.verbose:
                    print(f"  Added intra-chain SS bond: "
                          f"CYS2 {chain.chain_id}:{res1_seq}:SG - "
                          f"CYS2 {chain.chain_id}:{res2_seq}:SG")
        chain_top.bonds = sorted(set(chain_top.bonds))

        # Step 4: Apply terminal patches (CHARMM)
        if self.ff_type == 'charmm' and n_res > 0:
            self._apply_terminal_patches(chain, chain_top, rtp_names, atom_index_map)

        # Step 5: Build bond graph for angle/dihedral generation
        adj = defaultdict(set)
        for i, j in chain_top.bonds:
            adj[i].add(j)
            adj[j].add(i)

        # Step 6: Generate angles
        chain_top.angles = self._generate_angles(adj)

        # Step 7: Generate proper dihedrals
        chain_top.dihedrals = self._generate_dihedrals(
            adj, rtp_names, atom_index_map, n_res
        )

        # Step 8: Generate 1-4 pairs from dihedrals
        chain_top.pairs = self._generate_pairs(chain_top.dihedrals)

        # Step 9: Resolve impropers from RTP
        chain_top.impropers = self._resolve_impropers(rtp_names, atom_index_map, n_res)

        # Step 10: Resolve CMAP (CHARMM)
        chain_top.cmap = self._resolve_cmap(rtp_names, atom_index_map, n_res)

        # Step 11: Renumber atoms to be contiguous (patches may create gaps)
        self._renumber_atoms(chain_top)

        return chain_top

    def _resolve_atom_ref(self, ref, res_i, atom_index_map, n_res):
        """Resolve atom reference like '-C', '+N', or 'CA' to global index."""
        if ref.startswith('-'):
            target_res = res_i - 1
            atom_name = ref[1:]
        elif ref.startswith('+'):
            target_res = res_i + 1
            atom_name = ref[1:]
        else:
            target_res = res_i
            atom_name = ref

        if target_res < 0 or target_res >= n_res:
            return None

        return atom_index_map.get((target_res, atom_name))

    def _renumber_atoms(self, chain_top):
        """Renumber all atom indices to be contiguous 1..N based on list order."""
        if not chain_top.atoms:
            return
        # Build remap based on current position in the atoms list
        remap = {}
        for new_idx, atom in enumerate(chain_top.atoms, 1):
            remap[atom.index] = new_idx

        for atom in chain_top.atoms:
            atom.index = remap[atom.index]

        def remap_tuple(t):
            return tuple(remap.get(x, x) if isinstance(x, int) else x for x in t)

        chain_top.bonds = [remap_tuple(b) for b in chain_top.bonds]
        chain_top.pairs = [remap_tuple(p) for p in chain_top.pairs]
        chain_top.angles = [remap_tuple(a) for a in chain_top.angles]
        chain_top.dihedrals = [remap_tuple(d) for d in chain_top.dihedrals]
        chain_top.impropers = [remap_tuple(i) for i in chain_top.impropers]
        chain_top.cmap = [remap_tuple(c) for c in chain_top.cmap]

    def _apply_terminal_patches(self, chain, chain_top, rtp_names, atom_index_map):
        """Apply CHARMM terminal patches (TDB). Only for protein chains."""
        n_res = len(chain.residues)

        # Only apply patches to protein residues (not CGenFF small molecules)
        first_res = chain.residues[0].resname
        if first_res not in STANDARD_AA and first_res not in PDB_TO_GMX:
            return
        if first_res == 'GLY':
            n_patch_name = 'GLY-NH3+'
        elif first_res == 'PRO':
            n_patch_name = 'PRO-NH2+'
        else:
            n_patch_name = 'NH3+'

        # Determine C-terminal patch
        c_patch_name = 'COO-'

        # Apply N-terminal patch
        if n_patch_name in self.n_patches:
            self._apply_patch(self.n_patches[n_patch_name], 0, chain_top, atom_index_map)

        # Apply C-terminal patch
        if c_patch_name in self.c_patches:
            self._apply_patch(self.c_patches[c_patch_name], n_res - 1, chain_top, atom_index_map)

    def _apply_patch(self, patch, res_i, chain_top, atom_index_map):
        """Apply a single terminal patch to the chain topology."""
        # Save coordinates of atoms that will be deleted (for reuse by added atoms)
        deleted_coords = {}
        for atom_name in patch.delete:
            idx = atom_index_map.get((res_i, atom_name))
            if idx is not None:
                for atom in chain_top.atoms:
                    if atom.index == idx:
                        deleted_coords[atom_name] = (atom.x, atom.y, atom.z)
                        break

        # Delete atoms
        delete_indices = set()
        for atom_name in patch.delete:
            idx = atom_index_map.get((res_i, atom_name))
            if idx is not None:
                delete_indices.add(idx)

        if delete_indices:
            chain_top.atoms = [a for a in chain_top.atoms if a.index not in delete_indices]
            chain_top.bonds = [(i, j) for i, j in chain_top.bonds
                               if i not in delete_indices and j not in delete_indices]

        # Replace atoms (change type, mass, charge)
        for name, new_type, mass, charge in patch.replace:
            idx = atom_index_map.get((res_i, name))
            if idx is not None:
                for atom in chain_top.atoms:
                    if atom.index == idx:
                        atom.atom_type = new_type
                        atom.mass = mass
                        atom.charge = charge
                        break

        # Add atoms — insert after reference atom in the residue
        if patch.add:
            max_idx = max(a.index for a in chain_top.atoms) if chain_top.atoms else 0
            resnr = None
            resname = 'UNK'
            for atom in chain_top.atoms:
                if atom_index_map.get((res_i, atom.atomname)) == atom.index:
                    resnr = atom.resnr
                    resname = atom.resname
                    break

            for add_entry in patch.add:
                count = add_entry['count']
                name_base = add_entry['name']
                atype = add_entry['type']
                mass = add_entry['mass']
                charge = add_entry['charge']
                cgnr_val = add_entry['cgnr']

                # Find insertion position: after the first reference atom
                ref_atoms = add_entry['ref_atoms']
                insert_after_idx = None
                if ref_atoms:
                    insert_after_idx = atom_index_map.get((res_i, ref_atoms[0]))

                # Find position in atoms list to insert + get reference coords
                insert_pos = len(chain_top.atoms)
                ref_x, ref_y, ref_z = 0.0, 0.0, 0.0
                ref_chain_id = ' '
                if insert_after_idx is not None:
                    for pos, atom in enumerate(chain_top.atoms):
                        if atom.index == insert_after_idx:
                            insert_pos = pos + 1
                            ref_x, ref_y, ref_z = atom.x, atom.y, atom.z
                            ref_chain_id = atom.chain_id
                            break

                # Build coordinate list for added atoms:
                # - Use deleted atom coords when available (e.g. OT1 gets O's coords)
                # - Otherwise use reference atom coords with small offsets to avoid overlap
                add_coords = []
                for k in range(count):
                    atom_name = f"{name_base}{k + 1}" if count > 1 else name_base
                    # Try to find a matching deleted atom's coordinates
                    # COO-: deleted O → use for OT1; offset for OT2
                    # NH3+: deleted HN → use for H1; offset for H2, H3
                    coord_found = False
                    if k == 0:
                        # First added atom: try deleted atoms that share the same
                        # element (O→OT1, HN→H1)
                        for dname, dcoord in deleted_coords.items():
                            # Match by element: H* deleted → first H added
                            if dname[0] == name_base[0]:
                                add_coords.append(dcoord)
                                coord_found = True
                                break
                    if not coord_found:
                        # Offset from reference atom to avoid overlap
                        offset = 0.1 * (k + 1)  # 1 Angstrom increments
                        add_coords.append((ref_x + offset, ref_y + offset, ref_z))

                new_atoms = []
                for k in range(count):
                    max_idx += 1
                    atom_name = f"{name_base}{k + 1}" if count > 1 else name_base
                    actual_cgnr = cgnr_val if cgnr_val > 0 else (
                        chain_top.atoms[insert_pos - 1].cgnr if insert_pos > 0 else 1
                    )
                    ax, ay, az = add_coords[k]
                    new_atom = AtomEntry(
                        index=max_idx,
                        atom_type=atype,
                        resnr=resnr or (res_i + 1),
                        resname=resname,
                        atomname=atom_name,
                        cgnr=actual_cgnr,
                        charge=charge,
                        mass=mass,
                        x=ax, y=ay, z=az,
                        chain_id=ref_chain_id,
                        orig_resseq=chain_top.atoms[insert_pos - 1].orig_resseq if insert_pos > 0 else 0,
                        orig_resname=chain_top.atoms[insert_pos - 1].orig_resname if insert_pos > 0 else '',
                    )
                    new_atoms.append(new_atom)
                    atom_index_map[(res_i, atom_name)] = max_idx

                    # Add bond to reference atom
                    if insert_after_idx is not None:
                        bond = (min(insert_after_idx, max_idx), max(insert_after_idx, max_idx))
                        chain_top.bonds.append(bond)

                # Insert at the right position
                for offset, new_atom in enumerate(new_atoms):
                    chain_top.atoms.insert(insert_pos + offset, new_atom)

        # Add impropers from patch
        for imp in patch.impropers:
            indices = []
            for atom_name in imp:
                idx = atom_index_map.get((res_i, atom_name))
                if idx is not None:
                    indices.append(idx)
            if len(indices) == 4:
                chain_top.impropers.append(tuple(indices))

    def _generate_angles(self, adj):
        """Generate all angles from bond graph."""
        angles = set()
        for j in sorted(adj.keys()):
            neighbors = sorted(adj[j])
            for idx_a, i in enumerate(neighbors):
                for k in neighbors[idx_a + 1:]:
                    angles.add((i, j, k))
        return sorted(angles)

    def _generate_dihedrals(self, adj, rtp_names, atom_index_map, n_res):
        """Generate proper dihedrals from bond graph + RTP explicit dihedrals."""
        # Generate all possible dihedrals from connectivity
        generated = set()
        for j in sorted(adj.keys()):
            for k in sorted(adj[j]):
                if k <= j:
                    continue
                for i in sorted(adj[j]):
                    if i == k:
                        continue
                    for l in sorted(adj[k]):
                        if l == j:
                            continue
                        generated.add((i, j, k, l))

        dihedrals = sorted(generated)

        # Add explicit RTP dihedrals (AMBER ILDN corrections etc.)
        explicit = []
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for dih in rtp_res.dihedrals:
                indices = []
                for ref in dih[:4]:
                    idx = self._resolve_atom_ref(ref, res_i, atom_index_map, n_res)
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    if len(dih) == 5:
                        explicit.append((*indices, dih[4]))
                    else:
                        explicit.append(tuple(indices))

        # Merge: explicit dihedrals with type names go at the end
        result = [(i, j, k, l) for i, j, k, l in dihedrals]
        result.extend(explicit)
        return result

    def _generate_pairs(self, dihedrals):
        """Generate 1-4 pairs from proper dihedrals."""
        pairs = set()
        for dih in dihedrals:
            i, l = dih[0], dih[3]
            pair = (min(i, l), max(i, l))
            pairs.add(pair)
        return sorted(pairs)

    def _resolve_impropers(self, rtp_names, atom_index_map, n_res):
        """Resolve all improper dihedrals from RTP entries."""
        impropers = []
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for imp in rtp_res.impropers:
                indices = []
                for ref in imp:
                    idx = self._resolve_atom_ref(ref, res_i, atom_index_map, n_res)
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    impropers.append(tuple(indices))
        return impropers

    def _resolve_cmap(self, rtp_names, atom_index_map, n_res):
        """Resolve CMAP entries from RTP."""
        cmaps = []
        for res_i, rtp_name in enumerate(rtp_names):
            rtp_res = self.residues[rtp_name]
            for cm in rtp_res.cmap:
                indices = []
                for ref in cm:
                    idx = self._resolve_atom_ref(ref, res_i, atom_index_map, n_res)
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 5:
                    cmaps.append(tuple(indices))
        return cmaps

    def build_glycan_chain(self, tree, link_atoms, all_chains, glycan_links):
        """Build topology for a glycan tree (one moleculetype).

        tree: list of (chain_id, resseq) in topological order
        link_atoms: dict (chain_id, resseq) -> set of O atoms that are linked
        all_chains: list of PDBChain objects
        glycan_links: list of (don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom)
        """
        # Build residue lookup from PDB chains
        res_lookup = {}  # (chain, resseq) -> PDBResidue
        for chain in all_chains:
            for res in chain.residues:
                res_lookup[(chain.chain_id, res.resseq)] = res

        bt = self.carb_bonded_types or self.bonded_types
        first_ch, first_rs = tree[0]
        first_res = res_lookup.get((first_ch, first_rs))
        chain_name = f"Glycan_{first_ch.strip() or 'X'}_{first_rs}"

        chain_top = ChainTopology(
            name=chain_name,
            nrexcl=bt.nrexcl,
        )

        atom_index_map = {}  # (tree_idx, rtp_atom_name) -> global atom index
        global_idx = 0

        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                print(f"WARNING: Sugar {ch}:{rs} not found in PDB", file=sys.stderr)
                continue

            # Map PDB name -> CHARMM name (with auto-detect for BGALNA etc.)
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                print(f"WARNING: No CHARMM RTP for {res.resname} ({ch}:{rs})",
                      file=sys.stderr)
                continue

            rtp_res = self.residues[rtp_name]
            pdb_atom_names = {a[0] for a in res.atoms}
            rtp_atom_names = [a[0] for a in rtp_res.atoms]

            # Determine which HO atoms to remove at linked positions
            # and redistribute their charge to the bonded O atom
            linked_os = link_atoms.get((ch, rs), set())
            remove_ho = {}  # ho_name -> o_name (for charge transfer)
            rtp_charges = {a[0]: a[2] for a in rtp_res.atoms}
            for o_name in linked_os:
                ho_name = 'HO' + o_name[1:]
                if ho_name in rtp_charges:
                    remove_ho[ho_name] = o_name

            # Build RTP->PDB mapping using sugar-specific atom map
            # Invert the PDB->CHARMM map to get CHARMM->PDB
            carb_rtp_to_pdb = {}
            pdb_resname = res.resname
            if pdb_resname in CARB_ATOM_MAP:
                for pdb_aname, charmm_aname in CARB_ATOM_MAP[pdb_resname].items():
                    carb_rtp_to_pdb[charmm_aname] = pdb_aname

            # Also add ARN-based renames
            for rtp_aname in rtp_atom_names:
                if rtp_aname in carb_rtp_to_pdb:
                    continue
                for (resn, ff_name), gmx_name in self.arn_reverse.items():
                    if resn == '*' and ff_name == rtp_aname:
                        carb_rtp_to_pdb[rtp_aname] = gmx_name
                        break

            rtp_to_pdb = _match_atom_names(rtp_atom_names, pdb_atom_names,
                                           carb_rtp_to_pdb)
            pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in res.atoms}

            # Build charge adjustments: O gets HO's charge when HO is removed
            charge_adjust = {}  # o_name -> extra charge
            # When glycosidic bond forms, O type changes from hydroxyl to ether
            type_change = {}  # o_name -> new_type
            for ho_name, o_name in remove_ho.items():
                pdb_o_name = rtp_to_pdb.get(o_name, o_name)
                if pdb_o_name not in pdb_coords and o_name in linked_os:
                    # O not in PDB (removed by CHARMM-GUI at linkage):
                    # redistribute O + HO combined charge to anomeric carbon
                    c_name = 'C2' if rtp_name in ('ANE5AC', 'BNE5AC') else 'C1'
                    charge_adjust[c_name] = charge_adjust.get(c_name, 0.0) + \
                        rtp_charges.get(o_name, 0.0) + rtp_charges[ho_name]
                else:
                    charge_adjust[o_name] = charge_adjust.get(o_name, 0.0) + rtp_charges[ho_name]
                    # OC311 (hydroxyl) -> OC3C61 (ether) for linked O
                    type_change[o_name] = 'OC3C61'

            for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
                # Skip HO atoms at linked positions
                if atom_name in remove_ho:
                    continue

                # Apply charge redistribution and type change for linked O atoms
                adj_charge = charge + charge_adjust.get(atom_name, 0.0)
                adj_type = type_change.get(atom_name, atom_type)

                pdb_name = rtp_to_pdb.get(atom_name, atom_name)
                # For carbs, skip atoms not in PDB: H atoms, linked O atoms,
                # and any O1/O2 not present (glycosidic bond sites where
                # CHARMM-GUI removes the bridging O)
                if pdb_name not in pdb_coords:
                    if (atom_name.startswith('H') or atom_name in linked_os
                            or atom_name in ('O1', 'O2')):
                        if self.verbose:
                            print(f"    Skipping {rtp_name}:{atom_name} "
                                  f"(not in PDB {ch}:{rs})")
                        continue

                global_idx += 1
                mass = self.atom_masses.get(adj_type, 0.0)
                x, y, z = pdb_coords.get(pdb_name, (0.0, 0.0, 0.0))
                chain_top.atoms.append(AtomEntry(
                    index=global_idx,
                    atom_type=adj_type,
                    resnr=tree_idx + 1,
                    resname=rtp_name,
                    atomname=pdb_name if pdb_name in pdb_coords else atom_name,
                    cgnr=cgnr,
                    charge=adj_charge,
                    mass=mass,
                    x=x, y=y, z=z,
                    chain_id=ch,
                    orig_resseq=rs,
                    orig_resname=res.resname,
                ))
                atom_index_map[(tree_idx, atom_name)] = global_idx

            # Intra-residue bonds from RTP (skip bonds involving removed atoms)
            for a1, a2 in rtp_res.bonds:
                if a1 in remove_ho or a2 in remove_ho:
                    continue
                idx1 = atom_index_map.get((tree_idx, a1))
                idx2 = atom_index_map.get((tree_idx, a2))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Add inter-residue glycosidic bonds
        tree_pos = {(ch, rs): i for i, (ch, rs) in enumerate(tree)}
        for don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom in glycan_links:
            don_key = (don_ch, don_rs)
            acc_key = (acc_ch, acc_rs)
            if don_key in tree_pos and acc_key in tree_pos:
                don_idx = tree_pos[don_key]
                acc_idx = tree_pos[acc_key]
                idx1 = atom_index_map.get((don_idx, don_atom))
                idx2 = atom_index_map.get((acc_idx, acc_atom))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Deduplicate bonds
        chain_top.bonds = sorted(set(chain_top.bonds))

        # Build bond graph and generate angles/dihedrals/pairs
        adj = defaultdict(set)
        for i, j in chain_top.bonds:
            adj[i].add(j)
            adj[j].add(i)

        # Generate angles — all use default ftype from [ bondedtypes ]
        chain_top.angles = list(self._generate_angles(adj))

        # For carbs, generate all dihedrals from connectivity
        generated = set()
        for j in sorted(adj.keys()):
            for k in sorted(adj[j]):
                if k <= j:
                    continue
                for i in sorted(adj[j]):
                    if i == k:
                        continue
                    for l in sorted(adj[k]):
                        if l == j:
                            continue
                        generated.add((i, j, k, l))
        chain_top.dihedrals = sorted(generated)

        chain_top.pairs = self._generate_pairs(chain_top.dihedrals)

        # Resolve impropers from RTP
        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                continue
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                continue
            rtp_res = self.residues[rtp_name]
            for imp in rtp_res.impropers:
                indices = []
                for ref in imp:
                    idx = atom_index_map.get((tree_idx, ref))
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    chain_top.impropers.append(tuple(indices))

        # Renumber
        self._renumber_atoms(chain_top)

        return chain_top

    def build_glycolipid_chain(self, ceramide_res, tree, link_atoms,
                               all_chains, glycan_links, ceramide_link):
        """Build topology for a glycolipid (ceramide + sugar tree as one moleculetype).

        ceramide_res: PDBResidue for the ceramide
        tree: list of (chain_id, resseq) for sugars in topological order
        link_atoms: dict (chain_id, resseq) -> set of linked O atoms
        all_chains: list of PDBChain objects
        glycan_links: list of sugar-sugar links
        ceramide_link: (cer_chain, cer_resseq, cer_atom, sugar_chain, sugar_resseq, sugar_atom)
        """
        # Build residue lookup
        res_lookup = {}
        for chain in all_chains:
            for res in chain.residues:
                res_lookup[(chain.chain_id, res.resseq)] = res

        bt = self.carb_bonded_types or self.bonded_types
        cer_ch = ceramide_link[0]
        cer_rs = ceramide_link[1]
        chain_name = f"Glycolipid_{cer_ch.strip() or 'X'}_{cer_rs}"

        chain_top = ChainTopology(
            name=chain_name,
            nrexcl=bt.nrexcl,
        )

        atom_index_map = {}  # (resnr, rtp_atom_name) -> global atom index
        global_idx = 0

        # --- Step 1: Build ceramide residue ---
        cer_rtp_name = PDB_TO_LIPID.get(ceramide_res.resname, ceramide_res.resname)
        if cer_rtp_name not in self.residues:
            print(f"WARNING: No RTP entry for ceramide {ceramide_res.resname} "
                  f"(tried {cer_rtp_name})", file=sys.stderr)
            return None

        rtp_res = self.residues[cer_rtp_name]
        pdb_atom_names = {a[0] for a in ceramide_res.atoms}
        pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in ceramide_res.atoms}
        rtp_charges = {a[0]: a[2] for a in rtp_res.atoms}

        # At linkage: remove ceramide HO1 and O1 (not in PDB — CHARMM-GUI
        # already removed them). Redistribute their combined charge to C1S.
        cer_remove_ho = {}
        cer_charge_adjust = {}
        cer_type_change = {}
        if 'HO1' in rtp_charges:
            cer_remove_ho['HO1'] = 'O1'
            # O1 not in PDB: redistribute O1+HO1 combined charge to C1S
            if 'O1' not in pdb_atom_names:
                o1_charge = rtp_charges.get('O1', 0.0)
                ho1_charge = rtp_charges['HO1']
                cer_charge_adjust['C1S'] = o1_charge + ho1_charge
            else:
                # O1 is in PDB: standard redistribution
                cer_charge_adjust['O1'] = rtp_charges['HO1']
                cer_type_change['O1'] = 'OC301'

        resnr = 1  # ceramide is residue 1
        for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
            if atom_name in cer_remove_ho:
                continue

            adj_charge = charge + cer_charge_adjust.get(atom_name, 0.0)
            adj_type = cer_type_change.get(atom_name, atom_type)

            # Skip atoms not in PDB (H atoms, and O1 if already removed)
            if atom_name not in pdb_coords:
                if atom_name.startswith('H') or atom_name in ('O1', 'HO1'):
                    if self.verbose:
                        print(f"    Skipping {cer_rtp_name}:{atom_name} (not in PDB)")
                    continue

            global_idx += 1
            mass = self.atom_masses.get(adj_type, 0.0)
            x, y, z = pdb_coords.get(atom_name, (0.0, 0.0, 0.0))
            chain_top.atoms.append(AtomEntry(
                index=global_idx,
                atom_type=adj_type,
                resnr=resnr,
                resname=cer_rtp_name,
                atomname=atom_name,
                cgnr=cgnr,
                charge=adj_charge,
                mass=mass,
                x=x, y=y, z=z,
                chain_id=cer_ch,
                orig_resseq=cer_rs,
                orig_resname=ceramide_res.resname,
            ))
            atom_index_map[(0, atom_name)] = global_idx

        # Intra-residue bonds for ceramide
        for a1, a2 in rtp_res.bonds:
            if a1 in cer_remove_ho or a2 in cer_remove_ho:
                continue
            idx1 = atom_index_map.get((0, a1))
            idx2 = atom_index_map.get((0, a2))
            if idx1 is not None and idx2 is not None:
                chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Ceramide impropers
        for imp in rtp_res.impropers:
            indices = []
            for ref in imp:
                idx = atom_index_map.get((0, ref))
                if idx is not None:
                    indices.append(idx)
            if len(indices) == 4:
                chain_top.impropers.append(tuple(indices))

        # --- Step 2: Build sugar tree (reusing glycan chain logic) ---
        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                print(f"WARNING: Sugar {ch}:{rs} not found in PDB", file=sys.stderr)
                continue

            # Map PDB name -> CHARMM name (with auto-detect for BGALNA etc.)
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                print(f"WARNING: No CHARMM RTP for {res.resname} ({ch}:{rs})",
                      file=sys.stderr)
                continue

            rtp_res = self.residues[rtp_name]
            pdb_atom_names_set = {a[0] for a in res.atoms}
            rtp_atom_names = [a[0] for a in rtp_res.atoms]

            # Determine linked O atoms
            linked_os = link_atoms.get((ch, rs), set())
            remove_ho = {}
            rtp_charges = {a[0]: a[2] for a in rtp_res.atoms}
            for o_name in linked_os:
                ho_name = 'HO' + o_name[1:]
                if ho_name in rtp_charges:
                    remove_ho[ho_name] = o_name

            # Build RTP->PDB atom name mapping
            carb_rtp_to_pdb = {}
            pdb_resname = res.resname
            if pdb_resname in CARB_ATOM_MAP:
                for pdb_aname, charmm_aname in CARB_ATOM_MAP[pdb_resname].items():
                    carb_rtp_to_pdb[charmm_aname] = pdb_aname

            for rtp_aname in rtp_atom_names:
                if rtp_aname in carb_rtp_to_pdb:
                    continue
                for (resn, ff_name), gmx_name in self.arn_reverse.items():
                    if resn == '*' and ff_name == rtp_aname:
                        carb_rtp_to_pdb[rtp_aname] = gmx_name
                        break

            rtp_to_pdb = _match_atom_names(rtp_atom_names, pdb_atom_names_set,
                                           carb_rtp_to_pdb)
            pdb_coords = {a[0]: (a[1], a[2], a[3]) for a in res.atoms}

            charge_adjust = {}
            type_change = {}
            # Determine which O atom links to ceramide (if any)
            cer_linked_o = ceramide_link[5] if (ch, rs) == (ceramide_link[3], ceramide_link[4]) else None
            for ho_name, o_name in remove_ho.items():
                pdb_o_name = rtp_to_pdb.get(o_name, o_name)
                if pdb_o_name not in pdb_coords and o_name in linked_os:
                    # O not in PDB (removed at linkage): redistribute to anomeric C
                    c_name = 'C2' if rtp_name in ('ANE5AC', 'BNE5AC') else 'C1'
                    charge_adjust[c_name] = charge_adjust.get(c_name, 0.0) + \
                        rtp_charges.get(o_name, 0.0) + rtp_charges[ho_name]
                else:
                    charge_adjust[o_name] = charge_adjust.get(o_name, 0.0) + rtp_charges[ho_name]
                    # Ceramide-linked O becomes OC301 (linear ether),
                    # sugar-sugar linked O becomes OC3C61 (cyclic ether)
                    if o_name == cer_linked_o:
                        type_change[o_name] = 'OC301'
                    else:
                        type_change[o_name] = 'OC3C61'

            resnr_sugar = tree_idx + 2  # ceramide is resnr 1
            for atom_name, atom_type, charge, cgnr in rtp_res.atoms:
                if atom_name in remove_ho:
                    continue

                adj_charge = charge + charge_adjust.get(atom_name, 0.0)
                adj_type = type_change.get(atom_name, atom_type)

                pdb_name = rtp_to_pdb.get(atom_name, atom_name)
                if pdb_name not in pdb_coords:
                    if (atom_name.startswith('H') or atom_name in linked_os
                            or atom_name in ('O1', 'O2')):
                        if self.verbose:
                            print(f"    Skipping {rtp_name}:{atom_name} "
                                  f"(not in PDB {ch}:{rs})")
                        continue

                global_idx += 1
                mass = self.atom_masses.get(adj_type, 0.0)
                x, y, z = pdb_coords.get(pdb_name, (0.0, 0.0, 0.0))
                chain_top.atoms.append(AtomEntry(
                    index=global_idx,
                    atom_type=adj_type,
                    resnr=resnr_sugar,
                    resname=rtp_name,
                    atomname=pdb_name if pdb_name in pdb_coords else atom_name,
                    cgnr=cgnr,
                    charge=adj_charge,
                    mass=mass,
                    x=x, y=y, z=z,
                    chain_id=ch,
                    orig_resseq=rs,
                    orig_resname=res.resname,
                ))
                atom_index_map[(tree_idx + 1, atom_name)] = global_idx

            # Intra-residue bonds
            for a1, a2 in rtp_res.bonds:
                if a1 in remove_ho or a2 in remove_ho:
                    continue
                idx1 = atom_index_map.get((tree_idx + 1, a1))
                idx2 = atom_index_map.get((tree_idx + 1, a2))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # --- Step 3: Add inter-residue bonds ---
        # Ceramide C1S — sugar O1 bond (sugar O1 bridges ceramide C1S to sugar C1)
        cer_atom = ceramide_link[2]  # e.g. 'C1S'
        sugar_atom = ceramide_link[5]  # e.g. 'O1'
        cer_idx = atom_index_map.get((0, cer_atom))
        root_sugar_idx = atom_index_map.get((1, sugar_atom))
        if cer_idx is not None and root_sugar_idx is not None:
            chain_top.bonds.append((min(cer_idx, root_sugar_idx),
                                    max(cer_idx, root_sugar_idx)))

        # Sugar-sugar glycosidic bonds
        tree_pos = {(ch, rs): i + 1 for i, (ch, rs) in enumerate(tree)}
        for don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom in glycan_links:
            don_key = (don_ch, don_rs)
            acc_key = (acc_ch, acc_rs)
            if don_key in tree_pos and acc_key in tree_pos:
                don_tidx = tree_pos[don_key]
                acc_tidx = tree_pos[acc_key]
                idx1 = atom_index_map.get((don_tidx, don_atom))
                idx2 = atom_index_map.get((acc_tidx, acc_atom))
                if idx1 is not None and idx2 is not None:
                    chain_top.bonds.append((min(idx1, idx2), max(idx1, idx2)))

        # Deduplicate bonds
        chain_top.bonds = sorted(set(chain_top.bonds))

        # --- Step 4: Build bond graph and enumerate angles/dihedrals/pairs ---
        adj = defaultdict(set)
        for i, j in chain_top.bonds:
            adj[i].add(j)
            adj[j].add(i)

        chain_top.angles = list(self._generate_angles(adj))

        generated = set()
        for j in sorted(adj.keys()):
            for k in sorted(adj[j]):
                if k <= j:
                    continue
                for i in sorted(adj[j]):
                    if i == k:
                        continue
                    for l in sorted(adj[k]):
                        if l == j:
                            continue
                        generated.add((i, j, k, l))
        chain_top.dihedrals = sorted(generated)

        chain_top.pairs = self._generate_pairs(chain_top.dihedrals)

        # Sugar impropers
        for tree_idx, (ch, rs) in enumerate(tree):
            res = res_lookup.get((ch, rs))
            if res is None:
                continue
            rtp_name = _resolve_sugar_rtp(res.resname, res.atoms, self.residues)
            if rtp_name is None:
                continue
            rtp_res = self.residues[rtp_name]
            for imp in rtp_res.impropers:
                indices = []
                for ref in imp:
                    idx = atom_index_map.get((tree_idx + 1, ref))
                    if idx is not None:
                        indices.append(idx)
                if len(indices) == 4:
                    chain_top.impropers.append(tuple(indices))

        # Renumber
        self._renumber_atoms(chain_top)

        return chain_top


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------
def _write_moleculetype(f, chain_top, bonded_types):
    """Write a chain moleculetype section to an open file handle."""
    bt = bonded_types

    # [ moleculetype ]
    f.write("[ moleculetype ]\n")
    f.write("; Name            nrexcl\n")
    f.write(f"{chain_top.name:<18s}{chain_top.nrexcl}\n\n")

    # [ atoms ]
    f.write("[ atoms ]\n")
    f.write(";   nr       type  resnr residue  atom   cgnr     charge       mass\n")
    for a in chain_top.atoms:
        f.write(f"{a.index:6d} {a.atom_type:>10s} {a.resnr:6d} {a.resname:>6s} "
                f"{a.atomname:>6s} {a.cgnr:6d} {a.charge:11.4f} {a.mass:11.4f}\n")
    f.write("\n")

    # [ bonds ]
    if chain_top.bonds:
        f.write("[ bonds ]\n")
        f.write(";  ai    aj funct\n")
        for i, j in chain_top.bonds:
            f.write(f"{i:5d} {j:5d} {bt.bond_type:5d}\n")
        f.write("\n")

    # [ pairs ]
    if chain_top.pairs:
        f.write("[ pairs ]\n")
        f.write(";  ai    aj funct\n")
        for i, j in chain_top.pairs:
            f.write(f"{i:5d} {j:5d}     1\n")
        f.write("\n")

    # [ angles ]
    if chain_top.angles:
        f.write("[ angles ]\n")
        f.write(";  ai    aj    ak funct\n")
        for angle in chain_top.angles:
            if len(angle) == 4:
                i, j, k, ftype = angle
                f.write(f"{i:5d} {j:5d} {k:5d} {ftype:5d}\n")
            else:
                i, j, k = angle
                f.write(f"{i:5d} {j:5d} {k:5d} {bt.angle_type:5d}\n")
        f.write("\n")

    # [ dihedrals ] — proper
    if chain_top.dihedrals:
        f.write("[ dihedrals ]\n")
        f.write(";  ai    aj    ak    al funct\n")
        for dih in chain_top.dihedrals:
            if len(dih) == 5:
                # Named dihedral type (AMBER ILDN)
                i, j, k, l, dtype = dih
                f.write(f"{i:5d} {j:5d} {k:5d} {l:5d} {bt.dihedral_type:5d}  ; {dtype}\n")
            else:
                i, j, k, l = dih
                f.write(f"{i:5d} {j:5d} {k:5d} {l:5d} {bt.dihedral_type:5d}\n")
        f.write("\n")

    # [ dihedrals ] — improper
    if chain_top.impropers:
        f.write("[ dihedrals ] ; impropers\n")
        f.write(";  ai    aj    ak    al funct\n")
        for i, j, k, l in chain_top.impropers:
            f.write(f"{i:5d} {j:5d} {k:5d} {l:5d} {bt.improper_type:5d}\n")
        f.write("\n")

    # [ cmap ]
    if chain_top.cmap:
        f.write("[ cmap ]\n")
        f.write(";  ai    aj    ak    al    am funct\n")
        for i, j, k, l, m in chain_top.cmap:
            f.write(f"{i:5d} {j:5d} {k:5d} {l:5d} {m:5d}     1\n")
        f.write("\n")

    # Position restraints include (stays as separate file with #ifdef)
    posre_name = f"posre_{chain_top.name}.itp"
    f.write("; Include Position restraint file\n")
    f.write("#ifdef POSRES\n")
    f.write(f'#include "{posre_name}"\n')
    f.write("#endif\n\n")


def write_pdb(chain_tops, path, extra_pdb_lines=None, cryst1=None):
    """Write a PDB file with atom names matching the topology.

    extra_pdb_lines: list of raw PDB ATOM/HETATM lines to append
    (for ions/BUF particles not built as chain topologies).
    """
    serial = 0
    with open(path, 'w') as f:
        f.write("REMARK    Generated by dvbfixer top\n")
        if cryst1:
            f.write(cryst1 if cryst1.endswith('\n') else cryst1 + '\n')
        for ct in chain_tops:
            for atom in ct.atoms:
                serial += 1
                # PDB columns: 13-16 atom name, 17 altLoc, 18-20 resName
                name = atom.atomname
                if len(name) < 4:
                    name_field = f" {name:<3s}"
                else:
                    name_field = f"{name:<4s}"

                # Determine element from atom name
                elem = name[0] if name[0].isalpha() else name[1]

                # Use orig_resname (PDB name, includes protonation renames)
                # atom.resname is the RTP name (e.g. BGLCNA for glycans)
                resname = atom.orig_resname
                # PDB cols: 17=altLoc, 18-20=resName, 21=space, 22=chainID
                # 4-char resnames (ASPP, GLUP): cols 18-21, no space before chainID
                if len(resname) <= 3:
                    res_chain = f" {resname:>3s} {atom.chain_id}"  # " ASP A"
                else:
                    res_chain = f" {resname:<4s}{atom.chain_id}"   # " ASPP A" (no gap)

                f.write(
                    f"ATOM  {serial % 100000:5d} {name_field}"
                    f"{res_chain}"
                    f"{atom.orig_resseq:4d}"
                    f"    "
                    f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
                    f"  1.00  0.00"
                    f"          "
                    f"{elem:>2s}\n"
                )
            # TER record after each chain
            if ct.atoms:
                last = ct.atoms[-1]
                serial += 1
                resname = last.orig_resname
                if len(resname) <= 3:
                    res_chain = f" {resname:>3s} {last.chain_id}"
                else:
                    res_chain = f" {resname:<4s}{last.chain_id}"
                f.write(
                    f"TER   {serial % 100000:5d}      "
                    f"{res_chain}"
                    f"{last.orig_resseq:4d}\n"
                )
        # Append extra molecules (ions/water/small molecules) with renumbered serials
        if extra_pdb_lines:
            resseq = 0
            prev_resid = None
            for line in extra_pdb_lines:
                serial += 1
                # Track residue changes for resseq renumbering
                orig_resid = line[17:26]  # resname + chain + resseq
                if orig_resid != prev_resid:
                    resseq += 1
                    prev_resid = orig_resid
                # Rewrite serial (cols 7-11), clear chain ID (col 22),
                # and renumber resseq (cols 23-26)
                f.write(f"{line[:6]}{serial % 100000:5d}{line[11:21]} "
                        f"{resseq % 10000:4d}{line[26:]}")
        f.write("END\n")


def write_posre(chain_top, path, fc=1000.0):
    """Write position restraint file for heavy atoms."""
    with open(path, 'w') as f:
        f.write("; Position restraints for heavy atoms\n")
        f.write("; Generated by dvbfixer top\n\n")
        f.write("[ position_restraints ]\n")
        f.write(";  ai  funct    fcx      fcy      fcz\n")
        for atom in chain_top.atoms:
            # Restrain heavy atoms only (not hydrogen)
            if not atom.atomname.startswith('H') and atom.atomname not in ('HN',):
                f.write(f"{atom.index:5d}     1  {fc:.1f}  {fc:.1f}  {fc:.1f}\n")


def _read_ff_content(path, keep_posres_ifdef=False):
    """Read a force field file, stripping preprocessor directives.

    Removes #include, #define, #ifdef, #ifndef, #else, #endif lines
    so the content can be inlined directly into the .top file.

    If keep_posres_ifdef=True, preserves #ifdef POSRES_* blocks
    (used for ions.itp which may have position restraint sections).
    """
    lines = []
    in_posres_ifdef = False
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if keep_posres_ifdef:
                if stripped.startswith('#ifdef POSRES'):
                    in_posres_ifdef = True
                    lines.append(line)
                    continue
                if in_posres_ifdef and stripped == '#endif':
                    in_posres_ifdef = False
                    lines.append(line)
                    continue
                if in_posres_ifdef:
                    lines.append(line)
                    continue
            if stripped.startswith(('#include', '#define', '#ifdef', '#ifndef',
                                    '#else', '#endif', '#error')):
                continue
            lines.append(line)
    return ''.join(lines)


def _parse_defaults(ff_dir):
    """Extract [ defaults ] section from forcefield.itp."""
    ff_itp = ff_dir / 'forcefield.itp'
    in_defaults = False
    defaults_lines = []
    with open(ff_itp) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('[ defaults ]'):
                in_defaults = True
                defaults_lines.append(line)
                continue
            if in_defaults:
                if stripped.startswith('[') or stripped.startswith('#include'):
                    break
                defaults_lines.append(line)
    return ''.join(defaults_lines)


def _dedup_atomtypes(content):
    """Remove duplicate atomtype entries, keeping the last definition.

    CHARMM ffnonbonded.itp has duplicate entries (e.g. HT, OT) from
    #ifdef HEAVY_H blocks — one with heavy mass, one with real mass.
    After stripping preprocessor directives, both remain. We keep the
    last definition for each atomtype name.
    """
    lines = content.split('\n')
    in_atomtypes = False
    atomtype_lines = []   # (line_idx, atomtype_name, line)
    other_lines = []      # (line_idx, line)

    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('[ atomtypes ]'):
            in_atomtypes = True
            other_lines.append((idx, line))
            continue
        if in_atomtypes and stripped.startswith('['):
            in_atomtypes = False
            other_lines.append((idx, line))
            continue
        if in_atomtypes and stripped and not stripped.startswith(';'):
            parts = stripped.split()
            if parts:
                atomtype_lines.append((idx, parts[0], line))
            continue
        other_lines.append((idx, line))

    # Keep only last occurrence of each atomtype
    seen = {}
    for idx, name, line in atomtype_lines:
        seen[name] = (idx, line)

    # Reconstruct
    all_entries = other_lines + [(idx, line) for idx, line in seen.values()]
    all_entries.sort(key=lambda x: x[0])
    return '\n'.join(line for _, line in all_entries)


# Extra bonded parameters for CHARMM glycosidic linkage sites.
# When OC311->OC3C61 at linkages, some atom type combos (CC321-OC3C61,
# CC3161-OC3C61-CC3162, etc.) are not in the standard FF distribution.
# Parameters by analogy with existing CC321D/CC321C/CC311D variants.
_GLYCAN_LINKAGE_PARAMS = """\
; ======================================================================
; Extra parameters for glycosidic linkage sites (by analogy)
; ======================================================================

; --- Extra bondtypes ---
[ bondtypes ]
; i       j     func    b0          kb
  CC321   OC3C61     1   0.14150000    301248.00 ; from CC321D OC3C61

; --- Extra angletypes ---
[ angletypes ]
; i       j       k     func    theta0      ktheta      rub         kub
  HCA2    CC321   OC3C61     5   109.500000   376.560000   0.00000000         0.00 ; from HCA2 CC321D OC3C61
  CC3163  CC321   OC3C61     5   111.500000   376.560000   0.00000000         0.00 ; from OC3C61 CC321D CC311C
  CC3161  OC3C61  CC3162     5   109.700000   794.960000   0.00000000         0.00 ; from CC3163 OC3C61 CC3162
  CC3162  CC3161  OC3C61     5   106.000000   376.560000   0.00000000         0.00 ; from CC3161 CC3162 OC3C61
  CC321   OC3C61  CC3162     5   109.700000   794.960000   0.00000000         0.00 ; from CC321D OC3C61 CC321C
  OC3C61  CC3162  OC3C61     5   112.000000   753.120000   0.00000000         0.00 ; from OC301 CC3162 OC3C61

; --- Extra dihedraltypes ---
[ dihedraltypes ]
; i       j       k       l     func    phi0        kphi        mult
; C-C-O-C glycosidic torsions (from CC3161 CC3162 OC3C61 CC3163)
  CC3161  CC3163  CC321   OC3C61     9     0.000000     0.836800     3 ; from par27 X CT1 CT2 X
  CC3161  OC3C61  CC3162  HCA1       9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  CC3161  OC3C61  CC3162  OC3C61     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  CC3162  CC3161  CC3161  OC3C61     9   180.000000     1.297040     3 ; from CC3161 CC3161 CC3162 OC3C61
  CC3162  CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  CC3163  CC321   OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC321 CC3163 OC3C61 CC3162
  CC321   CC3163  CC3161  OC3C61     9     0.000000     0.836800     3 ; from par27 X CT1 CT2 X
  CC321   OC3C61  CC3162  CC3161     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3162
  CC321   OC3C61  CC3162  HCA1       9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  CC321   OC3C61  CC3162  OC3C61     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  HCA1    CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  HCA1    CC3162  CC3161  OC3C61     9     0.000000     0.836800     3 ; from HCA1 CC3161 CC3162 OC3C61
  HCA1    CC3163  CC321   OC3C61     9     0.000000     0.836800     3 ; from HCA2 CC321 CC3163 OC3C61
  HCA2    CC321   OC3C61  CC3162     9     0.000000     0.836800     3 ; from par27 X CT2 OC30A X
  OC3C61  CC3161  CC3161  CC3163     9   180.000000     1.297040     3 ; from CC3161 CC3161 CC3162 OC3C61
  OC3C61  CC3162  CC3161  OC3C61     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  OC3C61  CC3162  OC3C61  CC3163     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  OC3C61  CC3163  CC321   OC3C61     9     0.000000     0.836800     3 ; from par27 X CT2 CT1 X
; Additional C-C-O-C and C-O-C-C linkage torsions
  CC3161  CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC3161 CC3162 OC3C61 CC3163
  CC3161  OC3C61  CC3162  CC3161     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3162
  CC3163  CC3161  OC3C61  CC3162     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3162

; ======================================================================
; Extra parameters for sialic acid (C2) glycosidic linkage sites
; ANE5AC C2 (type CC3062) links to parent sugar O (type OC3C61)
; ======================================================================

; --- Extra angletypes (sialic acid linkage) ---
[ angletypes ]
; i       j       k     func    theta0      ktheta      rub         kub
  CC3161  OC3C61  CC3062     5   109.700000   794.960000   0.00000000         0.00 ; from CC3163 OC3C61 CC3062 (sugar C-O-sialic C2)
  OC3C61  CC3062  OC3C61     5   112.000000   753.120000   0.00000000         0.00 ; from OC301 CC3062 OC3C61 (ring O - C2 - linked O)

; --- Extra dihedraltypes (sialic acid linkage) ---
[ dihedraltypes ]
; i       j       k       l     func    phi0        kphi        mult
; Torsions around CC3161-OC3C61-CC3062 sialic acid linkage
  CC3161  CC3161  OC3C61  CC3062     9     0.000000     0.836800     3 ; from CC3161 CC3163 OC3C61 CC3062
  CC3161  OC3C61  CC3062  CC2O2      9     0.000000     0.836800     3 ; from CC3163 OC3C61 CC3062 CC2O2
  CC3161  OC3C61  CC3062  OC3C61     9     0.000000     0.836800     3 ; from CC3163 OC3C61 CC3062 OC3C61
  CC3161  OC3C61  CC3062  CC3261     9     0.000000     0.836800     3 ; from CC3163 OC3C61 CC3062 CC3261
  HCA1    CC3161  OC3C61  CC3062     9     0.000000     1.188256     3 ; from HCA1 CC3163 OC3C61 CC3062
  OC3C61  CC3161  CC3161  OC3C61     9     0.000000     0.836800     3 ; from OC3C61 CC3162 CC3161 OC3C61
  OC3C61  CC3062  OC3C61  CC3163     9     0.000000     0.836800     3 ; from OC3C61 CC3162 OC3C61 CC3163
  OC3C61  CC3161  OC3C61  CC3062     9     0.000000     0.836800     3 ; from OC3C61 CC3161 OC3C61 CC3162
  NC2D1   CC3161  CC3161  OC3C61     9     0.000000     0.836800     3 ; from OC3C61 CC3261 CC3161 NC2D1

"""


def write_top(chain_tops, path, ff_dir, ff_type, bonded_types_list,
              water_model='tip3p', system_name='Protein',
              has_interchain_ss=False, extra_molecules=None,
              small_mol_itps=None):
    """Write topology: FF params in ffparams.itp, rest in topol.top.

    extra_molecules: list of (name, count) for ions/BUF/small molecules/water.
    small_mol_itps: list of small molecule .itp filenames to include.
    """
    ff_dir = Path(ff_dir)
    ff_name = ff_dir.name
    out_dir = Path(path).parent

    water_files = {
        'tip3p': 'tip3p.itp',
        'spc': 'spc.itp',
        'spce': 'spce.itp',
        'tip4p': 'tip4p.itp',
    }

    # --- Write ffparams.itp with all FF parameters ---
    ffparams_path = out_dir / 'ffparams.itp'
    with open(ffparams_path, 'w') as f:
        f.write("; Force field parameters — generated by dvbfixer top\n")
        f.write(f"; Force field: {ff_name}\n\n")

        # [ defaults ]
        f.write(_parse_defaults(ff_dir))
        f.write("\n")

        # [ atomtypes ] from ffnonbonded.itp (deduplicated)
        f.write("; Non-bonded parameters (from ffnonbonded.itp)\n")
        nb_content = _read_ff_content(ff_dir / 'ffnonbonded.itp')
        f.write(_dedup_atomtypes(nb_content))
        f.write("\n")

        # [ bondtypes ], [ angletypes ], [ dihedraltypes ] from ffbonded.itp
        f.write("; Bonded parameters (from ffbonded.itp)\n")
        f.write(_read_ff_content(ff_dir / 'ffbonded.itp'))
        f.write("\n")

        # CHARMM extras: cmap.itp, nbfix.itp
        if ff_type == 'charmm':
            cmap_path = ff_dir / 'cmap.itp'
            if cmap_path.exists():
                f.write("; CMAP parameters (from cmap.itp)\n")
                f.write(_read_ff_content(cmap_path))
                f.write("\n")

            nbfix_path = ff_dir / 'nbfix.itp'
            if nbfix_path.exists():
                f.write("; NBFIX parameters (from nbfix.itp)\n")
                f.write(_read_ff_content(nbfix_path))
                f.write("\n")

            # Extra parameters for glycosidic linkage sites.
            # At linkage sites OC311->OC3C61, creating atom type combos
            # not in the standard FF. Parameters by analogy with existing
            # glycosidic linkage params (CC321D/CC311D/CC321C variants).
            f.write(_GLYCAN_LINKAGE_PARAMS)

    # --- Write topol.top ---
    with open(path, 'w') as f:
        f.write("; Generated by dvbfixer top\n")
        f.write(f"; Force field: {ff_name}\n")
        f.write(f"; Water model: {water_model}\n")
        doc_path = ff_dir / 'forcefield.doc'
        if doc_path.exists():
            with open(doc_path) as doc:
                for doc_line in doc:
                    doc_line = doc_line.strip()
                    if doc_line and not all(c == '*' for c in doc_line):
                        f.write(f"; {doc_line}\n")
                        break
        f.write("\n")

        # Include FF parameters
        f.write('#include "ffparams.itp"\n\n')

        # Chain moleculetypes (separate .itp files)
        for ct, bt in zip(chain_tops, bonded_types_list):
            chain_itp = out_dir / f"{ct.name}.itp"
            with open(chain_itp, 'w') as cf:
                cf.write(f"; Moleculetype: {ct.name}\n")
                cf.write(f"; Generated by dvbfixer top\n\n")
                _write_moleculetype(cf, ct, bt)
            f.write(f'#include "{ct.name}.itp"\n')

        # Small molecule .itp includes (ACET, ACEH, etc.)
        if small_mol_itps:
            for mol_name in small_mol_itps:
                f.write(f'#include "{mol_name}.itp"\n')

        f.write("\n")

        # Water moleculetype (separate .itp)
        water_itp = water_files.get(water_model, 'tip3p.itp')
        water_src = ff_dir / water_itp
        if water_src.exists():
            water_out = out_dir / 'water.itp'
            with open(water_out, 'w') as wf:
                wf.write("; Water topology\n")
                wf.write("; Generated by dvbfixer top\n\n")
                _write_water_topology(wf, water_src)
            f.write('#include "water.itp"\n')

        # Ion moleculetypes (separate .itp)
        ions_path = ff_dir / 'ions.itp'
        if ions_path.exists():
            ions_out = out_dir / 'ions.itp'
            with open(ions_out, 'w') as ionf:
                ionf.write("; Ion topology\n")
                ionf.write("; Generated by dvbfixer top\n\n")
                ionf.write(_read_ff_content(ions_path, keep_posres_ifdef=True))
            f.write('#include "ions.itp"\n')
        f.write("\n")

        # [ system ]
        f.write("[ system ]\n")
        f.write(f"; Name\n")
        f.write(f"{system_name}\n\n")

        # [ molecules ]
        f.write("[ molecules ]\n")
        f.write("; Compound        #mols\n")
        for ct in chain_tops:
            f.write(f"{ct.name:<18s}1\n")
        if extra_molecules:
            for mol_name, mol_count in extra_molecules:
                f.write(f"{mol_name:<18s}{mol_count}\n")

        # Inter-chain SS bonds (must come after [ molecules ])
        if has_interchain_ss:
            f.write('\n; WARNING: The interchain_ss.itp include MUST remain at the end of this file,\n')
            f.write('; after [ molecules ] and after any SOL/ion entries added by gmx solvate/genion.\n')
            f.write('; If you add solvent, move this line below the SOL and ion molecule entries.\n')
            f.write('#include "interchain_ss.itp"\n')


def _write_water_topology(f, water_path):
    """Write water moleculetype from .itp file, using rigid (settles) version.

    Water .itp files have #ifndef FLEXIBLE / #else blocks. We extract only
    the rigid (settles) section since that's the standard for MD.
    """
    lines = []
    in_flexible = False
    skip_block = False

    with open(water_path) as wf:
        for line in wf:
            stripped = line.strip()

            # Track preprocessor blocks
            if stripped == '#ifndef FLEXIBLE' or stripped == '#ifdef FLEXIBLE':
                # #ifndef FLEXIBLE → rigid section follows (keep it)
                # #ifdef FLEXIBLE → flexible section follows (skip it)
                in_flexible = stripped == '#ifdef FLEXIBLE'
                skip_block = in_flexible
                continue
            elif stripped == '#else':
                # Toggle: if we were in rigid, now flexible (skip); vice versa
                skip_block = not skip_block
                continue
            elif stripped == '#endif':
                skip_block = False
                in_flexible = False
                continue
            elif stripped.startswith('#'):
                continue

            if not skip_block:
                lines.append(line)

    f.write(''.join(lines))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='dvbfixer top',
        description='Generate GROMACS topology files from PDB',
    )
    parser.add_argument('input', help='Input PDB file')
    parser.add_argument('-o', '--output', help='Output .top file (default: topol.top)')
    parser.add_argument('--ff', choices=['amber', 'charmm'], default='amber',
                        help='Force field (default: amber)')
    parser.add_argument('--ff-dir', help='Custom force field directory')
    parser.add_argument('--water', default='tip3p',
                        choices=['tip3p', 'spc', 'spce', 'tip4p'],
                        help='Water model (default: tip3p)')
    parser.add_argument('--ignh', action='store_true',
                        help='Ignore hydrogens in input PDB')
    parser.add_argument('--ss', action='append', default=[],
                        help='Disulfide bond: CHAIN1:NUM1:CHAIN2:NUM2 (repeatable)')
    parser.add_argument('--his', action='append', default=[],
                        help='HIS protonation: CHAIN:NUM:STATE (HIE/HID/HIP, repeatable)')
    parser.add_argument('--protonate', default=None,
                        help='Protonate residues. "all" protonates every ASP->ASPP, '
                             'GLU->GLUP, HIS->HSP. Comma-separated list protonates '
                             'specific residues: CHAIN:NUM[:STATE],... '
                             '(e.g. --protonate all, --protonate H:66,K:50:GLUP).')
    parser.add_argument('--merge', action='store_true',
                        help='Merge all chains into single moleculetype')
    parser.add_argument('--pdb', help='Output PDB file with topology-matched atom names')
    parser.add_argument('--acpype', action='store_true',
                        help='Use ACPYPE pipeline (AMBER14+GLYCAM -> ParmEd -> GROMACS). '
                             'Handles mixed 1-4 scaling via [ pairs_nb ].')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)

    if not input_path.exists():
        print(f"Error: {input_path} not found", file=sys.stderr)
        sys.exit(1)

    # Convert GRO to temp PDB if needed
    tmp_pdb = None
    orig_input_path = input_path  # preserve for output path defaults
    if input_path.suffix.lower() == '.gro':
        print(f"Converting GRO to PDB via MDAnalysis...")
        tmp_pdb = _gro_to_pdb(input_path)
        input_path = tmp_pdb

    # ACPYPE mode: OpenMM + ParmEd + ACPYPE pipeline
    if args.acpype:
        from dvbfixer.acpype_export import export_gromacs

        if args.ff == 'charmm':
            print("WARNING: --acpype always uses AMBER14+GLYCAM, ignoring --ff charmm",
                  file=sys.stderr)

        # Parse --ss flags into extra_ss set
        extra_ss = set()
        for ss_spec in args.ss:
            parts = ss_spec.split(':')
            if len(parts) == 4:
                extra_ss.add((parts[0], int(parts[1])))
                extra_ss.add((parts[2], int(parts[3])))

        if args.output:
            out_dir = Path(args.output).parent or Path('.')
            basename = Path(args.output).stem
        else:
            out_dir = input_path.parent or Path('.')
            basename = input_path.stem

        export_gromacs(input_path, out_dir, basename=basename,
                       extra_ss=extra_ss or None, verbose=args.verbose)
        return

    # Determine FF directory
    if args.ff_dir:
        ff_dir = Path(args.ff_dir)
    else:
        ff_name = FF_CHOICES[args.ff]
        ff_dir = FF_DIR / ff_name
        if not ff_dir.exists():
            print(f"Error: Force field directory not found: {ff_dir}", file=sys.stderr)
            print(f"Use --ff-dir to specify the path", file=sys.stderr)
            sys.exit(1)

    ff_name = ff_dir.name
    print(f"Using force field: {ff_name}")

    # Read PDB
    chains = read_pdb_chains(input_path)
    if not chains:
        print("Error: No chains found in PDB", file=sys.stderr)
        sys.exit(1)

    # Build topology builder first (need its residue dict for chain filtering)
    builder = TopologyBuilder(ff_dir, args.ff, args.verbose)

    # Strip water and ion residues (handled separately via counting)
    ion_names_for_filter = set()
    ions_path_check = ff_dir / 'ions.itp'
    if ions_path_check.exists():
        ion_names_for_filter = _parse_ion_names(ions_path_check)
    skip_resnames = _WATER_RESNAMES | ion_names_for_filter

    for chain in chains:
        chain.residues = [r for r in chain.residues
                          if r.resname not in skip_resnames]

    # Remove empty chains (were all water/ions)
    chains = [c for c in chains if c.residues]

    # Split chains by nearest-atom distance (detect separate molecules in
    # GROMACS PDB output where multiple molecules share a chain ID)
    split_chains_list = []
    for chain in chains:
        split_chains_list.extend(_split_chain_by_distance(chain, gap_cutoff=4.0))
    chains = split_chains_list

    if chains:
        resnames = set()
        for c in chains:
            for r in c.residues:
                resnames.add(r.resname)
        print(f"Found {len(chains)} chain(s), residue types: {', '.join(sorted(resnames))}")

    protein_chains = []
    other_chains = []
    small_mol_counts = {}  # resname -> count (for single-residue small molecules)
    sugar_names = set(PDB_TO_CARB.keys())
    charmm_sugar_names_filter = set(PDB_TO_CARB.values())
    for chain in chains:
        has_protein = any(
            r.resname in STANDARD_AA or r.resname in PDB_TO_GMX
            for r in chain.residues
        )
        if has_protein:
            # Keep only protein residues in protein chains
            chain.residues = [
                r for r in chain.residues
                if r.resname in STANDARD_AA or r.resname in PDB_TO_GMX
            ]
            protein_chains.append(chain)
        elif len(chain.residues) == 1 and chain.residues[0].resname in builder.residues:
            # Single-residue molecule (detected by distance split) — count it
            rn = chain.residues[0].resname
            small_mol_counts[rn] = small_mol_counts.get(rn, 0) + 1
        else:
            # Check if chain has any FF-recognized residues
            has_known = any(
                r.resname in builder.residues or r.resname in sugar_names
                or _is_ceramide(r.resname)
                for r in chain.residues
            )
            if has_known:
                other_chains.append(chain)

    if not protein_chains and not other_chains and not small_mol_counts:
        print("Error: No recognized chains found", file=sys.stderr)
        sys.exit(1)

    # Detect SS bonds
    ss_bonds = read_ssbonds(input_path)

    # Parse explicit --ss flags
    for ss_spec in args.ss:
        parts = ss_spec.split(':')
        if len(parts) == 4:
            ss_bonds.append((parts[0], int(parts[1]), parts[2], int(parts[3])))

    # Build per-chain SS residue sets and intra-chain SS pairs
    chain_ss = defaultdict(set)
    intrachain_ss = defaultdict(list)  # chain_id -> [(resseq1, resseq2)]
    for ch1, res1, ch2, res2 in ss_bonds:
        chain_ss[ch1].add(res1)
        chain_ss[ch2].add(res2)
        if ch1 == ch2:
            intrachain_ss[ch1].append((res1, res2))

    if ss_bonds and args.verbose:
        print(f"Disulfide bonds: {len(ss_bonds)}")
        for ch1, res1, ch2, res2 in ss_bonds:
            print(f"  {ch1}:{res1} - {ch2}:{res2}")

    # Apply protonation overrides (--protonate and --his)
    # --protonate without args: all ASP->ASPP, GLU->GLUP, HIS->HSP
    # --protonate with args: specific CHAIN:NUM[:STATE] overrides
    # Default protonated forms per FF
    if args.ff == 'charmm':
        _PROT_DEFAULTS = {'ASP': 'ASPP', 'GLU': 'GLUP', 'HIS': 'HSP',
                          'HIE': 'HSP', 'HID': 'HSP', 'HSE': 'HSP', 'HSD': 'HSP'}
    else:
        _PROT_DEFAULTS = {'ASP': 'ASH', 'GLU': 'GLH', 'HIS': 'HIP',
                          'HIE': 'HIP', 'HID': 'HIP', 'HSE': 'HIP', 'HSD': 'HIP'}

    protonate_all = args.protonate == 'all'
    prot_overrides = {}
    if args.protonate and args.protonate != 'all':
        # Validate: must contain ':' (CHAIN:NUM format), not a filename
        if ':' not in args.protonate:
            print(f"ERROR: Invalid --protonate value '{args.protonate}'. "
                  f"Use --protonate (all) or --protonate CHAIN:NUM[:STATE],... "
                  f"Note: place --protonate after the input file.",
                  file=sys.stderr)
            sys.exit(1)
        for spec in args.protonate.split(','):
            parts = spec.split(':')
            if len(parts) == 3:
                prot_overrides[(parts[0], int(parts[1]))] = parts[2]
            elif len(parts) == 2:
                # CHAIN:NUM without STATE — use default protonated form
                prot_overrides[(parts[0], int(parts[1]))] = None  # resolve later

    his_overrides = {}
    for his_spec in args.his:
        parts = his_spec.split(':')
        if len(parts) == 3:
            his_overrides[(parts[0], int(parts[1]))] = parts[2]

    # Apply overrides or auto-detect from PDB atoms
    for chain in protein_chains:
        for res in chain.residues:
            key = (res.chain_id, res.resseq)
            if key in prot_overrides:
                state = prot_overrides[key]
                if state is None:
                    state = _PROT_DEFAULTS.get(res.resname)
                if state:
                    old_name = res.resname
                    res.resname = state
                    if args.verbose:
                        print(f"  {old_name} {res.chain_id}:{res.resseq} -> "
                              f"{res.resname} (--protonate)")
            elif protonate_all and res.resname in _PROT_DEFAULTS:
                old_name = res.resname
                res.resname = _PROT_DEFAULTS[old_name]
                if args.verbose:
                    print(f"  {old_name} {res.chain_id}:{res.resseq} -> "
                          f"{res.resname} (--protonate all)")
            elif key in his_overrides:
                res.resname = his_overrides[key]
            elif res.resname == 'HIS':
                # Auto-detect protonation from H atoms present
                atom_names = {a[0] for a in res.atoms}
                has_hd1 = 'HD1' in atom_names
                has_he2 = 'HE2' in atom_names
                if has_hd1 and has_he2:
                    res.resname = 'HIP'  # doubly protonated
                elif has_hd1:
                    res.resname = 'HID'  # delta protonated
                else:
                    res.resname = 'HIE'  # epsilon protonated (default)
                if args.verbose and res.resname != 'HIE':
                    print(f"  HIS {res.chain_id}:{res.resseq} -> {res.resname} "
                          f"(auto-detected from H atoms)")

    # Add missing protonation H atoms with proper geometry via PDBFixer
    if args.protonate or args.his:
        _add_protonation_hydrogens(protein_chains, args.input, args.ff,
                                   verbose=args.verbose)

    # Ignore hydrogens if requested
    if args.ignh:
        for chain in protein_chains:
            for res in chain.residues:
                res.atoms = [(n, x, y, z) for n, x, y, z in res.atoms
                             if not n.startswith('H') and n not in ('1H', '2H', '3H')]

    # Build topologies (builder already created above for chain filtering)
    chain_tops = []

    for chain in protein_chains:
        if args.verbose:
            print(f"\nBuilding topology for chain {chain.chain_id} "
                  f"({len(chain.residues)} residues)")

        ct = builder.build_chain(chain, chain_ss.get(chain.chain_id, set()),
                                 intrachain_ss.get(chain.chain_id, []))
        if ct is not None:
            chain_tops.append(ct)
            n_bonds = len(ct.bonds)
            n_angles = len(ct.angles)
            n_dihedrals = len(ct.dihedrals)
            print(f"Chain {chain.chain_id}: {len(ct.atoms)} atoms, {n_bonds} bonds, "
                  f"{n_angles} angles, {n_dihedrals} dihedrals")
        else:
            print(f"WARNING: Failed to build topology for chain {chain.chain_id}",
                  file=sys.stderr)

    # Detect glycan/glycolipid links early to know which residues are handled
    glycan_links = detect_glycan_links(chains)
    glycolipid_ceramide_resseqs = set()  # (chain_id, resseq) of ceramides in glycolipids
    glycolipid_sugar_resseqs = set()     # sugars in glycolipid trees

    if glycan_links:
        trees = build_glycan_trees(glycan_links, chains)
        for tree, link_atoms, prot_links, cer_links in trees:
            if cer_links:
                for cl in cer_links:
                    glycolipid_ceramide_resseqs.add((cl[0], cl[1]))
                for ch, rs in tree:
                    glycolipid_sugar_resseqs.add((ch, rs))

    charmm_sugar_names = set(PDB_TO_CARB.values())

    # Build non-protein chains (no terminal patches, no SS bonds)
    # Skip chains that are entirely sugars or glycolipids — handled later
    for chain in other_chains:
        all_sugar_or_lipid = all(
            r.resname in sugar_names or r.resname in charmm_sugar_names
            or (chain.chain_id, r.resseq) in glycolipid_ceramide_resseqs
            for r in chain.residues
        )
        if all_sugar_or_lipid:
            if args.verbose:
                print(f"\nSkipping chain {chain.chain_id} (sugar/glycolipid residues, "
                      f"handled by glycan/glycolipid detection)")
            continue

        if args.verbose:
            print(f"\nBuilding topology for non-protein chain {chain.chain_id} "
                  f"({len(chain.residues)} residues)")

        ct = builder.build_chain(chain)
        if ct is not None:
            # Rename from default Protein_ to Other_
            ct.name = ct.name.replace('Protein_', 'Other_')
            chain_tops.append(ct)
            print(f"Chain {chain.chain_id}: {len(ct.atoms)} atoms, "
                  f"{len(ct.bonds)} bonds (non-protein)")
        else:
            print(f"WARNING: Failed to build topology for chain {chain.chain_id}",
                  file=sys.stderr)

    if not chain_tops and not glycan_links and not small_mol_counts:
        print("Error: No topologies built", file=sys.stderr)
        sys.exit(1)

    # Detect and build glycan/glycolipid chains
    glycan_tops = []
    protein_sugar_links = []  # for intermolecular interactions
    if glycan_links:
        # Build residue lookup for glycolipid ceramide
        res_lookup = {}
        for chain in chains:
            for res in chain.residues:
                res_lookup[(chain.chain_id, res.resseq)] = res

        for tree, link_atoms, prot_links, cer_links in trees:
            if not tree:
                continue
            # Filter glycan_links relevant to this tree
            tree_set = set(tree)
            tree_links = [
                gl for gl in glycan_links
                if (gl[0], gl[1]) in tree_set or (gl[3], gl[4]) in tree_set
            ]
            # Only sugar-sugar links for building the glycan molecule
            sugar_sugar_links = [
                gl for gl in tree_links
                if (gl[0], gl[1]) in tree_set and (gl[3], gl[4]) in tree_set
            ]

            if cer_links:
                # Glycolipid: build ceramide + sugar tree as one moleculetype
                cl = cer_links[0]  # use first ceramide link
                cer_res = res_lookup.get((cl[0], cl[1]))
                if cer_res is not None:
                    ct = builder.build_glycolipid_chain(
                        cer_res, tree, link_atoms, chains,
                        sugar_sugar_links, cl)
                    if ct is not None:
                        glycan_tops.append(ct)
                        print(f"Glycolipid {ct.name}: {len(ct.atoms)} atoms, "
                              f"{len(ct.bonds)} bonds")
                    else:
                        print(f"WARNING: Failed to build glycolipid topology",
                              file=sys.stderr)
            else:
                # Pure glycan tree
                ct = builder.build_glycan_chain(tree, link_atoms, chains,
                                                sugar_sugar_links)
                if ct is not None:
                    glycan_tops.append(ct)
                    print(f"Glycan {ct.name}: {len(ct.atoms)} atoms, "
                          f"{len(ct.bonds)} bonds")

            # Collect protein-sugar links relevant to this tree
            for pl in prot_links:
                if (pl[3], pl[4]) in tree_set and pl not in protein_sugar_links:
                    protein_sugar_links.append(pl)

    chain_tops.extend(glycan_tops)

    # Build small molecule topologies (single-residue CGenFF molecules)
    small_mol_tops = []  # list of (ChainTopology, count)
    if small_mol_counts:
        for resname, count in small_mol_counts.items():
            # Build a single-residue chain for this molecule type
            dummy_chain = PDBChain(chain_id='X')
            # Find a representative residue from the original parsed chains
            rep_res = None
            for ch in read_pdb_chains(input_path):
                for r in ch.residues:
                    if r.resname == resname:
                        rep_res = r
                        break
                if rep_res is not None:
                    break
            if rep_res is None:
                print(f"WARNING: Cannot find {resname} residue for topology",
                      file=sys.stderr)
                continue
            rep_res.resseq = 1  # reset to 1 (GRO may have global numbering)
            dummy_chain.residues = [rep_res]
            ct = builder.build_chain(dummy_chain)
            if ct is not None:
                ct.name = resname
                small_mol_tops.append((ct, count))
                print(f"Small molecule {resname}: {len(ct.atoms)} atoms, "
                      f"{count} copies")

    if not chain_tops and not small_mol_tops:
        print("Error: No topologies built", file=sys.stderr)
        sys.exit(1)

    # Merge chains if requested
    if args.merge and len(chain_tops) > 1:
        merged = _merge_chains(chain_tops)
        chain_tops = [merged]
        print(f"Merged into single moleculetype: {merged.name}")

    # Output paths
    if args.output:
        top_path = Path(args.output)
    else:
        top_path = orig_input_path.parent / 'topol.top'

    out_dir = top_path.parent

    # Collect inter-chain SS bonds (both chains different)
    interchain_ss = []
    for ch1, res1, ch2, res2 in ss_bonds:
        if ch1 != ch2:
            interchain_ss.append((ch1, res1, ch2, res2))

    # Build per-chain bonded_types list for write_top
    bonded_types_list = []
    for ct in chain_tops:
        bt = (builder.carb_bonded_types
              if (ct.name.startswith('Glycan_') or ct.name.startswith('Glycolipid_'))
              and builder.carb_bonded_types
              else builder.bonded_types)
        bonded_types_list.append(bt)

    # Write position restraint files (still separate, used with #ifdef POSRES)
    for ct in chain_tops:
        posre_path = out_dir / f"posre_{ct.name}.itp"
        write_posre(ct, posre_path)
        print(f"Wrote {posre_path}")

    # Position restraints for small molecules too
    for ct, count in small_mol_tops:
        posre_path = out_dir / f"posre_{ct.name}.itp"
        write_posre(ct, posre_path)
        print(f"Wrote {posre_path}")

    # Write small molecule .itp files
    for ct, count in small_mol_tops:
        itp_path = out_dir / f"{ct.name}.itp"
        bt = builder.bonded_types
        if builder.carb_bonded_types:
            bt = builder.carb_bonded_types
        with open(itp_path, 'w') as f:
            f.write(f"; Moleculetype: {ct.name}\n")
            f.write(f"; Generated by dvbfixer top\n\n")
            _write_moleculetype(f, ct, bt)
        print(f"Wrote {itp_path}")

    # Write inter-chain bond file if needed (SS bonds + protein-glycan bonds)
    has_interchain_bonds = interchain_ss or protein_sugar_links
    if has_interchain_bonds:
        ss_path = out_dir / "interchain_ss.itp"
        _write_interchain_ss(interchain_ss, chain_tops, protein_chains, ss_path,
                             args.ff, protein_sugar_links)
        n_ss = len(interchain_ss)
        n_pg = len(protein_sugar_links)
        parts = []
        if n_ss:
            parts.append(f"{n_ss} inter-chain SS bond(s)")
        if n_pg:
            parts.append(f"{n_pg} protein-glycan bond(s)")
        print(f"Wrote {ss_path} ({', '.join(parts)})")
        print("WARNING: interchain_ss.itp must stay at the end of topol.top, after [ molecules ].")
        print("         After gmx solvate/genion, move the #include line below SOL/ion entries.")

    # Detect ions/BUF, small molecules, and water in PDB (preserving PDB order)
    ions_path = ff_dir / 'ions.itp'
    ion_names = set()
    if ions_path.exists():
        ion_names = _parse_ion_names(ions_path)
    small_mol_names_set = {ct.name for ct, _ in small_mol_tops}
    # Count all extra molecules preserving PDB order
    countable_names = ion_names | small_mol_names_set | _WATER_RESNAMES
    extra_molecules = _count_molecules(input_path, countable_names)
    # Fix water count: _count_molecules uses (chain, resseq) dedup which breaks
    # for large systems where PDB wraps resseq at 9999. Use atom count / 3 instead.
    water_count = _count_water(input_path)
    extra_molecules = [
        ('SOL', water_count) if name in _WATER_RESNAMES else (name, count)
        for name, count in extra_molecules
    ]
    if extra_molecules:
        for mol_name, mol_count in extra_molecules:
            print(f"Found {mol_count} {mol_name} molecule(s) in PDB")

    # Write TOP file with modular .itp includes
    small_mol_names = [ct.name for ct, _ in small_mol_tops]
    system_name = orig_input_path.stem
    write_top(chain_tops, top_path, ff_dir, args.ff, bonded_types_list,
              args.water, system_name,
              has_interchain_ss=has_interchain_bonds,
              extra_molecules=extra_molecules,
              small_mol_itps=small_mol_names)
    print(f"Wrote {out_dir / 'ffparams.itp'}")
    for ct in chain_tops:
        print(f"Wrote {out_dir / ct.name}.itp")
    print(f"Wrote {out_dir / 'water.itp'}")
    print(f"Wrote {out_dir / 'ions.itp'}")
    print(f"Wrote {top_path}")

    # Write output PDB with topology-matched atom names
    if args.pdb:
        pdb_path = Path(args.pdb)
    else:
        pdb_path = out_dir / 'conf.pdb'
    # Extract CRYST1 (box vectors) from input PDB
    cryst1_line = None
    with open(input_path) as f:
        for line in f:
            if line.startswith('CRYST1'):
                cryst1_line = line
                break

    # Collect ion/BUF/water PDB lines for output
    extra_pdb_lines = []
    if extra_molecules:
        extra_mol_names = ion_names | _WATER_RESNAMES | {ct.name for ct, _ in small_mol_tops}
        extra_pdb_lines = _extract_molecule_lines(input_path, extra_mol_names)
    write_pdb(chain_tops, pdb_path, extra_pdb_lines=extra_pdb_lines,
              cryst1=cryst1_line)
    print(f"Wrote {pdb_path}")

    # Clean up temp PDB from GRO conversion
    if tmp_pdb is not None:
        tmp_pdb.unlink(missing_ok=True)


def _write_interchain_ss(ss_bonds, chain_tops, protein_chains, path, ff_type,
                         protein_sugar_links=None):
    """Write inter-chain bond topology (SS bonds + protein-glycan bonds).

    Uses [ intermolecular_interactions ] which must be #included in .top
    after [ molecules ].
    """
    # Compute global atom offsets per chain topology (in [ molecules ] order)
    chain_offsets = {}  # ct.name -> offset
    offset = 0
    for ct in chain_tops:
        chain_offsets[ct.name] = offset
        offset += len(ct.atoms)

    def find_atom(chain_id, resseq, atomname):
        """Find global atom index by chain_id, resseq, atomname."""
        for ct in chain_tops:
            ct_offset = chain_offsets[ct.name]
            for atom in ct.atoms:
                if (atom.chain_id == chain_id and
                    atom.orig_resseq == resseq and
                    atom.atomname == atomname):
                    return ct_offset + atom.index
        return None

    with open(path, 'w') as f:
        f.write("; Inter-chain bonds (SS + protein-glycan)\n")
        f.write("; Generated by dvbfixer top\n\n")
        f.write("[ intermolecular_interactions ]\n\n")
        f.write("[ bonds ]\n")
        f.write(";  ai    aj  funct   r0 (nm)   k (kJ/mol/nm^2)\n")

        # SS bonds
        for ch1, res1, ch2, res2 in ss_bonds:
            idx1 = find_atom(ch1, res1, 'SG')
            idx2 = find_atom(ch2, res2, 'SG')
            if idx1 is not None and idx2 is not None:
                f.write(f"{idx1:5d} {idx2:5d}     6    0.204   250000\n")
                f.write(f"; {ch1}:{res1}:SG - {ch2}:{res2}:SG\n")

        # Protein-glycan bonds (ASN ND2 - NAG C1)
        if protein_sugar_links:
            for prot_ch, prot_rs, prot_atom, sug_ch, sug_rs in protein_sugar_links:
                idx1 = find_atom(prot_ch, prot_rs, prot_atom)
                idx2 = find_atom(sug_ch, sug_rs, 'C1')
                if idx1 is not None and idx2 is not None:
                    f.write(f"{idx1:5d} {idx2:5d}     6    0.1430   250000\n")
                    f.write(f"; {prot_ch}:{prot_rs}:{prot_atom} - "
                            f"{sug_ch}:{sug_rs}:C1\n")


def _merge_chains(chain_tops):
    """Merge multiple chain topologies into one."""
    merged = ChainTopology(
        name="Protein",
        nrexcl=chain_tops[0].nrexcl,
    )

    offset = 0
    for ct in chain_tops:
        for atom in ct.atoms:
            new_atom = AtomEntry(
                index=atom.index + offset,
                atom_type=atom.atom_type,
                resnr=atom.resnr,
                resname=atom.resname,
                atomname=atom.atomname,
                cgnr=atom.cgnr,
                charge=atom.charge,
                mass=atom.mass,
            )
            merged.atoms.append(new_atom)

        for i, j in ct.bonds:
            merged.bonds.append((i + offset, j + offset))
        for i, j in ct.pairs:
            merged.pairs.append((i + offset, j + offset))
        for i, j, k in ct.angles:
            merged.angles.append((i + offset, j + offset, k + offset))
        for dih in ct.dihedrals:
            if len(dih) == 5:
                i, j, k, l, t = dih
                merged.dihedrals.append((i + offset, j + offset, k + offset, l + offset, t))
            else:
                i, j, k, l = dih
                merged.dihedrals.append((i + offset, j + offset, k + offset, l + offset))
        for i, j, k, l in ct.impropers:
            merged.impropers.append((i + offset, j + offset, k + offset, l + offset))
        for cm in ct.cmap:
            merged.cmap.append(tuple(x + offset for x in cm))

        max_idx = max(a.index for a in ct.atoms) if ct.atoms else 0
        offset += max_idx

    return merged


if __name__ == '__main__':
    main()
