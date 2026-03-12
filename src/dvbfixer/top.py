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
}


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
    chain_lines = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            chain_id = line[21]
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


def detect_glycan_links(chains):
    """Auto-detect glycosidic bonds from C1-O distances using parsed chain data.

    Detects:
    - Sugar-sugar links: C1 of one sugar within 2.0 A of Ox of another
    - Protein-sugar links: C1 of NAG within 2.0 A of ASN ND2

    Returns list of GlycanLink tuples.
    """
    sugar_names = set(PDB_TO_CARB.keys())
    links = []

    # Collect C1 atoms from sugars and ND2/O* from sugars+ASN
    c1_atoms = []
    target_atoms = []

    for chain in chains:
        ch = chain.chain_id
        for res in chain.residues:
            for aname, x, y, z in res.atoms:
                if res.resname in sugar_names and aname == 'C1':
                    c1_atoms.append((ch, res.resname, res.resseq, x, y, z))

                if res.resname in sugar_names and aname.startswith('O') and aname != 'O5':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

                if res.resname == 'ASN' and aname == 'ND2':
                    target_atoms.append((ch, res.resname, res.resseq, aname, x, y, z))

    # Find C1-target pairs within 2.0 A
    for ch1, rn1, rs1, x1, y1, z1 in c1_atoms:
        for ch2, rn2, rs2, aname2, x2, y2, z2 in target_atoms:
            if ch1 == ch2 and rs1 == rs2:
                continue
            d = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
            if d < 2.0:
                # Link: acceptor (has C1) bonds to donor (has Ox/ND2)
                links.append((ch2, rs2, aname2, ch1, rs1, 'C1'))

    return links


def build_glycan_trees(glycan_links, chains):
    """Build glycan trees from linkage list.

    Returns list of glycan trees, each a list of (chain_id, resseq, resname)
    in topological order (root first).
    """
    # Build adjacency: sugar -> list of child sugars
    sugar_residues = set()
    for chain in chains:
        for res in chain.residues:
            if res.resname in PDB_TO_CARB:
                sugar_residues.add((chain.chain_id, res.resseq))

    # Find root sugars (connected to protein ASN, not to another sugar)
    children = defaultdict(list)  # parent -> [child]
    parent_of = {}  # child -> parent
    protein_links = []  # [(protein_chain, protein_resseq, protein_atom, sugar_chain, sugar_resseq)]

    for don_ch, don_rs, don_atom, acc_ch, acc_rs, acc_atom in glycan_links:
        acc_key = (acc_ch, acc_rs)
        don_key = (don_ch, don_rs)
        if don_key in sugar_residues and acc_key in sugar_residues:
            # sugar-sugar link: donor is parent, acceptor is child
            children[don_key].append((acc_key, don_atom))
            parent_of[acc_key] = (don_key, don_atom)
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
            for child, donor_atom in children.get(node, []):
                # Record which O on parent is linked
                if node not in link_atoms:
                    link_atoms[node] = set()
                link_atoms[node].add(donor_atom)
                # Child's C1-O1 is linked (remove HO1)
                if child not in link_atoms:
                    link_atoms[child] = set()
                link_atoms[child].add('O1')
                queue.append(child)

        # Also mark root's O1 as linked (bonds to protein or is reducing end)
        for pl in protein_links:
            if (pl[3], pl[4]) == root:
                if root not in link_atoms:
                    link_atoms[root] = set()
                link_atoms[root].add('O1')

        trees.append((tree, link_atoms, protein_links))

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

    def build_chain(self, chain, ss_residues=None):
        """Build topology for a single chain.

        ss_residues: set of resseq involved in SS bonds in this chain.
        """
        if ss_residues is None:
            ss_residues = set()

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
            name=f"Protein_chain_{chain.chain_id}",
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

            # Add missing protonation H atoms with estimated coordinates.
            # Only for ASP/GLU/HIS protonation variants (ASPP HD2, GLUP HE2,
            # HSP/HSD/HSE HD1/HE2).
            _PROT_RTP_NAMES = {'ASPP', 'GLUP', 'HSP', 'HSD', 'HSE',
                               'ASH', 'GLH', 'HIP', 'HID', 'HIE',
                               'ASPH', 'GLUH', 'HISH', 'HISD', 'HISE'}
            if rtp_name.upper() in _PROT_RTP_NAMES:
                rtp_bond_graph = {}
                for b1, b2 in rtp_res.bonds:
                    if b1.startswith(('-', '+')) or b2.startswith(('-', '+')):
                        continue
                    rtp_bond_graph.setdefault(b1, set()).add(b2)
                    rtp_bond_graph.setdefault(b2, set()).add(b1)

                for atom_name, _, _, _ in rtp_res.atoms:
                    if rtp_to_pdb.get(atom_name) is not None:
                        continue
                    if not atom_name.startswith('H'):
                        continue
                    for neighbor in rtp_bond_graph.get(atom_name, []):
                        nb_pdb = rtp_to_pdb.get(neighbor)
                        if nb_pdb and nb_pdb in pdb_coords:
                            hx, hy, hz = pdb_coords[nb_pdb]
                            pdb_coords[atom_name] = (hx + 1.0, hy, hz)
                            rtp_to_pdb[atom_name] = atom_name
                            if self.verbose:
                                print(f"    Adding {rtp_name}:{atom_name} near "
                                      f"{nb_pdb} ({res.chain_id}:{res.resseq})")
                            break

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
        """Apply CHARMM terminal patches (TDB)."""
        n_res = len(chain.residues)

        # Determine N-terminal patch
        first_res = chain.residues[0].resname
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
        chain_name = f"Glycan_{first_ch}_{first_rs}"

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

            # Map PDB name -> CHARMM name
            rtp_name = PDB_TO_CARB.get(res.resname)
            if rtp_name is None or rtp_name not in self.residues:
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
                # For carbs, skip atoms not in PDB: H atoms (may be missing)
                # and linked O atoms (O1 at glycosidic bond sites)
                if pdb_name not in pdb_coords:
                    if atom_name.startswith('H') or atom_name in linked_os:
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
            rtp_name = PDB_TO_CARB.get(res.resname)
            if rtp_name is None or rtp_name not in self.residues:
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


def write_pdb(chain_tops, path, extra_pdb_lines=None):
    """Write a PDB file with atom names matching the topology.

    extra_pdb_lines: list of raw PDB ATOM/HETATM lines to append
    (for ions/BUF particles not built as chain topologies).
    """
    serial = 0
    with open(path, 'w') as f:
        f.write("REMARK    Generated by dvbfixer top\n")
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
                    f"ATOM  {serial:5d} {name_field}"
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
                    f"TER   {serial:5d}      "
                    f"{res_chain}"
                    f"{last.orig_resseq:4d}\n"
                )
        # Append ion/BUF particles with renumbered serials
        if extra_pdb_lines:
            for line in extra_pdb_lines:
                serial += 1
                # Rewrite serial number (cols 7-11)
                f.write(f"{line[:6]}{serial:5d}{line[11:]}")
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

"""


def write_top(chain_tops, path, ff_dir, ff_type, bonded_types_list,
              water_model='tip3p', system_name='Protein',
              has_interchain_ss=False, extra_molecules=None):
    """Write topology: FF params in ffparams.itp, rest in topol.top.

    extra_molecules: list of (name, count) for ions/BUF/etc found in PDB.
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

    print(f"Found {len(chains)} chain(s): {', '.join(c.chain_id for c in chains)}")

    # Build topology builder first (need its residue dict for chain filtering)
    builder = TopologyBuilder(ff_dir, args.ff, args.verbose)

    # Classify chains into protein vs non-protein
    protein_chains = []
    other_chains = []
    sugar_names = set(PDB_TO_CARB.keys())
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
        else:
            # Check if chain has any FF-recognized residues
            has_known = any(
                r.resname in builder.residues or r.resname in sugar_names
                for r in chain.residues
            )
            if has_known:
                other_chains.append(chain)

    if not protein_chains and not other_chains:
        print("Error: No recognized chains found", file=sys.stderr)
        sys.exit(1)

    # Detect SS bonds
    ss_bonds = read_ssbonds(input_path)

    # Parse explicit --ss flags
    for ss_spec in args.ss:
        parts = ss_spec.split(':')
        if len(parts) == 4:
            ss_bonds.append((parts[0], int(parts[1]), parts[2], int(parts[3])))

    # Build per-chain SS residue sets
    chain_ss = defaultdict(set)
    for ch1, res1, ch2, res2 in ss_bonds:
        chain_ss[ch1].add(res1)
        chain_ss[ch2].add(res2)

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

        ct = builder.build_chain(chain, chain_ss.get(chain.chain_id, set()))
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

    # Build non-protein chains (no terminal patches, no SS bonds)
    # Skip chains that are entirely sugars — they're handled by glycan detection
    for chain in other_chains:
        all_sugar = all(
            r.resname in sugar_names or r.resname in PDB_TO_CARB.values()
            for r in chain.residues
        )
        if all_sugar:
            if args.verbose:
                print(f"\nSkipping chain {chain.chain_id} (all sugar residues, "
                      f"handled by glycan detection)")
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

    if not chain_tops:
        print("Error: No topologies built", file=sys.stderr)
        sys.exit(1)

    # Detect and build glycan chains
    glycan_tops = []
    glycan_links = detect_glycan_links(chains)
    if glycan_links:
        trees = build_glycan_trees(glycan_links, chains)
        protein_sugar_links = []  # for intermolecular interactions
        for tree, link_atoms, prot_links in trees:
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

            ct = builder.build_glycan_chain(tree, link_atoms, chains, sugar_sugar_links)
            if ct is not None:
                glycan_tops.append(ct)
                print(f"Glycan {ct.name}: {len(ct.atoms)} atoms, "
                      f"{len(ct.bonds)} bonds")

            # Collect protein-sugar links for intermolecular interactions
            for pl in prot_links:
                protein_sugar_links.append(pl)

    chain_tops.extend(glycan_tops)

    # Merge chains if requested
    if args.merge and len(chain_tops) > 1:
        merged = _merge_chains(chain_tops)
        chain_tops = [merged]
        print(f"Merged into single moleculetype: {merged.name}")

    # Output paths
    if args.output:
        top_path = Path(args.output)
    else:
        top_path = input_path.parent / 'topol.top'

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
              if ct.name.startswith('Glycan_') and builder.carb_bonded_types
              else builder.bonded_types)
        bonded_types_list.append(bt)

    # Write position restraint files (still separate, used with #ifdef POSRES)
    for ct in chain_tops:
        posre_path = out_dir / f"posre_{ct.name}.itp"
        write_posre(ct, posre_path)
        print(f"Wrote {posre_path}")

    # Write inter-chain SS bond file if needed
    if interchain_ss:
        ss_path = out_dir / "interchain_ss.itp"
        _write_interchain_ss(interchain_ss, chain_tops, protein_chains, ss_path,
                             args.ff)
        print(f"Wrote {ss_path} ({len(interchain_ss)} inter-chain SS bond(s))")
        print("WARNING: interchain_ss.itp must stay at the end of topol.top, after [ molecules ].")
        print("         After gmx solvate/genion, move the #include line below SOL/ion entries.")

    # Detect ions/BUF particles in PDB
    ions_path = ff_dir / 'ions.itp'
    extra_molecules = []
    if ions_path.exists():
        ion_names = _parse_ion_names(ions_path)
        extra_molecules = _count_molecules(input_path, ion_names)
        if extra_molecules:
            for mol_name, mol_count in extra_molecules:
                print(f"Found {mol_count} {mol_name} molecule(s) in PDB")

    # Write TOP file with modular .itp includes
    system_name = input_path.stem
    write_top(chain_tops, top_path, ff_dir, args.ff, bonded_types_list,
              args.water, system_name,
              has_interchain_ss=bool(interchain_ss),
              extra_molecules=extra_molecules)
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
    # Collect ion/BUF PDB lines for output
    extra_pdb_lines = []
    if extra_molecules:
        extra_pdb_lines = _extract_molecule_lines(input_path, ion_names)
    write_pdb(chain_tops, pdb_path, extra_pdb_lines=extra_pdb_lines)
    print(f"Wrote {pdb_path}")


def _write_interchain_ss(ss_bonds, chain_tops, protein_chains, path, ff_type):
    """Write inter-chain disulfide bond topology to a separate .itp file.

    Uses [ intermolecular_interactions ] which must be #included in .top
    after [ molecules ].
    """
    # Build mapping: (chain_id) -> (topology, chain)
    chain_id_to_top = {}
    for ct, chain in zip(chain_tops, protein_chains):
        chain_id_to_top[chain.chain_id] = (ct, chain)

    # Compute chain offsets in global numbering
    chain_offsets = {}
    offset = 0
    for ct, chain in zip(chain_tops, protein_chains):
        chain_offsets[chain.chain_id] = offset
        offset += len(ct.atoms)

    def find_sg(chain_id, resseq):
        if chain_id not in chain_id_to_top:
            return None
        ct, chain = chain_id_to_top[chain_id]
        for res_i, res in enumerate(chain.residues):
            if res.resseq == resseq:
                resnr = res_i + 1
                for atom in ct.atoms:
                    if atom.resnr == resnr and atom.atomname == 'SG':
                        return chain_offsets[chain_id] + atom.index
                break
        return None

    with open(path, 'w') as f:
        f.write("; Inter-chain disulfide bonds\n")
        f.write("; Generated by dvbfixer top\n\n")
        f.write("[ intermolecular_interactions ]\n\n")
        f.write("[ bonds ]\n")
        f.write(";  ai    aj  funct   r0 (nm)   k (kJ/mol/nm^2)\n")
        for ch1, res1, ch2, res2 in ss_bonds:
            idx1 = find_sg(ch1, res1)
            idx2 = find_sg(ch2, res2)
            if idx1 is not None and idx2 is not None:
                f.write(f"{idx1:5d} {idx2:5d}     6    0.204   250000\n")
                f.write(f"; {ch1}:{res1}:SG - {ch2}:{res2}:SG\n")


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
