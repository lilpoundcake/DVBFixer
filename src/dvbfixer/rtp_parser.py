"""Parsers for GROMACS force field topology files (RTP, R2B, ARN, TDB, ATP)."""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class BondedTypes:
    """Default function types from [ bondedtypes ] header in RTP."""
    bond_type: int = 1
    angle_type: int = 1
    dihedral_type: int = 9
    improper_type: int = 4
    all_dihedrals: int = 1
    nrexcl: int = 3
    HH14: int = 1
    remove_dih: int = 0


@dataclass
class RTPResidue:
    """A residue entry from an RTP file."""
    name: str
    atoms: list = field(default_factory=list)      # [(name, type, charge, cgnr)]
    bonds: list = field(default_factory=list)       # [(atom1, atom2)]
    impropers: list = field(default_factory=list)   # [(a1, a2, a3, a4)]
    dihedrals: list = field(default_factory=list)   # [(a1, a2, a3, a4)] or [(a1,a2,a3,a4,type_name)]
    cmap: list = field(default_factory=list)        # [(a1, a2, a3, a4, a5)]


@dataclass
class TerminalPatch:
    """A terminal patch from a TDB file."""
    name: str
    delete: list = field(default_factory=list)      # [atom_name, ...]
    replace: list = field(default_factory=list)     # [(name, type, mass, charge)]
    add: list = field(default_factory=list)          # [(count, method, name, *ref_atoms, type, mass, charge, cgnr)]
    impropers: list = field(default_factory=list)   # [(a1, a2, a3, a4)]


def _strip_comment(line):
    """Remove ; comment from line."""
    idx = line.find(';')
    return line[:idx] if idx >= 0 else line


def parse_rtp(path):
    """Parse an RTP file. Returns (BondedTypes, dict[resname -> RTPResidue])."""
    path = Path(path)
    bonded_types = BondedTypes()
    residues = {}
    current_res = None
    current_section = None

    with open(path) as f:
        for line in f:
            stripped = _strip_comment(line).strip()
            if not stripped:
                continue

            # Check for section headers
            if stripped.startswith('[') and stripped.endswith(']'):
                section_name = stripped[1:-1].strip()

                # Top-level sections (residue names) vs sub-sections
                if section_name == 'bondedtypes':
                    current_res = None
                    current_section = 'bondedtypes'
                    continue
                elif section_name in ('atoms', 'bonds', 'impropers', 'dihedrals', 'cmap'):
                    current_section = section_name
                    continue
                else:
                    # New residue
                    current_res = RTPResidue(name=section_name)
                    residues[section_name] = current_res
                    current_section = None
                    continue

            parts = stripped.split()

            if current_section == 'bondedtypes':
                if len(parts) >= 8:
                    bonded_types = BondedTypes(
                        bond_type=int(parts[0]),
                        angle_type=int(parts[1]),
                        dihedral_type=int(parts[2]),
                        improper_type=int(parts[3]),
                        all_dihedrals=int(parts[4]),
                        nrexcl=int(parts[5]),
                        HH14=int(parts[6]),
                        remove_dih=int(parts[7]),
                    )
                current_section = None
                continue

            if current_res is None:
                continue

            if current_section == 'atoms' and len(parts) >= 4:
                current_res.atoms.append((
                    parts[0],           # name
                    parts[1],           # type
                    float(parts[2]),    # charge
                    int(parts[3]),      # charge group
                ))
            elif current_section == 'bonds' and len(parts) >= 2:
                current_res.bonds.append((parts[0], parts[1]))
            elif current_section == 'impropers' and len(parts) >= 4:
                current_res.impropers.append(tuple(parts[:4]))
            elif current_section == 'dihedrals' and len(parts) >= 4:
                # AMBER ILDN has named dihedral types as 5th column
                current_res.dihedrals.append(tuple(parts[:5]) if len(parts) >= 5 else tuple(parts[:4]))
            elif current_section == 'cmap' and len(parts) >= 5:
                current_res.cmap.append(tuple(parts[:5]))

    return bonded_types, residues


def parse_r2b(path):
    """Parse residue-to-building-block file.

    Returns dict[gmx_name -> (main, nter, cter, twter)].
    '-' means no entry available.
    """
    path = Path(path)
    mapping = {}
    with open(path) as f:
        for line in f:
            stripped = _strip_comment(line).strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) >= 5:
                gmx_name = parts[0]
                mapping[gmx_name] = (parts[1], parts[2], parts[3], parts[4])
            elif len(parts) >= 2:
                # Minimal: just main block
                mapping[parts[0]] = (parts[1], '-', '-', '-')
    return mapping


def parse_arn(path):
    """Parse atom renaming file.

    Returns dict[(resname, gromacs_name) -> ff_name].
    resname '*' means wildcard (all residues).
    """
    path = Path(path)
    mapping = {}
    with open(path) as f:
        for line in f:
            stripped = _strip_comment(line).strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) >= 3:
                mapping[(parts[0], parts[1])] = parts[2]
    return mapping


def parse_atomtypes(path):
    """Parse atomtypes.atp. Returns dict[type_name -> mass]."""
    path = Path(path)
    types = {}
    with open(path) as f:
        for line in f:
            stripped = _strip_comment(line).strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    types[parts[0]] = float(parts[1])
                except ValueError:
                    continue
    return types


def parse_tdb(path):
    """Parse a terminal database file (.n.tdb or .c.tdb).

    Returns dict[patch_name -> TerminalPatch].
    """
    path = Path(path)
    patches = {}
    current_patch = None
    current_section = None

    with open(path) as f:
        for line in f:
            stripped = _strip_comment(line).strip()
            if not stripped:
                continue

            if stripped.startswith('[') and stripped.endswith(']'):
                section_name = stripped[1:-1].strip()
                if section_name in ('delete', 'replace', 'add', 'impropers'):
                    current_section = section_name
                else:
                    current_patch = TerminalPatch(name=section_name)
                    patches[section_name] = current_patch
                    current_section = None
                continue

            if current_patch is None:
                continue

            parts = stripped.split()

            if current_section == 'delete' and len(parts) >= 1:
                current_patch.delete.append(parts[0])
            elif current_section == 'replace' and len(parts) >= 4:
                current_patch.replace.append((
                    parts[0],           # atom name
                    parts[1],           # new type
                    float(parts[2]),    # mass
                    float(parts[3]),    # charge
                ))
            elif current_section == 'add' and len(parts) >= 2:
                # Add lines come in pairs:
                # Line 1: count method name ref_atoms...
                # Line 2: type mass charge cgnr
                # But they can also be on a single line depending on format
                # The format is: count method name ref1 ref2 ref3 [ref4]
                #                type mass charge cgnr
                # We detect by first token being a number
                if parts[0].isdigit():
                    # This is a geometry/reference line
                    current_patch._pending_add = parts
                elif hasattr(current_patch, '_pending_add') and current_patch._pending_add:
                    # This is the type/mass/charge line
                    ref = current_patch._pending_add
                    count = int(ref[0])
                    method = int(ref[1])
                    name = ref[2]
                    ref_atoms = ref[3:]
                    atype = parts[0]
                    mass = float(parts[1])
                    charge = float(parts[2])
                    cgnr = int(parts[3]) if len(parts) > 3 else -1
                    current_patch.add.append({
                        'count': count,
                        'method': method,
                        'name': name,
                        'ref_atoms': ref_atoms,
                        'type': atype,
                        'mass': mass,
                        'charge': charge,
                        'cgnr': cgnr,
                    })
                    current_patch._pending_add = None
            elif current_section == 'impropers' and len(parts) >= 4:
                current_patch.impropers.append(tuple(parts[:4]))

    # Clean up _pending_add attributes
    for patch in patches.values():
        if hasattr(patch, '_pending_add'):
            del patch._pending_add

    return patches
