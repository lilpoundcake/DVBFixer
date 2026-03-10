"""Split PDB chains empirically using C-N distance, residue numbering, and nearest-atom gaps."""

import argparse
import sys
from pathlib import Path

import numpy as np

CHAIN_IDS = (
    list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    + list("abcdefghijklmnopqrstuvwxyz")
    + [str(i) for i in range(10)]
)

# Max C(i)->N(i+1) distance in angstroms. Normal peptide bond ~1.33 A; in MD
# trajectories all real bonds stay < 2.0 A, while chain breaks are > 3.5 A.
DEFAULT_DISTANCE_CUTOFF = 2.5

# For residues lacking C/N backbone (sugars, ligands, etc.) use the minimum
# distance between any atom of the last residue and any atom of the next.
DEFAULT_GAP_CUTOFF = 15.0

WATER_RESNAMES = {'HOH', 'WAT', 'TIP3', 'TIP', 'SOL', 'T3P', 'T4P', 'T5P'}

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Common protonation/tautomer variants
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM",
    # Terminal patches
    "ACE", "NME", "NH2",
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer split",
        description="Empirically split chains in a PDB file. Detects chain breaks by "
        "residue number resets and/or C-N inter-residue distance, assigns unique "
        "chain IDs, and inserts TER records."
    )
    p.add_argument("input", help="Input PDB file")
    p.add_argument("-o", "--output", help="Output PDB file (default: <input>_split.pdb)")
    p.add_argument(
        "-d", "--distance-cutoff", type=float, default=DEFAULT_DISTANCE_CUTOFF,
        help=f"C->N peptide bond cutoff in angstroms (default: {DEFAULT_DISTANCE_CUTOFF})"
    )
    p.add_argument(
        "-g", "--gap-cutoff", type=float, default=DEFAULT_GAP_CUTOFF,
        help=f"Min nearest-atom distance between consecutive residues to call a break "
             f"when C/N atoms are missing (default: {DEFAULT_GAP_CUTOFF} A)"
    )
    p.add_argument(
        "--no-distance", action="store_true",
        help="Disable all distance-based detection, use only residue numbering"
    )
    p.add_argument(
        "--no-renumber", action="store_true",
        help="Keep original residue numbers (default: renumber per chain starting from 1)"
    )
    p.add_argument(
        "--keep-water", action="store_true",
        help="Keep water molecules (HOH, WAT, TIP3, SOL) in output (default: remove)"
    )
    p.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print detected chain info"
    )
    return p.parse_args(argv)


def is_atom_line(line):
    return line.startswith(("ATOM  ", "HETATM"))


def get_coords(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def get_resid(line):
    """Return (resSeq, iCode) tuple — the full residue identifier."""
    return (int(line[22:26].strip()), line[26] if len(line) > 26 else " ")


def get_resseq(line):
    """Return just the integer residue sequence number."""
    return int(line[22:26].strip())


def get_atom_name(line):
    return line[12:16].strip()


def get_resname(line):
    return line[17:20].strip()


def set_chain_id(line, chain_id):
    return line[:21] + chain_id + line[22:]


def set_resid(line, resid, icode=" "):
    """Set residue sequence number and insertion code (col 22-26 + col 26)."""
    return line[:22] + f"{resid:4d}" + icode + line[27:]


def make_ter(serial, resname, chain_id, resid):
    return f"TER   {serial:>5d}      {resname:<3s} {chain_id}{resid:>4d}\n"


def find_chain_breaks(atom_lines, distance_cutoff, gap_cutoff, use_distance):
    """Return sorted list of atom_lines indices where new chains start.

    Detection criteria (applied in order):
      1. Residue sequence number goes backwards — always a chain break.
      2. C->N peptide bond distance exceeds distance_cutoff — chain break
         (only for residues that have both C and N backbone atoms).
      3. Nearest-atom distance between consecutive residues exceeds gap_cutoff —
         fallback for non-protein residues (sugars, ligands, ions) that lack
         C/N backbone atoms.
    """
    breaks = {0}

    # Group consecutive atoms into residues.
    # Each entry: (first_atom_idx, resseq, resname, C_coord, N_coord, all_coords)
    residues = []
    cur_resid = None  # (resSeq, iCode) tuple
    cur_start = 0
    cur_resname = None
    c_coord = None
    n_coord = None
    all_coords = []

    def flush():
        if cur_resid is not None:
            coords_arr = np.array(all_coords) if all_coords else None
            residues.append((cur_start, cur_resid[0], cur_resname, c_coord, n_coord, coords_arr))

    for i, line in enumerate(atom_lines):
        rid = get_resid(line)  # (resSeq, iCode)
        aname = get_atom_name(line)

        if rid != cur_resid:
            flush()
            cur_resid = rid
            cur_start = i
            cur_resname = get_resname(line)
            c_coord = None
            n_coord = None
            all_coords = []

        all_coords.append(get_coords(line))

        if aname == "C":
            c_coord = get_coords(line)
        elif aname == "N":
            n_coord = get_coords(line)

    flush()

    for i in range(1, len(residues)):
        _, resseq_prev, rn_prev, c_prev, _, coords_prev = residues[i - 1]
        first_idx_cur, resseq_cur, rn_cur, _, n_cur, coords_cur = residues[i]

        # Criterion 1: residue sequence number goes strictly backwards
        # (equal resSeq with different iCode = insertion, not a break)
        if resseq_cur < resseq_prev:
            breaks.add(first_idx_cur)
            continue

        if not use_distance:
            continue

        both_protein = rn_prev in STANDARD_RESIDUES and rn_cur in STANDARD_RESIDUES

        # Criterion 2: C->N peptide bond distance (protein residues only)
        if both_protein and c_prev is not None and n_cur is not None:
            if np.linalg.norm(c_prev - n_cur) > distance_cutoff:
                breaks.add(first_idx_cur)
            continue

        # Criterion 3: nearest-atom gap (non-protein or missing backbone atoms)
        if coords_prev is not None and coords_cur is not None:
            diff = coords_prev[:, None, :] - coords_cur[None, :, :]
            min_dist = np.sqrt((diff ** 2).sum(axis=2)).min()
            if min_dist > gap_cutoff:
                breaks.add(first_idx_cur)

    return sorted(breaks)


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_split")
    use_distance = not args.no_distance
    renumber = not args.no_renumber

    with open(input_path) as f:
        lines = f.readlines()

    atom_lines = []
    atom_orig_indices = []
    for i, line in enumerate(lines):
        if is_atom_line(line):
            atom_lines.append(line)
            atom_orig_indices.append(i)

    if not atom_lines:
        print("No ATOM/HETATM records found.", file=sys.stderr)
        sys.exit(1)

    breaks = find_chain_breaks(atom_lines, args.distance_cutoff, args.gap_cutoff, use_distance)
    n_chains = len(breaks)

    if n_chains > len(CHAIN_IDS):
        print(f"Too many chains ({n_chains}) for available chain IDs.", file=sys.stderr)
        sys.exit(1)

    # Map each atom index to its chain index
    chain_for_atom = np.zeros(len(atom_lines), dtype=int)
    for ci in range(n_chains):
        start = breaks[ci]
        end = breaks[ci + 1] if ci + 1 < n_chains else len(atom_lines)
        chain_for_atom[start:end] = ci

    # Build per-chain residue renumbering maps
    resid_maps = {}
    if renumber:
        for ci in range(n_chains):
            start = breaks[ci]
            end = breaks[ci + 1] if ci + 1 < n_chains else len(atom_lines)
            seen = {}
            counter = 1
            for j in range(start, end):
                rid = get_resid(atom_lines[j])  # (resSeq, iCode) tuple
                if rid not in seen:
                    seen[rid] = counter
                    counter += 1
            resid_maps[ci] = seen

    if args.verbose:
        print(f"Detected {n_chains} chain(s):")
        for ci in range(n_chains):
            start = breaks[ci]
            end = (breaks[ci + 1] if ci + 1 < n_chains else len(atom_lines)) - 1
            n_atoms = end - start + 1
            n_res = len(resid_maps[ci]) if renumber else len(set(get_resid(atom_lines[j]) for j in range(start, end + 1)))
            print(
                f"  Chain {CHAIN_IDS[ci]}: {n_atoms:>5d} atoms, {n_res:>4d} residues, "
                f"{get_resname(atom_lines[start])} {get_resseq(atom_lines[start])} - "
                f"{get_resname(atom_lines[end])} {get_resseq(atom_lines[end])}"
            )

    # Set of break atom indices for quick TER insertion lookup
    break_set = set(breaks)

    # Build output
    output_lines = []
    ai = 0  # atom_lines index
    serial = 1

    for line in lines:
        if is_atom_line(line) and ai < len(atom_lines):
            ci = int(chain_for_atom[ai])
            chain_id = CHAIN_IDS[ci]

            new_line = set_chain_id(line, chain_id)
            if renumber:
                rid = get_resid(line)
                new_line = set_resid(new_line, resid_maps[ci][rid], icode=" ")
            new_line = new_line[:6] + f"{serial:5d}" + new_line[11:]
            output_lines.append(new_line)

            # Insert TER after last atom of a chain (before next chain starts)
            next_ai = ai + 1
            if next_ai in break_set:
                serial += 1
                resname = get_resname(new_line)
                resid = int(new_line[22:26].strip())
                output_lines.append(make_ter(serial, resname, chain_id, resid))

            serial += 1
            ai += 1
        elif line.startswith("TER"):
            continue  # skip original TER records
        else:
            output_lines.append(line)

    if not args.keep_water:
        filtered = []
        for line in output_lines:
            if is_atom_line(line) and get_resname(line) in WATER_RESNAMES:
                continue
            if line.startswith("TER") and len(line) > 20 and line[17:20].strip() in WATER_RESNAMES:
                continue
            filtered.append(line)
        output_lines = filtered

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    print(f"Wrote {output_path} with {n_chains} chain(s)")
