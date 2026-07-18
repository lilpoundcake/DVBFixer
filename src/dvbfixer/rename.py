"""Rename non-canonical residues to standard PDB names.

Converts AMBER protonation names (HIE/HID/HIP, ASH, GLH, CYX, CYM, LYN),
CHARMM names (HSD/HSE/HSP), and other common variants back to their canonical
three-letter codes. Text-based — does not modify coordinates or atoms.
"""

import argparse
import sys
from pathlib import Path

# Non-canonical -> canonical residue name mapping
CANONICAL_MAP = {
    # AMBER protonation variants
    'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
    'ASH': 'ASP',
    'GLH': 'GLU',
    'CYX': 'CYS', 'CYM': 'CYS',
    'LYN': 'LYS',
    # CHARMM variants
    'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
    # Selenomethionine
    'MSE': 'MET',
    # Caps (sometimes present)
    'ACE': 'ACE',  # keep as-is (not a standard AA)
    'NME': 'NME',
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer rename",
        description="Rename non-canonical residue names to standard PDB names. "
        "Converts AMBER (HIE/HID/HIP, ASH, GLH, CYX, CYM, LYN), "
        "CHARMM (HSD/HSE/HSP), and MSE to their canonical forms.",
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file")
    io.add_argument("-o", "--output", help="Output PDB file (default: <input>_canon.pdb)")

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument("-v", "--verbose", action="store_true",
                      help="Print each rename")

    return p.parse_args(argv)


def canonicalize_pdb(input_path, output_path, verbose=False):
    """Rename non-canonical residues in a PDB file. Returns number of renames."""
    with open(input_path) as f:
        lines = f.readlines()

    renamed = {}  # (chain, resid, old_name) -> new_name for dedup reporting
    out_lines = []

    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM', 'TER   ')):
            resname = line[17:20].strip()
            if resname in CANONICAL_MAP:
                new_name = CANONICAL_MAP[resname]
                if resname != new_name:
                    chain = line[21]
                    resid = line[22:26].strip()
                    key = (chain, resid, resname)
                    if key not in renamed:
                        renamed[key] = new_name
                        if verbose:
                            print(f"  {chain}:{resname}{resid} -> {new_name}")
                    line = line[:17] + f"{new_name:>3s}" + line[20:]
        out_lines.append(line)

    with open(output_path, 'w') as f:
        f.writelines(out_lines)

    return len(renamed)


def canonicalize_in_place(input_path, verbose=False):
    """Canonicalize a PDB file, overwriting it. Returns number of renames."""
    import tempfile
    tmp = Path(tempfile.mktemp(suffix='.pdb'))
    n = canonicalize_pdb(input_path, tmp, verbose)
    if n > 0:
        import shutil
        shutil.move(str(tmp), str(input_path))
    elif tmp.exists():
        tmp.unlink()
    return n


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_renamed")

    n = canonicalize_pdb(input_path, output_path, args.verbose)
    if n > 0:
        print(f"Renamed {n} non-canonical residue(s)")
    else:
        print("No non-canonical residues found")
    print(f"Wrote {output_path}")
