"""Renumber PDB residues using SEQRES alignment to remove insertion codes
and preserve correct gap positions. Updates all PDB sections that reference
residue numbers: ATOM, HETATM, TER, HELIX, SHEET, SSBOND, LINK, CISPEP,
HET, DBREF, SEQADV, CONECT, and relevant REMARKs (465, 500, 610)."""

import argparse
import sys
from pathlib import Path

WATER_RESNAMES = {'HOH', 'WAT', 'TIP3', 'TIP', 'SOL', 'T3P', 'T4P', 'T5P'}

# Three-letter → one-letter AA code (used by antibody numbering)
_AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # Common protonation / non-canonical variants — map to parent.
    'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'HSD': 'H', 'HSE': 'H', 'HSP': 'H',
    'CYX': 'C', 'CYM': 'C', 'ASH': 'D', 'GLH': 'E', 'GLUP': 'E', 'ASPP': 'D',
    'LYN': 'K', 'MSE': 'M',
    # GLYCAM glycoprotein residues map to the parent AA.
    'NLN': 'N', 'OLS': 'S', 'OLT': 'T',
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer renumber",
        description="Read SEQRES from a PDB file, align ATOM residues to the full "
        "sequence, and renumber to remove insertion codes while preserving "
        "gap positions. Updates all PDB sections referencing residue numbers."
    )
    p.add_argument("input", help="Input PDB file")
    p.add_argument("-o", "--output", help="Output PDB file (default: <input>_renum.pdb)")
    p.add_argument(
        "--keep-water", action="store_true",
        help="Keep water molecules (HOH, WAT, TIP3, SOL) in output (default: remove)"
    )
    p.add_argument(
        "--rename", action="store_true",
        help="Rename non-canonical residues (AMBER/CHARMM) to standard names before processing"
    )
    p.add_argument(
        "--scheme", choices=["seqres", "kabat", "chothia", "imgt", "martin", "eu", "aho"],
        default="seqres",
        help="Antibody numbering scheme. Default 'seqres' uses SEQRES-based "
             "sequential numbering (the original behaviour). The four V-domain "
             "schemes (kabat/chothia/imgt/martin/aho) are produced by ANARCI; "
             "'eu' uses Kabat for V-domains and EU positions for C-domains. "
             "Constant-region numbering always uses EU (Edelman 1969) regardless "
             "of the V-scheme — Kabat/Chothia/Martin don't define C-domain "
             "positions. Non-antibody chains fall back to SEQRES."
    )
    p.add_argument(
        "--chain-scheme", action="append", default=[],
        metavar="CHAIN:SCHEME",
        help="Per-chain scheme override (e.g. H:kabat). Repeatable. Wins over --scheme."
    )
    p.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print alignment details and gap positions"
    )
    return p.parse_args(argv)


def parse_seqres(lines):
    """Return dict: chain_id -> [resname, ...] from SEQRES records."""
    seqres = {}
    for line in lines:
        if not line.startswith("SEQRES"):
            continue
        chain = line[11]
        seqres.setdefault(chain, []).extend(line[19:].split())
    return seqres


def get_atom_residues(lines, chain):
    """Return ordered list of unique (resSeq, iCode, resname) for a chain."""
    residues = []
    seen = set()
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[21] != chain:
            continue
        key = (int(line[22:26].strip()), line[26])
        if key not in seen:
            seen.add(key)
            residues.append((key[0], key[1], line[17:20].strip()))
    return residues


def align_to_seqres(atom_residues, seqres):
    """Align ATOM residues to SEQRES via subsequence matching.
    Returns dict: (old_resSeq, old_iCode) -> new_resSeq.
    Non-SEQRES residues (waters, ligands) get None initially.
    """
    mapping = {}
    j = 0
    for resseq, icode, resname in atom_residues:
        search_j = j
        while search_j < len(seqres) and seqres[search_j] != resname:
            search_j += 1
        if search_j < len(seqres):
            mapping[(resseq, icode)] = search_j + 1
            j = search_j + 1
        else:
            mapping[(resseq, icode)] = None
    return mapping


# ---------------------------------------------------------------------------
# Helpers to remap a (chain, resSeq, iCode) reference in a PDB line
# ---------------------------------------------------------------------------

def remap_resid(line, chain_col, seq_cols, ic_col, chain_maps):
    """Remap a single residue reference in a line.
    chain_col: 0-based index of chain character
    seq_cols: (start, end) 0-based slice for resSeq (4 chars)
    ic_col: 0-based index of insertion code character
    Returns modified line, or original if no mapping found.

    `chain_maps[chain]` values can be either an int (legacy SEQRES path) or
    a (new_resseq, new_icode) tuple (antibody path with insertion codes).
    """
    if len(line) <= ic_col:
        return line
    chain = line[chain_col]
    seq_str = line[seq_cols[0]:seq_cols[1]].strip()
    if not seq_str or not seq_str.lstrip('-').isdigit():
        return line
    old_seq = int(seq_str)
    old_ic = line[ic_col]
    old_key = (old_seq, old_ic)
    if chain in chain_maps and old_key in chain_maps[chain]:
        val = chain_maps[chain][old_key]
        if isinstance(val, tuple):
            new_seq, new_ic = val
            if not new_ic:
                new_ic = " "
        else:
            new_seq, new_ic = val, " "
        line = line[:seq_cols[0]] + f"{new_seq:4d}" + line[seq_cols[1]:]
        line = line[:ic_col] + new_ic + line[ic_col + 1:]
    return line


def pad_line(line, min_len):
    """Ensure line is at least min_len characters (pad with spaces)."""
    if len(line.rstrip('\n')) < min_len:
        stripped = line.rstrip('\n')
        line = stripped + ' ' * (min_len - len(stripped)) + '\n'
    return line


# ---------------------------------------------------------------------------
# Section-specific updaters
# ---------------------------------------------------------------------------

def update_helix(line, chain_maps):
    """HELIX: init residue at (19, 21-24, 25), end residue at (31, 33-36, 37)."""
    line = pad_line(line, 38)
    line = remap_resid(line, 19, (21, 25), 25, chain_maps)
    line = remap_resid(line, 31, (33, 37), 37, chain_maps)
    return line


def update_sheet(line, chain_maps):
    """SHEET: init (21, 22-25, 26), end (32, 33-36, 37),
    cur H-bond (49, 50-53, 54), prev H-bond (64, 65-68, 69)."""
    line = pad_line(line, 38)
    line = remap_resid(line, 21, (22, 26), 26, chain_maps)
    line = remap_resid(line, 32, (33, 37), 37, chain_maps)
    if len(line.rstrip('\n')) > 54:
        line = pad_line(line, 70)
        line = remap_resid(line, 49, (50, 54), 54, chain_maps)
        line = remap_resid(line, 64, (65, 69), 69, chain_maps)
    return line


def update_ssbond(line, chain_maps):
    """SSBOND: res1 (15, 17-20, 21), res2 (29, 31-34, 35)."""
    line = pad_line(line, 36)
    line = remap_resid(line, 15, (17, 21), 21, chain_maps)
    line = remap_resid(line, 29, (31, 35), 35, chain_maps)
    return line


def update_link(line, chain_maps):
    """LINK: res1 (21, 22-25, 26), res2 (51, 52-55, 56)."""
    line = pad_line(line, 57)
    line = remap_resid(line, 21, (22, 26), 26, chain_maps)
    line = remap_resid(line, 51, (52, 56), 56, chain_maps)
    return line


def update_cispep(line, chain_maps):
    """CISPEP: res1 (15, 17-20, 21), res2 (29, 31-34, 35)."""
    line = pad_line(line, 36)
    line = remap_resid(line, 15, (17, 21), 21, chain_maps)
    line = remap_resid(line, 29, (31, 35), 35, chain_maps)
    return line


def update_het(line, chain_maps):
    """HET: (12, 13-16, 17)."""
    line = pad_line(line, 18)
    line = remap_resid(line, 12, (13, 17), 17, chain_maps)
    return line


def update_dbref(line, chain_maps):
    """DBREF: chain at 12, seqBegin (14-17, 18), seqEnd (20-23, 24)."""
    line = pad_line(line, 25)
    line = remap_resid(line, 12, (14, 18), 18, chain_maps)
    line = remap_resid(line, 12, (20, 24), 24, chain_maps)
    return line


def update_seqadv(line, chain_maps):
    """SEQADV: (16, 18-21, 22)."""
    line = pad_line(line, 23)
    line = remap_resid(line, 16, (18, 22), 22, chain_maps)
    return line


def update_remark_465_610(line, chain_maps):
    """REMARK 465/610 data lines: resName(15-17) chain(19) seqNum(21-24) iCode(25)."""
    if len(line.rstrip('\n')) < 25:
        return line
    chain = line[19]
    seq_str = line[21:25].strip()
    if chain not in chain_maps or not seq_str or not seq_str.lstrip('-').isdigit():
        return line
    line = pad_line(line, 26)
    line = remap_resid(line, 19, (21, 25), 25, chain_maps)
    return line


def is_remark500_residue_data(line):
    """Check if a REMARK 500 line has M RES CSSEQI data format at cols 13-23."""
    if len(line) < 24:
        return False
    resname = line[14:17]
    if not resname.strip() or not all(c.isalpha() or c == ' ' for c in resname):
        return False
    chain = line[18]
    if chain not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz':
        return False
    seq_str = line[19:23].strip()
    if not seq_str or not seq_str.lstrip('-').isdigit():
        return False
    return True


def update_remark_500(line, chain_maps):
    """REMARK 500 data lines with M RES CSSEQI format: chain(18) seqNum(19-22) iCode(23)."""
    if not is_remark500_residue_data(line):
        return line
    line = pad_line(line, 24)
    line = remap_resid(line, 18, (19, 23), 23, chain_maps)
    return line


def update_conect(line, serial_map):
    """CONECT: remap atom serial numbers. Each serial is 5 chars starting at col 6."""
    parts = []
    parts.append(line[:6])  # "CONECT"
    i = 6
    while i + 5 <= len(line.rstrip('\n')):
        s = line[i:i + 5].strip()
        if s and s.isdigit():
            old_serial = int(s)
            new_serial = serial_map.get(old_serial, old_serial)
            parts.append(f"{new_serial:5d}")
        else:
            parts.append(line[i:i + 5])
        i += 5
    return ''.join(parts) + '\n'


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_stem(input_path.stem + "_renum")

    if args.rename:
        from dvbfixer.rename import canonicalize_pdb
        import tempfile as _tf
        _tmp = Path(_tf.mktemp(suffix='.pdb'))
        n = canonicalize_pdb(input_path, _tmp, args.verbose)
        if n > 0:
            print(f"Canonicalized {n} non-canonical residue(s)")
            input_path = _tmp
        elif _tmp.exists():
            _tmp.unlink()

    with open(input_path) as f:
        lines = f.readlines()

    seqres = parse_seqres(lines)
    if not seqres:
        print("No SEQRES records — renumbering sequentially (resolving insertion codes)")

    # Collect chain IDs from ATOM records (preserving order)
    atom_chains = []
    seen_chains = set()
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            c = line[21]
            if c not in seen_chains:
                seen_chains.add(c)
                atom_chains.append(c)

    # Parse per-chain scheme overrides (--chain-scheme H:kabat,L:chothia)
    chain_scheme_override = {}
    for spec in args.chain_scheme:
        if ':' not in spec:
            print(f"Error: invalid --chain-scheme format '{spec}' "
                  f"(expected CHAIN:SCHEME)", file=sys.stderr)
            sys.exit(1)
        ch, sch = spec.split(':', 1)
        sch = sch.lower()
        if sch not in {"seqres", "kabat", "chothia", "imgt", "martin", "eu", "aho"}:
            print(f"Error: unknown scheme '{sch}' for chain {ch}", file=sys.stderr)
            sys.exit(1)
        chain_scheme_override[ch] = sch

    # Build renumbering map per chain: {(old_resSeq, old_iCode): new_resSeq}
    # or {(old_resSeq, old_iCode): (new_resSeq, new_iCode)} for antibody path.
    chain_maps = {}

    for chain in atom_chains:
        atom_res = get_atom_residues(lines, chain)

        # Pick scheme for this chain: per-chain override > global > seqres
        scheme = chain_scheme_override.get(chain, args.scheme)

        # Try antibody numbering when an antibody scheme is selected
        if scheme != "seqres":
            from dvbfixer.antibody import number_chain_to_mapping
            if args.verbose:
                print(f"Chain {chain}: scheme={scheme}")
            ab_map, info = number_chain_to_mapping(atom_res, scheme, verbose=args.verbose)
            if info["is_antibody"] and ab_map:
                # Merge with SEQRES fallback for any residues not placed
                placed_keys = set(ab_map.keys())
                fallback = {}
                if chain in seqres:
                    seqres_map = align_to_seqres(atom_res, seqres[chain])
                    next_num = len(seqres[chain]) + 1
                    for key in seqres_map:
                        if seqres_map[key] is None:
                            seqres_map[key] = next_num
                            next_num += 1
                    fallback = seqres_map
                else:
                    for i, (rs, ic, _) in enumerate(atom_res, 1):
                        fallback[(rs, ic)] = i
                # Antibody wins where it placed something; SEQRES fills the rest
                merged = {}
                for key in fallback:
                    if key in placed_keys:
                        merged[key] = ab_map[key]
                    else:
                        merged[key] = fallback[key]
                chain_maps[chain] = merged
                if args.verbose and info["warnings"]:
                    for w in info["warnings"]:
                        print(f"  warning: {w}")
                continue
            else:
                if args.verbose:
                    print(f"  no antibody domain detected — falling back to SEQRES")

        if chain in seqres:
            mapping = align_to_seqres(atom_res, seqres[chain])

            next_num = len(seqres[chain]) + 1
            for key in mapping:
                if mapping[key] is None:
                    mapping[key] = next_num
                    next_num += 1

            chain_maps[chain] = mapping

            if args.verbose:
                matched = [(k, v) for k, v in mapping.items()
                           if v <= len(seqres[chain])]
                non_seq = [(k, v) for k, v in mapping.items()
                           if v > len(seqres[chain])]
                print(f"Chain {chain}: {len(seqres[chain])} SEQRES residues, "
                      f"{len(matched)} ATOM matched, {len(non_seq)} non-SEQRES")

                insertions = [(k, v) for k, v in mapping.items()
                              if k[1] != ' ' and v <= len(seqres[chain])]
                if insertions:
                    print(f"  Resolved insertion codes:")
                    for (old_seq, old_ic), new_seq in insertions:
                        rn = next(rn for s, ic, rn in atom_res if s == old_seq and ic == old_ic)
                        print(f"    {rn} {old_seq}{old_ic.strip()} -> {new_seq}")

                occupied = set(v for v in mapping.values() if v <= len(seqres[chain]))
                gaps = []
                in_gap = False
                gap_start = None
                for pos in range(1, len(seqres[chain]) + 1):
                    if pos not in occupied:
                        if not in_gap:
                            gap_start = pos
                            in_gap = True
                    else:
                        if in_gap:
                            gaps.append((gap_start, pos - 1))
                            in_gap = False
                if in_gap:
                    gaps.append((gap_start, len(seqres[chain])))
                if gaps:
                    print(f"  Gaps (missing ATOM residues):")
                    for gs, ge in gaps:
                        resnames = seqres[chain][gs - 1:ge]
                        print(f"    positions {gs}-{ge} ({ge - gs + 1} residues): "
                              f"{' '.join(resnames)}")
        else:
            mapping = {}
            for i, (resseq, icode, _) in enumerate(atom_res, 1):
                mapping[(resseq, icode)] = i
            chain_maps[chain] = mapping
            if args.verbose:
                print(f"Chain {chain}: no SEQRES, renumbered {len(atom_res)} residues sequentially")

    # Also populate chain_maps with REMARK 465 missing residues (they appear in
    # non-ATOM sections but need remapping too)
    for line in lines:
        for prefix in ('REMARK 465', 'REMARK 610'):
            if line.startswith(prefix) and len(line) > 25:
                chain = line[19]
                seq_str = line[21:25].strip()
                if chain in chain_maps and seq_str and seq_str.lstrip('-').isdigit():
                    old_key = (int(seq_str), line[25] if len(line) > 25 else ' ')
                    if old_key not in chain_maps[chain]:
                        if chain in seqres:
                            resname = line[15:18].strip()
                            for si, sr in enumerate(seqres[chain]):
                                if sr == resname and (si + 1) not in set(chain_maps[chain].values()):
                                    chain_maps[chain][old_key] = si + 1
                                    break

    # ------------------------------------------------------------------
    # Pass 1: rewrite ATOM/HETATM/TER, build old->new serial map
    # ------------------------------------------------------------------
    output_lines = []
    serial = 0
    serial_map = {}  # old_serial -> new_serial

    for line in lines:
        rec = line[:6].strip()

        if rec in ("ATOM", "HETATM"):
            old_serial = int(line[6:11].strip()) if line[6:11].strip().isdigit() else None
            chain = line[21]
            old_key = (int(line[22:26].strip()), line[26])
            serial += 1
            if old_serial is not None:
                serial_map[old_serial] = serial
            if chain in chain_maps and old_key in chain_maps[chain]:
                val = chain_maps[chain][old_key]
                if isinstance(val, tuple):
                    new_resseq, new_icode = val
                    if not new_icode:
                        new_icode = " "
                else:
                    new_resseq, new_icode = val, " "
                new_line = (line[:6] + f"{serial:5d}" + line[11:22]
                            + f"{new_resseq:4d}" + new_icode + line[27:])
            else:
                new_line = line[:6] + f"{serial:5d}" + line[11:]
            output_lines.append(new_line)

        elif rec == "TER":
            serial += 1
            old_serial = int(line[6:11].strip()) if len(line) > 10 and line[6:11].strip().isdigit() else None
            if old_serial is not None:
                serial_map[old_serial] = serial
            if len(line) > 26:
                chain = line[21]
                seq_str = line[22:26].strip()
                icode = line[26] if len(line) > 26 else ' '
                if seq_str and seq_str.isdigit() and chain in chain_maps:
                    old_key = (int(seq_str), icode)
                    if old_key in chain_maps[chain]:
                        val = chain_maps[chain][old_key]
                        if isinstance(val, tuple):
                            new_seq, new_icode = val
                            if not new_icode:
                                new_icode = " "
                        else:
                            new_seq, new_icode = val, " "
                        resname = line[17:20]
                        output_lines.append(
                            f"TER   {serial:>5d}      {resname} {chain}{new_seq:>4d}{new_icode}\n"
                        )
                        continue
            output_lines.append(f"TER   {serial:>5d}" + line[11:])

        elif rec == "HELIX":
            output_lines.append(update_helix(line, chain_maps))
        elif rec == "SHEET":
            output_lines.append(update_sheet(line, chain_maps))
        elif rec == "SSBOND":
            output_lines.append(update_ssbond(line, chain_maps))
        elif rec == "LINK":
            output_lines.append(update_link(line, chain_maps))
        elif rec == "CISPEP":
            output_lines.append(update_cispep(line, chain_maps))
        elif rec == "HET":
            output_lines.append(update_het(line, chain_maps))
        elif rec == "DBREF":
            output_lines.append(update_dbref(line, chain_maps))
        elif rec == "SEQADV":
            output_lines.append(update_seqadv(line, chain_maps))
        elif rec == "CONECT":
            output_lines.append(update_conect(line, serial_map))
        elif line.startswith("REMARK 465") or line.startswith("REMARK 610"):
            output_lines.append(update_remark_465_610(line, chain_maps))
        elif line.startswith("REMARK 500"):
            output_lines.append(update_remark_500(line, chain_maps))
        else:
            output_lines.append(line)

    if not args.keep_water:
        filtered = []
        for line in output_lines:
            if line.startswith(("ATOM  ", "HETATM", "ANISOU")):
                resname = line[17:20].strip()
                if resname in WATER_RESNAMES:
                    continue
            elif line.startswith("TER") and len(line) > 20:
                resname = line[17:20].strip()
                if resname in WATER_RESNAMES:
                    continue
            filtered.append(line)
        output_lines = filtered

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    print(f"Wrote {output_path}")
