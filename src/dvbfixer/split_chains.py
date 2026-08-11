"""Split PDB/GRO chains empirically using C-N distance, residue numbering, and nearest-atom gaps."""

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
ION_RESNAMES = {'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU',
                'NA+', 'CL-', 'K+', 'MG2+', 'CA2+', 'ZN2+', 'BUF', 'BUFF'}
SOLVENT_IONS = WATER_RESNAMES | ION_RESNAMES

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Common protonation/tautomer variants
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM",
    "ASH", "GLH", "LYN", "ASPP", "GLUP",
    # GLYCAM glycoprotein residues (still part of the protein backbone)
    "NLN", "OLS", "OLT",
    # Terminal patches
    "ACE", "NME", "NH2",
    # Selenomethionine and other commonly-treated-as-protein residues
    "MSE",
}

# Default threshold: when the total number of detected chains exceeds this,
# small-molecule chains (ions, ligands, lipids, etc.) no longer get unique
# chain IDs — only protein chains do. Overridable via --max-chains.
DEFAULT_MAX_CHAIN_THRESHOLD = 26


def parse_gro(path):
    """Convert a GROMACS .gro file to PDB lines via MDAnalysis.

    Preserves all residue names exactly as-is (GLUP, ASPP, etc.).
    """
    import tempfile
    import warnings

    import MDAnalysis as mda

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(path))
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            tmp_path = tmp.name
        u.atoms.write(tmp_path)

    with open(tmp_path) as f:
        lines = f.readlines()
    Path(tmp_path).unlink()
    return lines


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer split",
        description="Split chains empirically, or extract PDB biological assemblies "
        "from REMARK 350 BIOMT records."
    )
    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB or GRO file")
    io.add_argument("-o", "--output", help="Output PDB file (default: <input>_split.pdb)")
    io.add_argument(
        "--assembly", metavar="ID|all",
        help="Extract one REMARK 350 biological assembly, or all assemblies. "
             "PDB input only; empirical splitting remains the default."
    )

    detection = p.add_argument_group("Chain-break detection")
    detection.add_argument(
        "-d", "--distance-cutoff", type=float, default=DEFAULT_DISTANCE_CUTOFF,
        help=f"C->N peptide bond cutoff in angstroms (default: {DEFAULT_DISTANCE_CUTOFF})"
    )
    detection.add_argument(
        "-g", "--gap-cutoff", type=float, default=DEFAULT_GAP_CUTOFF,
        help=f"Min nearest-atom distance between consecutive residues to call a break "
             f"when C/N atoms are missing (default: {DEFAULT_GAP_CUTOFF} A)"
    )
    detection.add_argument(
        "--no-distance", action="store_true",
        help="Disable all distance-based detection, use only residue numbering"
    )
    detection.add_argument(
        "--max-chains", type=int, default=DEFAULT_MAX_CHAIN_THRESHOLD,
        help=f"When more than this many chains are detected, do NOT assign "
             f"chain IDs to small-molecule chains (ions, ligands, lipids, "
             f"single-residue HETATMs) — they stay with blank chain ID. "
             f"Protein chains always get chain IDs. "
             f"Default: {DEFAULT_MAX_CHAIN_THRESHOLD}. Set higher to keep "
             f"the original assign-all behaviour (max {len(CHAIN_IDS)})."
    )

    content = p.add_argument_group("Content / renumbering")
    numbering = content.add_mutually_exclusive_group()
    numbering.add_argument(
        "--renumber", action="store_true", default=None,
        help="Renumber residues per chain (default in empirical mode)"
    )
    numbering.add_argument(
        "--no-renumber", action="store_false", dest="renumber",
        help="Keep original residue numbers (default in assembly mode)"
    )
    content.add_argument(
        "--keep-water", action="store_true",
        help="Keep water molecules (HOH, WAT, TIP3, SOL) in output (default: remove)"
    )

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print detected chain info"
    )

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p, batch=True)
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
    # Cols 17-19 standard; col 20 may hold 4th char (e.g. GLUP, ASPP from GRO)
    name = line[17:20].strip()
    if len(line) > 20 and line[20] not in (' ', '\n'):
        name = line[17:21].strip()
    return name


def set_chain_id(line, chain_id):
    return line[:21] + chain_id + line[22:]


def set_resid(line, resid, icode=" "):
    """Set residue sequence number and insertion code (col 22-26 + col 26)."""
    return line[:22] + f"{resid:4d}" + icode + line[27:]


def make_ter(serial, resname, chain_id, resid):
    return f"TER   {serial:>5d}      {resname:<3s} {chain_id}{resid:>4d}\n"


def find_model_blocks(lines):
    """Return a list of (model_label, atom_line_idx_range) for each MODEL.

    Each entry is ((model_start_line, model_end_line), model_label_or_None).
    `model_start_line` is the index of MODEL (or 0 if none), `model_end_line`
    is the index of ENDMDL (or len(lines) if none). Returns None if no MODEL
    records are present (single-block case).
    """
    blocks = []
    current_start = None
    current_label = None
    for i, line in enumerate(lines):
        if line.startswith("MODEL"):
            current_start = i
            current_label = line[10:14].strip() or str(len(blocks) + 1)
        elif line.startswith("ENDMDL"):
            if current_start is not None:
                blocks.append((current_start, i, current_label))
                current_start = None
                current_label = None
    return blocks if blocks else None


def _classify_chains(atom_lines, breaks):
    """Return per-chain classification: 'protein' or 'small_molecule'.

    A chain is 'protein' when its residue set is dominated (>= 50%) by
    standard amino acids. Glycan trees, single-ion chains, lipid chains,
    and arbitrary ligands fall into 'small_molecule'.
    """
    kinds = []
    n = len(atom_lines)
    for ci in range(len(breaks)):
        s = breaks[ci]
        e = breaks[ci + 1] if ci + 1 < len(breaks) else n
        residues = {}
        for j in range(s, e):
            rid = get_resid(atom_lines[j])
            if rid not in residues:
                residues[rid] = get_resname(atom_lines[j])
        if not residues:
            kinds.append('small_molecule')
            continue
        n_std = sum(1 for rn in residues.values() if rn in STANDARD_RESIDUES)
        if n_std / len(residues) >= 0.5:
            kinds.append('protein')
        else:
            kinds.append('small_molecule')
    return kinds


def _assign_chain_ids(kinds, max_chain_threshold):
    """Return a list of chain IDs per chain (one entry per `kinds`).

    When `len(kinds) > max_chain_threshold`, protein chains are assigned
    chain IDs from CHAIN_IDS in order and small-molecule chains keep an
    empty/blank chain ID (' ' written into PDB column 22). Otherwise every
    chain gets a chain ID as before.
    """
    if len(kinds) <= max_chain_threshold:
        return [CHAIN_IDS[i] for i in range(len(kinds))]
    assigned = []
    cursor = 0
    for k in kinds:
        if k == 'protein':
            if cursor >= len(CHAIN_IDS):
                # Run out even of protein chain IDs — emit blank and warn
                # downstream; caller already validates protein count.
                assigned.append(' ')
            else:
                assigned.append(CHAIN_IDS[cursor])
                cursor += 1
        else:
            assigned.append(' ')
    return assigned


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

        # Criterion 2: C->N peptide bond distance (any residues with backbone C/N)
        if c_prev is not None and n_cur is not None:
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

    is_gro = input_path.suffix.lower() == '.gro'
    if args.assembly:
        if is_gro:
            print("--assembly requires a PDB input with REMARK 350 records.", file=sys.stderr)
            sys.exit(1)
        from dvbfixer.biological_assembly import (
            AssemblyError,
            assembly_output_path,
            parse_biological_assemblies,
            render_biological_assembly,
        )

        try:
            lines = input_path.read_text().splitlines(keepends=True)
            assemblies = parse_biological_assemblies(lines)
            all_mode = args.assembly.lower() == "all"
            identifiers = list(assemblies) if all_mode else [args.assembly]
            unknown = [identifier for identifier in identifiers if identifier not in assemblies]
            if unknown:
                available = ", ".join(assemblies)
                raise AssemblyError(
                    f"unknown biological assembly {unknown[0]!r}; available: {available}"
                )
            rendered = []
            for identifier in identifiers:
                output_path = assembly_output_path(
                    input_path, args.output, identifier, all_mode
                )
                output_lines = render_biological_assembly(
                    lines, assemblies[identifier], CHAIN_IDS,
                    renumber=bool(args.renumber),
                )
                rendered.append((output_path, output_lines))
        except (AssemblyError, OSError, ValueError) as exc:
            print(f"Biological assembly error: {exc}", file=sys.stderr)
            sys.exit(1)
        for output_path, output_lines in rendered:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text("".join(output_lines))
            print(f"Wrote biological assembly to {output_path}")
        return

    if args.output:
        output_path = Path(args.output)
    else:
        # Always output PDB (GRO has no chain IDs)
        output_path = input_path.with_stem(input_path.stem + "_split").with_suffix(".pdb")
    use_distance = not args.no_distance
    renumber = True if args.renumber is None else args.renumber

    if is_gro:
        lines = parse_gro(input_path)
    else:
        with open(input_path) as f:
            lines = f.readlines()

    model_blocks = find_model_blocks(lines)
    if model_blocks and len(model_blocks) > 1:
        _process_multi_model(
            lines, model_blocks, output_path, args, use_distance, renumber, is_gro
        )
        return

    # Single-block path (no MODEL records, or just one MODEL).
    # If there's exactly one MODEL, restrict atom scanning to inside it.
    block_range = None
    if model_blocks and len(model_blocks) == 1:
        m_start, m_end, _label = model_blocks[0]
        block_range = (m_start + 1, m_end)

    atom_lines, atom_orig_indices, solvent_lines = _collect_atoms(lines, block_range)

    if not atom_lines:
        print("No ATOM/HETATM records found.", file=sys.stderr)
        sys.exit(1)

    breaks = find_chain_breaks(atom_lines, args.distance_cutoff, args.gap_cutoff, use_distance)
    n_chains = len(breaks)

    kinds = _classify_chains(atom_lines, breaks)
    n_protein = kinds.count('protein')
    if n_protein > len(CHAIN_IDS):
        print(f"Too many protein chains ({n_protein}) for available chain IDs "
              f"({len(CHAIN_IDS)}).", file=sys.stderr)
        sys.exit(1)

    chain_id_assignment = _assign_chain_ids(kinds, args.max_chains)
    if args.verbose and n_chains > args.max_chains:
        n_sm_blank = sum(1 for cid in chain_id_assignment if cid == ' ')
        print(f"Detected {n_chains} chain(s) > --max-chains={args.max_chains}: "
              f"{n_protein} protein chain(s) get IDs; "
              f"{n_sm_blank} small-molecule chain(s) keep blank chain ID.")

    chain_for_atom = _build_chain_for_atom(breaks, len(atom_lines))
    resid_maps = _build_resid_maps(atom_lines, breaks, n_chains) if renumber else {}

    if args.verbose:
        _print_chain_summary(atom_lines, breaks, n_chains, resid_maps, renumber,
                              chain_ids=chain_id_assignment)

    break_set = set(breaks)
    output_lines = _render_block(
        lines, atom_lines, atom_orig_indices, chain_for_atom, break_set,
        resid_maps, renumber, start_serial=1,
        chain_id_assignment=chain_id_assignment,
    )

    if args.keep_water and solvent_lines:
        end_idx = len(output_lines)
        for i, line in enumerate(output_lines):
            if line.startswith("END"):
                end_idx = i
                break
        for sl in solvent_lines:
            output_lines.insert(end_idx, sl)
            end_idx += 1

    if is_gro and (not output_lines or not output_lines[-1].startswith("END")):
        output_lines.append("END\n")

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    print(f"Wrote {output_path} with {n_chains} chain(s)")


def _collect_atoms(lines, block_range=None):
    """Pull out ATOM/HETATM lines from a region of `lines`, separating solvent.

    Returns (atom_lines, atom_orig_indices, solvent_lines). When `block_range`
    is given, only lines in [start, end) are scanned; otherwise all lines.
    """
    start, end = (0, len(lines)) if block_range is None else block_range
    atom_lines = []
    atom_orig_indices = []
    solvent_lines = []
    for i in range(start, end):
        line = lines[i]
        if is_atom_line(line):
            if get_resname(line) in SOLVENT_IONS:
                solvent_lines.append(line)
            else:
                atom_lines.append(line)
                atom_orig_indices.append(i)
    return atom_lines, atom_orig_indices, solvent_lines


def _build_chain_for_atom(breaks, n_atoms):
    """Map each atom-list index to its chain index."""
    chain_for_atom = np.zeros(n_atoms, dtype=int)
    for ci in range(len(breaks)):
        s = breaks[ci]
        e = breaks[ci + 1] if ci + 1 < len(breaks) else n_atoms
        chain_for_atom[s:e] = ci
    return chain_for_atom


def _build_resid_maps(atom_lines, breaks, n_chains):
    """Build per-chain (resSeq, iCode) → new-resnum maps starting at 1 per chain."""
    resid_maps = {}
    for ci in range(n_chains):
        start = breaks[ci]
        end = breaks[ci + 1] if ci + 1 < n_chains else len(atom_lines)
        seen = {}
        counter = 1
        for j in range(start, end):
            rid = get_resid(atom_lines[j])
            if rid not in seen:
                seen[rid] = counter
                counter += 1
        resid_maps[ci] = seen
    return resid_maps


def _chain_signature(atom_lines, breaks):
    """Compact structural fingerprint of a model's chain layout for comparison
    across MODELs. Returns a tuple of (n_atoms, n_residues, first_resname,
    last_resname) per chain. Two MODELs are considered the same complex when
    their signatures are equal.
    """
    sig = []
    n = len(atom_lines)
    for ci in range(len(breaks)):
        start = breaks[ci]
        end = breaks[ci + 1] if ci + 1 < len(breaks) else n
        residues = set()
        for j in range(start, end):
            residues.add(get_resid(atom_lines[j]))
        sig.append((
            end - start,
            len(residues),
            get_resname(atom_lines[start]),
            get_resname(atom_lines[end - 1]),
        ))
    return tuple(sig)


def _print_chain_summary(atom_lines, breaks, n_chains, resid_maps, renumber,
                          chain_ids=CHAIN_IDS, prefix=""):
    print(f"{prefix}Detected {n_chains} chain(s):")
    for ci in range(n_chains):
        start = breaks[ci]
        end = (breaks[ci + 1] if ci + 1 < n_chains else len(atom_lines)) - 1
        n_atoms = end - start + 1
        if renumber and ci in resid_maps:
            n_res = len(resid_maps[ci])
        else:
            n_res = len({get_resid(atom_lines[j]) for j in range(start, end + 1)})
        print(
            f"{prefix}  Chain {chain_ids[ci]}: {n_atoms:>5d} atoms, {n_res:>4d} residues, "
            f"{get_resname(atom_lines[start])} {get_resseq(atom_lines[start])} - "
            f"{get_resname(atom_lines[end])} {get_resseq(atom_lines[end])}"
        )


def _render_block(lines, atom_lines, atom_orig_indices, chain_for_atom,
                   break_set, resid_maps, renumber, start_serial,
                   chain_id_assignment, line_range=None):
    """Emit output lines for a single block (whole file or one MODEL).

    `chain_id_assignment` is a sequence indexed by chain index; for single-block
    operation pass CHAIN_IDS itself. For multi-MODEL operation, pass the
    consensus chain ID list (so every MODEL uses the same IDs).
    `line_range` restricts which slice of `lines` to iterate.
    """
    start, end = (0, len(lines)) if line_range is None else line_range
    atom_orig_set = set(atom_orig_indices)
    # Map original line index → atom_lines index
    orig_to_ai = {idx: ai for ai, idx in enumerate(atom_orig_indices)}

    output_lines = []
    serial = start_serial
    for i in range(start, end):
        line = lines[i]
        if is_atom_line(line) and get_resname(line) in SOLVENT_IONS:
            continue
        if is_atom_line(line) and i in atom_orig_set:
            ai = orig_to_ai[i]
            ci = int(chain_for_atom[ai])
            chain_id = chain_id_assignment[ci]
            new_line = set_chain_id(line, chain_id)
            if renumber:
                rid = get_resid(line)
                new_line = set_resid(new_line, resid_maps[ci][rid], icode=" ")
            new_line = new_line[:6] + f"{serial:5d}" + new_line[11:]
            output_lines.append(new_line)
            next_ai = ai + 1
            # `break_set` only ever holds chain-START indices, so it never
            # contains `len(atom_orig_indices)` — the position right after
            # the very LAST atom of this block. Without this explicit
            # check, the final chain in the output never got a closing
            # TER record.
            if next_ai in break_set or next_ai == len(atom_orig_indices):
                serial += 1
                resname = get_resname(new_line)
                resid = int(new_line[22:26].strip())
                output_lines.append(make_ter(serial, resname, chain_id, resid))
            serial += 1
        elif line.startswith("TER"):
            continue
        else:
            output_lines.append(line)
    return output_lines


def _process_multi_model(lines, model_blocks, output_path, args, use_distance,
                          renumber, is_gro):
    """Process an input with multiple MODEL records.

    Each MODEL is split independently (so chain breaks within each state are
    detected correctly), then the per-MODEL chain signatures are compared. If
    all MODELs share the same signature (same complex sampled at different
    states), every MODEL is rewritten with the SAME chain IDs (A, B, C, …),
    matching the expected behaviour of an MD trajectory. If MODELs differ
    structurally, each one falls back to its own independent chain IDs and a
    warning is emitted.
    """
    n_models = len(model_blocks)
    per_model = []  # one entry per MODEL: dict of analysis results
    for (m_start, m_end, label) in model_blocks:
        atom_lines, atom_orig_indices, solvent_lines = _collect_atoms(
            lines, block_range=(m_start + 1, m_end)
        )
        if not atom_lines:
            print(f"MODEL {label}: no ATOM/HETATM records.", file=sys.stderr)
            sys.exit(1)
        breaks = find_chain_breaks(
            atom_lines, args.distance_cutoff, args.gap_cutoff, use_distance
        )
        per_model.append({
            "label": label,
            "m_start": m_start,
            "m_end": m_end,
            "atom_lines": atom_lines,
            "atom_orig_indices": atom_orig_indices,
            "solvent_lines": solvent_lines,
            "breaks": breaks,
            "signature": _chain_signature(atom_lines, breaks),
        })

    signatures = [m["signature"] for m in per_model]
    all_same = all(s == signatures[0] for s in signatures)

    if all_same:
        n_chains = len(per_model[0]["breaks"])
        # Classify chains in MODEL 1 (signatures match, so other MODELs share
        # the same protein/small-molecule layout) and pick chain IDs honouring
        # --max-chains.
        kinds = _classify_chains(
            per_model[0]["atom_lines"], per_model[0]["breaks"]
        )
        n_protein = kinds.count('protein')
        if n_protein > len(CHAIN_IDS):
            print(f"Too many protein chains ({n_protein}) for available chain "
                  f"IDs ({len(CHAIN_IDS)}).", file=sys.stderr)
            sys.exit(1)
        assignment = _assign_chain_ids(kinds, args.max_chains)
        chain_id_lists = [assignment] * n_models
        if args.verbose:
            print(f"Detected {n_models} MODEL records, all with identical "
                  f"chain layout — reusing chain IDs across MODELs.")
            if n_chains > args.max_chains:
                n_sm_blank = sum(1 for cid in assignment if cid == ' ')
                print(f"  {n_chains} chain(s) > --max-chains={args.max_chains}: "
                      f"{n_protein} protein chain(s) get IDs; "
                      f"{n_sm_blank} small-molecule chain(s) keep blank chain ID.")
            _print_chain_summary(
                per_model[0]["atom_lines"], per_model[0]["breaks"], n_chains,
                _build_resid_maps(per_model[0]["atom_lines"],
                                   per_model[0]["breaks"], n_chains) if renumber else {},
                renumber, chain_ids=assignment,
            )
    else:
        print(f"Warning: {n_models} MODEL records have differing chain layouts; "
              f"assigning independent chain IDs per MODEL.", file=sys.stderr)
        chain_id_lists = []
        offset = 0
        for m in per_model:
            kinds_m = _classify_chains(m["atom_lines"], m["breaks"])
            n_c = len(kinds_m)
            n_prot_m = kinds_m.count('protein')
            if n_c > args.max_chains:
                # Only assign IDs to protein chains; small molecules get blank.
                assignment = []
                cursor = 0
                for k in kinds_m:
                    if k == 'protein':
                        if offset + cursor >= len(CHAIN_IDS):
                            assignment.append(' ')
                        else:
                            assignment.append(CHAIN_IDS[offset + cursor])
                            cursor += 1
                    else:
                        assignment.append(' ')
                if offset + n_prot_m > len(CHAIN_IDS):
                    print("Too many protein chains across MODELs for "
                          "available IDs.", file=sys.stderr)
                    sys.exit(1)
                offset += n_prot_m
                chain_id_lists.append(assignment)
            else:
                ids = CHAIN_IDS[offset:offset + n_c]
                if len(ids) < n_c:
                    print("Too many total chains across MODELs for available IDs.",
                          file=sys.stderr)
                    sys.exit(1)
                chain_id_lists.append(ids)
                offset += n_c

    # Render output. Walk the original `lines` and emit each segment:
    # - lines before MODEL 1 (header / SEQRES / etc.)
    # - for each MODEL: MODEL header + atom block (with chain IDs from
    #   `chain_id_lists[i]`) + ENDMDL
    # - lines after final ENDMDL (footer)
    output_lines = []
    prev_end = 0
    all_solvent = []
    for i, m in enumerate(per_model):
        # Pre-MODEL lines (header for MODEL 1, otherwise inter-MODEL lines)
        output_lines.extend(lines[prev_end:m["m_start"] + 1])  # include MODEL line itself
        n_c = len(m["breaks"])
        chain_for_atom = _build_chain_for_atom(m["breaks"], len(m["atom_lines"]))
        resid_maps = _build_resid_maps(
            m["atom_lines"], m["breaks"], n_c
        ) if renumber else {}
        break_set = set(m["breaks"])
        block_out = _render_block(
            lines, m["atom_lines"], m["atom_orig_indices"], chain_for_atom,
            break_set, resid_maps, renumber, start_serial=1,
            chain_id_assignment=chain_id_lists[i],
            line_range=(m["m_start"] + 1, m["m_end"]),
        )
        output_lines.extend(block_out)
        output_lines.append(lines[m["m_end"]])  # ENDMDL line
        prev_end = m["m_end"] + 1
        all_solvent.append(m["solvent_lines"])

    # Trailing lines after last ENDMDL
    output_lines.extend(lines[prev_end:])

    # Re-append solvent only if --keep-water (one copy per MODEL would change
    # atom counts; for simplicity insert each MODEL's solvent inside its block
    # — but this requires re-walking. Skipping that complexity: append all
    # solvent before END once, matching single-block behaviour.)
    if args.keep_water and any(all_solvent):
        end_idx = len(output_lines)
        for i, line in enumerate(output_lines):
            if line.startswith("END") and not line.startswith("ENDMDL"):
                end_idx = i
                break
        # Use solvent from the first MODEL only (all MODELs should have the
        # same solvent in a same-complex trajectory).
        for sl in all_solvent[0]:
            output_lines.insert(end_idx, sl)
            end_idx += 1

    if is_gro and (not output_lines or not output_lines[-1].startswith("END")):
        output_lines.append("END\n")

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    n_chains_first = len(per_model[0]["breaks"])
    if all_same:
        print(f"Wrote {output_path} with {n_models} MODEL(s), "
              f"{n_chains_first} chain(s) per MODEL (consistent across MODELs)")
    else:
        total = sum(len(m["breaks"]) for m in per_model)
        print(f"Wrote {output_path} with {n_models} MODEL(s), "
              f"{total} total chain(s)")
