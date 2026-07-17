"""Sequence-aware Kabsch superposition of one PDB onto a reference.

Used inside `dvbfixer zbs` after every pipeline step (default ON, opt-out
via `--no-align-to-input`) so interim outputs stay in the same Cartesian
frame as the user's original input — no accumulated rigid-body drift from
successive OpenMM minimizations.

Only exposed as an internal helper; there's no standalone
`dvbfixer align` subcommand.

**v2 (2026-07)** — atom correspondences are established by GLOBAL
protein-sequence alignment per chain (via `Bio.Align.PairwiseAligner`).
The previous version keyed on `(chain, resseq, icode, atomname)` and
broke the moment `renumber` rewrote residue numbers: nearly every key
mismatched, Kabsch fit to random-atom pairs, and RMSD floored around
~6 Å on `test/gaff_test/1VCU.pdb` even for coord-identical steps.

Sequence-based pairing (what PyMOL's `align` and MDAnalysis's
`fasta2select` do) is dependency-free (biopython + MDAnalysis already in
env) and robust to `renumber`, model loop insertion, protonation renames
(HIS→HIE/HID/HIP, etc.), and altLoc removal. The Kabsch fit itself is a
plain SVD on the matched-atom coord pairs — that's the structural part.
"""

from __future__ import annotations

from pathlib import Path

_BACKBONE_ATOMS = ('N', 'CA', 'C', 'O')

# AMBER + CHARMM protonation-variant → canonical parent, so
# `.residues.sequence()` returns matching letters across renames.
_CANONICAL = {
    'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
    'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
    'ASH': 'ASP', 'ASPP': 'ASP',
    'GLH': 'GLU', 'GLUP': 'GLU',
    'CYX': 'CYS', 'CYM': 'CYS',
    'LYN': 'LYS', 'LSN': 'LYS',
    'MSE': 'MET',
}


def _canonical_sequence(mda_residues):
    """Return one-letter sequence + list of matching MDA Residue objects.

    Non-protein residues are skipped. AMBER/CHARMM protonation variants
    are folded to their parent so a HIS in the reference matches a HIE in
    the mobile.
    """
    from Bio.SeqUtils import IUPACData
    three_to_one = {k.upper(): v for k, v in
                    IUPACData.protein_letters_3to1_extended.items()}
    seq = []
    res_list = []
    for r in mda_residues:
        rn = _CANONICAL.get(r.resname, r.resname)
        letter = three_to_one.get(rn)
        if letter is None:
            continue
        seq.append(letter)
        res_list.append(r)
    return ''.join(seq), res_list


def _kabsch(P, Q):
    """Return (R, t) that best superposes P onto Q. P, Q: (N, 3) arrays."""
    import numpy as np
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    H = (P - Pc).T @ (Q - Qc)
    U, _S, Vt = np.linalg.svd(H)
    d = float(1.0 if np.linalg.det(Vt.T @ U.T) > 0 else -1.0)
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Qc - Pc @ R.T
    return R, t


def _fit_atoms_for_selection(res_pairs, selection):
    """From matched residue pairs, harvest coord pairs for the Kabsch fit.

    Returns two (N, 3) numpy arrays P, Q (mobile then reference).
    """
    import numpy as np
    sel = selection.lower()
    if sel == 'ca':
        want = {'CA'}
    elif sel == 'backbone':
        want = set(_BACKBONE_ATOMS)
    elif sel == 'heavy':
        want = None  # keep every non-H atom whose name matches on both sides
    elif sel == 'all':
        want = None  # keep every atom whose name matches on both sides
    else:
        raise ValueError(f"unknown selection '{selection}'")

    P, Q = [], []
    for m_res, r_res in res_pairs:
        m_atoms = {a.name: a for a in m_res.atoms}
        r_atoms = {a.name: a for a in r_res.atoms}
        common = set(m_atoms) & set(r_atoms)
        if want is not None:
            common &= want
        elif sel == 'heavy':
            common = {n for n in common if not n.startswith('H')}
        # deterministic ordering — sorted() over the atom-name intersection
        for name in sorted(common):
            P.append(m_atoms[name].position)
            Q.append(r_atoms[name].position)
    if not P:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.asarray(P), np.asarray(Q)


def kabsch_align_pdb(input_pdb, reference_pdb, output_pdb, *,
                     selection='backbone', verbose=False):
    """Kabsch-superpose `input_pdb` onto `reference_pdb`, write `output_pdb`.

    Correspondences are established by per-chain protein-sequence
    alignment (`Bio.Align.PairwiseAligner`, global mode). The Kabsch
    rotation + translation is computed on the matched-residue atoms
    (backbone/ca/heavy/all — controlled by `selection`) and applied to
    EVERY atom in the input universe (protein + HETATMs + water etc.).

    Returns:
        (rmsd_before, rmsd_after, n_matched_atoms) — floats/int. On any
        failure (unreadable PDB, no protein overlap, biopython missing)
        the function copies the input to the output unchanged and returns
        (None, None, 0).
    """
    try:
        import MDAnalysis as mda
        import numpy as np
        from Bio.Align import PairwiseAligner
    except ImportError as e:
        if verbose:
            print(f"  [align] optional dependency missing ({e}); passing "
                  f"input through unchanged")
        _copy_file(input_pdb, output_pdb)
        return None, None, 0

    try:
        u_mobile = mda.Universe(str(input_pdb))
        u_ref = mda.Universe(str(reference_pdb))
    except Exception as e:
        if verbose:
            print(f"  [align] could not load PDBs ({e}); passing through")
        _copy_file(input_pdb, output_pdb)
        return None, None, 0

    aligner = PairwiseAligner(mode='global',
                              match_score=2, mismatch_score=-1,
                              open_gap_score=-10, extend_gap_score=-1)

    # Per-chain sequence pairing → flat list of matched (mobile_res, ref_res).
    matched_pairs = []
    common_chains = sorted(set(u_mobile.atoms.chainIDs)
                           & set(u_ref.atoms.chainIDs))
    for ch in common_chains:
        mob_ca = u_mobile.select_atoms(
            f"protein and name CA and chainID {ch}")
        ref_ca = u_ref.select_atoms(
            f"protein and name CA and chainID {ch}")
        if len(mob_ca) < 3 or len(ref_ca) < 3:
            continue
        mob_seq, mob_residues = _canonical_sequence(mob_ca.residues)
        ref_seq, ref_residues = _canonical_sequence(ref_ca.residues)
        if not mob_seq or not ref_seq:
            continue
        try:
            aln = aligner.align(mob_seq, ref_seq)[0]
        except Exception:
            continue
        # aln.aligned is (mobile_blocks, reference_blocks); each is an
        # (N, 2) int array of [start, end) index pairs.
        mob_blocks, ref_blocks = aln.aligned
        for (m0, m1), (r0, r1) in zip(mob_blocks, ref_blocks):
            for k in range(int(m1 - m0)):
                matched_pairs.append(
                    (mob_residues[int(m0) + k], ref_residues[int(r0) + k]))

    if len(matched_pairs) < 3:
        if verbose:
            print(f"  [align] not enough protein overlap "
                  f"({len(matched_pairs)} matched residues); passing "
                  f"through")
        _copy_file(input_pdb, output_pdb)
        return None, None, 0

    P, Q = _fit_atoms_for_selection(matched_pairs, selection)
    if len(P) < 3:
        if verbose:
            print(f"  [align] no atoms matched for selection "
                  f"'{selection}' across {len(matched_pairs)} residue "
                  f"pairs; passing through")
        _copy_file(input_pdb, output_pdb)
        return None, None, 0

    rmsd_before = float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))
    R, t = _kabsch(P, Q)

    # Post-fit RMSD on the same atom pairs.
    P_new = P @ R.T + t
    rmsd_after = float(np.sqrt(np.mean(np.sum((P_new - Q) ** 2, axis=1))))

    # Apply (R, t) to every ATOM/HETATM line in the INPUT FILE and pass
    # every other record (SEQRES, HELIX, SHEET, SSBOND, LINK, CISPEP,
    # HET, DBREF, SEQADV, CONECT, REMARK, HEADER, TITLE, MODEL/ENDMDL,
    # TER, END, CRYST1) through unchanged. MDAnalysis's PDB writer would
    # strip all headers, which breaks downstream steps that need SEQRES
    # (model) or CONECT (prepare/minimize).
    _apply_transform_preserving_headers(input_pdb, output_pdb, R, t)

    if verbose:
        print(f"  [align] {selection} RMSD: {rmsd_before:.3f} → "
              f"{rmsd_after:.3f} Å ({len(P)} atoms, "
              f"{len(matched_pairs)} residues)")

    return rmsd_before, rmsd_after, len(P)


def _copy_file(src, dst):
    if str(Path(src).resolve()) == str(Path(dst).resolve()):
        return
    with open(src) as f_in, open(dst, 'w') as f_out:
        f_out.writelines(f_in)


def _apply_transform_preserving_headers(input_pdb, output_pdb, R, t):
    """Apply the Kabsch (R, t) to every ATOM/HETATM line in `input_pdb`;
    write to `output_pdb`. Every other line (SEQRES, HELIX, SHEET, SSBOND,
    LINK, CISPEP, HET, DBREF, SEQADV, CONECT, REMARK, HEADER, TITLE,
    MODEL/ENDMDL, TER, END, CRYST1) is passed through byte-identical.

    Read-then-write so input_pdb == output_pdb is safe.
    """
    with open(input_pdb) as f:
        lines = f.readlines()
    out = []
    r00, r01, r02 = R[0, 0], R[0, 1], R[0, 2]
    r10, r11, r12 = R[1, 0], R[1, 1], R[1, 2]
    r20, r21, r22 = R[2, 0], R[2, 1], R[2, 2]
    t0, t1, t2 = t[0], t[1], t[2]
    for line in lines:
        if (line.startswith(('ATOM  ', 'HETATM')) and len(line) >= 54):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                out.append(line)
                continue
            xn = r00 * x + r01 * y + r02 * z + t0
            yn = r10 * x + r11 * y + r12 * z + t1
            zn = r20 * x + r21 * y + r22 * z + t2
            new_coord = f"{xn:8.3f}{yn:8.3f}{zn:8.3f}"
            out.append(line[:30] + new_coord + line[54:])
        else:
            out.append(line)
    with open(output_pdb, 'w') as f:
        f.writelines(out)
