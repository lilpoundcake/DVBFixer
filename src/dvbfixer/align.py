"""Kabsch superposition of one PDB onto a reference.

Used inside `dvbfixer zbs` after every pipeline step (default ON, opt-out
via `--no-align-to-input`) so that interim outputs stay in the same
Cartesian frame as the user's original input — no accumulated drift from
successive OpenMM minimizations.

Only exposed as an internal helper. No standalone CLI subcommand.
"""

from __future__ import annotations

from pathlib import Path


_BACKBONE_ATOMS = ('N', 'CA', 'C', 'O')


def _load_atom_records(pdb_path):
    """Parse ATOM+HETATM lines. Returns list of (line, key, x, y, z, is_hetatm).

    key is (chain, resseq_str, icode_str, atomname). resseq is kept as a
    STRING to preserve insertion codes and PDB negative-number quirks.
    """
    out = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                out.append((line, None, None, None, None, False))
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                out.append((line, None, None, None, None, False))
                continue
            chain = line[21]
            resseq = line[22:26].strip()
            icode = line[26].strip() if len(line) > 26 else ''
            atomname = line[12:16].strip()
            key = (chain, resseq, icode, atomname)
            is_het = line.startswith('HETATM')
            out.append((line, key, x, y, z, is_het))
    return out


def _select_indices(records, selection, is_reference=False):
    """Return list of dict indices whose atom matches the selection."""
    sel = selection.lower()
    picked = []
    for i, (_line, key, _x, _y, _z, is_het) in enumerate(records):
        if key is None:
            continue
        atomname = key[3]
        if sel == 'ca':
            if atomname == 'CA' and not is_het:
                picked.append(i)
        elif sel == 'backbone':
            if atomname in _BACKBONE_ATOMS and not is_het:
                picked.append(i)
        elif sel == 'heavy':
            if not atomname.startswith('H'):
                picked.append(i)
        elif sel == 'all':
            picked.append(i)
        else:
            raise ValueError(f"unknown selection '{selection}'")
    return picked


def _kabsch(P, Q):
    """Return (R, t) such that (P @ R.T + t) best superposes P onto Q.

    P, Q: (N, 3) numpy arrays. Both must have the same length.
    """
    import numpy as np
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    H = P0.T @ Q0
    U, _S, Vt = np.linalg.svd(H)
    # Correct for reflection so det(R) == +1.
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Qc - Pc @ R.T
    return R, t


def kabsch_align_pdb(input_pdb, reference_pdb, output_pdb, *,
                     selection='backbone', verbose=False):
    """Kabsch-superpose `input_pdb` onto `reference_pdb`, write `output_pdb`.

    Atoms are matched by (chain, resseq, icode, atomname). The Kabsch
    rotation + translation is computed on the intersection filtered by
    `selection` and applied to EVERY atom line in the input (ATOM +
    HETATM). Non-coordinate lines pass through unchanged.

    Args:
        input_pdb: source PDB to be aligned.
        reference_pdb: reference PDB defining the target frame.
        output_pdb: destination PDB (may equal input_pdb).
        selection: `'backbone'` (default; N/CA/C/O of standard AAs),
            `'ca'` (CA only), `'heavy'` (all non-H atoms of any residue),
            or `'all'` (every atom).
        verbose: print RMSD + atom counts.

    Returns:
        (rmsd_before, rmsd_after, n_matched) — floats/int for tests and logging.
        rmsd_before / rmsd_after are None if the alignment atom set was empty.
    """
    import numpy as np

    inp_recs = _load_atom_records(input_pdb)
    ref_recs = _load_atom_records(reference_pdb)

    inp_sel_idx = _select_indices(inp_recs, selection)
    ref_map = {}
    for i in _select_indices(ref_recs, selection, is_reference=True):
        key = ref_recs[i][1]
        ref_map[key] = i

    P_list = []
    Q_list = []
    for i in inp_sel_idx:
        key = inp_recs[i][1]
        j = ref_map.get(key)
        if j is None:
            continue
        _, _, xp, yp, zp, _ = inp_recs[i]
        _, _, xq, yq, zq, _ = ref_recs[j]
        P_list.append((xp, yp, zp))
        Q_list.append((xq, yq, zq))

    n_matched = len(P_list)
    if n_matched < 3:
        if verbose:
            print(f"  [align] not enough matching atoms in selection "
                  f"'{selection}' (found {n_matched}); writing input "
                  f"unchanged")
        _copy_file(input_pdb, output_pdb)
        return None, None, n_matched

    P = np.asarray(P_list)
    Q = np.asarray(Q_list)
    rmsd_before = float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))

    R, t = _kabsch(P, Q)

    # Apply to every input coordinate line.
    out_lines = []
    all_new_coords = []
    for rec in inp_recs:
        line, key, x, y, z, _is_het = rec
        if key is None:
            out_lines.append(line)
            continue
        v = np.array([x, y, z])
        v_new = v @ R.T + t
        all_new_coords.append(v_new)
        # Replace cols 30-54 with new coords; preserve everything else.
        new_coord = f"{v_new[0]:8.3f}{v_new[1]:8.3f}{v_new[2]:8.3f}"
        out_lines.append(line[:30] + new_coord + line[54:])

    # Post-alignment RMSD on the selection.
    P_new = P @ R.T + t
    rmsd_after = float(np.sqrt(np.mean(np.sum((P_new - Q) ** 2, axis=1))))

    with open(output_pdb, 'w') as f:
        f.writelines(out_lines)

    if verbose:
        print(f"  [align] {selection} RMSD: {rmsd_before:.3f} → "
              f"{rmsd_after:.3f} Å ({n_matched} atoms)")

    return rmsd_before, rmsd_after, n_matched


def _copy_file(src, dst):
    if str(Path(src).resolve()) == str(Path(dst).resolve()):
        return
    with open(src) as f_in, open(dst, 'w') as f_out:
        f_out.writelines(f_in)
