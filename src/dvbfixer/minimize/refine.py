"""xtb and OpenBabel refinement passes for ``dvbfixer minimize``.

Split out of the flat ``minimize.py`` in Phase 2.1 of the revision plan.
Contains the two BioLuminate-style refinement engines that run AFTER
OpenMM minimization to polish arbitrary heterogens without needing a
per-ligand template:

- :func:`refine_with_xtb` — invokes the ``xtb`` binary (GFN-FF universal
  force field) on a temp XYZ file.
- :func:`refine_with_obminimize` — invokes OpenBabel's ``obminimize`` on a
  temp PDB (via CLI when the whole system is free, via pybel when frozen
  atoms are needed since the CLI has no freeze flag).

Both engines can restrict work to just the heterogen residues via
:func:`_extract_heterogen_subsystem`; when they do, the anchor residue's
backbone + the heterogen-side linkage atom stay frozen so the
OpenMM-AMBER interface geometry survives.

Also exports the Kabsch glycan-tracking helper
:func:`_rigid_track_glycan_trees` — used by the strip-and-splice fallback
inside ``pipeline.minimize`` to follow the post-minimize protein anchor
without re-parametrising the glycan tree.
"""

from __future__ import annotations

# Note: each refinement engine imports its own OpenMM bits inside the
# function body (nanometer conversion, Modeller for subsystem extract).
# Kept lazy so ``from dvbfixer.minimize.refine import _find_binary``
# doesn't pull OpenMM when the caller only needs the CLI helper.


def _find_binary(name):
    """Locate a binary, checking PATH and the current Python env's bin dir."""
    import os as _os
    import shutil as _sh
    import sys as _sys
    found = _sh.which(name)
    if found:
        return found
    # Fall back to Python env's bin directory (handles direct-executable env)
    py_bin = _os.path.dirname(_sys.executable)
    candidate = _os.path.join(py_bin, name)
    if _os.access(candidate, _os.X_OK):
        return candidate
    return None


def _write_xyz(path, topology, positions, comment=""):
    """Write topology + positions to XYZ format (used by xtb)."""
    from openmm.unit import nanometer
    atoms = list(topology.atoms())
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for atom in atoms:
            p = positions[atom.index].value_in_unit(nanometer)
            x = float(p[0]) * 10.0  # nm → Å
            y = float(p[1]) * 10.0
            z = float(p[2]) * 10.0
            sym = atom.element.symbol
            f.write(f"{sym:<3s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def _read_xyz_coords(path):
    """Read coordinates from XYZ file. Returns list of Vec3 in nm."""
    from openmm import Vec3
    coords = []
    with open(path) as f:
        n = int(f.readline().strip())
        f.readline()  # comment line
        for _ in range(n):
            parts = f.readline().split()
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            # Å → nm
            coords.append(Vec3(x / 10.0, y / 10.0, z / 10.0))
    return coords


def _build_frozen_atom_list(topology, heterogen_only):
    """Return list of atom indices (1-based) to freeze during xtb opt.

    If heterogen_only=True, freezes protein atoms; otherwise freezes nothing.
    """
    if not heterogen_only:
        return []
    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    frozen = []
    for atom in topology.atoms():
        if atom.residue.name in known:
            frozen.append(atom.index + 1)  # 1-based for xtb
    return frozen


def refine_with_xtb(topology, positions, cycles=200, heterogens_only=False,
                    verbose=False):
    """Refine geometry with xtb GFN-FF universal force field.

    GFN-FF auto-parametrizes any organic molecule from connectivity rules.
    For whole-protein systems use heterogens_only=True (extract sub-system),
    otherwise xtb tries to handle the full system which is slow.
    """
    from openmm import Vec3
    from openmm.unit import Quantity, nanometer

    xtb_bin = _find_binary('xtb')
    if xtb_bin is None:
        print("WARNING: xtb binary not found in PATH — skipping xtb refinement")
        return positions

    if heterogens_only:
        sub_top, sub_pos, idx_map, anchor_indices = _extract_heterogen_subsystem(
            topology, positions
        )
        n_sub = sum(1 for _ in sub_top.atoms())
        if n_sub == 0:
            return positions
        print(f"\n=== xtb GFN-FF refinement ({cycles} cycles, "
              f"heterogens-only — {n_sub} atoms, "
              f"{len(anchor_indices)} protein anchors frozen) ===")
        new_sub_pos = _run_xtb(sub_top, sub_pos, xtb_bin, cycles, verbose,
                                frozen_indices=anchor_indices)
        if new_sub_pos is None:
            return positions
        coords = []
        for p in positions:
            v = p.value_in_unit(nanometer)
            coords.append(Vec3(float(v[0]), float(v[1]), float(v[2])))
        for full_idx, sub_idx in idx_map.items():
            if sub_idx in anchor_indices:
                continue  # anchors didn't move; skip
            p = new_sub_pos[sub_idx].value_in_unit(nanometer)
            coords[full_idx] = Vec3(float(p[0]), float(p[1]), float(p[2]))
        return Quantity(coords, nanometer)
    else:
        n = sum(1 for _ in topology.atoms())
        if n > 5000:
            print(f"\nINFO: full-system xtb on {n} atoms takes hours — "
                  f"auto-switching to --refine-heterogens-only (use that flag "
                  f"explicitly to silence this notice)")
            return refine_with_xtb(
                topology, positions, cycles=cycles,
                heterogens_only=True, verbose=verbose,
            )
        print(f"\n=== xtb GFN-FF refinement ({cycles} cycles, "
              f"whole system — {n} atoms) ===")
        new_pos = _run_xtb(topology, positions, xtb_bin, cycles, verbose)
        return new_pos if new_pos is not None else positions


def _run_xtb(topology, positions, xtb_bin, cycles, verbose, frozen_indices=None):
    """Internal: write XYZ, run xtb --opt --gfnff, return refined Quantity or None.

    frozen_indices: set of sub-topology atom indices (0-based) to freeze.
    Written to xcontrol with $fix atoms (1-based).
    """
    import os as _os
    import shutil as _sh
    import subprocess as _sp
    import tempfile as _tf

    from openmm.unit import Quantity, nanometer

    workdir = _tf.mkdtemp(prefix='dvbfixer_xtb_')
    old_cwd = _os.getcwd()
    try:
        input_xyz = _os.path.join(workdir, 'input.xyz')
        _write_xyz(input_xyz, topology, positions, "dvbfixer xtb input")

        cmd = [xtb_bin, 'input.xyz', '--opt', '--gfnff',
               '--cycles', str(cycles), '--norestart']

        # Freeze protein anchor atoms via xcontrol $fix block (1-based indices)
        if frozen_indices:
            xc_path = _os.path.join(workdir, 'xcontrol')
            sorted_idx = sorted(i + 1 for i in frozen_indices)
            # Group consecutive indices into ranges for compactness
            ranges = []
            start = prev = sorted_idx[0]
            for i in sorted_idx[1:]:
                if i == prev + 1:
                    prev = i
                else:
                    ranges.append((start, prev))
                    start = prev = i
            ranges.append((start, prev))
            range_str = ",".join(
                f"{s}" if s == e else f"{s}-{e}" for s, e in ranges
            )
            with open(xc_path, 'w') as f:
                f.write("$fix\n")
                f.write(f"  atoms: {range_str}\n")
                f.write("$end\n")
            cmd.extend(['--input', 'xcontrol'])
        _os.chdir(workdir)
        if verbose:
            print(f"  Running: {' '.join(cmd)}")
        try:
            result = _sp.run(cmd, capture_output=True, text=True, timeout=3600)
        except _sp.TimeoutExpired:
            print("WARNING: xtb timeout")
            return None
        if result.returncode != 0:
            print(f"WARNING: xtb failed (returncode {result.returncode})")
            if verbose:
                print("--- stderr ---")
                print((result.stderr or "")[-3000:])
                print("--- stdout (last 30 lines) ---")
                print('\n'.join(result.stdout.splitlines()[-30:]))
            return None

        opt_xyz = _os.path.join(workdir, 'xtbopt.xyz')
        if not _os.path.exists(opt_xyz):
            print("WARNING: xtb produced no xtbopt.xyz")
            return None

        new_coords = _read_xyz_coords(opt_xyz)
        for line in result.stdout.splitlines():
            if 'TOTAL ENERGY' in line:
                print(f"  xtb {line.strip()}")
                break
        print(f"  xtb refined {len(new_coords)} atoms")
        return Quantity(new_coords, nanometer)
    finally:
        _os.chdir(old_cwd)
        _sh.rmtree(workdir, ignore_errors=True)


def _rigid_track_glycan_trees(in_top, in_pos, result_array, min_pos_map,
                               verbose=False):
    """Rigid-body transform each glycan tree to follow its protein anchor.

    The legacy strip-and-splice path minimizes protein-only, then restores
    HETATM coords verbatim. But the protein anchor (ASN/SER/THR side chain)
    moves during minimization, leaving the glycan in the wrong place — the
    ND2-C1 bond stretches and the glycan can clash with the moved protein.

    Fix: for each protein-heterogen bond, compute the Kabsch (rotation +
    translation) transform from prep anchor heavy atoms → post-min anchor
    heavy atoms, then apply it to every atom in the BFS-connected glycan
    tree. Result: glycan keeps its relative orientation to the anchor amide;
    bond length and stereochemistry are preserved.

    in_top/in_pos: original (prep) topology + positions.
    result_array: numpy (n_atoms, 3) — protein already overwritten with
                  minimized positions, HETATM still at prep positions.
    min_pos_map: dict (chain, resid, atomname) → minimized position (nm).
    """
    import numpy as np
    from openmm.unit import nanometer as nm_unit

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS
    atoms = list(in_top.atoms())

    # Find glycosylated protein residues + which heterogen atom they connect to.
    het_atom_set = {a.index for a in atoms if a.residue.name not in known}
    anchor_to_het_atoms = {}  # (chain, rid) → set of bonded heterogen atom indices
    for b in in_top.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_set and bi not in het_atom_set:
            r = b[1].residue
            anchor_to_het_atoms.setdefault((r.chain.id, r.id), set()).add(ai)
        elif bi in het_atom_set and ai not in het_atom_set:
            r = b[0].residue
            anchor_to_het_atoms.setdefault((r.chain.id, r.id), set()).add(bi)

    if not anchor_to_het_atoms:
        return result_array

    # Build heterogen adjacency for BFS through the glycan tree.
    het_adj = {idx: set() for idx in het_atom_set}
    for b in in_top.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_set and bi in het_atom_set:
            het_adj[ai].add(bi)
            het_adj[bi].add(ai)

    # Geometry recipe per (residue, anchor_atom). For each recipe, we
    # compute the IDEAL linkage atom (e.g. NAG C1) position from the post-min
    # protein side-chain geometry — bypassing the AMBER-placed HD21/HD22 which
    # may have been put on the wrong side.
    #
    # N-linked (ASN/NLN ND2): trans amide. C1 lies in the amide plane (defined
    # by CG, OD1, ND2), at 120° from the CG-ND2 axis, on the OPPOSITE side of
    # the plane from OD1. This is the canonical E,Z amide configuration.
    #
    # O-linked (SER/THR): sp3 tetrahedral. C1 substitutes for the hydroxyl H,
    # at the position the H currently occupies (anti-periplanar to CB).
    _IDEAL_GEOMETRY = {
        ('ASN', 'ND2'): {'bond_nm': 0.145, 'type': 'amide_trans',
                          'axis': 'CG', 'plane_neighbor': 'OD1', 'angle_deg': 120.0},
        ('NLN', 'ND2'): {'bond_nm': 0.145, 'type': 'amide_trans',
                          'axis': 'CG', 'plane_neighbor': 'OD1', 'angle_deg': 120.0},
        ('SER', 'OG'):  {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG'},
        ('OLS', 'OG'):  {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG'},
        ('THR', 'OG1'): {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG1'},
        ('OLT', 'OG1'): {'bond_nm': 0.143, 'type': 'replace_H', 'h_atom': 'HG1'},
    }

    def _ideal_linkage_pos(anchor_atom_name, anchor_res_name, anchor_ch, anchor_rid):
        """Compute the IDEAL linkage atom position from the post-min anchor
        side-chain geometry. Returns ideal_pos (numpy, nm) or None on missing
        reference atoms.
        """
        geom = _IDEAL_GEOMETRY.get((anchor_res_name, anchor_atom_name))
        if geom is None:
            return None
        akey = (anchor_ch, anchor_rid, anchor_atom_name)
        if akey not in min_pos_map:
            return None
        a_pos = np.array(min_pos_map[akey])

        if geom['type'] == 'replace_H':
            hkey = (anchor_ch, anchor_rid, geom['h_atom'])
            if hkey not in min_pos_map:
                return None
            h_pos = np.array(min_pos_map[hkey])
            v = h_pos - a_pos
            n = np.linalg.norm(v)
            if n < 1e-6:
                return None
            return a_pos + geom['bond_nm'] * (v / n)

        if geom['type'] == 'amide_trans':
            ax_key = (anchor_ch, anchor_rid, geom['axis'])
            pn_key = (anchor_ch, anchor_rid, geom['plane_neighbor'])
            if ax_key not in min_pos_map or pn_key not in min_pos_map:
                return None
            ax_pos = np.array(min_pos_map[ax_key])         # e.g. CG
            pn_pos = np.array(min_pos_map[pn_key])         # e.g. OD1
            # u_axis = unit vector from anchor (ND2) toward axis atom (CG)
            u_axis = ax_pos - a_pos
            n = np.linalg.norm(u_axis)
            if n < 1e-6:
                return None
            u_axis = u_axis / n
            # OD1 relative to anchor (ND2)
            v_pn = pn_pos - a_pos
            # In-plane perpendicular: component of OD1-ND2 perpendicular to
            # CG-ND2 axis. This points TOWARD OD1 in the amide plane.
            v_perp = v_pn - np.dot(v_pn, u_axis) * u_axis
            n_perp = np.linalg.norm(v_perp)
            if n_perp < 1e-6:
                return None
            u_perp_toward_pn = v_perp / n_perp
            # Trans-amide C1 direction:
            #   angle between u_axis and u_C1 = 120° (CG-ND2-C1)
            #   u_C1 in amide plane, on OPPOSITE side from OD1 across CG-ND2.
            # In-plane decomposition:
            #   u_C1 = cos(angle) * (-u_axis_from_anchor_to_axis) (so angle from u_axis = 120°)
            # equivalently: u_C1 = cos(180° - angle) * u_axis + sin(angle) * (-u_perp_toward_pn)
            # but simpler:
            theta = np.radians(geom['angle_deg'])
            # u_C1 is at angle theta from u_axis (measured at ND2), in the amide
            # plane, on the side AWAY from OD1.
            #   u_C1 = cos(theta) * u_axis + sin(theta) * (-u_perp_toward_pn)
            u_c1 = np.cos(theta) * u_axis + np.sin(theta) * (-u_perp_toward_pn)
            return a_pos + geom['bond_nm'] * u_c1

        return None

    n_trees = 0
    n_atoms_moved = 0
    for (anchor_ch, anchor_rid), root_het_atoms in anchor_to_het_atoms.items():
        # Anchor residue in the input topology
        anchor_res = None
        for r in in_top.residues():
            if r.chain.id == anchor_ch and r.id == anchor_rid:
                anchor_res = r
                break
        if anchor_res is None:
            continue

        # Collect heavy atoms of the anchor side chain + linkage neighbor for
        # the Kabsch fit. We exclude H atoms (their positions are dominated by
        # AMBER addHydrogens, not by minimization itself).
        prep_pts = []
        post_pts = []
        for a in anchor_res.atoms():
            if a.element.symbol == 'H':
                continue
            key = (a.residue.chain.id, a.residue.id, a.name)
            if key not in min_pos_map:
                continue
            pp = in_pos[a.index].value_in_unit(nm_unit)
            mp = min_pos_map[key]
            prep_pts.append([float(pp[0]), float(pp[1]), float(pp[2])])
            post_pts.append([float(mp[0]), float(mp[1]), float(mp[2])])

        if len(prep_pts) < 3:
            continue

        # Inject the IDEAL linkage atom position as an extra Kabsch reference.
        # Without this, the bad prep-side glycosidic geometry (e.g. CG-ND2-C1
        # ≈ 167° instead of 122° after Modeller loop modeling) is faithfully
        # reproduced at the post-min anchor — the Kabsch is dominated by the
        # anchor-side coordinates and preserves the relative orientation.
        # Injecting the ideal C1 (from post-min amide trigonal-planar geometry)
        # pulls the fit toward correct stereochemistry.
        anchor_atoms_by_name = {a.name: a for a in anchor_res.atoms()}
        # Use the first bonded heterogen atom as the linkage atom (e.g. NAG C1)
        link_het_idx = next(iter(root_het_atoms))
        link_het_atom = atoms[link_het_idx]
        # Find the protein anchor atom name (e.g. ND2 for ASN)
        anchor_link_atom_name = None
        for b in in_top.bonds():
            ai, bi = b[0].index, b[1].index
            if ai == link_het_idx and bi in {a.index for a in anchor_res.atoms()}:
                anchor_link_atom_name = atoms[bi].name
                break
            if bi == link_het_idx and ai in {a.index for a in anchor_res.atoms()}:
                anchor_link_atom_name = atoms[ai].name
                break
        # Resolve the anchor residue name. The legacy path renames NLN→ASN
        # in topology BEFORE AMBER addHydrogens, so anchor_res.name may now
        # be 'ASN' even though the prep input had 'NLN'. Both keys point to
        # the same trigonal-planar geometry recipe (CG+HD21), so either name
        # works — we just need a registered entry.
        post_anchor_name = anchor_res.name
        ideal_pos = None
        if anchor_link_atom_name:
            for try_name in (anchor_res.name,
                             {'ASN': 'NLN', 'SER': 'OLS', 'THR': 'OLT',
                              'NLN': 'ASN', 'OLS': 'SER', 'OLT': 'THR'}.get(anchor_res.name)):
                if try_name is None:
                    continue
                if (try_name, anchor_link_atom_name) in _IDEAL_GEOMETRY:
                    post_anchor_name = try_name
                    break
            ideal_pos = _ideal_linkage_pos(anchor_link_atom_name,
                                            post_anchor_name,
                                            anchor_ch, anchor_rid)
        if ideal_pos is not None:
            prep_link = in_pos[link_het_idx].value_in_unit(nm_unit)
            prep_pts.append([float(prep_link[0]), float(prep_link[1]),
                              float(prep_link[2])])
            post_pts.append(ideal_pos.tolist())

        prep_arr = np.array(prep_pts)
        post_arr = np.array(post_pts)
        # Kabsch
        prep_c = prep_arr.mean(axis=0)
        post_c = post_arr.mean(axis=0)
        P = prep_arr - prep_c
        Q = post_arr - post_c
        H = P.T @ Q
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        t = post_c - R @ prep_c

        # BFS to collect glycan tree atoms (heterogen atoms reachable from
        # the anchor's bonded heterogen atoms).
        seen = set()
        queue = list(root_het_atoms)
        while queue:
            cur = queue.pop()
            if cur in seen:
                continue
            seen.add(cur)
            for nb in het_adj.get(cur, ()):
                if nb not in seen:
                    queue.append(nb)

        # Apply transform to each glycan atom
        for idx in seen:
            p = in_pos[idx].value_in_unit(nm_unit)
            v = np.array([float(p[0]), float(p[1]), float(p[2])])
            v_new = R @ v + t
            result_array[idx] = v_new
            n_atoms_moved += 1

        # Post-Kabsch correction: snap the linkage atom exactly to the ideal
        # post-min trigonal-planar position by translating the whole glycan
        # tree by (ideal - tracked). The Kabsch fit only approximately satisfies
        # the ideal C1 target (it gets averaged with the 5 ASN reference
        # points), but for correct stereochemistry the C1 needs to land on
        # the exact amide plane direction. Translation preserves intra-glycan
        # geometry; the next obminimize pass relaxes any residual strain.
        if anchor_link_atom_name and 'ideal_pos' in locals() and ideal_pos is not None:
            tracked_link = result_array[link_het_idx]
            correction = ideal_pos - tracked_link
            for idx in seen:
                result_array[idx] = result_array[idx] + correction

            # Also snap the amide H (HD21 for ASN/NLN) to its canonical
            # trigonal-planar position — cis to OD1, in the amide plane.
            # AMBER places HD21 at some valid sp2 position but not always at
            # the canonical cis-OD1 spot; this enforces CG-ND2-HD21 ≈ 120°
            # and the H in the amide plane.
            if (post_anchor_name, anchor_link_atom_name) == ('ASN', 'ND2') or \
               (post_anchor_name, anchor_link_atom_name) == ('NLN', 'ND2'):
                # Find HD21 atom in in_top (prep NLN topology has HD21)
                hd21_in_idx = None
                for a in anchor_res.atoms():
                    if a.name == 'HD21':
                        hd21_in_idx = a.index
                        break
                if hd21_in_idx is not None:
                    # Compute ideal HD21 position (cis to OD1, in amide plane).
                    # Same formula as for C1 but with opposite perp sign.
                    nd2_key = (anchor_ch, anchor_rid, 'ND2')
                    cg_key = (anchor_ch, anchor_rid, 'CG')
                    od1_key = (anchor_ch, anchor_rid, 'OD1')
                    if all(k in min_pos_map for k in (nd2_key, cg_key, od1_key)):
                        nd2_p = np.array(min_pos_map[nd2_key])
                        cg_p = np.array(min_pos_map[cg_key])
                        od1_p = np.array(min_pos_map[od1_key])
                        u_ax = cg_p - nd2_p
                        u_ax = u_ax / np.linalg.norm(u_ax)
                        v_pn = od1_p - nd2_p
                        v_perp = v_pn - np.dot(v_pn, u_ax) * u_ax
                        n_perp = np.linalg.norm(v_perp)
                        if n_perp > 1e-6:
                            u_perp_to_od1 = v_perp / n_perp
                            theta = np.radians(120.0)
                            # HD21: trigonal-planar, cis to OD1 → same perp sign as OD1
                            u_h = np.cos(theta) * u_ax + np.sin(theta) * u_perp_to_od1
                            result_array[hd21_in_idx] = nd2_p + 0.101 * u_h  # N-H 1.01 Å
        n_trees += 1

    if verbose and n_trees:
        print(f"  Rigid-tracked {n_trees} glycan tree(s) "
              f"to follow protein anchor ({n_atoms_moved} atoms)")
    return result_array


def _restore_glycosylated_h(out_top, out_pos, in_top, in_pos):
    """Restore side-chain positions of glycosylated ASN/SER/THR (NLN/OLS/OLT)
    to their input values to preserve the protein-glycan stereochemistry.

    After legacy strip-and-splice minimize, the AMBER ASN template treats the
    glycosylated residue as a normal amide. The minimization can:
    - Flip HD21/HD22 to the wrong side of ND2 (collision with linked C1)
    - Rotate the side chain so the amide plane no longer aligns with the sugar

    Restoring the entire side chain (CB, HB*, CG, OD1, ND2, HD21 for ASN;
    CB, OG, HG for SER; CB, CG2, OG1, HG1 for THR) snaps the protein-glycan
    interface back to the prepared geometry. Backbone atoms remain refined.
    """
    from openmm import Vec3
    from openmm.unit import Quantity, nanometer

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS

    glycosylated_residues = set()
    het_atom_set = {a.index for a in in_top.atoms() if a.residue.name not in known}
    for b in in_top.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_set and bi not in het_atom_set:
            r = b[1].residue
            glycosylated_residues.add((r.chain.id, r.id))
        elif bi in het_atom_set and ai not in het_atom_set:
            r = b[0].residue
            glycosylated_residues.add((r.chain.id, r.id))
    for r in in_top.residues():
        if r.name in ('NLN', 'OLS', 'OLT'):
            glycosylated_residues.add((r.chain.id, r.id))

    if not glycosylated_residues:
        return out_top, out_pos

    # Atoms to restore: ONLY the amide/hydroxyl group at the glycosidic linkage
    # (NOT HB2/HB3 or other tetrahedral H — those need AMBER-quality minimization
    # for proper sp3 angles). Restoring CG-OD1-ND2-HD21 preserves amide planarity
    # and HD21 orientation w.r.t. the linked sugar. CB and HB* are NOT touched —
    # they get the minimized sp3 geometry.
    _SIDECHAIN = {
        'ASN': {'CG', 'OD1', 'ND2', 'HD21', 'HD22'},
        'NLN': {'CG', 'OD1', 'ND2', 'HD21'},
        'SER': {'OG', 'HG', 'HG1'},
        'OLS': {'OG'},
        'THR': {'OG1', 'HG1'},
        'OLT': {'OG1'},
    }

    input_positions = {}
    for atom in in_top.atoms():
        key = (atom.residue.chain.id, atom.residue.id, atom.name)
        p = in_pos[atom.index].value_in_unit(nanometer)
        input_positions[key] = (float(p[0]), float(p[1]), float(p[2]))

    n_restored = 0
    out_pos_list = []
    for atom in out_top.atoms():
        p = out_pos[atom.index].value_in_unit(nanometer)
        v = Vec3(float(p[0]), float(p[1]), float(p[2]))
        res_key = (atom.residue.chain.id, atom.residue.id)
        if res_key in glycosylated_residues:
            # Determine which side-chain atom set to use. Output may have
            # ASN (renamed from NLN); input may have NLN — try both.
            sc_atoms = (_SIDECHAIN.get(atom.residue.name, set())
                        | _SIDECHAIN.get('NLN', set())
                        | _SIDECHAIN.get('OLS', set())
                        | _SIDECHAIN.get('OLT', set()))
            if atom.name in sc_atoms:
                ipos = input_positions.get(
                    (atom.residue.chain.id, atom.residue.id, atom.name)
                )
                if ipos is not None:
                    v = Vec3(*ipos)
                    n_restored += 1
        out_pos_list.append(v)
    if n_restored:
        print(f"  Restored {n_restored} side-chain atoms on glycosylated residues")
    return out_top, Quantity(out_pos_list, nanometer)


def _extract_heterogen_subsystem(topology, positions, padding_residues=False):
    """Build a sub-topology of heterogen residues PLUS protein anchor atoms.

    For each protein-heterogen bond (e.g. ASN ND2 → NAG C1), the protein
    atom is included as a single "anchor" atom (without its full residue).
    This preserves the glycosidic bond constraint during refinement; downstream
    the anchor atoms are frozen so the protein doesn't move.

    Returns (sub_topology, sub_positions, atom_index_map, anchor_sub_indices).
    atom_index_map[full_atom_index] = sub_atom_index for all included atoms.
    anchor_sub_indices = set of sub_topology atom indices that are protein
                        anchors (should be frozen during refinement).
    """
    from openmm import Vec3
    from openmm.app import Topology
    from openmm.unit import Quantity, nanometer

    from dvbfixer.ffutils import PROTEIN_RESIDUES, SOLVENT_IONS

    known = PROTEIN_RESIDUES | SOLVENT_IONS

    # Find protein atoms bonded to heterogens, then include their ENTIRE
    # residues as anchors. Including just the single linkage atom (e.g. ASN
    # ND2 alone) leaves xtb confused about bond chemistry — it perceives a
    # lone N and pulls the C-N bond to ~1 Å. Including the full ASN gives
    # xtb the proper sp3 N context.
    het_atom_indices = set()
    for res in topology.residues():
        if res.name not in known:
            for atom in res.atoms():
                het_atom_indices.add(atom.index)
    anchor_residues = set()  # residues whose atoms anchor a heterogen
    linkage_het_atoms = set()  # heterogen atoms at the cross-residue bond
    for b in topology.bonds():
        ai, bi = b[0].index, b[1].index
        if ai in het_atom_indices and bi not in het_atom_indices:
            anchor_residues.add(b[1].residue.index)
            linkage_het_atoms.add(ai)
        elif bi in het_atom_indices and ai not in het_atom_indices:
            anchor_residues.add(b[0].residue.index)
            linkage_het_atoms.add(bi)
    anchor_protein_atoms = set()
    for res in topology.residues():
        if res.index in anchor_residues:
            for atom in res.atoms():
                anchor_protein_atoms.add(atom.index)
    # ALSO freeze the heterogen linkage atom itself (e.g. NAG C1 bonded to
    # ASN ND2). Reason: UFF/MMFF don't have explicit amide planarity terms,
    # so during refinement the linkage C can rotate ~15° out of the amide
    # plane. Freezing it preserves the OpenMM-AMBER geometry of the
    # protein-glycan interface; the rest of the sugar refines freely.

    sub_top = Topology()
    sub_pos_list = []
    atom_index_map = {}
    anchor_sub_indices = set()

    # First pass: add heterogen residues
    for chain in topology.chains():
        new_chain = None
        for res in chain.residues():
            if res.name in known:
                continue
            if new_chain is None:
                new_chain = sub_top.addChain(chain.id)
            new_res = sub_top.addResidue(res.name, new_chain, res.id,
                                          res.insertionCode)
            for atom in res.atoms():
                new_atom = sub_top.addAtom(atom.name, atom.element, new_res)
                atom_index_map[atom.index] = new_atom.index
                # Linkage heterogen atom (e.g. NAG C1) is INTENTIONALLY LEFT
                # FREE so UFF can minimize the protein-glycan bond length and
                # angle. UFF's sp2 N type (N_2) has equilibrium angle 120° and
                # an improper torsion that keeps the amide near-planar, so the
                # geometry settles naturally instead of being snapped rigid.
                p = positions[atom.index].value_in_unit(nanometer)
                sub_pos_list.append(Vec3(float(p[0]), float(p[1]), float(p[2])))

    # Second pass: add full anchor residues with their original names/atoms.
    # Use a separate "Z" chain so they don't merge with heterogen residues.
    if anchor_protein_atoms:
        anchor_chain = sub_top.addChain("Z")
        from collections import defaultdict
        anchors_by_res = defaultdict(list)
        full_atoms = list(topology.atoms())
        for ai in anchor_protein_atoms:
            anchors_by_res[full_atoms[ai].residue].append(ai)
        for orig_res, atom_indices in anchors_by_res.items():
            anc_res = sub_top.addResidue(orig_res.name, anchor_chain,
                                         str(orig_res.id),
                                         orig_res.insertionCode)
            # Freeze only the BACKBONE + first sidechain carbon of the anchor
            # residue. The amide group (CG, OD1, ND2, HD21, HD22 for ASN; OG/HG
            # for SER; OG1/HG1 for THR) and the linkage atom are LEFT FREE so
            # UFF/MMFF can actually minimize the protein-glycan bond length,
            # angle, and amide planarity. This is the answer to "why doesn't
            # the tool minimize the ASN-glycan bond" — previously the whole
            # anchor was frozen and the bond was only rigid-tracked, not
            # minimized by any FF.
            _FREEZE_ANCHOR_ATOMS = {
                'ASN': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'NLN': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'SER': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'OLS': {'N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'C', 'O'},
                'THR': {'N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG2',
                        'HG21', 'HG22', 'HG23', 'C', 'O'},
                'OLT': {'N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG2',
                        'HG21', 'HG22', 'HG23', 'C', 'O'},
            }
            freeze_names = _FREEZE_ANCHOR_ATOMS.get(
                orig_res.name,
                # fallback: freeze everything for unfamiliar anchors
                {a.name for a in [full_atoms[ai] for ai in atom_indices]}
            )
            for ai in atom_indices:
                a = full_atoms[ai]
                new_atom = sub_top.addAtom(a.name, a.element, anc_res)
                atom_index_map[ai] = new_atom.index
                if a.name in freeze_names:
                    anchor_sub_indices.add(new_atom.index)
                p = positions[ai].value_in_unit(nanometer)
                sub_pos_list.append(Vec3(float(p[0]), float(p[1]), float(p[2])))

    # Carry all bonds where both atoms made it into the sub-topology
    sub_atoms = list(sub_top.atoms())
    for b in topology.bonds():
        a1 = atom_index_map.get(b[0].index)
        a2 = atom_index_map.get(b[1].index)
        if a1 is not None and a2 is not None:
            sub_top.addBond(sub_atoms[a1], sub_atoms[a2])

    sub_pos = Quantity(sub_pos_list, nanometer)
    return sub_top, sub_pos, atom_index_map, anchor_sub_indices


def refine_with_obminimize(topology, positions, ff='MMFF94s', steps=500,
                            heterogens_only=False, verbose=False):
    """Refine geometry with OpenBabel obminimize (MMFF94/UFF/GAFF).

    OpenBabel auto-types any organic molecule via SMARTS rules.
    Returns refined positions (Quantity).
    """
    from openmm import Vec3
    from openmm.unit import Quantity, nanometer

    if _find_binary('obminimize') is None:
        print("WARNING: obminimize binary not found in PATH — skipping refinement")
        return positions

    if heterogens_only:
        sub_top, sub_pos, idx_map, anchor_indices = _extract_heterogen_subsystem(
            topology, positions
        )
        n_sub = sum(1 for _ in sub_top.atoms())
        if n_sub == 0:
            return positions
        print(f"\n=== OpenBabel obminimize refinement ({ff}, {steps} steps, "
              f"heterogens-only — {n_sub} atoms, "
              f"{len(anchor_indices)} protein anchors frozen) ===")
        new_sub_pos = _run_obminimize(sub_top, sub_pos, ff, steps, verbose,
                                       frozen_indices=anchor_indices)
        if new_sub_pos is None:
            return positions
        coords = []
        for p in positions:
            v = p.value_in_unit(nanometer)
            coords.append(Vec3(float(v[0]), float(v[1]), float(v[2])))
        for full_idx, sub_idx in idx_map.items():
            if sub_idx in anchor_indices:
                continue  # anchors didn't move
            p = new_sub_pos[sub_idx].value_in_unit(nanometer)
            coords[full_idx] = Vec3(float(p[0]), float(p[1]), float(p[2]))
        return Quantity(coords, nanometer)
    else:
        n_atoms = sum(1 for _ in topology.atoms())
        if n_atoms > 5000:
            print(f"\nINFO: full-system obminimize on {n_atoms} atoms is "
                  f"memory-intensive — auto-switching to --refine-heterogens-only "
                  f"(use that flag explicitly to silence this notice)")
            return refine_with_obminimize(
                topology, positions, ff=ff, steps=steps,
                heterogens_only=True, verbose=verbose,
            )
        print(f"\n=== OpenBabel obminimize refinement ({ff}, {steps} steps, "
              f"whole system — {n_atoms} atoms) ===")
        new_pos = _run_obminimize(topology, positions, ff, steps, verbose)
        return new_pos if new_pos is not None else positions


def _run_obminimize_pybel(topology, positions, ff, steps, frozen_indices, verbose):
    """Run OpenBabel via Python API with atom freezing.

    OBFFConstraints.AddAtomConstraint(idx) freezes atoms by 1-based index.
    Used when frozen_indices is non-empty (CLI obminimize has no freeze flag).
    """
    import os as _os
    import tempfile as _tf

    from openmm import Vec3
    from openmm.app import PDBFile
    from openmm.unit import Quantity, nanometer

    try:
        from openbabel import openbabel as ob
        from openbabel import pybel
    except ImportError:
        print("WARNING: openbabel Python bindings missing — skipping refinement")
        return None

    workdir = _tf.mkdtemp(prefix='dvbfixer_obmin_')
    try:
        in_pdb = _os.path.join(workdir, 'in.pdb')
        with open(in_pdb, 'w') as f:
            PDBFile.writeFile(topology, positions, f, keepIds=True)

        # Try requested FF, fall back to UFF on setup failure
        ffs_to_try = [ff] if ff == 'UFF' else [ff, 'UFF']
        for try_ff in ffs_to_try:
            mol = next(pybel.readfile('pdb', in_pdb))
            obff = ob.OBForceField.FindForceField(try_ff)
            if obff is None:
                continue
            constraints = ob.OBFFConstraints()
            for idx in sorted(frozen_indices):
                constraints.AddAtomConstraint(idx + 1)  # 1-based
            if not obff.Setup(mol.OBMol, constraints):
                if verbose:
                    print(f"  {try_ff} Setup failed; trying next FF")
                continue
            obff.SteepestDescent(steps)
            obff.GetCoordinates(mol.OBMol)
            used_ff = try_ff
            if try_ff != ff:
                print(f"  ({ff} setup failed; refined with UFF instead)")
            print(f"  obminimize refined {mol.OBMol.NumAtoms()} atoms "
                  f"({used_ff}, {len(frozen_indices)} frozen)")
            # Extract coords
            n_atoms = sum(1 for _ in topology.atoms())
            if mol.OBMol.NumAtoms() != n_atoms:
                print("WARNING: atom count mismatch after obminimize")
                return None
            coords = []
            for ai in range(1, mol.OBMol.NumAtoms() + 1):
                ob_a = mol.OBMol.GetAtom(ai)
                # OpenBabel stores Å; convert to nm
                coords.append(Vec3(ob_a.GetX() / 10.0,
                                   ob_a.GetY() / 10.0,
                                   ob_a.GetZ() / 10.0))
            return Quantity(coords, nanometer)
        print("WARNING: obminimize Python API: no FF could be set up")
        return None
    finally:
        import shutil as _sh
        _sh.rmtree(workdir, ignore_errors=True)


def _run_obminimize(topology, positions, ff, steps, verbose, frozen_indices=None):
    """Internal: run obminimize via CLI (no freeze) or pybel API (with freeze).

    When the requested FF (MMFF94/MMFF94s) lacks parameters for some atoms,
    automatically retries with UFF (universal FF, covers any element).
    """
    if frozen_indices:
        return _run_obminimize_pybel(
            topology, positions, ff, steps, frozen_indices, verbose
        )
    import os as _os
    import shutil as _sh
    import subprocess as _sp
    import tempfile as _tf

    from openmm import Vec3
    from openmm.app import PDBFile
    from openmm.unit import Quantity, nanometer

    workdir = _tf.mkdtemp(prefix='dvbfixer_obmin_')
    try:
        in_pdb = _os.path.join(workdir, 'in.pdb')
        out_pdb = _os.path.join(workdir, 'out.pdb')
        with open(in_pdb, 'w') as f:
            PDBFile.writeFile(topology, positions, f, keepIds=True)
        obmin_bin = _find_binary('obminimize')

        # Try requested FF, fall back to UFF if param lookup fails
        ffs_to_try = [ff]
        if ff != 'UFF':
            ffs_to_try.append('UFF')

        result = None
        used_ff = None
        for try_ff in ffs_to_try:
            cmd = [obmin_bin, '-ff', try_ff, '-n', str(steps), in_pdb]
            try:
                result = _sp.run(cmd, capture_output=True, text=True, timeout=3600)
            except _sp.TimeoutExpired:
                print("WARNING: obminimize timeout")
                return None
            stderr = result.stderr or ""
            stdout = result.stdout or ""
            param_missing = (
                'could not setup force field' in stderr.lower()
                or 'could not find van der waals' in stdout.lower()
                or 'could not find van der waals' in stderr.lower()
            )
            if result.returncode == 0 and stdout.strip():
                used_ff = try_ff
                if try_ff != ff:
                    print(f"  ({ff} lacked parameters; refined with UFF instead)")
                break
            if param_missing and try_ff != ffs_to_try[-1]:
                print(f"  {try_ff} parameters incomplete — retrying with UFF...")
                continue
            # Hard failure or unexpected output — give up
            print(f"WARNING: obminimize failed with {try_ff} "
                  f"(rc={result.returncode}). Output left unrefined.")
            if verbose:
                print((stderr or stdout)[-1500:])
            return None

        if used_ff is None:
            return None
        with open(out_pdb, 'w') as f:
            f.write(result.stdout)
        new_pdb = PDBFile(out_pdb)
        if sum(1 for _ in new_pdb.topology.atoms()) != sum(1 for _ in topology.atoms()):
            print("WARNING: obminimize output atom count mismatch")
            return None
        coords = []
        for p in new_pdb.positions:
            v = p.value_in_unit(nanometer)
            coords.append(Vec3(float(v[0]), float(v[1]), float(v[2])))
        print(f"  obminimize refined {len(coords)} atoms ({used_ff})")
        return Quantity(coords, nanometer)
    finally:
        _sh.rmtree(workdir, ignore_errors=True)


