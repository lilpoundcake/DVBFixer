"""Cluster glycan conformations from MD trajectories using torsion angle RMSD.

Implements the GFDB approach (Glycan Fragment Database):
1. Auto-detect glycosidic linkages from topology/first frame
2. Extract phi/psi/omega torsion angles across all frames
3. Build pairwise circular-RMSD distance matrix
4. GROMOS-style clustering (default 30 degree cutoff)
5. Output: torsion CSV, cluster assignments, summary, representative PDBs

Two clustering modes:
  --mode global     : cluster on all torsions simultaneously (default)
  --mode per-linkage: cluster each linkage independently, combine states

Torsion angle definitions (crystallographic convention):
  phi   = O5-C1-Ox-C'x          (ring O - anomeric C - glycosidic O - parent Cx)
  psi   = C1-Ox-C'x-C'(x-1)    (anomeric C - glycosidic O - parent Cx - parent C(x-1))
  omega = Ox-C'6-C'5-O'5        (only for 1->6 linkages)

For sialic acid (Neu5Ac): anomeric C = C2, ring O = O6.
"""

import argparse
import json
import sys
from collections import namedtuple
from pathlib import Path

import numpy as np

from dvbfixer.cli_types import nonnegative_int, positive_float, positive_int
from dvbfixer.residue_registry import (
    CHARMM_SUGAR_RESNAMES,
    PDB_SUGAR_RESNAMES,
    SIALIC_RESNAMES,
    is_glycam_sugar,
    is_sugar_resname,
)

# ---------------------------------------------------------------------------
# Sugar residue recognition
# ---------------------------------------------------------------------------

_PDB_SUGARS = PDB_SUGAR_RESNAMES
_CHARMM_SUGARS = CHARMM_SUGAR_RESNAMES
SIALIC_RESIDUES = SIALIC_RESNAMES
SUGAR_RESIDUES = _PDB_SUGARS | _CHARMM_SUGARS | SIALIC_RESIDUES


def _is_glycam_name(resname):
    return is_glycam_sugar(resname)


def _is_sugar(resname):
    return is_sugar_resname(resname)


def _is_sialic(resname):
    if resname in SIALIC_RESIDUES:
        return True
    if _is_glycam_name(resname) and resname[1] == 'S':
        return True
    return False


# ---------------------------------------------------------------------------
# Linkage detection
# ---------------------------------------------------------------------------

LinkageInfo = namedtuple('LinkageInfo', [
    'child_resindex', 'parent_resindex',
    'child_resname', 'parent_resname',
    'child_resid', 'parent_resid',
    'position',
    'phi_indices',
    'psi_indices',
    'omega_indices',
    'label',
])


def _find_atom_index(residue, name):
    try:
        return residue.atoms.select_atoms(f"name {name}")[0].index
    except (IndexError, Exception):
        return None


def detect_linkages(universe, select=None, verbose=False):
    """Detect glycosidic linkages from the first frame."""
    universe.trajectory[0]

    sugar_residues = []
    for res in universe.residues:
        if _is_sugar(res.resname):
            sugar_residues.append(res)

    if not sugar_residues:
        return []

    if verbose:
        resnames = sorted(set(r.resname for r in sugar_residues))
        print(f"Found {len(sugar_residues)} sugar residues: {', '.join(resnames)}")

    anomeric_atoms = []
    target_atoms = []

    for res in sugar_residues:
        sialic = _is_sialic(res.resname)
        anom_name = 'C2' if sialic else 'C1'
        anom_idx = _find_atom_index(res, anom_name)
        if anom_idx is not None:
            anomeric_atoms.append((universe.atoms[anom_idx], res, sialic))

        for aname in ['O1', 'O2', 'O3', 'O4', 'O6']:
            idx = _find_atom_index(res, aname)
            if idx is not None:
                target_atoms.append((universe.atoms[idx], res))

    for res in universe.residues:
        if res.resname in ('ASN', 'NLN'):
            idx = _find_atom_index(res, 'ND2')
            if idx is not None:
                target_atoms.append((universe.atoms[idx], res))
        elif res.resname in ('SER', 'OLS'):
            idx = _find_atom_index(res, 'OG')
            if idx is not None:
                target_atoms.append((universe.atoms[idx], res))
        elif res.resname in ('THR', 'OLT'):
            idx = _find_atom_index(res, 'OG1')
            if idx is not None:
                target_atoms.append((universe.atoms[idx], res))

    linkages = []
    for anom_atom, child_res, sialic in anomeric_atoms:
        anom_pos = anom_atom.position
        for tgt_atom, parent_res in target_atoms:
            if parent_res.resindex == child_res.resindex:
                continue
            d = np.linalg.norm(tgt_atom.position - anom_pos)
            if d < 2.0:
                lk = _resolve_linkage_atoms(
                    universe, child_res, parent_res,
                    anom_atom, tgt_atom, sialic, verbose)
                if lk is not None:
                    linkages.append(lk)

    seen = set()
    unique = []
    for lk in linkages:
        key = (lk.child_resindex, lk.parent_resindex)
        if key not in seen:
            seen.add(key)
            unique.append(lk)

    return unique


def _resolve_linkage_atoms(universe, child_res, parent_res,
                           anom_atom, tgt_atom, is_sialic, verbose):
    tgt_name = tgt_atom.name
    if tgt_name.startswith('O') and tgt_name[1:].isdigit():
        position = int(tgt_name[1:])
    elif tgt_name == 'ND2':
        position = -1
    elif tgt_name in ('OG', 'OG1'):
        position = -2
    else:
        return None

    ring_o_name = 'O6' if is_sialic else 'O5'
    ring_o_idx = _find_atom_index(child_res, ring_o_name)
    anom_c_idx = anom_atom.index
    glyco_o_idx = tgt_atom.index

    if ring_o_idx is None:
        if verbose:
            print(f"  WARNING: {child_res.resname}:{child_res.resid} missing {ring_o_name}")
        return None

    if position > 0:
        parent_cx_name = f'C{position}'
        parent_cx_idx = _find_atom_index(parent_res, parent_cx_name)
        if parent_cx_idx is None:
            if verbose:
                print(f"  WARNING: {parent_res.resname}:{parent_res.resid} "
                      f"missing {parent_cx_name}")
            return None
        parent_cx_m1_name = 'O5' if position == 1 else f'C{position - 1}'
        parent_cx_m1_idx = _find_atom_index(parent_res, parent_cx_m1_name)
        if parent_cx_m1_idx is None:
            if verbose:
                print(f"  WARNING: {parent_res.resname}:{parent_res.resid} "
                      f"missing {parent_cx_m1_name}")
            return None
    elif position == -1:
        parent_cx_idx = _find_atom_index(parent_res, 'CG')
        parent_cx_m1_idx = _find_atom_index(parent_res, 'CB')
        if parent_cx_idx is None or parent_cx_m1_idx is None:
            return None
    elif position == -2:
        parent_cx_idx = _find_atom_index(parent_res, 'CB')
        parent_cx_m1_idx = _find_atom_index(parent_res, 'CA')
        if parent_cx_idx is None or parent_cx_m1_idx is None:
            return None
    else:
        return None

    phi_indices = (ring_o_idx, anom_c_idx, glyco_o_idx, parent_cx_idx)
    psi_indices = (anom_c_idx, glyco_o_idx, parent_cx_idx, parent_cx_m1_idx)

    omega_indices = None
    if position == 6:
        parent_c5_idx = _find_atom_index(parent_res, 'C5')
        parent_o5_idx = _find_atom_index(parent_res, 'O5')
        if parent_c5_idx is not None and parent_o5_idx is not None:
            omega_indices = (glyco_o_idx, parent_cx_idx, parent_c5_idx, parent_o5_idx)

    pos_str = 'N-link' if position == -1 else 'O-link' if position == -2 else f'1->{position}'
    label = (f"{child_res.resname}:{child_res.resid}-"
             f"[{pos_str}]-"
             f"{parent_res.resname}:{parent_res.resid}")

    return LinkageInfo(
        child_resindex=child_res.resindex,
        parent_resindex=parent_res.resindex,
        child_resname=child_res.resname,
        parent_resname=parent_res.resname,
        child_resid=child_res.resid,
        parent_resid=parent_res.resid,
        position=position,
        phi_indices=phi_indices,
        psi_indices=psi_indices,
        omega_indices=omega_indices,
        label=label,
    )


# ---------------------------------------------------------------------------
# Torsion angle extraction
# ---------------------------------------------------------------------------

def extract_torsions(universe, linkages, stride=1, begin=None, end=None,
                     verbose=False):
    """Compute phi/psi/omega for all linkages across trajectory frames.

    Returns:
        torsion_data: ndarray (n_frames, n_torsions) in degrees
        torsion_labels: list of str labels
        frame_indices: list of actual frame indices used
        linkage_columns: list of (linkage_index, [col_indices]) mapping
    """
    from MDAnalysis.lib.distances import calc_dihedrals

    all_quads = []
    torsion_labels = []
    linkage_columns = []  # (linkage_idx, [col_indices])
    for li, lk in enumerate(linkages):
        cols = []
        cols.append(len(all_quads))
        all_quads.append(lk.phi_indices)
        torsion_labels.append(f"{lk.label} phi")
        cols.append(len(all_quads))
        all_quads.append(lk.psi_indices)
        torsion_labels.append(f"{lk.label} psi")
        if lk.omega_indices is not None:
            cols.append(len(all_quads))
            all_quads.append(lk.omega_indices)
            torsion_labels.append(f"{lk.label} omega")
        linkage_columns.append((li, cols))

    n_torsions = len(all_quads)
    if n_torsions == 0:
        return np.empty((0, 0)), [], [], []

    idx_a = np.array([q[0] for q in all_quads])
    idx_b = np.array([q[1] for q in all_quads])
    idx_c = np.array([q[2] for q in all_quads])
    idx_d = np.array([q[3] for q in all_quads])

    start = begin if begin is not None else 0
    stop = end if end is not None else len(universe.trajectory)
    frames_to_read = range(start, stop, stride)
    n_frames = len(frames_to_read)

    torsion_data = np.empty((n_frames, n_torsions), dtype=np.float64)
    frame_indices = []

    if verbose:
        print(f"Extracting {n_torsions} torsion angles from {n_frames} frames...")

    for fi, ts_idx in enumerate(frames_to_read):
        universe.trajectory[ts_idx]
        pos = universe.atoms.positions
        angles = calc_dihedrals(
            pos[idx_a], pos[idx_b], pos[idx_c], pos[idx_d])
        torsion_data[fi] = np.degrees(angles)
        frame_indices.append(ts_idx)

    return torsion_data, torsion_labels, frame_indices, linkage_columns


# ---------------------------------------------------------------------------
# Circular distance (arctan2 method — numerically robust)
# ---------------------------------------------------------------------------

def _circ_diff_rad(a_deg, b_deg):
    """Circular difference in radians (robust arctan2 method)."""
    a = np.radians(a_deg)
    b = np.radians(b_deg)
    return np.arctan2(np.sin(a - b), np.cos(a - b))


def compute_distance_matrix(torsion_data, chunk_size=5000, verbose=False):
    """Compute pairwise circular-RMSD distance matrix in degrees.

    Uses arctan2(sin, cos) for robust circular difference.
    """
    n = len(torsion_data)
    if n == 0:
        return np.empty((0, 0), dtype=np.float32)

    if verbose:
        print(f"Computing {n}x{n} distance matrix...")

    data_rad = np.radians(torsion_data)
    dist = np.zeros((n, n), dtype=np.float32)

    for i in range(n):
        if i < n - 1:
            diff = data_rad[i] - data_rad[i + 1:]
            circ_diff = np.arctan2(np.sin(diff), np.cos(diff))
            rmsd = np.degrees(np.sqrt(np.mean(circ_diff ** 2, axis=1)))
            dist[i, i + 1:] = rmsd
            dist[i + 1:, i] = rmsd

        if verbose and (i + 1) % 2000 == 0:
            print(f"  {i + 1}/{n} rows...")

    return dist


# ---------------------------------------------------------------------------
# GROMOS clustering (vectorized)
# ---------------------------------------------------------------------------

def gromos_cluster(dist_matrix, cutoff=30.0, max_clusters=100, verbose=False):
    """GROMOS-style clustering (vectorized numpy implementation).

    Returns:
        cluster_ids: ndarray (N,) cluster ID per frame (0-based)
        centers: list of frame indices (GROMOS centers)
        populations: list of cluster sizes
    """
    n = len(dist_matrix)
    available = np.ones(n, dtype=bool)
    cluster_ids = np.full(n, -1, dtype=np.int32)
    centers = []
    populations = []
    cluster_id = 0

    while np.any(available) and cluster_id < max_clusters:
        avail_idx = np.where(available)[0]
        sub_dist = dist_matrix[np.ix_(avail_idx, avail_idx)]
        neighbor_counts = np.sum(sub_dist <= cutoff, axis=1)

        best_local = np.argmax(neighbor_counts)
        neighbor_mask = sub_dist[best_local] <= cutoff
        members = avail_idx[neighbor_mask]

        cluster_ids[members] = cluster_id
        centers.append(int(avail_idx[best_local]))
        populations.append(len(members))
        available[members] = False

        if verbose:
            print(f"  Cluster {cluster_id}: center={avail_idx[best_local]}, "
                  f"members={len(members)}, remaining={np.sum(available)}")

        cluster_id += 1

    return cluster_ids, centers, populations


# ---------------------------------------------------------------------------
# Per-linkage clustering
# ---------------------------------------------------------------------------

def per_linkage_cluster(torsion_data, linkage_columns, linkages, cutoff=30.0,
                        verbose=False):
    """Cluster each linkage independently, then combine into compound states.

    Returns same format as gromos_cluster.
    """
    n_frames = len(torsion_data)
    n_linkages = len(linkage_columns)

    # Cluster each linkage independently
    per_lk_ids = np.zeros((n_frames, n_linkages), dtype=np.int32)

    for li, (lk_idx, cols) in enumerate(linkage_columns):
        lk_data = torsion_data[:, cols]  # just this linkage's torsions
        if verbose:
            print(f"  Clustering linkage {linkages[lk_idx].label} "
                  f"({len(cols)} torsions)...")
        dm = compute_distance_matrix(lk_data)
        cids, _, _ = gromos_cluster(dm, cutoff=cutoff)
        per_lk_ids[:, li] = cids

    # Combine: each unique combination of per-linkage cluster IDs = compound state
    # Convert each row to a tuple for grouping
    state_map = {}
    compound_ids = np.empty(n_frames, dtype=np.int32)
    next_id = 0

    for i in range(n_frames):
        state = tuple(per_lk_ids[i])
        if state not in state_map:
            state_map[state] = next_id
            next_id += 1
        compound_ids[i] = state_map[state]

    # Sort by population (largest first)
    unique_ids, counts = np.unique(compound_ids, return_counts=True)
    sorted_order = np.argsort(-counts)

    id_remap = {}
    for new_id, old_pos in enumerate(sorted_order):
        id_remap[unique_ids[old_pos]] = new_id

    cluster_ids = np.array([id_remap[c] for c in compound_ids], dtype=np.int32)
    n_clusters = len(unique_ids)

    # Find centers (frame closest to centroid per cluster)
    centers = []
    populations = []
    for cid in range(n_clusters):
        mask = cluster_ids == cid
        pop = np.sum(mask)
        populations.append(int(pop))

        cluster_data = torsion_data[mask]
        cluster_indices = np.where(mask)[0]
        mean_angles = _circular_mean(cluster_data)
        diffs = _circ_diff_rad(cluster_data, mean_angles[np.newaxis, :])
        rmsds = np.sqrt(np.mean(diffs ** 2, axis=1))
        centers.append(int(cluster_indices[np.argmin(rmsds)]))

    if verbose:
        print(f"  Per-linkage states: {[len(set(per_lk_ids[:, i])) for i in range(n_linkages)]}")
        print(f"  Compound states: {n_clusters}")

        # Show per-linkage breakdown
        inv_map = {v: k for k, v in state_map.items()}
        inv_remap = {v: k for k, v in id_remap.items()}
        for cid in range(min(n_clusters, 10)):
            old_id = inv_remap[cid]
            state = inv_map[old_id]
            state_str = ', '.join(
                f"{linkages[lk_idx].label}={s}"
                for (lk_idx, _), s in zip(linkage_columns, state))
            print(f"    State {cid} ({populations[cid]} frames): {state_str}")

    return cluster_ids, centers, populations


# ---------------------------------------------------------------------------
# Representatives and helpers
# ---------------------------------------------------------------------------

def _circular_mean(angles_deg):
    rad = np.radians(angles_deg)
    return np.degrees(np.arctan2(
        np.mean(np.sin(rad), axis=0),
        np.mean(np.cos(rad), axis=0)))


def _circular_std(angles_deg):
    rad = np.radians(angles_deg)
    R = np.sqrt(np.mean(np.sin(rad), axis=0)**2 +
                np.mean(np.cos(rad), axis=0)**2)
    return np.degrees(np.sqrt(-2.0 * np.log(np.clip(R, 1e-10, 1.0))))


def find_representatives(torsion_data, cluster_ids):
    """Find the medoid (frame closest to circular mean) for each cluster."""
    n_clusters = cluster_ids.max() + 1
    representatives = []
    for cid in range(n_clusters):
        mask = cluster_ids == cid
        cluster_data = torsion_data[mask]
        cluster_indices = np.where(mask)[0]
        mean_angles = _circular_mean(cluster_data)
        diffs = _circ_diff_rad(cluster_data, mean_angles[np.newaxis, :])
        rmsds = np.sqrt(np.mean(diffs ** 2, axis=1))
        representatives.append(int(cluster_indices[np.argmin(rmsds)]))
    return representatives


# ---------------------------------------------------------------------------
# Output functions
# ---------------------------------------------------------------------------

def write_torsions_csv(prefix, torsion_data, torsion_labels, frame_indices):
    path = f"{prefix}_torsions.csv"
    with open(path, 'w') as f:
        f.write('frame,' + ','.join(torsion_labels) + '\n')
        for i, fi in enumerate(frame_indices):
            vals = ','.join(f'{v:.2f}' for v in torsion_data[i])
            f.write(f'{fi},{vals}\n')
    return path


def write_clusters_csv(prefix, cluster_ids, centers, representatives,
                       frame_indices):
    path = f"{prefix}_clusters.csv"
    center_set = set(centers)
    rep_set = set(representatives)
    with open(path, 'w') as f:
        f.write('frame,cluster,is_gromos_center,is_representative\n')
        for i, fi in enumerate(frame_indices):
            f.write(f'{fi},{cluster_ids[i]},{i in center_set},{i in rep_set}\n')
    return path


def write_json_summary(prefix, linkages, torsion_data, torsion_labels,
                       cluster_ids, centers, populations, representatives,
                       frame_indices, args):
    path = f"{prefix}_summary.json"
    n_frames = len(torsion_data)
    n_clusters = len(centers)

    clusters_info = []
    for cid in range(n_clusters):
        mask = cluster_ids == cid
        cluster_data = torsion_data[mask]
        means = _circular_mean(cluster_data)
        stds = _circular_std(cluster_data)
        torsions = {}
        for ti, label in enumerate(torsion_labels):
            torsions[label] = {'mean': round(float(means[ti]), 1),
                               'std': round(float(stds[ti]), 1)}
        clusters_info.append({
            'id': cid,
            'population': int(populations[cid]),
            'fraction': round(populations[cid] / n_frames, 4),
            'center_frame': int(frame_indices[centers[cid]]),
            'representative_frame': int(frame_indices[representatives[cid]]),
            'torsion_angles': torsions,
        })

    summary = {
        'topology': args.topology,
        'trajectory': args.trajectory,
        'n_frames': n_frames,
        'n_linkages': len(linkages),
        'n_torsions': len(torsion_labels),
        'cutoff_degrees': args.cutoff,
        'mode': args.mode,
        'n_clusters': n_clusters,
        'linkages': [lk.label for lk in linkages],
        'clusters': clusters_info,
    }

    with open(path, 'w') as f:
        json.dump(summary, f, indent=2)
    return path


def write_summary(prefix, linkages, torsion_data, torsion_labels,
                  cluster_ids, centers, populations, representatives,
                  frame_indices, args):
    path = f"{prefix}_summary.txt"
    n_frames = len(torsion_data)
    n_clusters = len(centers)

    with open(path, 'w') as f:
        f.write("Glycan Conformational Clustering Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Topology:    {args.topology}\n")
        f.write(f"Trajectory:  {args.trajectory}\n")
        f.write(f"Frames:      {n_frames}")
        if args.stride > 1:
            f.write(f" (stride {args.stride})")
        f.write("\n")
        f.write(f"Mode:        {args.mode}\n")
        f.write(f"Linkages:    {len(linkages)}\n")
        f.write(f"Torsions:    {len(torsion_labels)}\n")
        f.write(f"Cutoff:      {args.cutoff:.1f} degrees\n")
        f.write(f"Clusters:    {n_clusters}\n\n")

        f.write("Detected linkages:\n")
        for lk in linkages:
            omega_str = " + omega" if lk.omega_indices else ""
            f.write(f"  {lk.label}  (phi + psi{omega_str})\n")
        f.write("\n")

        f.write(f"{'Cluster':<10s} {'Population':<12s} {'Fraction':<10s} "
                f"{'Representative':<15s}\n")
        f.write("-" * 47 + "\n")
        for cid in range(n_clusters):
            pop = populations[cid]
            frac = pop / n_frames
            rep = frame_indices[representatives[cid]]
            f.write(f"{cid:<10d} {pop:<12d} {frac:<10.3f} {rep:<15d}\n")
        f.write("\n")

        f.write("Average torsion angles per cluster (circular mean +/- std):\n\n")
        for cid in range(n_clusters):
            mask = cluster_ids == cid
            cluster_data = torsion_data[mask]
            f.write(f"Cluster {cid} ({populations[cid]} frames):\n")
            means = _circular_mean(cluster_data)
            stds = _circular_std(cluster_data)
            for ti, label in enumerate(torsion_labels):
                f.write(f"  {label:<40s}  {means[ti]:7.1f} +/- {stds[ti]:5.1f}\n")
            f.write("\n")

    return path


def _detect_align_resid(linkages, universe):
    """Auto-detect the best residue for alignment.

    Priority:
      1. Protein attachment residue (ASN/SER/THR) if protein-linked
      2. Root sugar (the sugar that is a parent but never a child) otherwise
    """
    # Check for protein-linked glycan
    for lk in linkages:
        if lk.position == -1:  # N-linked to ASN
            return lk.parent_resid, lk.parent_resname
        if lk.position == -2:  # O-linked to SER/THR
            return lk.parent_resid, lk.parent_resname

    # No protein link — find root sugar: appears as parent but never as child
    parent_resids = {lk.parent_resid for lk in linkages}
    child_resids = {lk.child_resid for lk in linkages}
    roots = parent_resids - child_resids
    if roots:
        root_resid = min(roots)  # pick lowest resid if multiple
        for lk in linkages:
            if lk.parent_resid == root_resid:
                return root_resid, lk.parent_resname

    # Fallback: first parent in first linkage
    if linkages:
        return linkages[0].parent_resid, linkages[0].parent_resname

    return None, None


def _align_to_reference(universe, ref_positions, align_select):
    """Align current frame to reference positions using Kabsch superposition."""
    mobile_atoms = universe.select_atoms(align_select)
    if len(mobile_atoms) == 0 or len(ref_positions) == 0:
        return
    if len(mobile_atoms) != len(ref_positions):
        return

    mobile_pos = mobile_atoms.positions.copy()
    ref_pos = ref_positions.copy()

    # Center both
    mobile_center = mobile_pos.mean(axis=0)
    ref_center = ref_pos.mean(axis=0)
    mobile_pos -= mobile_center
    ref_pos -= ref_center

    # Kabsch: find optimal rotation
    H = mobile_pos.T @ ref_pos
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T

    # Apply: translate to origin, rotate, translate to ref center
    all_pos = universe.atoms.positions
    all_pos -= mobile_center
    all_pos = all_pos @ R.T
    all_pos += ref_center
    universe.atoms.positions = all_pos


def write_representative_pdbs(universe, prefix, representatives,
                              frame_indices, linkages=None,
                              select=None, separate=False,
                              align_resid=None):
    if select:
        ag = universe.select_atoms(select)
    else:
        ag = universe.atoms

    # Determine alignment reference
    if align_resid is not None:
        ref_resid = align_resid
        ref_resname = '?'
    elif linkages:
        ref_resid, ref_resname = _detect_align_resid(linkages, universe)
    else:
        ref_resid = None
        ref_resname = None

    align_select = None
    ref_positions = None
    if ref_resid is not None:
        align_select = f"resid {ref_resid} and not name H*"
        # Load first representative as reference
        universe.trajectory[frame_indices[representatives[0]]]
        ref_atoms = universe.select_atoms(align_select)
        if len(ref_atoms) > 0:
            ref_positions = ref_atoms.positions.copy()

    def _write_atom_lines(f_out, atoms):
        for atom in atoms:
            r = atom.residue
            aname = atom.name
            if len(aname) < 4:
                aname_fmt = f" {aname:<3s}"
            else:
                aname_fmt = f"{aname:<4s}"
            resname = r.resname[:4]
            chain = r.segid[0] if r.segid else ' '
            x, y, z = atom.position
            f_out.write(
                f"ATOM  {atom.index + 1:5d} {aname_fmt}"
                f"{resname:>4s}{chain:1s}{r.resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00\n")

    if separate:
        paths = []
        for cid, rep_idx in enumerate(representatives):
            universe.trajectory[frame_indices[rep_idx]]
            if ref_positions is not None and cid > 0:
                _align_to_reference(universe, ref_positions, align_select)
            path = f"{prefix}_cluster_{cid}.pdb"
            ag.write(path)
            paths.append(path)
        return paths
    else:
        path = f"{prefix}_representatives.pdb"
        with open(path, 'w') as f:
            for cid, rep_idx in enumerate(representatives):
                frame_idx = frame_indices[rep_idx]
                universe.trajectory[frame_idx]
                if ref_positions is not None and cid > 0:
                    _align_to_reference(universe, ref_positions, align_select)
                f.write(f"MODEL     {cid + 1:>4d}\n")
                f.write(f"REMARK   Cluster {cid}, frame {frame_idx}\n")
                if ref_resid is not None:
                    f.write(f"REMARK   Aligned on resid {ref_resid}\n")
                _write_atom_lines(f, ag)
                f.write("ENDMDL\n")
            f.write("END\n")
        return [path]


# ---------------------------------------------------------------------------
# Plotting (plotly interactive HTML)
# ---------------------------------------------------------------------------

def _get_plotly():
    try:
        import plotly.graph_objects as go
        import plotly.subplots as sp
        return go, sp
    except ImportError:
        print("WARNING: plotly not available, skipping plots. "
              "Install with: pip install plotly", file=sys.stderr)
        return None, None


# Plotly qualitative color palette (tab10 equivalent)
_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
]


def plot_ramachandran(prefix, torsion_data, torsion_labels, cluster_ids,
                      linkages, linkage_columns, frame_indices):
    """Ramachandran scatter + free energy surface per linkage (interactive HTML)."""
    go, sp = _get_plotly()
    if go is None:
        return []

    n_clusters = int(cluster_ids.max() + 1)
    paths = []
    frames = np.array(frame_indices)

    for li, lk in enumerate(linkages):
        _, cols = linkage_columns[li]
        phi_col = cols[0]
        psi_col = cols[1]
        phi = torsion_data[:, phi_col]
        psi = torsion_data[:, psi_col]

        fig = sp.make_subplots(
            rows=1, cols=2,
            subplot_titles=[f'{lk.label} — by cluster',
                            f'{lk.label} — free energy surface'],
            horizontal_spacing=0.08)

        # Left: scatter colored by cluster
        for cid in range(n_clusters):
            mask = cluster_ids == cid
            if not np.any(mask):
                continue
            color = _COLORS[cid % len(_COLORS)]
            fig.add_trace(go.Scattergl(
                x=phi[mask], y=psi[mask],
                mode='markers',
                marker=dict(size=3, color=color, opacity=0.4),
                name=f'Cl {cid} ({np.sum(mask)})',
                text=[f'Frame {f}' for f in frames[mask]],
                hovertemplate='φ=%{x:.1f}° ψ=%{y:.1f}°<br>%{text}',
            ), row=1, col=1)

        fig.update_xaxes(title_text='φ (deg)', range=[-180, 180],
                         row=1, col=1)
        fig.update_yaxes(title_text='ψ (deg)', range=[-180, 180],
                         scaleanchor='x', scaleratio=1, row=1, col=1)

        # Right: free energy surface as heatmap
        h, xedges, yedges = np.histogram2d(
            phi, psi, bins=72, range=[[-180, 180], [-180, 180]])
        h = h.T
        with np.errstate(divide='ignore'):
            fe = -np.log(h / h.sum())
        fe[np.isinf(fe)] = np.nan
        vmax = float(np.nanpercentile(fe, 95))

        fig.add_trace(go.Heatmap(
            z=fe, x0=-180, dx=360/72, y0=-180, dy=360/72,
            colorscale='RdYlBu_r', zmin=0, zmax=vmax,
            colorbar=dict(title='-ln P (a.u.)', x=1.0),
            hovertemplate='φ=%{x:.0f}° ψ=%{y:.0f}°<br>-ln P=%{z:.2f}',
        ), row=1, col=2)

        fig.update_xaxes(title_text='φ (deg)', range=[-180, 180],
                         row=1, col=2)
        fig.update_yaxes(title_text='ψ (deg)', range=[-180, 180],
                         scaleanchor='x2', scaleratio=1, row=1, col=2)

        fig.update_layout(
            height=550, width=1200,
            title_text=lk.label,
            showlegend=True,
            legend=dict(font_size=10),
        )

        safe_label = lk.label.replace(':', '_').replace('->', '_')
        path = f"{prefix}_rama_{safe_label}.html"
        fig.write_html(path)
        paths.append(path)

    return paths


def plot_timeseries(prefix, torsion_data, torsion_labels, cluster_ids,
                    frame_indices):
    """Time series of each dihedral colored by cluster (interactive HTML)."""
    go, sp = _get_plotly()
    if go is None:
        return None

    n_dihedrals = len(torsion_labels)
    n_clusters = int(cluster_ids.max() + 1)
    frames = np.array(frame_indices)

    fig = sp.make_subplots(
        rows=n_dihedrals, cols=1,
        shared_xaxes=True,
        subplot_titles=torsion_labels,
        vertical_spacing=0.02)

    # Only add legend entries once (first subplot)
    legend_added = set()

    for j in range(n_dihedrals):
        for cid in range(n_clusters):
            mask = cluster_ids == cid
            if not np.any(mask):
                continue
            color = _COLORS[cid % len(_COLORS)]
            show_legend = cid not in legend_added
            fig.add_trace(go.Scattergl(
                x=frames[mask], y=torsion_data[mask, j],
                mode='markers',
                marker=dict(size=2, color=color, opacity=0.4),
                name=f'Cl {cid}',
                showlegend=show_legend,
                legendgroup=f'cl{cid}',
                hovertemplate=f'{torsion_labels[j]}<br>'
                              'Frame %{x}<br>%{y:.1f}°',
            ), row=j + 1, col=1)
            legend_added.add(cid)

        fig.update_yaxes(title_text='deg', range=[-180, 180],
                         row=j + 1, col=1)

    fig.update_xaxes(title_text='Frame', row=n_dihedrals, col=1)
    fig.update_layout(
        height=200 * n_dihedrals + 100,
        width=1200,
        title_text='Torsion angle time series',
        showlegend=True,
    )

    path = f"{prefix}_timeseries.html"
    fig.write_html(path)
    return path


def plot_populations(prefix, populations):
    """Bar chart of cluster populations with cumulative line (interactive HTML)."""
    go, _ = _get_plotly()
    if go is None:
        return None

    n = len(populations)
    total = sum(populations)
    fracs = [p / total * 100 for p in populations]
    cum = list(np.cumsum(fracs))
    x = list(range(n))

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=x, y=fracs,
        marker_color='steelblue',
        text=[str(p) for p in populations],
        textposition='outside',
        name='Population',
        hovertemplate='Cluster %{x}<br>%{y:.1f}% (%{text} frames)',
    ))

    fig.add_trace(go.Scatter(
        x=x, y=cum,
        mode='lines+markers',
        marker=dict(color='red', size=6),
        line=dict(color='red'),
        name='Cumulative',
        yaxis='y2',
        hovertemplate='Cluster %{x}<br>Cumulative: %{y:.1f}%',
    ))

    fig.update_layout(
        title='Cluster Populations',
        xaxis_title='Cluster ID',
        yaxis_title='Population (%)',
        yaxis2=dict(title='Cumulative (%)', overlaying='y', side='right',
                    range=[0, 105]),
        width=max(600, n * 40 + 200),
        height=500,
        showlegend=True,
    )

    path = f"{prefix}_populations.html"
    fig.write_html(path)
    return path


def plot_cutoff_scan(prefix, dist_matrix, selected_cutoff):
    """Scan cutoff values and plot #clusters vs cutoff (interactive HTML)."""
    go, _ = _get_plotly()
    if go is None:
        return None

    cutoffs = list(np.arange(5, 91, 5, dtype=float))
    n_clusters_list = []
    for c in cutoffs:
        cids, _, _ = gromos_cluster(dist_matrix, cutoff=c, max_clusters=200)
        n_clusters_list.append(int(cids.max() + 1) if len(cids) > 0 else 0)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=cutoffs, y=n_clusters_list,
        mode='lines+markers',
        marker=dict(color='steelblue', size=8),
        line=dict(color='steelblue'),
        hovertemplate='Cutoff: %{x}°<br>Clusters: %{y}',
    ))

    fig.add_vline(x=selected_cutoff, line_dash='dash', line_color='red',
                  annotation_text=f'Selected: {selected_cutoff}°')

    fig.update_layout(
        title='Cutoff scan — #clusters vs cutoff',
        xaxis_title='Cutoff (degrees)',
        yaxis_title='Number of clusters',
        width=700, height=450,
    )

    path = f"{prefix}_cutoff_scan.html"
    fig.write_html(path)
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog='dvbfixer cluster',
        description='Cluster glycan conformations from MD trajectory '
                    'using glycosidic torsion angle RMSD (GFDB method).',
    )
    io = p.add_argument_group('Input / output')
    io.add_argument('topology',
                    help='Topology file (.tpr, .pdb, .gro)')
    io.add_argument('trajectory',
                    help='Trajectory file (.xtc, .trr, .dcd)')
    io.add_argument('-o', '--output', default=None,
                    help='Output prefix (default: trajectory stem)')

    frames = p.add_argument_group('Frame selection')
    frames.add_argument('--stride', type=positive_int, default=1,
                        help='Read every Nth frame (default: 1)')
    frames.add_argument('--begin', type=nonnegative_int, default=None,
                        help='First frame (0-based)')
    frames.add_argument('--end', type=nonnegative_int, default=None,
                        help='Last frame (exclusive)')
    frames.add_argument('--select', default=None,
                        help='MDAnalysis selection for output PDB atoms')

    clustering = p.add_argument_group('Clustering')
    clustering.add_argument('--cutoff', type=positive_float, default=30.0,
                            help='RMSD cutoff in degrees (default: 30.0)')
    clustering.add_argument('--mode', choices=['global', 'per-linkage'],
                            default='per-linkage',
                            help='Clustering mode: global (all torsions at once) or '
                                 'per-linkage (each linkage independently, then combine)')

    representative = p.add_argument_group('Representative PDBs')
    representative.add_argument('--align-resid', type=int, default=None,
                                help='Residue ID to align representative PDBs on '
                                     '(default: auto-detect protein attachment or root sugar)')
    representative.add_argument('--no-align', action='store_true',
                                help='Disable alignment of representative PDBs')
    representative.add_argument('--separate-pdb', action='store_true',
                                help='Write each cluster as separate PDB '
                                     '(default: multi-MODEL PDB)')

    diag = p.add_argument_group('Diagnostics')
    diag.add_argument('--plot', action='store_true',
                      help='Generate interactive HTML plots (requires plotly)')
    diag.add_argument('-v', '--verbose', action='store_true',
                      help='Verbose output')

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    prefix = args.output or Path(args.trajectory).stem

    # 1. Load universe
    import MDAnalysis as mda
    if args.verbose:
        print(f"Loading {args.topology} + {args.trajectory}...")
    u = mda.Universe(args.topology, args.trajectory)
    if args.verbose:
        print(f"  {len(u.atoms)} atoms, {len(u.trajectory)} frames")

    # 2. Detect linkages
    linkages = detect_linkages(u, select=args.select, verbose=args.verbose)
    if not linkages:
        resnames = sorted(set(r.resname for r in u.residues))
        print(f"ERROR: No glycosidic linkages detected.\n"
              f"  Residue names: {', '.join(resnames[:30])}",
              file=sys.stderr)
        sys.exit(1)

    print(f"Detected {len(linkages)} glycosidic linkage(s):")
    for lk in linkages:
        omega_str = " + omega" if lk.omega_indices else ""
        print(f"  {lk.label}  (phi + psi{omega_str})")

    # 3. Extract torsions
    torsion_data, torsion_labels, frame_indices, linkage_columns = \
        extract_torsions(u, linkages,
                         stride=args.stride, begin=args.begin, end=args.end,
                         verbose=args.verbose)

    n_frames, n_torsions = torsion_data.shape
    print(f"Extracted {n_torsions} torsion angles from {n_frames} frames")

    if n_frames < 2:
        print("Only 1 frame — writing torsion angles, skipping clustering.")
        write_torsions_csv(prefix, torsion_data, torsion_labels, frame_indices)
        return

    # 4. Cluster
    print(f"Clustering mode: {args.mode}, cutoff: {args.cutoff}°")
    if args.mode == 'per-linkage':
        cluster_ids, centers, populations = per_linkage_cluster(
            torsion_data, linkage_columns, linkages,
            cutoff=args.cutoff, verbose=args.verbose)
    else:
        dist_matrix = compute_distance_matrix(
            torsion_data, verbose=args.verbose)
        cluster_ids, centers, populations = gromos_cluster(
            dist_matrix, cutoff=args.cutoff, verbose=args.verbose)
        # Save distance matrix
        np.save(f"{prefix}_distmatrix.npy", dist_matrix)

    n_clusters = len(centers)
    print(f"Found {n_clusters} cluster(s):")
    for cid in range(min(n_clusters, 20)):
        frac = populations[cid] / n_frames * 100
        print(f"  Cluster {cid}: {populations[cid]} frames ({frac:.1f}%)")
    if n_clusters > 20:
        print(f"  ... and {n_clusters - 20} more")

    # 5. Representatives
    representatives = find_representatives(torsion_data, cluster_ids)

    # 6. Write outputs
    p = write_torsions_csv(prefix, torsion_data, torsion_labels, frame_indices)
    print(f"Wrote {p}")

    p = write_clusters_csv(prefix, cluster_ids, centers, representatives,
                           frame_indices)
    print(f"Wrote {p}")

    p = write_summary(prefix, linkages, torsion_data, torsion_labels,
                      cluster_ids, centers, populations, representatives,
                      frame_indices, args)
    print(f"Wrote {p}")

    p = write_json_summary(prefix, linkages, torsion_data, torsion_labels,
                           cluster_ids, centers, populations, representatives,
                           frame_indices, args)
    print(f"Wrote {p}")

    paths = write_representative_pdbs(
        u, prefix, representatives, frame_indices,
        linkages=None if args.no_align else linkages,
        select=args.select, separate=args.separate_pdb,
        align_resid=None if args.no_align else args.align_resid)
    for p in paths:
        print(f"Wrote {p}")

    # 7. Plots
    if args.plot:
        for p in plot_ramachandran(prefix, torsion_data, torsion_labels,
                                   cluster_ids, linkages, linkage_columns,
                                   frame_indices):
            print(f"Wrote {p}")

        p = plot_timeseries(prefix, torsion_data, torsion_labels,
                            cluster_ids, frame_indices)
        if p:
            print(f"Wrote {p}")

        p = plot_populations(prefix, populations)
        if p:
            print(f"Wrote {p}")

        if args.mode == 'global':
            p = plot_cutoff_scan(prefix, dist_matrix, args.cutoff)
            if p:
                print(f"Wrote {p}")


if __name__ == "__main__":
    main()
