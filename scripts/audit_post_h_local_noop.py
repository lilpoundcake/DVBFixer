"""Independent NumPy audit of local policy outputs, without classifier calls."""
from __future__ import annotations

from dataclasses import asdict

from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def audit(value, proposal):
    import numpy as np

    # Input validity is established independently of proposal generation as well.
    atoms = value.atoms
    p = np.asarray(value.positions_nm, dtype=float)
    q = np.asarray(proposal.result.positions_nm, dtype=float)
    if p.shape != q.shape or not np.all(np.isfinite(q)):
        raise ValueError('audit invalid output inventory/coordinates')
    ids = [a.identity for a in atoms]
    if len(ids) != len(set(ids)) or proposal.original_new_hydrogens != value.new_hydrogens:
        raise ValueError('audit identity/ownership mismatch')
    adjacent = [[] for _ in atoms]
    for x, y in value.bonds:
        adjacent[x].append(y)
        adjacent[y].append(x)
    supported, unsupported, interventions, warned = [], [], [], []
    wanted = p.copy()
    minimum_contact, radius_error = None, 0.

    def volume(origin, neighbors):
        v = p[neighbors]-origin
        divisor = float(np.prod(np.linalg.norm(v, axis=1)))
        if divisor <= 1e-30:
            raise ValueError('audit degenerate vectors')
        return float(np.linalg.det(v))/divisor

    for h, atom in enumerate(atoms):
        if atom.identity not in value.new_hydrogens:
            continue
        if atom.element != 'H' or len(adjacent[h]) != 1:
            raise ValueError('audit invalid original new H')
        parent = adjacent[h][0]
        center = atoms[parent]
        heavy = [i for i in adjacent[parent] if atoms[i].element not in ('H', 'D')]
        hs = [i for i in adjacent[parent] if atoms[i].element in ('H', 'D')]
        inventory = CANONICAL_HEAVY_ATOMS.get(VARIANT_TO_PARENT.get(center.residue_name, center.residue_name), ())
        if (center.element != 'C' or len(heavy) != 3 or hs != [h] or center.identity.name not in inventory
                or any(atoms[i].identity.component != center.identity.component for i in [parent, h]+heavy)
                or any(atoms[i].identity.name not in inventory for i in heavy)):
            unsupported.append(atom.identity)
            continue
        supported.append(atom.identity)
        heavy.sort(key=lambda i: atoms[i].identity)
        hv, fv = volume(p[parent], heavy), volume(p[h], heavy)
        if abs(hv) <= 1e-4 or abs(fv) <= 1e-4:
            raise ValueError('audit unclassifiable supported center')
        radius = float(np.linalg.norm(p[h]-p[parent]))
        distances = np.linalg.norm(p[heavy]-p[parent], axis=1)
        admitted = 0.05 <= radius <= 0.15 and bool(np.all((distances >= 0.10) & (distances <= 0.20)))
        if hv*fv > 0:
            if not admitted:
                warned.append((atom.identity, center.residue_name, radius, distances))
            continue
        if not admitted:
            raise ValueError('audit moved center outside strict admission')
        interventions.append(h)
        direction = np.sum(p[parent]-p[heavy], axis=0)
        length = float(np.linalg.norm(direction))
        if length <= 1e-12:
            raise ValueError('audit zero heavy sum')
        wanted[h] = p[parent]+radius*direction/length
        cv = volume(q[h], heavy)
        if abs(cv) <= 1e-4 or hv*cv <= 0:
            raise ValueError('audit incompatible proposed sign')
        radius_error = max(radius_error, abs(float(np.linalg.norm(q[h]-q[parent]))-radius))
    changed = np.where(np.any(p != q, axis=1))[0].tolist()
    if changed != interventions or proposal.result.moved != tuple(ids[i] for i in interventions):
        raise ValueError('audit exact changed-atom mask mismatch')
    if proposal.intervention_hydrogens != frozenset(ids[i] for i in interventions):
        raise ValueError('audit intervention mask mismatch')
    if proposal.result.supported != tuple(supported) or proposal.result.unsupported != tuple(unsupported):
        raise ValueError('audit support classification mismatch')
    if float(np.max(np.abs(q-wanted))) > 1e-14 or radius_error > 1e-14:
        raise ValueError('audit strict formula/radius mismatch')
    if len(proposal.warnings) != len(warned):
        raise ValueError('audit warning count mismatch')
    for w, (identity, residue, radius, distances) in zip(proposal.warnings, warned, strict=True):
        if (w.atom != identity or w.residue_name != residue
                or abs(w.radius_nm-radius) > 1e-14
                or not np.allclose(w.heavy_bond_lengths_nm, distances, atol=1e-14, rtol=0.)
                or w.reason != 'compatible_noop_outside_distance_admission_not_geometry_certified'):
            raise ValueError('audit warning identity/distance mismatch')
    for h in interventions:
        parent = adjacent[h][0]
        others = [i for i in range(len(atoms)) if i not in (h, parent)]
        contact = float(np.min(np.linalg.norm(q[others]-q[h], axis=1)))
        if contact < 0.05:
            raise ValueError('audit final nonparent contact')
        minimum_contact = contact if minimum_contact is None else min(contact, minimum_contact)
    return dict(passed=True, original_new_h_count=len(value.new_hydrogens), supported=len(supported),
        unsupported=len(unsupported), moved=[asdict(ids[i]) for i in interventions],
        compatible_distance_warnings=len(warned), maximum_radius_error_nm=radius_error,
        minimum_moved_nonparent_distance_nm=minimum_contact,
        observations_and_noops_exact=True)
