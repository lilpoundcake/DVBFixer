"""Chemical-exclusion geometry metrics matching the preceding 9CQD comparison.

No distance-only interchain bonds are ever added to the clash exclusion graph.
"""

from __future__ import annotations

import io
from dataclasses import asdict

import numpy as np
from openmm.app import PDBFile
from openmm.unit import angstrom

from dvbfixer.diagnose.geometry import VDW_RADII_ANGSTROM
from dvbfixer.diagnose.steric import _expand_exclusions
from dvbfixer.model.diffusion.contract import AtomIdentity


def segment_distance(p0, p1, q0, q1):
    """Minimum distance on closed segments, including zero-length segments."""
    p0, p1, q0, q1 = (np.asarray(x, dtype=float) for x in (p0, p1, q0, q1))
    u, v, w = p1 - p0, q1 - q0, p0 - q0
    a, b, c, d, e = float(u @ u), float(u @ v), float(v @ v), float(u @ w), float(v @ w)
    if a == 0 and c == 0:
        return float(np.linalg.norm(w))
    if a == 0:
        return float(np.linalg.norm(w - np.clip(e / c, 0, 1) * v))
    if c == 0:
        return float(np.linalg.norm(w + np.clip(-d / a, 0, 1) * u))
    candidates = [
        (0, np.clip(e / c, 0, 1)),
        (1, np.clip((e + b) / c, 0, 1)),
        (np.clip(-d / a, 0, 1), 0),
        (np.clip((b - d) / a, 0, 1), 1),
    ]
    den = a * c - b * b
    if den > np.finfo(float).eps * a * c * 16:
        s, t = (b * e - c * d) / den, (a * e - b * d) / den
        if 0 <= s <= 1 and 0 <= t <= 1:
            candidates.append((s, t))
    return min(float(np.linalg.norm(w + s * u - t * v)) for s, t in candidates)


def atom_id(atom):
    r = atom.residue
    return AtomIdentity(r.chain.id, r.id, r.insertionCode.strip(), atom.name)


def score(pdb_text, source_text, request):
    pdb = PDBFile(io.StringIO(pdb_text))
    source = PDBFile(io.StringIO(source_text))
    atoms = list(pdb.topology.atoms())
    axis = tuple(atom_id(a) for a in atoms)
    if len(set(axis)) != len(axis):
        raise ValueError("duplicate output identity")
    xyz = np.asarray(pdb.positions.value_in_unit(angstrom))
    source_xyz = {
        atom_id(a): np.asarray(p.value_in_unit(angstrom))
        for a, p in zip(source.topology.atoms(), source.positions)
    }
    indexes = {a: i for i, a in enumerate(axis)}
    generated, fixed = set(request.generated_atoms), set(request.fixed_atoms)
    if set(axis) != fixed | generated:
        raise ValueError("output identity mismatch")
    source_bonds = {frozenset((atom_id(a), atom_id(b))) for a, b in source.topology.bonds()}
    explicit = {frozenset((link.atom1, link.atom2)) for link in request.retained_explicit_links}
    positions = {}
    for p in request.sequence_placements:
        positions.update(
            {
                (r.chain, r.residue_number, r.insertion_code): i
                for i, r in zip(p.observed_target_indices, p.observed_residues)
            }
        )
    for gap in request.gaps:
        positions.update(
            {
                (r.chain, r.residue_number, r.insertion_code): i
                for i, r in zip(
                    range(gap.target_interval.start, gap.target_interval.stop),
                    gap.generated_residues,
                )
            }
        )
    edges = set()
    for a, b in pdb.topology.bonds():
        aa, bb = atom_id(a), atom_id(b)
        same = (aa.chain, aa.residue_number, aa.insertion_code) == (
            bb.chain,
            bb.residue_number,
            bb.insertion_code,
        )
        pair = frozenset((aa, bb))
        peptide = False
        if aa.chain == bb.chain and {aa.atom_name, bb.atom_name} == {"C", "N"}:
            carbon, nitrogen = (aa, bb) if aa.atom_name == "C" else (bb, aa)
            ci = positions.get((carbon.chain, carbon.residue_number, carbon.insertion_code))
            ni = positions.get((nitrogen.chain, nitrogen.residue_number, nitrogen.insertion_code))
            peptide = ci is not None and ni == ci + 1
        valid_source = pair in source_bonds and (
            same or (peptide and float(np.linalg.norm(xyz[a.index] - xyz[b.index])) <= 1.8)
        )
        if (
            pair in explicit
            or valid_source
            or ((aa in generated or bb in generated) and (same or peptide))
        ):
            edges.add(frozenset((a.index, b.index)))
    for link in request.retained_explicit_links:
        if link.atom1 in indexes and link.atom2 in indexes:
            edges.add(frozenset((indexes[link.atom1], indexes[link.atom2])))
    neighbors = {}
    for edge in edges:
        i, j = tuple(edge)
        neighbors.setdefault(i, set()).add(j)
        neighbors.setdefault(j, set()).add(i)
    exclusions = _expand_exclusions(neighbors)
    gi = [indexes[a] for a in axis if a in generated]
    fi = [indexes[a] for a in axis if a in fixed]
    result = {}
    for label, selected in [("full_source", fi), ("X", [i for i in fi if axis[i].chain == "X"])]:
        distances = np.linalg.norm(xyz[gi, None, :] - xyz[None, selected, :], axis=-1)
        where = np.unravel_index(int(np.argmin(distances)), distances.shape)
        result[label] = {
            "minimum_distance_angstrom": float(distances[where]),
            "minimum_pair": [asdict(axis[gi[where[0]]]), asdict(axis[selected[where[1]]])],
            "pairs_below_1_9_angstrom": int(np.sum(distances < 1.9)),
        }
    radii = [VDW_RADII_ANGSTROM.get(a.element.symbol, 1.7) for a in atoms]
    severe = []
    for i in gi:
        for j in range(len(axis)):
            if i == j or (axis[j] in generated and j < i) or frozenset((i, j)) in exclusions:
                continue
            distance = float(np.linalg.norm(xyz[i] - xyz[j]))
            overlap = radii[i] + radii[j] - distance
            if overlap >= 0.9:
                severe.append(
                    {
                        "pair": [asdict(axis[i]), asdict(axis[j])],
                        "distance_angstrom": distance,
                        "overlap_angstrom": overlap,
                    }
                )
    backbone = {"N", "CA", "C"}
    bb_edges = [tuple(e) for e in edges if all(axis[i].atom_name in backbone for i in e)]
    gs = [e for e in bb_edges if any(axis[i] in generated for i in e)]
    fs = [e for e in bb_edges if all(axis[i] in fixed for i in e)]
    for label, subset in [
        ("full_source", fs),
        ("X", [e for e in fs if all(axis[i].chain == "X" for i in e)]),
    ]:
        minimum, pair = None, None
        for a, b in gs:
            for c, d in subset:
                if any(i == j or frozenset((i, j)) in exclusions for i in (a, b) for j in (c, d)):
                    continue
                distance = segment_distance(xyz[a], xyz[b], xyz[c], xyz[d])
                if minimum is None or distance < minimum:
                    minimum, pair = distance, [a, b, c, d]
        result[label].update(
            backbone_segment_minimum_angstrom=minimum,
            near_crossing_alarm=minimum is not None and minimum < 0.5,
            segment_pair=[asdict(axis[i]) for i in pair] if pair else None,
        )
    drift = np.asarray([xyz[indexes[a]] - source_xyz[a] for a in request.fixed_atoms])
    result.update(
        authoritative_severe_clash_count=len(severe),
        authoritative_severe_clashes=severe,
        fixed_heavy_rmsd_angstrom=float(np.sqrt(np.mean(np.sum(drift * drift, axis=1)))),
        fixed_heavy_max_displacement_angstrom=float(np.linalg.norm(drift, axis=1).max()),
        authoritative_edge_count=len(edges),
        generated_interchain_covalent_edges=sum(
            1
            for e in edges
            if len({axis[i].chain for i in e}) > 1 and any(axis[i] in generated for i in e)
        ),
    )
    return result
