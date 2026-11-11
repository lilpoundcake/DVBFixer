"""Private small-covalent-ring surface diagnostics, not topology qualification.

Ring membership and path bonds must come from explicit topology, not a geometric
guess or duplicated residue templates. A bounded centroid fan is still a chosen
surface: small nonplanar rings are reported undefined, and even planar piercings
do not establish a knot invariant or an unsaved optimization trajectory.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from dvbfixer.model.diffusion.backbone_intersections import BackboneSegment
from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.loop_surface_probe import probe_loop_surface


@dataclass(frozen=True)
class RingPiercing:
    bond: tuple[AtomIdentity, AtomIdentity]
    fraction: float
    point: tuple[float, float, float]
    surface_edge_ambiguous: bool
    ownership: str


@dataclass(frozen=True)
class CovalentRingProbe:
    ring: tuple[AtomIdentity, ...]
    piercings: tuple[RingPiercing, ...]
    planarity_max_angstrom: float
    undefined_reason: str | None
    path_bonds_tested: int
    method: str = 'experimental-topology-bound-small-ring-fan-v1'
    publication_gate: bool = False
    topology_or_trajectory_qualified: bool = False


def probe_covalent_ring(
    ring: tuple[AtomIdentity, ...],
    coordinates: Mapping[AtomIdentity, tuple[float, float, float]],
    bonds: frozenset[frozenset[AtomIdentity]],
    generated: frozenset[AtomIdentity],
    *,
    planarity_max_angstrom: float = .10,
    max_bonds: int = 50_000,
    max_tests: int = 200_000,
) -> CovalentRingProbe:
    """Probe explicit 3–8-member cycles against non-neighbor covalent bonds.

    Only geometric piercing is tested, not chemical energy or legitimate native
    threading. The engineering planarity limit is not an acceptance threshold.
    Missing identities, invalid topology and bounded-work excess fail closed.
    """
    if not 3 <= len(ring) <= 8 or len(set(ring)) != len(ring):
        raise ValueError('declare a unique ordered 3–8-member covalent ring')
    if not np.isfinite(planarity_max_angstrom) or planarity_max_angstrom <= 0:
        raise ValueError('invalid diagnostic planarity bound')
    if max_bonds < 1 or max_tests < 1 or len(bonds) > max_bonds:
        raise ValueError('covalent ring diagnostic budget exceeded')
    if any(len(bond) != 2 for bond in bonds):
        raise ValueError('topology edges must have two distinct atom identities')
    ring_edges = {frozenset((a, b)) for a, b in zip(ring, ring[1:]+ring[:1])}
    if not ring_edges <= bonds:
        raise ValueError('ring is not a closed cycle in the supplied covalent topology')
    if any(a not in coordinates for edge in bonds for a in edge):
        raise ValueError('covalent topology references missing coordinates')
    topology_atoms = {a for edge in bonds for a in edge}
    all_points = np.asarray([coordinates[a] for a in topology_atoms], dtype=float)
    if all_points.shape != (len(topology_atoms), 3) or not np.isfinite(all_points).all():
        raise ValueError('topology coordinates must be finite three-vectors')
    points = np.asarray([coordinates[a] for a in ring], dtype=float)
    if points.shape != (len(ring), 3) or not np.isfinite(points).all():
        raise ValueError('ring coordinates must be finite three-vectors')
    centered = points-points.mean(axis=0)
    _, singular, right = np.linalg.svd(centered)
    planarity = float(np.abs(centered @ right[-1]).max())
    if singular[1] < 1e-8:
        return CovalentRingProbe(ring, (), planarity, 'degenerate ring surface', 0)
    if planarity > planarity_max_angstrom:
        return CovalentRingProbe(ring, (), planarity, 'nonplanar ring exceeds declared surface diagnostic scope', 0)
    ring_set = set(ring)
    neighbors = {a for edge in bonds if edge & ring_set for a in edge} | ring_set
    paths = []
    pairs = []
    lower, upper = points.min(axis=0)-1e-5, points.max(axis=0)+1e-5
    for edge in sorted(bonds, key=lambda pair: tuple(sorted(pair))):
        if edge & neighbors:
            continue
        a, b = sorted(edge)
        endpoints = np.asarray((coordinates[a], coordinates[b]), dtype=float)
        if endpoints.shape != (2, 3) or not np.isfinite(endpoints).all():
            raise ValueError('bond coordinates must be finite three-vectors')
        if np.any(endpoints.max(axis=0) < lower) or np.any(endpoints.min(axis=0) > upper):
            continue
        pairs.append((a, b))
        # The existing segment primitive has no protein-template semantics;
        # these endpoints are a verified sidechain/other covalent edge too.
        paths.append(BackboneSegment(a, b, tuple(endpoints[0]), tuple(endpoints[1])))
    result = probe_loop_surface(points, tuple(paths), max_tests=max_tests)
    findings = []
    for hit in result.piercings:
        pair = pairs[hit.segment_index]
        ring_mobile, bond_mobile = bool(ring_set & generated), bool(set(pair) & generated)
        ownership = ('generated-generated' if ring_mobile and bond_mobile else
                     'generated-observed' if ring_mobile or bond_mobile else 'observed-observed')
        findings.append(RingPiercing(pair, hit.fraction, hit.point, hit.surface_edge_ambiguous, ownership))
    return CovalentRingProbe(ring, tuple(findings), planarity, None, len(paths))
