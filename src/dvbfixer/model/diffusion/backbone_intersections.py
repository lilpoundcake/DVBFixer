"""Experimental 3D segment alarms, not knot detection or a public fatal gate.

Callers supply explicit backbone bonds/continuity; this primitive never connects
nearby residues by proximity or infers topology from a drawing. Source-only
findings stay distinguishable. Open-path threading needs a separate declared
closure/surface convention and is NOT established by a segment-distance alarm.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from dvbfixer.model.diffusion.contract import AtomIdentity

METHOD = "experimental-backbone-segment-screen-v1"


@dataclass(frozen=True)
class BackboneSegment:
    first: AtomIdentity
    second: AtomIdentity
    start: tuple[float, float, float]
    stop: tuple[float, float, float]

    def __post_init__(self):
        points = np.asarray((self.start, self.stop), dtype=float)
        if self.first == self.second or points.shape != (2, 3) or not np.isfinite(points).all():
            raise ValueError("segment requires distinct identities and finite 3D endpoints")
        if np.linalg.norm(points[1] - points[0]) < 1e-8:
            raise ValueError("zero-length backbone segment")


@dataclass(frozen=True)
class SegmentAlarm:
    first_index: int
    second_index: int
    distance_angstrom: float
    first_fraction: float
    second_fraction: float
    ownership: str
    category: str
    topological_defect_established: bool = False


@dataclass(frozen=True)
class SegmentScreen:
    alarms: tuple[SegmentAlarm, ...]
    candidate_pairs: int
    threshold_angstrom: float
    intersection_tolerance_angstrom: float
    method: str = METHOD
    publication_gate: bool = False
    threading_tested: bool = False


def segment_distance(first: BackboneSegment, second: BackboneSegment):
    """Exact minimum on two closed line segments, including parallel endpoints."""
    a, b, c, d = (np.asarray(v, dtype=float) for v in (first.start, first.stop, second.start, second.stop))
    u, v, w = b - a, d - c, a - c
    uu, vv, uv = float(u @ u), float(v @ v), float(u @ v)
    candidates = []
    for s in (0., 1.):
        candidates.append((s, float(np.clip((a + s*u - c) @ v / vv, 0, 1))))
    for t in (0., 1.):
        candidates.append((float(np.clip((c + t*v - a) @ u / uu, 0, 1)), t))
    denominator = uu*vv - uv*uv
    if denominator > 1e-12*uu*vv:
        uw, vw = float(u @ w), float(v @ w)
        s, t = (uv*vw - vv*uw)/denominator, (uu*vw - uv*uw)/denominator
        if 0 <= s <= 1 and 0 <= t <= 1:
            candidates.append((s, t))
    return min((float(np.linalg.norm(a + s*u - c - t*v)), s, t) for s, t in candidates)


def screen_segments(segments: tuple[BackboneSegment, ...], generated: frozenset[AtomIdentity], *,
                    bonded_atoms: frozenset[frozenset[AtomIdentity]] = frozenset(),
                    threshold_angstrom: float = 0.75, intersection_tolerance_angstrom: float = 1e-5,
                    max_candidate_pairs: int = 200_000, max_findings: int = 10_000) -> SegmentScreen:
    """Spatially indexed diagnostic screen; budgets fail rather than report clean.

    Engineering distances are NOT calibrated protein-acceptance thresholds.
    A supplied disulfide/other bond participates in bonded-neighbor exclusion,
    but absence of segment alarms proves neither unthreaded topology nor a
    chain-passage-free minimization trajectory.
    """
    if not 0 < intersection_tolerance_angstrom <= threshold_angstrom or not np.isfinite(threshold_angstrom):
        raise ValueError("invalid diagnostic distances")
    if max_candidate_pairs < 1 or max_findings < 1:
        raise ValueError("diagnostic budgets must be positive")
    if not segments:
        return SegmentScreen((), 0, threshold_angstrom, intersection_tolerance_angstrom)
    edges = set(bonded_atoms) | {frozenset((s.first, s.second)) for s in segments}
    endpoints = [{s.first, s.second} for s in segments]
    coordinates = {}
    for segment in segments:
        for identity, xyz in ((segment.first, segment.start), (segment.second, segment.stop)):
            if identity in coordinates and not np.array_equal(coordinates[identity], xyz):
                raise ValueError("same atom has inconsistent segment coordinates")
            coordinates[identity] = xyz
    midpoints = np.array([(np.asarray(s.start) + s.stop)/2 for s in segments])
    radii = np.array([np.linalg.norm(np.asarray(s.stop) - s.start)/2 for s in segments])
    tree = cKDTree(midpoints)
    pairs = 0
    alarms = []
    for i, segment in enumerate(segments):
        for j in sorted(tree.query_ball_point(midpoints[i], radii[i] + radii.max() + threshold_angstrom)):
            if j <= i:
                continue
            pairs += 1
            if pairs > max_candidate_pairs:
                raise ValueError("segment screen candidate budget exceeded; no complete result")
            if endpoints[i] & endpoints[j] or any(frozenset((a, b)) in edges for a in endpoints[i] for b in endpoints[j]):
                continue
            distance, s, t = segment_distance(segment, segments[j])
            if distance > threshold_angstrom:
                continue
            mobile = bool(endpoints[i] & generated), bool(endpoints[j] & generated)
            ownership = "generated-generated" if all(mobile) else "generated-observed" if any(mobile) else "source-only"
            interior = 1e-6 < s < 1-1e-6 and 1e-6 < t < 1-1e-6
            category = "3d-interior-intersection" if distance <= intersection_tolerance_angstrom and interior else "3d-near-contact"
            alarms.append(SegmentAlarm(i, j, distance, s, t, ownership, category))
            if len(alarms) > max_findings:
                raise ValueError("segment screen finding budget exceeded; no complete result")
    return SegmentScreen(tuple(alarms), pairs, threshold_angstrom, intersection_tolerance_angstrom)
