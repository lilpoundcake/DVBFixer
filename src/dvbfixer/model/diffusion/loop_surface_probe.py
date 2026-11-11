"""Private straight-chord/fan surface probe, NOT a protein knot classifier.

The observed/generated loop is artificially closed by a straight chord. A
backbone path piercing this chosen surface can be legitimate, depend on closure,
or disappear when an open terminus is moved. Every hit is therefore an
experimental suspected-threading diagnostic, never a proven defect or gate.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from dvbfixer.model.diffusion.backbone_intersections import BackboneSegment


@dataclass(frozen=True)
class SurfacePiercing:
    segment_index: int
    fraction: float
    point: tuple[float, float, float]
    orientation: int
    surface_edge_ambiguous: bool
    topological_defect_established: bool = False


@dataclass(frozen=True)
class LoopSurfaceProbe:
    piercings: tuple[SurfacePiercing, ...]
    degenerate_triangles: int
    parallel_or_coplanar_pairs: int
    closure: str = "straight chord from last vertex to first; centroid triangle fan"
    method: str = "experimental-straight-chord-fan-probe-v1"
    closure_ambiguous: bool = True
    open_path_end_sensitive: bool = True
    disulfide_topology_tested: bool = False
    publication_gate: bool = False


def probe_loop_surface(vertices, paths: tuple[BackboneSegment, ...], *, tolerance_angstrom: float = 1e-5,
                       max_tests: int = 200_000) -> LoopSurfaceProbe:
    """Intersect open paths with one explicit artificial spanning surface.

    Internal fan-edge duplicates are merged per path segment. Parallel/tangent
    pairs remain visibly unresolved; absence of hits is never an unthreaded pass.
    The caller must omit the loop's own residues and preserve source ownership.
    """
    polygon = np.asarray(vertices, dtype=float)
    if polygon.ndim != 2 or polygon.shape[1] != 3 or len(polygon) < 3 or not np.isfinite(polygon).all():
        raise ValueError("surface requires at least three finite 3D loop vertices")
    if not np.isfinite(tolerance_angstrom) or tolerance_angstrom <= 0:
        raise ValueError("invalid surface tolerance")
    if len(paths)*len(polygon) > max_tests:
        raise ValueError("surface probe budget exceeded; no complete result")
    if any(np.linalg.norm(a-b) < tolerance_angstrom for a, b in zip(polygon, np.roll(polygon, -1, axis=0))):
        raise ValueError("degenerate loop boundary")
    center = polygon.mean(axis=0)
    hits = []
    degenerate = 0
    parallel = 0
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        edge1, edge2 = a-center, b-center
        normal = np.cross(edge1, edge2)
        area = np.linalg.norm(normal)
        if area < tolerance_angstrom**2:
            degenerate += 1
            continue
        normal /= area
        for index, path in enumerate(paths):
            start, stop = np.asarray(path.start), np.asarray(path.stop)
            direction = stop-start
            denominator = float(direction @ normal)
            if abs(denominator) <= tolerance_angstrom:
                parallel += 1
                continue
            fraction = float((center-start) @ normal / denominator)
            if not 0 < fraction < 1:
                continue
            point = start+fraction*direction
            # Barycentric triangle membership, with dimensionless edge tolerance.
            dot11, dot12, dot22 = edge1@edge1, edge1@edge2, edge2@edge2
            delta = point-center
            denominator_bary = dot11*dot22-dot12*dot12
            u = float((dot22*(delta@edge1)-dot12*(delta@edge2))/denominator_bary)
            v = float((dot11*(delta@edge2)-dot12*(delta@edge1))/denominator_bary)
            epsilon = tolerance_angstrom/max(np.linalg.norm(edge1), np.linalg.norm(edge2))
            if u < -epsilon or v < -epsilon or u+v > 1+epsilon:
                continue
            ambiguous = bool(min(abs(u), abs(v), abs(1-u-v)) <= epsilon)
            hit = SurfacePiercing(index, fraction, tuple(float(x) for x in point), 1 if denominator > 0 else -1, ambiguous)
            duplicate = next((i for i, old in enumerate(hits) if old.segment_index == index and
                              np.linalg.norm(np.asarray(old.point)-point) <= tolerance_angstrom), None)
            if duplicate is None:
                hits.append(hit)
    return LoopSurfaceProbe(tuple(sorted(hits, key=lambda h: (h.segment_index, h.fraction))), degenerate, parallel)
