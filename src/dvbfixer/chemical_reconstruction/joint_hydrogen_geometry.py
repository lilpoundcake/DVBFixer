"""Joint proper orientation of generated H groups; no observed-coordinate fit output."""
from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations

from dvbfixer.domain.structure_identity import ExactAtomRef

from .authority_geometry import geometry_view
from .contracts import ContractError, ReconstructionRequest, _identity, _point, digest_record
from .measurement import Point, _cross, _dot, _sub, _unit
from .metric_geometry import GeneratedMetricGeometryPlan
from .proposals import GeometryPlan, PlacementKind


@dataclass(frozen=True)
class DirectionGroup:
    parent: ExactAtomRef
    heavy_neighbours: tuple[ExactAtomRef, ...]
    hydrogens: tuple[ExactAtomRef, ...]

    def __post_init__(self) -> None:
        _identity(self.parent)
        if (not isinstance(self.heavy_neighbours, tuple) or not isinstance(self.hydrogens, tuple)
                or (len(self.heavy_neighbours), len(self.hydrogens)) not in {(2, 2), (3, 1)}):
            raise ContractError("joint sp3 group needs exact 2-heavy/2-H or 3-heavy/1-H inventory")
        atoms = (self.parent, *self.heavy_neighbours, *self.hydrogens)
        for a in atoms:
            _identity(a)
        if len(set(atoms)) != len(atoms) or any(a.component != self.parent.component for a in atoms):
            raise ContractError("joint group needs distinct same-component identities")
        if self.heavy_neighbours != tuple(sorted(self.heavy_neighbours)) or self.hydrogens != tuple(sorted(self.hydrogens)):
            raise ContractError("joint group identities must use exact canonical lexical order")


@dataclass(frozen=True, kw_only=True)
class JointMetricGeometryPlan(GeneratedMetricGeometryPlan):
    joint_groups: tuple[DirectionGroup, ...]

    def __post_init__(self) -> None:
        super().__post_init__()
        if (not isinstance(self.joint_groups, tuple) or not self.joint_groups or len(self.joint_groups) > 4096
                or any(not isinstance(g, DirectionGroup) for g in self.joint_groups)):
            raise ContractError("bounded explicit joint groups required")
        if len({g.parent for g in self.joint_groups}) != len(self.joint_groups):
            raise ContractError("duplicate joint-group parent")
        mapping = {h: g for g in self.joint_groups for h in g.hydrogens}
        if len(mapping) != sum(len(g.hydrogens) for g in self.joint_groups):
            raise ContractError("overlapping joint H groups")
        placements = [p for p in self.placements if p.kind == PlacementKind.JOINT]
        if len(placements) != len(mapping) or {p.atom for p in placements} != set(mapping):
            raise ContractError("joint placements must cover group H exactly once")
        for p in self.placements:
            if p.atom in mapping:
                g = mapping[p.atom]
                if p.kind != PlacementKind.JOINT or p.anchors != (g.parent, *g.heavy_neighbours[:2]):
                    raise ContractError("joint placement/group binding mismatch")

    @property
    def digest(self) -> str:
        return digest_record(("chemical-joint-H-geometry-plan-v1", self))


@dataclass(frozen=True)
class DirectionFit:
    rotation: tuple[Point, Point, Point]
    translation_norm: float
    directional_chord_rms: float

    def apply(self, vector: Point) -> Point:
        return _dot(self.rotation[0], vector), _dot(self.rotation[1], vector), _dot(self.rotation[2], vector)


def fit_joint_directions(source: tuple[Point, ...], target: tuple[Point, ...]) -> DirectionFit:
    """Origin-constrained SO(3) via the existing proper-fit owner on ± clouds."""
    import numpy as np
    from numpy.typing import NDArray

    from dvbfixer.model.diffusion.geometry import weighted_kabsch

    if not isinstance(source, tuple) or not isinstance(target, tuple) or len(source) != len(target) or len(source) not in {2, 3}:
        raise ValueError("two or three corresponding direction vectors required")
    for v in (*source, *target):
        _point(v)
    src, dst = tuple(_unit(v) for v in source), tuple(_unit(v) for v in target)
    for directions in (src, dst):
        if max(math.sqrt(_dot(_cross(a, b), _cross(a, b))) for a, b in combinations(directions, 2)) <= 1e-8:
            raise ValueError("collinear joint heavy directions")
    def cloud(vectors: tuple[Point, ...]) -> NDArray[np.float64]:
        return np.asarray(tuple(row for v in vectors for row in (v, (-v[0], -v[1], -v[2]))), dtype=np.float64)

    transform = weighted_kabsch(cloud(src), cloud(dst))
    rotation = transform.rotation
    translation = float(np.linalg.norm(transform.translation))
    if (not np.isfinite(rotation).all() or not np.allclose(rotation.T @ rotation, np.eye(3), rtol=0, atol=1e-10)
            or abs(float(np.linalg.det(rotation))-1.) > 1e-10 or not math.isfinite(translation) or translation > 1e-10):
        raise ValueError("joint direction fit is not proper origin-constrained SO(3)")
    rows = tuple(tuple(float(x) for x in row) for row in rotation)
    r: tuple[Point, Point, Point] = ((rows[0][0], rows[0][1], rows[0][2]),
                                    (rows[1][0], rows[1][1], rows[1][2]),
                                    (rows[2][0], rows[2][1], rows[2][2]))
    fit = DirectionFit(r, translation, 0.)
    residual = math.sqrt(sum(_dot(_sub(fit.apply(a), b), _sub(fit.apply(a), b)) for a, b in zip(src, dst))/len(src))
    return DirectionFit(r, translation, residual)


def joint_group_findings(request: ReconstructionRequest, plan: GeometryPlan) -> set[str]:
    if not isinstance(plan, JointMetricGeometryPlan):
        return {"joint-group-plan-required"} if any(p.kind == PlacementKind.JOINT for p in plan.placements) else set()
    errors = set()
    atoms = {ExactAtomRef(c.identity, a.name): a for c in request.components if c.resolved is not None for a in c.resolved.graph.atoms}
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in atoms}
    orders = {}
    for c in request.components:
        if c.resolved is None:
            continue
        for bond in c.resolved.graph.bonds:
            a, b = ExactAtomRef(c.identity, bond.atom1), ExactAtomRef(c.identity, bond.atom2)
            adjacent[a].add(b)
            adjacent[b].add(a)
            orders[frozenset((a, b))] = bond.order
    observed = {o.identity: o for o in request.observations}
    for g in plan.joint_groups:
        refs = {g.parent, *g.heavy_neighbours, *g.hydrogens}
        if not refs <= atoms.keys():
            errors.add("joint-group-unknown-authority-atom")
            continue
        p = atoms[g.parent]
        if p.element != "C" or p.formal_charge != 0 or p.aromatic:
            errors.add("joint-group-parent-chemistry-unqualified")
        if any(atoms[h].element != "H" for h in g.hydrogens) or any(atoms[h].element in {"H", "D"} for h in g.heavy_neighbours):
            errors.add("joint-group-element-mismatch")
        if (adjacent[g.parent] != {*g.heavy_neighbours, *g.hydrogens}
                or any(adjacent[h] != {g.parent} for h in g.hydrogens)
                or any(orders.get(frozenset((g.parent, n))) != 1 for n in adjacent[g.parent])):
            errors.add("joint-group-needs-complete-single-bond-sp3-inventory")
        if not set(g.hydrogens) <= request.ownership.generated:
            errors.add("joint-group-H-must-be-generated")
        active = request.ownership.generated | request.ownership.temporary | request.ownership.mobile
        if any(a not in observed or observed[a].occupancy != 1 or a in active for a in (g.parent, *g.heavy_neighbours)):
            errors.add("joint-group-needs-fixed-unit-occupancy-heavy-observations")
    return errors


def joint_group_points(request: ReconstructionRequest, group: DirectionGroup) -> dict[ExactAtomRef, Point]:
    ideals, bonds = geometry_view(request)
    observed = {o.identity: o.angstrom for o in request.observations}
    parent = group.parent
    fit = fit_joint_directions(tuple(_sub(ideals[n], ideals[parent]) for n in group.heavy_neighbours),
                               tuple(_sub(observed[n], observed[parent]) for n in group.heavy_neighbours))
    result = {}
    for h in group.hydrogens:
        direction = fit.apply(_unit(_sub(ideals[h], ideals[parent])))
        length = bonds[frozenset((parent, h))]
        result[h] = (observed[parent][0]+length*direction[0], observed[parent][1]+length*direction[1],
                     observed[parent][2]+length*direction[2])
    return result
