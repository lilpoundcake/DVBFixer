"""Explicit sourced scalar geometry, without changing template-plan semantics."""
from __future__ import annotations

from dataclasses import dataclass

from dvbfixer.domain.structure_identity import ExactAtomRef

from .authority_geometry import geometry_view
from .contracts import (
    ConstraintKind,
    ContractError,
    LinkKind,
    ReconstructionRequest,
    Unit,
    _identity,
    _number,
    _sha,
    digest_record,
)
from .proposals import GeometryPlan, PlacementKind, Strategy


def metric_key(atoms: tuple[ExactAtomRef, ...]) -> tuple[ExactAtomRef, ...]:
    return min(atoms, tuple(reversed(atoms)))


@dataclass(frozen=True)
class GeometryMetric:
    atoms: tuple[ExactAtomRef, ...]
    value: float
    unit: Unit
    source: str
    source_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.atoms, tuple) or len(self.atoms) not in {2, 3} or len(set(self.atoms)) != len(self.atoms):
            raise ContractError("metric needs distinct bond or ordered-angle identities")
        for atom in self.atoms:
            _identity(atom)
        _sha(self.source_sha256)
        if not isinstance(self.source, str) or not self.source.strip() or not _number(self.value) or self.value <= 0:
            raise ContractError("explicit positive metric and provenance required")
        if not isinstance(self.unit, Unit) or self.unit != (Unit.ANGSTROM if len(self.atoms) == 2 else Unit.DEGREE):
            raise ContractError("metric must use explicit angstrom or degree units")
        if len(self.atoms) == 3 and self.value >= 180:
            raise ContractError("metric angle must be below 180 degrees")


@dataclass(frozen=True)
class MetricGeometryPlan(GeometryPlan):
    metrics: tuple[GeometryMetric, ...] = ()

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.strategy != Strategy.INTERNAL:
            raise ContractError("metric-plan-strategy-unqualified")
        if (not isinstance(self.metrics, tuple) or not self.metrics or len(self.metrics) > 8192
                or any(not isinstance(m, GeometryMetric) for m in self.metrics)):
            raise ContractError("nonempty bounded immutable metric authority required")
        if len({metric_key(m.atoms) for m in self.metrics}) != len(self.metrics):
            raise ContractError("duplicate or reversed metric authority")

    @property
    def digest(self) -> str:
        return digest_record(("chemical-metric-geometry-plan-v1", self))


def explicit_metrics(plan: GeometryPlan) -> dict[tuple[ExactAtomRef, ...], float]:
    return {metric_key(m.atoms): m.value for m in plan.metrics} if isinstance(plan, MetricGeometryPlan) else {}


def metric_findings(request: ReconstructionRequest, plan: GeometryPlan) -> set[str]:
    if not isinstance(plan, MetricGeometryPlan):
        return {"planar-needs-explicit-metric-plan"} if any(p.kind == PlacementKind.PLANAR for p in plan.placements) else set()
    errors = set()
    ideals, bonds = geometry_view(request)
    available = set(ideals) | {o.identity for o in request.observations}
    connected = set(bonds) | {frozenset((c.link.component_atom, c.link.other_atom))
                              for c in request.connections if c.kind == LinkKind.COVALENT}
    metrics = explicit_metrics(plan)
    observed = {o.identity for o in request.observations}
    planar_atoms = {p.atom for p in plan.placements if p.kind == PlacementKind.PLANAR}
    for metric in plan.metrics:
        if not set(metric.atoms) & request.ownership.generated:
            errors.add("metric-cannot-override-only-observed-geometry")
        if not set(metric.atoms) <= available:
            errors.add("metric-unknown-atom")
        if not set(metric.atoms) & request.ownership.generated <= planar_atoms:
            errors.add("metric-generated-atom-strategy-unqualified")
        pairs = (metric.atoms,) if len(metric.atoms) == 2 else (metric.atoms[:2], metric.atoms[1:])
        if any(frozenset(pair) not in connected for pair in pairs):
            errors.add("metric-needs-authoritative-connectivity")
    for placement in plan.placements:
        if placement.kind != PlacementKind.PLANAR:
            continue
        a, b, c = placement.anchors
        d = placement.atom
        if not set(placement.anchors) <= observed:
            errors.add("planar-anchor-observations-required")
        if any(metric_key(atoms) not in metrics for atoms in ((a, d), (b, a, d), (c, a, d))):
            errors.add("planar-needs-explicit-length-and-both-angles")
        if not any(q.kind == ConstraintKind.PLANARITY and q.atoms == (*placement.anchors, d) for q in request.constraints):
            errors.add("planar-needs-explicit-plane-constraint")
    return errors
