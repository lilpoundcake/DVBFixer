"""Explicit sourced scalar geometry, without changing template-plan semantics."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from dvbfixer.domain.structure_identity import ExactAtomRef

from .authority_geometry import geometry_view
from .contracts import (
    ConstraintKind,
    ContractError,
    LinkKind,
    Operation,
    ReconstructionRequest,
    Unit,
    _identity,
    _number,
    _sha,
    canonical_graph_digest,
    digest_record,
)
from .measurement import angle_degrees
from .proposals import AuditPolicy, GeometryPlan, PlacementKind, Strategy


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


@dataclass(frozen=True, kw_only=True)
class GeneratedMetricGeometryPlan(MetricGeometryPlan):
    generated_policy: AuditPolicy

    def __post_init__(self) -> None:
        super().__post_init__()
        p = self.generated_policy
        if (not isinstance(p, AuditPolicy) or p.minimum_separation_angstrom != self.policy.minimum_separation_angstrom
                or p.max_anchor_rmsd_angstrom != self.policy.max_anchor_rmsd_angstrom
                or p.bond_tolerance_angstrom > self.policy.bond_tolerance_angstrom
                or p.angle_tolerance_degrees > self.policy.angle_tolerance_degrees):
            raise ContractError("generated bond/angle bounds must be explicit, no looser, with shared clash/anchor policy")

    @property
    def digest(self) -> str:
        return digest_record(("chemical-generated-metric-plan-v1", self))


@dataclass(frozen=True, kw_only=True)
class SourceMetricGeometryPlan(GeneratedMetricGeometryPlan):
    authority_graph_digest: str

    def __post_init__(self) -> None:
        super().__post_init__()
        _sha(self.authority_graph_digest)

    @property
    def digest(self) -> str:
        return digest_record(("chemical-source-only-metric-plan-v1", self))


def _source_metric_findings(request: ReconstructionRequest, plan: SourceMetricGeometryPlan) -> set[str]:
    targets = [c for c in request.components if c.identity in request.targets]
    if len(targets) != 1 or targets[0].resolved is None or request.operation != Operation.HYDROGENS:
        return {"source-metrics-need-one-resolved-H-target"}
    target = targets[0]
    assert target.resolved is not None
    graph = target.resolved.graph
    pin = canonical_graph_digest(graph)
    errors = set()
    if pin != plan.authority_graph_digest:
        errors.add("source-metric-graph-binding")
    generated = frozenset(ExactAtomRef(target.identity, a.name) for a in graph.atoms if a.element == "H")
    if (not generated or generated != request.ownership.generated or any(a.element == "D" for a in graph.atoms)):
        errors.add("source-metrics-need-exact-complete-protium-mask")
    ideals, _bonds = geometry_view(request)
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {ExactAtomRef(target.identity, a.name): set() for a in graph.atoms}
    expected = {}
    for bond in graph.bonds:
        a, b = ExactAtomRef(target.identity, bond.atom1), ExactAtomRef(target.identity, bond.atom2)
        adjacent[a].add(b)
        adjacent[b].add(a)
        if {a, b} & generated:
            expected[metric_key((a, b))] = bond.ideal_length
    try:
        for centre, neighbours in adjacent.items():
            for a, b in combinations(sorted(neighbours), 2):
                if {a, centre, b} & generated:
                    expected[metric_key((a, centre, b))] = angle_degrees(ideals[a], ideals[centre], ideals[b])
    except (KeyError, ValueError):
        return errors | {"source-metric-ideal-authority-unusable"}
    if set(expected) != {metric_key(m.atoms) for m in plan.metrics}:
        errors.add("source-metric-complete-coverage")
    for metric in plan.metrics:
        if metric.source_sha256 != pin:
            errors.add("source-metric-provenance-binding")
        if expected.get(metric_key(metric.atoms)) != metric.value:
            errors.add("source-metric-cannot-override-authority")
    return errors


def audit_policy_for(plan: GeometryPlan, atoms: tuple[ExactAtomRef, ...], generated: frozenset[ExactAtomRef]) -> AuditPolicy:
    if isinstance(plan, GeneratedMetricGeometryPlan) and set(atoms) & generated:
        return plan.generated_policy
    return plan.policy


def explicit_metrics(plan: GeometryPlan) -> dict[tuple[ExactAtomRef, ...], float]:
    return {metric_key(m.atoms): m.value for m in plan.metrics} if isinstance(plan, MetricGeometryPlan) else {}


def metric_findings(request: ReconstructionRequest, plan: GeometryPlan) -> set[str]:
    if not isinstance(plan, MetricGeometryPlan):
        return {"planar-needs-explicit-metric-plan"} if any(p.kind == PlacementKind.PLANAR for p in plan.placements) else set()
    errors = set()
    if isinstance(plan, SourceMetricGeometryPlan):
        errors.update(_source_metric_findings(request, plan))
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
        if not isinstance(plan, SourceMetricGeometryPlan) and not set(metric.atoms) & request.ownership.generated <= planar_atoms:
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
