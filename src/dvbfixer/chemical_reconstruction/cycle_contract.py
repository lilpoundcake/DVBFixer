"""Pure cycle-aware problem compilation, never an initializer or geometry parent."""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from itertools import combinations
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

from .admission import admit
from .authority_geometry import geometry_view
from .closure_diagnostics import is_cycle_edge
from .contracts import (
    Constraint,
    ConstraintKind,
    ContractError,
    Operation,
    Origin,
    PreparedContextRequest,
    ReconstructionRequest,
    Unit,
    _identity,
    _number,
    _sha,
    _text,
    _tuple,
    canonical_graph_digest,
    digest_record,
)
from .measurement import _cross, _sub, _unit, angle_degrees
from .projection_contracts import POLICY
from .proposals import AuditPolicy


class CycleContractError(ContractError):
    """A planning refusal, not a request for fallback coordinates."""


@dataclass(frozen=True)
class CycleMetric:
    kind: Literal["bond", "angle"]
    atoms: tuple[ExactAtomRef, ...]
    target: float
    tolerance: float
    role: Literal["fixed-only", "generated-involving"]

    def __post_init__(self) -> None:
        _tuple(self.atoms, ExactAtomRef, "cycle metric atoms", maximum=3)
        if (self.kind not in {"bond", "angle"} or len(self.atoms) != (2 if self.kind == "bond" else 3)
                or len(set(self.atoms)) != len(self.atoms) or not _number(self.target)
                or not _number(self.tolerance) or self.tolerance <= 0
                or self.role not in {"fixed-only", "generated-involving"}):
            raise CycleContractError("cycle-metric-shape")


@dataclass(frozen=True)
class GeneratedBlock:
    atoms: tuple[ExactAtomRef, ...]
    boundary_edges: tuple[tuple[ExactAtomRef, ExactAtomRef], ...]
    cyclic_atoms: tuple[ExactAtomRef, ...]

    def __post_init__(self) -> None:
        _tuple(self.atoms, ExactAtomRef, "generated cycle block", maximum=16)
        _tuple(self.boundary_edges, tuple, "observed boundary edges", maximum=128)
        _tuple(self.cyclic_atoms, ExactAtomRef, "cyclic generated atoms", maximum=16)
        if (not self.atoms or self.atoms != tuple(sorted(set(self.atoms)))
                or self.boundary_edges != tuple(sorted(set(self.boundary_edges)))
                or self.cyclic_atoms != tuple(sorted(set(self.cyclic_atoms)))
                or not set(self.cyclic_atoms) <= set(self.atoms)
                or any(len(edge) != 2 or edge[0] not in self.atoms or edge[1] in self.atoms for edge in self.boundary_edges)):
            raise CycleContractError("cycle-generated-block-shape")


@dataclass(frozen=True)
class CycleCoverage:
    request_digest: str
    observations_digest: str
    graph_digest: str
    variables: tuple[ExactAtomRef, ...]
    heavy_nodes: tuple[ExactAtomRef, ...]
    heavy_edges: tuple[tuple[ExactAtomRef, ExactAtomRef], ...]
    cycle_edges: tuple[tuple[ExactAtomRef, ExactAtomRef], ...]
    generated_cycle_edges: tuple[tuple[ExactAtomRef, ExactAtomRef], ...]
    metrics: tuple[CycleMetric, ...]
    constraint_roles: tuple[Literal["fixed-only", "generated-involving"], ...]
    constraint_scalar_count: int
    blocks: tuple[GeneratedBlock, ...]
    unassessed_context: tuple[ComponentInstanceRef, ...]

    def __post_init__(self) -> None:
        for pin in (self.request_digest, self.observations_digest, self.graph_digest):
            _sha(pin)
        _tuple(self.variables, ExactAtomRef, "cycle variables", maximum=16)
        _tuple(self.heavy_nodes, ExactAtomRef, "target heavy graph nodes", maximum=64)
        for edges in (self.heavy_edges, self.cycle_edges, self.generated_cycle_edges):
            _tuple(edges, tuple, "cycle graph edges", maximum=128)
            if edges != tuple(sorted(set(edges))) or any(len(e) != 2 or e[0] >= e[1] for e in edges):
                raise CycleContractError("cycle-normalized-edge-shape")
        _tuple(self.metrics, CycleMetric, "complete heavy metrics", maximum=512)
        _tuple(self.blocks, GeneratedBlock, "generated cycle blocks", maximum=16)
        _tuple(self.constraint_roles, str, "cycle constraint roles", maximum=512)
        if (type(self.constraint_scalar_count) is not int or self.constraint_scalar_count < 0
                or len(self.metrics)+self.constraint_scalar_count > 512):
            raise CycleContractError("cycle-metric-budget")


def constraint_request(source: ReconstructionRequest, extra: tuple[Constraint, ...]) -> ReconstructionRequest:
    """Distinct pure constraint provenance, NOT a projection or prepared output."""
    if type(source) not in {ReconstructionRequest, PreparedContextRequest}:
        raise CycleContractError("cycle-needs-original-source-request")
    _tuple(extra, Constraint, "additional cycle constraints", maximum=512)
    settings = source.origin.settings if source.origin else ()
    if any(key.startswith(("cycle_construction_", "metric_projection_")) for key, _ in settings):
        raise CycleContractError("cycle-nested-or-projection-source-refused")
    settings = (*settings, ("cycle_construction_source", source.digest), ("cycle_construction_constraints", digest_record(extra)))
    origin = (replace(source.origin, settings=settings) if source.origin else
              Origin("cycle-construction-request-v1", source.digest, settings))
    return replace(source, constraints=(*source.constraints, *extra), origin=origin)


def _coverage(source: ReconstructionRequest, extra: tuple[Constraint, ...], anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef],
                 policy: AuditPolicy) -> CycleCoverage:
    request = constraint_request(source, extra)
    if (digest_record(policy) != digest_record(POLICY) or source.operation != Operation.HEAVY or len(source.targets) != 1
            or source.ownership.mobile or source.ownership.temporary or not 1 <= len(source.ownership.generated) <= 16):
        raise CycleContractError("cycle-fixed-policy-single-target-heavy-envelope")
    admission = admit(request)
    findings = sorted(f.code for f in admission.findings if f.code != "no-declared-capability-matches")
    if findings:
        raise CycleContractError("cycle-source-admission:"+findings[0])
    target = next(iter(source.targets))
    component = next(c for c in source.components if c.identity == target)
    state = component.resolved
    if state is None:
        raise CycleContractError("cycle-resolved-target-required")
    elements = {ExactAtomRef(target, a.name): a.element for a in state.graph.atoms}
    heavy = {a for a, e in elements.items() if e.upper() not in {"H", "D"}}
    if len(heavy) > 64:
        raise CycleContractError("cycle-target-node-budget")
    observed = {o.identity: o for o in source.observations}
    selected = {a for a in observed if a.component == target}
    generated = source.ownership.generated
    if (selected | generated != heavy or selected & generated or not generated <= heavy
            or any(observed[a].occupancy != 1 or observed[a].unit != Unit.ANGSTROM for a in selected)):
        raise CycleContractError("cycle-complete-heavy-only-target-inventory")
    if any(generated & {c.link.component_atom, c.link.other_atom} for c in source.connections):
        raise CycleContractError("cycle-generated-external-link-unqualified")
    _tuple(anchors, ExactAtomRef, "cycle observed frame", maximum=3)
    if len(anchors) != 3 or len(set(anchors)) != 3 or not set(anchors) <= selected:
        raise CycleContractError("cycle-needs-three-distinct-fixed-target-frame-atoms")
    for atom in anchors:
        _identity(atom)
    first_point, second_point, third_point = (observed[x].angstrom for x in anchors)
    try:
        normal = _unit(_cross(_sub(second_point, first_point), _sub(third_point, first_point)))
        if not all(math.isfinite(v) for v in normal):
            raise ValueError("nonfinite observed frame arithmetic")
    except ValueError as exc:
        raise CycleContractError("cycle-degenerate-observed-frame") from exc
    ideals, bonds = geometry_view(source)
    edges = tuple(sorted((min(pair), max(pair)) for pair in bonds if pair <= heavy))
    if len(edges) > 128:
        raise CycleContractError("cycle-target-edge-budget")
    cycles = tuple(e for e in edges if is_cycle_edge(edges, e))
    active_cycles = tuple(e for e in cycles if set(e) & generated)
    if not active_cycles:
        raise CycleContractError("cycle-no-generated-involving-cycle")
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in heavy}
    for a, b in edges:
        adjacent[a].add(b)
        adjacent[b].add(a)
    metrics: list[CycleMetric] = []

    def add(kind: Literal["bond", "angle"], atoms: tuple[ExactAtomRef, ...], value: float) -> None:
        if len(metrics) >= 512:
            raise CycleContractError("cycle-metric-budget")
        role: Literal["fixed-only", "generated-involving"] = "generated-involving" if set(atoms) & generated else "fixed-only"
        tolerance = policy.bond_tolerance_angstrom if kind == "bond" else policy.angle_tolerance_degrees
        if role == "fixed-only":
            measured = math.dist(*(observed[x].angstrom for x in atoms)) if kind == "bond" else angle_degrees(*(observed[x].angstrom for x in atoms))
            if abs(measured-value) > tolerance:
                raise CycleContractError("cycle-fixed-authority-"+kind+"-conflict")
        metrics.append(CycleMetric(kind, atoms, value, tolerance, role))

    for a, b in edges:
        add("bond", (a, b), bonds[frozenset((a, b))])
    for centre in sorted(heavy):
        for first, last in combinations(sorted(adjacent[centre]), 2):
            add("angle", (first, centre, last), angle_degrees(ideals[first], ideals[centre], ideals[last]))
    roles: list[Literal["fixed-only", "generated-involving"]] = []
    scalar_count = 0
    for q in request.constraints:
        if q.kind == ConstraintKind.TORSION:
            raise CycleContractError("cycle-periodic-torsion-residual-not-declared")
        if q.kind in {ConstraintKind.DISTANCE, ConstraintKind.ANGLE, ConstraintKind.RING_CLOSURE}:
            if q.lower is None or q.upper is None or q.upper <= q.lower:
                raise CycleContractError("cycle-needs-positive-width-bounded-metric")
        if q.kind in {ConstraintKind.EXCLUSION, ConstraintKind.PLANARITY}:
            bound = q.lower if q.kind == ConstraintKind.EXCLUSION else q.upper
            if bound is None or bound <= 0:
                raise CycleContractError("cycle-needs-positive-exclusion-or-plane-bound")
        roles.append("generated-involving" if set(q.atoms) & generated else "fixed-only")
        scalar_count += len(q.atoms)-3 if q.kind == ConstraintKind.PLANARITY else 1
        if len(metrics)+scalar_count > 512:
            raise CycleContractError("cycle-metric-budget")
    cyclic_atoms = {a for e in cycles for a in e}
    pending = set(generated)
    blocks = []
    while pending:
        reached, stack = set(), [min(pending)]
        while stack:
            atom = stack.pop()
            if atom in reached:
                continue
            reached.add(atom)
            stack.extend(sorted((adjacent[atom] & generated)-reached, reverse=True))
        pending -= reached
        boundary = tuple(sorted((a, b) for a in reached for b in adjacent[a] & selected))
        if not boundary:
            raise CycleContractError("cycle-unanchored-generated-fragment")
        blocks.append(GeneratedBlock(tuple(sorted(reached)), boundary, tuple(sorted(reached & cyclic_atoms))))
    return CycleCoverage(request.digest, digest_record(source.observations), canonical_graph_digest(state.graph),
        tuple(sorted(generated)), tuple(sorted(heavy)), edges, cycles, active_cycles, tuple(metrics), tuple(roles),
        scalar_count, tuple(blocks), admission.unassessed_context)


@dataclass(frozen=True)
class CycleConstructionProblem:
    source: ReconstructionRequest
    extra_constraints: tuple[Constraint, ...]
    anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef]
    coverage: CycleCoverage
    policy: AuditPolicy = POLICY
    revision: Literal["cycle-construction-problem-v1-planning-only"] = "cycle-construction-problem-v1-planning-only"

    def __post_init__(self) -> None:
        expected = _coverage(self.source, self.extra_constraints, self.anchors, self.policy)
        if (self.revision != "cycle-construction-problem-v1-planning-only"
                or digest_record(self.coverage) != digest_record(expected)):
            raise CycleContractError("cycle-problem-exact-source-graph-metric-coverage-binding")

    @property
    def digest(self) -> str:
        return digest_record(("cycle-construction-problem-v1", self))

    @property
    def coordinates_produced(self) -> Literal[False]:
        return False

    @property
    def geometry_approved(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False

    @property
    def solver_defined(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False


def compile_cycle_problem(source: ReconstructionRequest, anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef], *,
                             extra_constraints: tuple[Constraint, ...] = ()) -> CycleConstructionProblem:
    coverage = _coverage(source, extra_constraints, anchors, POLICY)
    return CycleConstructionProblem(source, extra_constraints, anchors, coverage)


@dataclass(frozen=True)
class CyclePlanningAdmission:
    source_digest: str
    status: Literal["planned-contract-only", "refused"]
    reason: str
    problem: CycleConstructionProblem | None = None

    def __post_init__(self) -> None:
        _sha(self.source_digest)
        _text(self.reason, "cycle planning reason")
        if self.status == "planned-contract-only":
            if (not isinstance(self.problem, CycleConstructionProblem) or self.problem.source.digest != self.source_digest
                    or self.reason != "source-contract-only-no-initializer-solver-or-geometry"):
                raise CycleContractError("cycle-planning-status-binding")
        elif self.status != "refused" or self.problem is not None or not self.reason:
            raise CycleContractError("cycle-planning-refusal-binding")

    @property
    def coordinates_produced(self) -> Literal[False]:
        return False

    @property
    def geometry_approved(self) -> Literal[False]:
        return False

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False

    @property
    def solver_defined(self) -> Literal[False]:
        return False

    @property
    def md_ready(self) -> Literal[False]:
        return False


def plan_cycle_construction(source: ReconstructionRequest, anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef], *,
                                extra_constraints: tuple[Constraint, ...] = ()) -> CyclePlanningAdmission:
    if not isinstance(source, ReconstructionRequest):
        raise CycleContractError("cycle-typed-original-source-required")
    try:
        problem = compile_cycle_problem(source, anchors, extra_constraints=extra_constraints)
    except ValueError as exc:
        return CyclePlanningAdmission(source.digest, "refused", str(exc))
    return CyclePlanningAdmission(source.digest, "planned-contract-only", "source-contract-only-no-initializer-solver-or-geometry", problem)
