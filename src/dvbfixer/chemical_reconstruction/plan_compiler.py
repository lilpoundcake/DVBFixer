"""Deterministic chain-only plan compiler; no coordinates generated or scored."""
from __future__ import annotations

from dvbfixer.domain.structure_identity import ExactAtomRef

from .admission import admit
from .authority_geometry import geometry_view
from .contracts import ContractError, ReconstructionRequest
from .measurement import torsion_degrees
from .proposals import AuditPolicy, GeometryPlan, Placement, Strategy


class PlanCompilationError(ContractError):
    def __init__(self, code: str, blocked: tuple[ExactAtomRef, ...] = ()) -> None:
        self.code, self.blocked = code, blocked
        super().__init__(code)


def compile_chain_plan(request: ReconstructionRequest, policy: AuditPolicy, *,
                       initial_placements: tuple[Placement, ...] = ()) -> GeometryPlan:
    planning = admit(request)
    if any(f.code != "no-declared-capability-matches" for f in planning.findings):
        raise PlanCompilationError("request-not-admissible")
    if request.ownership.temporary or request.ownership.mobile:
        raise PlanCompilationError("unsupported-construction-ownership")
    if any(sum(o.identity.component == target for o in request.observations) < 3 for target in request.targets):
        raise PlanCompilationError("insufficient-observed-anchors")
    ideals, bonds = geometry_view(request)
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in ideals}
    for pair in bonds:
        a, b = sorted(pair)
        adjacent[a].add(b)
        adjacent[b].add(a)
    available = {o.identity for o in request.observations}
    pending = set(request.ownership.generated)
    GeometryPlan(request.digest, Strategy.INTERNAL, policy, initial_placements)  # Immutable shape check only.
    placements = list(initial_placements)
    for placement in initial_placements:
        if placement.atom not in pending or not set(placement.anchors) <= available:
            raise PlanCompilationError("initial-placement-ownership-or-dependency")
        pending.remove(placement.atom)
        available.add(placement.atom)
    while pending:
        chosen: Placement | None = None
        for d in sorted(pending):
            for c in sorted(adjacent[d] & available):
                for b in sorted((adjacent[c] & available)-{d}):
                    for a in sorted((adjacent[b] & available)-{c, d}):
                        if any(x.component != d.component for x in (a, b, c)):
                            continue
                        try:
                            torsion_degrees((ideals[a], ideals[b], ideals[c], ideals[d]))
                        except ValueError:
                            continue
                        chosen = Placement(d, (a, b, c))
                        break
                    if chosen is not None:
                        break
                if chosen is not None:
                    break
            if chosen is not None:
                break
        if chosen is None:
            raise PlanCompilationError("no-supported-chain-frame", tuple(sorted(pending)))
        placements.append(chosen)
        pending.remove(chosen.atom)
        available.add(chosen.atom)
    return GeometryPlan(request.digest, Strategy.INTERNAL, policy, tuple(placements))
