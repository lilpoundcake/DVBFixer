"""Full retained-context VDW overlap heuristic, not a force field or acceptance gate."""
from __future__ import annotations

import math
from dataclasses import dataclass

from dvbfixer.diagnose.geometry import VDW_RADII_ANGSTROM
from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import LinkKind, ReconstructionRequest, Unit, digest_record
from .proposals import Candidate, GeometryPlan
from .torsion_geometry import ConformerError

MAX_PAIR_EVALUATIONS = 1_000_000


def radii_digest() -> str:
    return digest_record(tuple(sorted(VDW_RADII_ANGSTROM.items())))


@dataclass(frozen=True)
class StericContext:
    request_digest: str
    radii_digest: str
    pairs: tuple[tuple[ExactAtomRef, ExactAtomRef, float], ...]


@dataclass(frozen=True)
class StericScore:
    sum_squared_overlap_angstrom2: float
    worst_overlap_angstrom: float
    pair_count: int
    radii_digest: str
    scope: str = "vdw-overlap-only-not-energy-or-acceptance-gate"


def compile_steric_context(request: ReconstructionRequest) -> StericContext:
    elements = {o.identity: o.element for o in request.observations}
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {}
    for c in request.components:
        if c.resolved is None:
            continue
        for atom in c.resolved.graph.atoms:
            ref = ExactAtomRef(c.identity, atom.name)
            if ref in request.ownership.generated:
                elements[ref] = atom.element
        for bond in c.resolved.graph.bonds:
            a, b = ExactAtomRef(c.identity, bond.atom1), ExactAtomRef(c.identity, bond.atom2)
            adjacent.setdefault(a, set()).add(b)
            adjacent.setdefault(b, set()).add(a)
    for connection in request.connections:
        if connection.kind == LinkKind.COVALENT:
            a, b = connection.link.component_atom, connection.link.other_atom
            adjacent.setdefault(a, set()).add(b)
            adjacent.setdefault(b, set()).add(a)
    if any(element.upper() not in VDW_RADII_ANGSTROM for element in elements.values()):
        raise ConformerError("conformer-unknown-context-radii")
    if not request.ownership.generated <= elements.keys():
        raise ConformerError("conformer-unknown-generated-element")
    checked, pairs = set(), []
    for a in sorted(request.ownership.generated):
        near = adjacent.get(a, set())
        excluded = {a} | near | {n for b in near for n in adjacent.get(b, set())}
        for b in sorted(elements):
            pair = tuple(sorted((a, b)))
            if b in excluded or pair in checked:
                continue
            checked.add(pair)
            if len(checked) > MAX_PAIR_EVALUATIONS:
                raise ConformerError("conformer-pair-budget-exceeded")
            radius = VDW_RADII_ANGSTROM[elements[a].upper()]+VDW_RADII_ANGSTROM[elements[b].upper()]
            pairs.append((pair[0], pair[1], radius))
    return StericContext(request.digest, radii_digest(), tuple(sorted(pairs)))


def score_sterics(request: ReconstructionRequest, plan: GeometryPlan, candidate: Candidate,
                   context: StericContext) -> StericScore:
    if (context.request_digest != request.digest or context.radii_digest != radii_digest()
            or candidate.request_digest != request.digest or plan.request_digest != request.digest
            or candidate.plan_digest != plan.digest or candidate.strategy != plan.strategy
            or {o.identity for o in candidate.additions} != request.ownership.generated
            or len(candidate.additions) != len(request.ownership.generated)
            or any(o.unit != Unit.ANGSTROM or o.occupancy != 1 for o in candidate.additions)):
        raise ConformerError("conformer-steric-score-binding")
    points = {o.identity: o.angstrom for o in request.observations}
    points.update({o.identity: o.position for o in candidate.additions})
    overlaps = [max(0., radius-math.dist(points[a], points[b])) for a, b, radius in context.pairs]
    return StericScore(math.fsum(x*x for x in overlaps), max(overlaps, default=0.), len(overlaps), context.radii_digest)
